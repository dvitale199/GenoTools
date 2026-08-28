#!/usr/bin/env python
# Copyright 2023 The GenoTools Authors. All Rights Reserved.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
# ==============================================================================

"""Compare two full ``--ancestry --all-sample --all-variant`` GenoTools runs.

``run_parity.py`` covers QC / full-pipeline / GWAS on a single flat cohort. This
script handles the *production* configuration, where ancestry prediction splits
the cohort and QC then runs per ancestry group. It compares the two runs' JSON
reports field by field and the per-ancestry output pfiles by genotype content.

Both CLIs emit the same JSON schema, so the comparison is symmetric.

Usage
-----
    python tests/scripts/compare_ancestry_run.py \
        --old  /path/work/old \
        --new  /path/work/new \
        --release-json /path/GP2_r12_final_post_genotools.json \
        --keep /path/subset10k_seed42.keep

``--old`` / ``--new`` are the ``--out`` prefixes handed to each CLI.

What is compared (old vs new)
-----------------------------
    ancestry labels     per-sample predicted label            exact
    ancestry counts     per-label totals                      exact
    model quality       test_accuracy, confusion matrix        exact
    QC metrics          pruned_count per (step, ancestry)      exact
    pruned samples      the actual pruned IDs, per step        exact
    related samples     related pairs per ancestry            exact
    step pass/fail      per-ancestry step status               exact
    output genotypes    {prefix}_{label}.pgen                  IDs + genotypes

``--release-json`` additionally cross-checks both runs against the released
full-cohort run: label concordance, and whether the samples that run pruned at
``callrate`` are pruned again here. Per-sample callrate is independent of cohort
composition, so those should reproduce exactly; ``sex``/``het`` and every variant
step depend on group composition and will legitimately differ on a subset.

Exit status is non-zero if any old-vs-new check diverges. The release
cross-check is informational and never fails the run.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from tests.regression.compare import (  # noqa: E402
    compare_genotypes,
    compare_pfiles,
    _resolve_plink2,
)

# JSON keys holding a per-ancestry step pass/fail map.
_PASS_FAIL_SUFFIX = "_pass_fail"


class Report:
    """Collects check outcomes and prints a summary."""

    def __init__(self) -> None:
        self.rows: list[tuple[str, bool, str]] = []

    def add(self, name: str, ok: bool, msg: str) -> None:
        self.rows.append((name, ok, msg))
        print(f"  [{'PASS' if ok else 'FAIL'}] {name}: {msg}")

    def failed(self) -> list[str]:
        return [n for n, ok, _ in self.rows if not ok]


def _load_json(prefix: Path) -> dict:
    # str-concat, not with_suffix: an --out prefix may legitimately contain a dot.
    path = Path(f"{prefix}.json")
    if not path.exists():
        raise SystemExit(f"ERROR: no JSON report at {path} — did that run finish?")
    with open(path) as f:
        return json.load(f)


def _frame(d: dict, key: str) -> pd.DataFrame | None:
    """Rebuild a DataFrame from the JSON's ``DataFrame.to_dict()`` payload."""
    if key not in d or d[key] is None:
        return None
    df = pd.DataFrame(d[key])
    if key == "QC":
        # Reports written by the refactor before the schema was restored named
        # this column "count" instead of "pruned_count", and emitted a
        # duplicate generic "outlier_count" row alongside each variant step's
        # own metric. Normalize both so old and new reports are comparable.
        if "count" in df.columns and "pruned_count" not in df.columns:
            df = df.rename(columns={"count": "pruned_count"})
        if {"level", "metric", "step"} <= set(df.columns):
            dupes = df[
                (df["level"] == "variant") & (df["metric"] == "outlier_count")
            ]
            if not dupes.empty:
                # Only drop the generic row where the step also reported its
                # own metric, so a step with nothing else is still represented.
                own = set(
                    df[
                        (df["level"] == "variant")
                        & (df["metric"] != "outlier_count")
                    ][["step", "ancestry"]].itertuples(index=False, name=None)
                )
                drop = dupes[
                    dupes[["step", "ancestry"]]
                    .apply(tuple, axis=1)
                    .isin(own)
                ].index
                df = df.drop(index=drop).reset_index(drop=True)
    return df


def _labels(d: dict) -> pd.Series | None:
    df = _frame(d, "ancestry_labels")
    if df is None or "IID" not in df.columns:
        return None
    return df.set_index(df["IID"].astype(str))["label"]


def check_labels(rep: Report, old: dict, new: dict) -> pd.DataFrame | None:
    """Per-sample predicted ancestry must agree exactly."""
    lo, ln = _labels(old), _labels(new)
    if lo is None or ln is None:
        rep.add("ancestry labels", False, "missing ancestry_labels in one report")
        return None

    if set(lo.index) != set(ln.index):
        rep.add(
            "ancestry labels",
            False,
            f"different sample sets (old {len(lo):,}, new {len(ln):,}, "
            f"old-only {len(set(lo.index) - set(ln.index)):,}, "
            f"new-only {len(set(ln.index) - set(lo.index)):,})",
        )
        return None

    joined = pd.DataFrame({"old": lo, "new": ln.reindex(lo.index)})
    disagree = joined[joined["old"] != joined["new"]]
    if disagree.empty:
        rep.add("ancestry labels", True, f"all {len(joined):,} samples got the same label")
    else:
        rep.add(
            "ancestry labels",
            False,
            f"{len(disagree):,}/{len(joined):,} samples labeled differently",
        )
        print("\n    old-vs-new label disagreement:")
        print(
            "    "
            + pd.crosstab(disagree["old"], disagree["new"]).to_string().replace("\n", "\n    ")
        )
    return joined


def check_counts(rep: Report, old: dict, new: dict) -> None:
    co, cn = _frame(old, "ancestry_counts"), _frame(new, "ancestry_counts")
    if co is None or cn is None:
        rep.add("ancestry counts", False, "missing ancestry_counts in one report")
        return
    o = co.set_index("label")["count"].sort_index()
    n = cn.set_index("label")["count"].sort_index()
    if o.equals(n):
        rep.add("ancestry counts", True, f"identical across {len(o)} labels")
    else:
        merged = pd.DataFrame({"old": o, "new": n}).fillna(0)
        rep.add("ancestry counts", False, "per-label counts differ")
        print("    " + merged[merged["old"] != merged["new"]].to_string().replace("\n", "\n    "))


def check_model_quality(rep: Report, old: dict, new: dict) -> None:
    ao, an = old.get("test_accuracy"), new.get("test_accuracy")
    if ao is None or an is None:
        rep.add("test accuracy", False, f"missing (old={ao}, new={an})")
    elif ao == an:
        rep.add("test accuracy", True, f"identical ({ao:.6f})")
    else:
        rep.add("test accuracy", False, f"old={ao:.6f} new={an:.6f} (diff {abs(ao - an):.2e})")

    mo, mn = _frame(old, "confusion_matrix"), _frame(new, "confusion_matrix")
    if mo is None or mn is None:
        rep.add("confusion matrix", False, "missing in one report")
        return
    mo, mn = mo.sort_index(axis=0).sort_index(axis=1), mn.sort_index(axis=0).sort_index(axis=1)
    if mo.equals(mn):
        rep.add("confusion matrix", True, f"identical ({mo.shape[0]}x{mo.shape[1]})")
    else:
        rep.add("confusion matrix", False, "differs")
        print("    old:\n    " + mo.to_string().replace("\n", "\n    "))
        print("    new:\n    " + mn.to_string().replace("\n", "\n    "))


def check_qc_metrics(rep: Report, old: dict, new: dict) -> None:
    """pruned_count must match for every (step, ancestry, metric)."""
    qo, qn = _frame(old, "QC"), _frame(new, "QC")
    if qo is None or qn is None:
        rep.add("QC metrics", False, "missing QC block in one report")
        return

    keys = ["step", "ancestry", "metric"]
    o = qo.set_index(keys)["pruned_count"].sort_index()
    n = qn.set_index(keys)["pruned_count"].sort_index()

    if set(o.index) != set(n.index):
        rep.add(
            "QC metrics",
            False,
            f"different (step, ancestry, metric) rows — old {len(o)}, new {len(n)}",
        )
        for label, missing in (("old-only", set(o.index) - set(n.index)),
                               ("new-only", set(n.index) - set(o.index))):
            if missing:
                print(f"    {label}: {sorted(missing)[:10]}")
        return

    merged = pd.DataFrame({"old": o, "new": n.reindex(o.index)})
    diff = merged[merged["old"] != merged["new"]]
    if diff.empty:
        rep.add("QC metrics", True, f"all {len(merged)} pruned_counts identical")
    else:
        rep.add("QC metrics", False, f"{len(diff)}/{len(merged)} pruned_counts differ")
        print("    " + diff.to_string().replace("\n", "\n    "))


def check_pruned_samples(rep: Report, old: dict, new: dict) -> None:
    """The actual pruned sample IDs must match, per step."""
    po, pn = _frame(old, "pruned_samples"), _frame(new, "pruned_samples")
    if po is None or pn is None:
        rep.add("pruned samples", False, "missing pruned_samples in one report")
        return

    def keyset(df: pd.DataFrame) -> dict[str, set[str]]:
        out: dict[str, set[str]] = {}
        for step, grp in df.groupby("step"):
            out[str(step)] = set(grp["IID"].astype(str))
        return out

    ko, kn = keyset(po), keyset(pn)
    steps = sorted(set(ko) | set(kn))
    problems = []
    for step in steps:
        so, sn = ko.get(step, set()), kn.get(step, set())
        if so != sn:
            problems.append(f"{step}: old-only {len(so - sn)}, new-only {len(sn - so)}")
    if problems:
        rep.add("pruned samples", False, "; ".join(problems))
    else:
        total = sum(len(v) for v in ko.values())
        rep.add("pruned samples", True, f"identical IDs across {len(steps)} step(s), {total:,} total")


def _scope_note(message: str) -> str:
    """Surface a narrowed genotype comparison in the PASS line.

    compare_genotypes falls back to autosomes when unset sex codes stop PLINK2
    from comparing chrX/chrY. A pass that covered less than everything has to
    say so, or a silently narrowed check reads as full coverage.
    """
    start = message.find("(autosomes only")
    return "" if start < 0 else " " + message[start:]


def check_related(rep: Report, old: dict, new: dict) -> None:
    ro, rn = _frame(old, "related_samples"), _frame(new, "related_samples")
    if ro is None and rn is None:
        rep.add("related samples", True, "neither run reported related pairs")
        return
    if ro is None or rn is None:
        rep.add("related samples", False, f"present in only one run (old={ro is not None}, new={rn is not None})")
        return

    def pairset(df: pd.DataFrame) -> set[tuple[str, str, str]]:
        c1 = "IID1" if "IID1" in df.columns else df.columns[1]
        c2 = "IID2" if "IID2" in df.columns else df.columns[3]
        anc = df["ancestry"] if "ancestry" in df.columns else pd.Series([""] * len(df))
        return {
            (str(a), *sorted((str(x), str(y))))  # unordered pair
            for a, x, y in zip(anc, df[c1], df[c2])
        }

    so, sn = pairset(ro), pairset(rn)
    if so == sn:
        rep.add("related samples", True, f"identical {len(so):,} pairs")
    else:
        rep.add("related samples", False, f"old-only {len(so - sn):,}, new-only {len(sn - so):,}")


def check_pass_fail(rep: Report, old: dict, new: dict) -> None:
    keys = sorted(k for k in set(old) | set(new) if k.endswith(_PASS_FAIL_SUFFIX))
    if not keys:
        rep.add("step pass/fail", False, "no *_pass_fail blocks found")
        return
    problems = []
    for key in keys:
        o, n = old.get(key, {}), new.get(key, {})
        for step in sorted(set(o) | set(n)):
            so = o.get(step, {}).get("status")
            sn = n.get(step, {}).get("status")
            if so != sn:
                problems.append(f"{key.replace(_PASS_FAIL_SUFFIX, '')}/{step}: old={so} new={sn}")
    if problems:
        shown = "; ".join(problems[:8])
        if len(problems) > 8:
            shown += f" (+{len(problems) - 8} more of {len(problems)} total)"
        rep.add("step pass/fail", False, shown)
    else:
        rep.add("step pass/fail", True, f"identical across {len(keys)} ancestry group(s)")


def check_output_pfiles(rep: Report, old_prefix: Path, new_prefix: Path,
                        labels: list[str], plink2: str) -> None:
    """Per-ancestry final pfiles must match on IDs and genotype content."""
    checked, missing = [], []
    for label in labels:
        o = old_prefix.parent / f"{old_prefix.name}_{label}"
        n = new_prefix.parent / f"{new_prefix.name}_{label}"
        if not o.with_suffix(".pgen").exists() or not n.with_suffix(".pgen").exists():
            missing.append(
                f"{label}(old={o.with_suffix('.pgen').exists()},new={n.with_suffix('.pgen').exists()})"
            )
            continue

        ids = compare_pfiles(o, n)
        if not ids.equal:
            rep.add(f"pfiles {label}", False, f"IDs differ: {ids.message}")
            continue
        # One group failing to compare should not abort the remaining groups
        # (or the release cross-check) the way an uncaught error does.
        try:
            geno = compare_genotypes(o, n, plink2_exec=plink2)
        except RuntimeError as exc:
            rep.add(
                f"pfiles {label}",
                False,
                f"genotype comparison could not run: {exc}\n"
                f"    IDs matched, so this is a comparison failure rather "
                f"than a parity failure. To reproduce:\n"
                f"    plink2 --pfile {o} --pgen-diff {n} include-missing "
                f"--out /tmp/pdiff_{label}",
            )
            continue
        if not geno.equal:
            rep.add(f"pfiles {label}", False, f"genotypes differ: {geno.message}")
            continue
        rep.add(
            f"pfiles {label}", True,
            f"identical IDs + genotypes{_scope_note(geno.message)}",
        )
        checked.append(label)

    if missing:
        rep.add("pfiles present", False, "output missing for: " + ", ".join(missing))
    elif checked:
        rep.add("pfiles present", True, f"compared {len(checked)} ancestry group(s)")


def cross_check_release(old: dict, new: dict, release_path: Path, keep: Path | None) -> None:
    """Informational: how the two runs line up with the released full-cohort run."""
    print(f"\n{'=' * 72}\nRELEASE CROSS-CHECK (informational — does not gate)\n{'=' * 72}")
    with open(release_path) as f:
        rel = json.load(f)

    rel_labels = pd.DataFrame(rel["ancestry_labels"])
    rel_map = rel_labels.set_index(rel_labels["IID"].astype(str))["label"]

    for name, d in (("old", old), ("new", new)):
        lab = _labels(d)
        if lab is None:
            print(f"  {name}: no ancestry_labels to compare")
            continue
        shared = lab.index.intersection(rel_map.index)
        agree = (lab.reindex(shared) == rel_map.reindex(shared)).sum()
        print(
            f"  {name}: {agree:,}/{len(shared):,} "
            f"({100 * agree / max(len(shared), 1):.2f}%) samples match the released label"
        )
        mism = lab.reindex(shared)[lab.reindex(shared) != rel_map.reindex(shared)]
        if len(mism):
            ct = pd.crosstab(rel_map.reindex(mism.index), mism)
            ct.index.name = "released"
            ct.columns.name = name
            print("    " + ct.to_string().replace("\n", "\n    "))

    # Per-sample callrate is cohort-composition independent, so the released
    # callrate outliers that are in this subset should be pruned again here.
    rel_pruned = pd.DataFrame(rel["pruned_samples"])
    rel_pruned["IID"] = rel_pruned["IID"].astype(str)
    subset_ids = None
    if keep is not None and keep.exists():
        subset_ids = set(pd.read_csv(keep, sep=r"\s+")["IID"].astype(str))

    rel_cr = set(rel_pruned[rel_pruned["step"] == "callrate"]["IID"])
    if subset_ids:
        rel_cr &= subset_ids
    print(f"\n  released callrate outliers inside this subset: {len(rel_cr):,}")
    for name, d in (("old", old), ("new", new)):
        p = _frame(d, "pruned_samples")
        if p is None:
            print(f"    {name}: no pruned_samples")
            continue
        got = set(p[p["step"] == "callrate"]["IID"].astype(str))
        print(
            f"    {name}: pruned {len(got):,} at callrate; "
            f"{len(rel_cr & got):,}/{len(rel_cr):,} of the released ones reproduced"
            + (f", {len(rel_cr - got):,} not re-pruned" if rel_cr - got else "")
        )


def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--old", required=True, type=Path, help="Old-CLI --out prefix")
    parser.add_argument("--new", required=True, type=Path, help="New-CLI --out prefix")
    parser.add_argument("--release-json", type=Path, default=None,
                        help="Released full-cohort JSON for an informational cross-check")
    parser.add_argument("--keep", type=Path, default=None,
                        help="Keep file used to build the subset (scopes the release cross-check)")
    parser.add_argument("--skip-genotypes", action="store_true",
                        help="Skip the output pfile check entirely - sample/variant "
                             "IDs and the plink2 --pgen-diff genotype comparison (faster)")
    args = parser.parse_args()

    old, new = _load_json(args.old), _load_json(args.new)

    plink2 = _resolve_plink2()
    if plink2 is None and not args.skip_genotypes:
        parser.error("plink2 not found (needed for the genotype diff); pass --skip-genotypes to skip")

    print("=" * 72)
    print("OLD vs NEW — ancestry + full QC")
    print("=" * 72)

    rep = Report()
    check_labels(rep, old, new)
    check_counts(rep, old, new)
    check_model_quality(rep, old, new)
    check_qc_metrics(rep, old, new)
    check_pruned_samples(rep, old, new)
    check_related(rep, old, new)
    check_pass_fail(rep, old, new)

    labels = sorted(
        k[: -len(_PASS_FAIL_SUFFIX)] for k in set(old) | set(new) if k.endswith(_PASS_FAIL_SUFFIX)
    )
    if args.skip_genotypes:
        print("  [SKIP] output pfiles: --skip-genotypes")
    elif labels:
        check_output_pfiles(rep, args.old, args.new, labels, plink2)
    else:
        rep.add("output pfiles", False, "could not infer ancestry labels from the reports")

    if args.release_json:
        cross_check_release(old, new, args.release_json, args.keep)

    print("\n" + "=" * 72)
    print("SUMMARY")
    print("=" * 72)
    failed = rep.failed()
    for name, ok, msg in rep.rows:
        print(f"  [{'PASS' if ok else 'FAIL'}] {name}: {msg.splitlines()[0]}")
    if failed:
        print(f"\n{len(failed)} check(s) FAILED: {', '.join(failed)}")
        return 1
    print(f"\nAll {len(rep.rows)} check(s) passed — old and new agree.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
