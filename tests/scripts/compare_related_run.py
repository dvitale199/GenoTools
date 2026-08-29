#!/usr/bin/env python3
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

"""Compare two flat ``--related`` runs, and prove the duplicate branch ran.

``duplicated_cutoff=0.354`` was never exercised: the 10k parity subset holds 52
related pairs and zero duplicates, so ``filter_relatedness``'s duplicate
``--king-cutoff`` call and the ``duplicate`` bin of its ``pd.cut`` never saw a
real positive. This compares old vs new on a subset built to contain known
duplicate pairs, and fails if the branch did *not* fire - a run where nothing
is a duplicate would otherwise pass every equality check while testing nothing.

Reuses ``compare_ancestry_run``'s Report and pfile comparators rather than
reimplementing them.

Usage:
    compare_related_run.py --old WORK/old --new WORK/new \\
        [--expect-duplicate-pairs FILE] [--skip-genotypes]
"""

from __future__ import annotations

import argparse
import importlib.util
import shutil
import sys
from pathlib import Path

import pandas as pd

_SIBLING = Path(__file__).resolve().parent / "compare_ancestry_run.py"
_spec = importlib.util.spec_from_file_location("compare_ancestry_run", _SIBLING)
assert _spec is not None and _spec.loader is not None
_cmp = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_cmp)

Report = _cmp.Report
_load_json = _cmp._load_json
_frame = _cmp._frame


def _scope_note(message: str) -> str:
    """Surface a narrowed genotype comparison in the PASS line.

    compare_genotypes falls back to autosomes when unset sex codes stop PLINK2
    from comparing chrX/chrY. A pass that covered less than everything has to
    say so, or a silently narrowed check reads as full coverage.
    """
    start = message.find("(autosomes only")
    return "" if start < 0 else " " + message[start:]


def _related_frame(d: dict) -> pd.DataFrame | None:
    return _frame(d, "related_samples")


def _pairs(df: pd.DataFrame) -> dict[tuple[str, str], str]:
    """Map each unordered pair to its REL classification."""
    c1 = "IID1" if "IID1" in df.columns else df.columns[1]
    c2 = "IID2" if "IID2" in df.columns else df.columns[3]
    rel = df["REL"] if "REL" in df.columns else pd.Series([""] * len(df))
    out: dict[tuple[str, str], str] = {}
    for x, y, r in zip(df[c1], df[c2], rel):
        key = tuple(sorted((str(x), str(y))))  # type: ignore[assignment]
        out[key] = str(r)
    return out


def check_related_metrics(rep: Report, old: dict, new: dict) -> None:
    """related_count and duplicated_count, per step row."""
    qo, qn = _frame(old, "QC"), _frame(new, "QC")
    if qo is None or qn is None:
        rep.add("QC metrics", False, f"missing QC block (old={qo is not None}, new={qn is not None})")
        return

    def rows(df: pd.DataFrame) -> dict[tuple[str, str], int]:
        return {
            (str(s), str(m)): int(c)
            for s, m, c in zip(df["step"], df["metric"], df["pruned_count"])
        }

    ro, rn = rows(qo), rows(qn)
    if set(ro) != set(rn):
        rep.add("QC metrics", False, f"different rows: old-only {sorted(set(ro) - set(rn))}, new-only {sorted(set(rn) - set(ro))}")
        return
    diffs = [f"{k}: old={ro[k]} new={rn[k]}" for k in sorted(ro) if ro[k] != rn[k]]
    if diffs:
        rep.add("QC metrics", False, "; ".join(diffs))
    else:
        rep.add("QC metrics", True, f"all {len(ro)} count(s) identical: " + ", ".join(f"{m}={ro[(s, m)]}" for s, m in sorted(ro)))


def check_duplicate_branch_fired(rep: Report, new: dict) -> None:
    """The point of the exercise: a zero here means nothing was tested."""
    qn = _frame(new, "QC")
    if qn is None:
        rep.add("duplicate branch fired", False, "no QC block")
        return
    dup = [
        int(c) for s, m, c in zip(qn["step"], qn["metric"], qn["pruned_count"])
        if str(m) == "duplicated_count"
    ]
    if not dup:
        rep.add("duplicate branch fired", False, "no duplicated_count row - was --related run?")
    elif sum(dup) == 0:
        rep.add(
            "duplicate branch fired", False,
            "duplicated_count=0: the subset holds no duplicates, so "
            "duplicated_cutoff is still unexercised",
        )
    else:
        rep.add("duplicate branch fired", True, f"duplicated_count={sum(dup)} samples removed at the 0.354 cutoff")


def check_rel_classification(rep: Report, old: dict, new: dict) -> None:
    """Every pair, and its REL bin - the 0.354 boundary decides 'duplicate'."""
    ro, rn = _related_frame(old), _related_frame(new)
    if ro is None or rn is None:
        rep.add("related pairs", False, f"present in only one run (old={ro is not None}, new={rn is not None})")
        return

    po, pn = _pairs(ro), _pairs(rn)
    if set(po) != set(pn):
        rep.add("related pairs", False, f"old-only {len(set(po) - set(pn))}, new-only {len(set(pn) - set(po))}")
        return

    mismatched = [f"{a}/{b}: old={po[(a, b)]} new={pn[(a, b)]}" for a, b in sorted(po) if po[(a, b)] != pn[(a, b)]]
    if mismatched:
        shown = "; ".join(mismatched[:6])
        if len(mismatched) > 6:
            shown += f" (+{len(mismatched) - 6} more)"
        rep.add("related pairs", False, f"REL differs: {shown}")
        return

    counts = pd.Series(list(pn.values())).value_counts().to_dict()
    breakdown = ", ".join(f"{k}={v}" for k, v in sorted(counts.items()))
    rep.add("related pairs", True, f"identical {len(pn):,} pairs and REL bins ({breakdown})")


def check_expected_duplicates(rep: Report, new: dict, expect_file: Path) -> None:
    """Pairs a prior release called duplicates must come back as duplicates."""
    want = set()
    for line in expect_file.read_text().splitlines():
        parts = line.split()
        if len(parts) >= 2 and not line.startswith("#"):
            want.add(tuple(sorted((parts[0], parts[1]))))
    if not want:
        rep.add("expected duplicates", False, f"no pairs read from {expect_file}")
        return

    rn = _related_frame(new)
    if rn is None:
        rep.add("expected duplicates", False, "new run reported no related pairs")
        return
    got = _pairs(rn)

    missing = sorted(p for p in want if p not in got)
    misclassified = sorted(
        f"{a}/{b}={got[(a, b)]}" for a, b in want
        if (a, b) in got and got[(a, b)] != "duplicate"
    )
    if missing or misclassified:
        detail = []
        if missing:
            detail.append(f"not reported at all: {missing[:5]}")
        if misclassified:
            detail.append(f"not classified duplicate: {misclassified[:5]}")
        rep.add("expected duplicates", False, "; ".join(detail))
    else:
        rep.add("expected duplicates", True, f"all {len(want)} known duplicate pairs classified 'duplicate'")


def check_pruned_samples(rep: Report, old: dict, new: dict) -> None:
    po, pn = _frame(old, "pruned_samples"), _frame(new, "pruned_samples")
    if po is None and pn is None:
        rep.add("pruned samples", True, "neither run pruned samples")
        return
    if po is None or pn is None:
        rep.add("pruned samples", False, f"present in only one run (old={po is not None}, new={pn is not None})")
        return

    def idset(df: pd.DataFrame) -> set[str]:
        col = "IID" if "IID" in df.columns else df.columns[1]
        return {str(v) for v in df[col]}

    so, sn = idset(po), idset(pn)
    if so == sn:
        rep.add("pruned samples", True, f"identical {len(so)} sample(s) removed")
    else:
        rep.add("pruned samples", False, f"old-only {sorted(so - sn)[:5]}, new-only {sorted(sn - so)[:5]}")


def check_output_pfiles(rep: Report, old_prefix: Path, new_prefix: Path, plink2: str | None) -> None:
    if not (old_prefix.with_suffix(".pgen").exists() and new_prefix.with_suffix(".pgen").exists()):
        rep.add("output pfiles", False, f"missing (old={old_prefix.with_suffix('.pgen').exists()}, new={new_prefix.with_suffix('.pgen').exists()})")
        return
    ids = _cmp.compare_pfiles(old_prefix, new_prefix)
    if not ids.equal:
        rep.add("output pfiles", False, f"IDs differ: {ids.message}")
        return
    if plink2 is None:
        rep.add("output pfiles", True, f"identical IDs ({ids.message}); genotypes skipped")
        return
    try:
        geno = _cmp.compare_genotypes(old_prefix, new_prefix, plink2_exec=plink2)
    except RuntimeError as exc:
        rep.add("output pfiles", False, f"genotype comparison could not run: {exc}")
        return
    rep.add(
        "output pfiles",
        geno.equal,
        f"identical IDs + genotypes{_scope_note(geno.message)}"
        if geno.equal
        else f"genotypes differ: {geno.message}",
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--old", type=Path, required=True, help="Output prefix of the 1.x run")
    parser.add_argument("--new", type=Path, required=True, help="Output prefix of the 2.0 run")
    parser.add_argument("--expect-duplicate-pairs", type=Path, default=None,
                        help="TSV of IID1 IID2 pairs a prior release called duplicates")
    parser.add_argument("--skip-genotypes", action="store_true",
                        help="Compare output pfile IDs only, skipping the plink2 --pgen-diff check")
    args = parser.parse_args()

    old, new = _load_json(args.old), _load_json(args.new)
    plink2 = None
    if not args.skip_genotypes:
        # Same order GenoTools itself resolves in: its own executables dir
        # first, PATH second, so the diff uses the plink2 that wrote the files.
        bundled = Path.home() / ".genotools/misc/executables/plink2"
        plink2 = str(bundled) if bundled.exists() else shutil.which("plink2")
    if not args.skip_genotypes and not plink2:
        parser.error("plink2 not found (needed for the genotype diff); pass --skip-genotypes to skip")

    print("=" * 72)
    print("OLD vs NEW - flat --related run (duplicate branch)")
    print("=" * 72)

    rep = Report()
    check_duplicate_branch_fired(rep, new)
    check_related_metrics(rep, old, new)
    check_rel_classification(rep, old, new)
    check_pruned_samples(rep, old, new)
    if args.expect_duplicate_pairs:
        check_expected_duplicates(rep, new, args.expect_duplicate_pairs)
    check_output_pfiles(rep, args.old, args.new, plink2)

    print("\n" + "=" * 72)
    print("SUMMARY")
    print("=" * 72)
    for name, ok, msg in rep.rows:
        print(f"  [{'PASS' if ok else 'FAIL'}] {name}: {msg.splitlines()[0]}")
    failed = rep.failed()
    if failed:
        print(f"\n{len(failed)} check(s) FAILED: {', '.join(failed)}")
        return 1
    print(f"\nAll {len(rep.rows)} check(s) passed - old and new agree, "
          f"and the duplicate branch ran.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
