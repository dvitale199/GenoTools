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

"""Run old-vs-new GenoTools parity on a real cohort.

The pytest parity suite (tests/regression/test_parity.py) runs on the small
synthetic dataset. This script is the same idea aimed at a *real* cohort: it runs
the pre-refactor CLI (from .venv-stable) and the new CLI on the same input and
compares their outputs.

Prerequisites
-------------
1. Build the pre-refactor baseline once:
       PYTHON=python3.11 bash tests/scripts/setup_stable_venv.sh origin/main
2. Run this script with the new code installed in the current environment
   (``pip install -e .``).

Usage
-----
    python tests/scripts/run_parity.py --geno /path/to/cohort --out /path/to/workdir
    python tests/scripts/run_parity.py --geno cohort --out work --scenarios qc gwas

``--geno`` is a pfile prefix (without .pgen/.pvar/.psam). The cohort is copied
into ``--out`` before running, so the source data dir is never written to.

Scenarios
---------
    qc             callrate, sex, het, geno, hwe        (IDs + genotype content)
    case_control   case_control                          (IDs + genotype content)
    haplotype      haplotype                             (IDs + genotype content)
    ld             ld 50 5 0.5                            (IDs + genotype content)
    full-pipeline  all_sample + all_variant              (IDs + genotype content)
    gwas           pca + gwas                            (tested-variant set + lambda)
    ancestry       ancestry + all_sample + all_variant   (full report + per-group genotypes)
    related        related                               (pair REL bins + both counts)

The ``ancestry`` scenario is the production configuration: ancestry prediction
splits the cohort and QC then runs per group, so the output is
``{out}_{LABEL}.pgen`` per label rather than a single fileset. It is compared with
the full check battery from ``compare_ancestry_run.py`` - labels, counts, model
quality, QC metrics, pruned samples, related pairs, per-group step pass/fail, and
per-group genotype content.

It is *not* in the default scenario set: it needs a reference panel and takes
roughly an hour per CLI on a 10k cohort, against seconds-to-minutes for the flat
scenarios. Select it explicitly (``--scenarios ancestry``) or via ``all``. Unlike
the flat scenarios it runs without ``--full-output``, because the per-ancestry
``*_hwe_tmp.bed`` intermediates run to ~10 GB per side and the comparison only
needs the promoted ``{out}_{LABEL}.*`` outputs.

The ``related`` scenario exists to exercise ``duplicated_cutoff`` (0.354),
which no parity run had ever reached: the 10k subset holds 52 related pairs and
zero duplicates, so ``filter_relatedness``'s duplicate ``--king-cutoff`` call
and the ``duplicate`` bin of its ``pd.cut`` never saw a positive. It compares
pair-level REL classifications and both counts (related_count and
duplicated_count) with ``compare_related_run.py``'s battery, and *fails* if
duplicated_count is zero - a run where nothing is a duplicate would otherwise
pass every equality check while testing nothing.

It is therefore selected by name only, never via ``all``: it needs a ``--geno``
subset built to contain pairs above the cutoff. To build one, take both members
of pairs a prior release classified ``duplicate`` (from its ``related_samples``
table) and pass the same pairs back via ``--expect-duplicate-pairs`` to assert
they are re-classified as duplicates.

GWAS parity is asserted at the tested-variant-set + genomic-inflation-lambda
level, NOT per-variant p-value: the new PCA pruning intentionally excludes
high-LD/MHC regions the old code left in (ratified "decision B" in
REFACTOR.md), so every GWAS p-value differs slightly by design. The
per-variant mismatch count is still reported, for information.

Exit status is non-zero if any selected scenario fails its parity check.
"""

from __future__ import annotations

import argparse
import platform
import shutil
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from tests.regression.compare import (  # noqa: E402
    compare_pfiles,
    compare_genotypes,
    compare_gwas_results,
    find_gwas_output,
    _lambda_gc,
    _resolve_plink2,
)

# token -> (old_flag, new_flag). Single-word steps are spelled the same.
_STEP_FLAGS = {
    "callrate": ("--callrate", "--callrate"),
    "sex": ("--sex", "--sex"),
    "het": ("--het", "--het"),
    "related": ("--related", "--related"),
    "geno": ("--geno", "--geno"),
    "case_control": ("--case_control", "--case-control"),
    "haplotype": ("--haplotype", "--haplotype"),
    "hwe": ("--hwe", "--hwe"),
    "ld": ("--ld", "--ld"),
    "ancestry": ("--ancestry", "--ancestry"),
    "all_sample": ("--all_sample", "--all-sample"),
    "all_variant": ("--all_variant", "--all-variant"),
    "pca": ("--pca", "--pca"),
    "gwas": ("--gwas", "--gwas"),
}
_STEP_ARGS = {"ld": ["50", "5", "0.5"]}  # old CLI requires all three ld values

SCENARIOS = {
    "qc": ["callrate", "sex", "het", "geno", "hwe"],
    "case_control": ["case_control"],
    "haplotype": ["haplotype"],
    "ld": ["ld"],
    "full-pipeline": ["all_sample", "all_variant"],
    "gwas": ["pca", "gwas"],
    "ancestry": ["ancestry", "all_sample", "all_variant"],
    "related": ["related"],
}
GWAS_SCENARIOS = {"gwas"}
ANCESTRY_SCENARIOS = {"ancestry"}
RELATED_SCENARIOS = {"related"}
DEFAULT_SCENARIOS = ["qc", "full-pipeline", "gwas"]

# "related" is not a default: it asserts that the duplicate branch actually
# fired, which needs a --geno subset containing pairs above duplicated_cutoff
# (0.354). On a subset without any it fails by design, rather than passing
# every equality check while exercising nothing.

# Where genotools caches the GP2 reference panel; overridable on the CLI.
_REF_DIR = Path.home() / ".genotools" / "ref" / "ref_panel"
DEFAULT_REF_PANEL = _REF_DIR / "ref_panel_gp2_prune_rm_underperform_pos_update"
DEFAULT_REF_LABELS = _REF_DIR / "ref_panel_ancestry_updated.txt"

# Ancestry runs an order of magnitude longer than the flat scenarios: ~50 min
# per CLI on 10k samples x 1.9M variants, which the 1-hour flat default would
# cut off mid-run.
FLAT_TIMEOUT = 3600
ANCESTRY_TIMEOUT = 21600


def _build_cmd(
    base: list[str],
    geno: Path,
    out: Path,
    tokens: list[str],
    new: bool,
    ref_panel: Path | None = None,
    ref_labels: Path | None = None,
    full_output: bool = True,
) -> list[str]:
    """Assemble one CLI invocation.

    Args:
        base: Executable prefix, e.g. ["/path/.venv-stable/bin/genotools"].
        geno: Input pfile prefix.
        out: --out prefix for this run.
        tokens: Step tokens from SCENARIOS.
        new: True for the 2.x CLI (hyphenated flags), False for pre-refactor.
        ref_panel: Reference panel prefix; required when "ancestry" is in tokens.
        ref_labels: Reference label file; required when "ancestry" is in tokens.
        full_output: Whether to keep intermediates. Off for ancestry, whose
            per-group intermediates run to ~10 GB per side.
    """
    cmd = list(base) + ["--pfile", str(geno), "--out", str(out)]
    for token in tokens:
        cmd.append(_STEP_FLAGS[token][1 if new else 0])
        cmd.extend(_STEP_ARGS.get(token, []))
    if "ancestry" in tokens:
        if ref_panel is None or ref_labels is None:
            raise ValueError("the ancestry scenario needs ref_panel and ref_labels")
        cmd += [
            "--ref-panel" if new else "--ref_panel", str(ref_panel),
            "--ref-labels" if new else "--ref_labels", str(ref_labels),
        ]
    if full_output:
        cmd.append("--full-output" if new else "--full_output")
    return cmd


def _run(cmd: list[str], cwd: str | None = None, timeout: int = 3600) -> subprocess.CompletedProcess:
    return subprocess.run(cmd, cwd=cwd, capture_output=True, text=True, timeout=timeout)


def _copy_pfileset(src_prefix: Path, dst_prefix: Path) -> None:
    dst_prefix.parent.mkdir(parents=True, exist_ok=True)
    for ext in (".pgen", ".pvar", ".psam"):
        shutil.copy2(src_prefix.with_suffix(ext), dst_prefix.with_suffix(ext))


def _check_king_hint(stable_bin: Path) -> None:
    """Warn if the old baseline will hang fetching KING at import on Linux."""
    if platform.system() != "Linux":
        return
    king = Path.home() / ".genotools" / "misc" / "executables" / "king"
    if not king.exists():
        print(
            "WARNING: on Linux the pre-refactor CLI downloads KING at import "
            "(kingrelatedness.com) and may hang. Re-run setup_stable_venv.sh "
            "(it pre-caches KING) or place a 'king' binary at "
            f"{king}.\n",
            file=sys.stderr,
        )


def _parity_qc(old_out: Path, new_out: Path, plink2: str) -> tuple[bool, str]:
    if not old_out.with_suffix(".pgen").exists():
        return False, "old CLI produced no output"
    if not new_out.with_suffix(".pgen").exists():
        return False, "new CLI produced no output"
    ids = compare_pfiles(old_out, new_out)
    if not ids.equal:
        return False, f"IDs differ: {ids.message}"
    geno = compare_genotypes(old_out, new_out, plink2_exec=plink2)
    if not geno.equal:
        return False, f"genotypes differ: {geno.message}"
    return True, "identical IDs + genotypes"


def _parity_gwas(old_out: Path, new_out: Path) -> tuple[bool, str]:
    old_glm = find_gwas_output(old_out)
    new_glm = find_gwas_output(new_out)
    if old_glm is None or new_glm is None:
        return False, f"missing GWAS output (old={old_glm}, new={new_glm})"

    # Full comparison for the informational per-variant mismatch count.
    detail = compare_gwas_results(old_glm, new_glm)

    import pandas as pd

    def _ids_and_lambda(path: Path) -> tuple[set, float]:
        df = pd.read_csv(path, sep=r"\s+", dtype={"#CHROM": str}, low_memory=False)
        if "TEST" in df.columns:
            df = df[df["TEST"] == "ADD"]
        return set(df["ID"]), _lambda_gc(df["P"])

    old_ids, lam_old = _ids_and_lambda(old_glm)
    new_ids, lam_new = _ids_and_lambda(new_glm)

    if old_ids != new_ids:
        return False, (
            f"different tested-variant sets (old-only {len(old_ids - new_ids)}, "
            f"new-only {len(new_ids - old_ids)})"
        )
    if abs(lam_old - lam_new) > 0.05:
        return False, f"lambda differs: old={lam_old:.5f}, new={lam_new:.5f}"

    return True, (
        f"same tested-variant set; lambda old={lam_old:.5f} new={lam_new:.5f}; "
        f"per-variant p-mismatches={len(detail.mismatched_variants)} "
        f"(expected under decision B: PCA excludes MHC/high-LD)"
    )


def _parity_ancestry(old_out: Path, new_out: Path, plink2: str,
                     release_json: Path | None = None,
                     keep: Path | None = None) -> tuple[bool, str]:
    """Compare two ancestry runs with compare_ancestry_run.py's check battery.

    The ancestry configuration has no single output fileset to diff: prediction
    splits the cohort and QC runs per group. Rather than reimplement that
    comparison, reuse the checks that already back compare_ancestry_run.py so
    the two entry points cannot drift apart.
    """
    # Imported here, not at module scope: it pulls in pandas and is needed only
    # when this scenario is selected.
    from tests.scripts.compare_ancestry_run import (
        _PASS_FAIL_SUFFIX,
        Report,
        _load_json,
        check_counts,
        check_labels,
        check_model_quality,
        check_output_pfiles,
        check_pass_fail,
        check_pruned_samples,
        check_qc_metrics,
        check_related,
        cross_check_release,
    )

    for prefix, side in ((old_out, "old"), (new_out, "new")):
        if not Path(f"{prefix}.json").exists():
            return False, f"{side} CLI produced no JSON report ({prefix}.json)"

    old, new = _load_json(old_out), _load_json(new_out)

    rep = Report()
    check_labels(rep, old, new)
    check_counts(rep, old, new)
    check_model_quality(rep, old, new)
    check_qc_metrics(rep, old, new)
    check_pruned_samples(rep, old, new)
    check_related(rep, old, new)
    check_pass_fail(rep, old, new)

    labels = sorted(
        k[: -len(_PASS_FAIL_SUFFIX)]
        for k in set(old) | set(new)
        if k.endswith(_PASS_FAIL_SUFFIX)
    )
    if labels:
        check_output_pfiles(rep, old_out, new_out, labels, plink2)
    else:
        rep.add("output pfiles", False, "could not infer ancestry labels from the reports")

    if release_json:
        cross_check_release(old, new, release_json, keep)

    failed = rep.failed()
    if failed:
        return False, f"{len(failed)} check(s) diverged: {', '.join(failed)}"
    return True, f"all {len(rep.rows)} checks identical across {len(labels)} ancestry group(s)"


def _parity_related(old_out: Path, new_out: Path, plink2: str,
                    expect_duplicate_pairs: Path | None = None) -> tuple[bool, str]:
    """Compare two --related runs with compare_related_run.py's check battery.

    Relatedness reports more than an output fileset: pair-level REL bins and
    two separate counts, one of which (duplicated_count) went unexercised for
    the whole refactor. Reuse the checks that back compare_related_run.py so
    the two entry points cannot drift apart.
    """
    from tests.scripts.compare_related_run import (
        _cmp,
        check_duplicate_branch_fired,
        check_expected_duplicates,
        check_output_pfiles,
        check_pruned_samples,
        check_rel_classification,
        check_related_metrics,
    )

    for prefix, side in ((old_out, "old"), (new_out, "new")):
        if not Path(f"{prefix}.json").exists():
            return False, f"{side} CLI produced no JSON report ({prefix}.json)"

    old, new = _cmp._load_json(old_out), _cmp._load_json(new_out)

    rep = _cmp.Report()
    check_duplicate_branch_fired(rep, new)
    check_related_metrics(rep, old, new)
    check_rel_classification(rep, old, new)
    check_pruned_samples(rep, old, new)
    if expect_duplicate_pairs:
        check_expected_duplicates(rep, new, expect_duplicate_pairs)
    check_output_pfiles(rep, old_out, new_out, plink2)

    failed = rep.failed()
    if failed:
        return False, f"{len(failed)} check(s) diverged: {', '.join(failed)}"
    return True, f"all {len(rep.rows)} checks identical, duplicate branch exercised"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--geno", required=True, type=Path, help="Real cohort pfile prefix (no extension)")
    parser.add_argument("--out", required=True, type=Path, help="Working directory for the two runs")
    parser.add_argument(
        "--stable-bin",
        type=Path,
        default=REPO_ROOT / ".venv-stable" / "bin" / "genotools",
        help="Path to the pre-refactor genotools (default: .venv-stable/bin/genotools)",
    )
    parser.add_argument(
        "--scenarios",
        nargs="+",
        choices=[*SCENARIOS.keys(), "all"],
        default=DEFAULT_SCENARIOS,
        help=f"Scenarios to run (default: {' '.join(DEFAULT_SCENARIOS)})",
    )
    parser.add_argument(
        "--ref-panel", type=Path, default=DEFAULT_REF_PANEL,
        help=f"Ancestry reference panel prefix (default: {DEFAULT_REF_PANEL})",
    )
    parser.add_argument(
        "--ref-labels", type=Path, default=DEFAULT_REF_LABELS,
        help=f"Ancestry reference labels (default: {DEFAULT_REF_LABELS})",
    )
    parser.add_argument(
        "--expect-duplicate-pairs", type=Path, default=None,
        help="TSV of IID1 IID2 pairs a prior release called duplicates "
             "(related scenario only)",
    )
    parser.add_argument(
        "--release-json", type=Path, default=None,
        help="Released full-cohort JSON for an informational cross-check (ancestry only)",
    )
    parser.add_argument(
        "--keep", type=Path, default=None,
        help="Keep file used to build the subset, scoping the release cross-check",
    )
    parser.add_argument(
        "--timeout", type=int, default=None,
        help=f"Per-CLI timeout in seconds (default: {FLAT_TIMEOUT}, "
             f"or {ANCESTRY_TIMEOUT} when an ancestry scenario is selected)",
    )
    args = parser.parse_args()

    # "all" does not pull in "related": that scenario asserts the duplicate
    # branch fired, which needs input built to contain duplicate pairs. Under
    # "all" on an ordinary cohort it would fail for the wrong reason, so it is
    # selected by name only.
    if "all" in args.scenarios:
        scenarios = [n for n in SCENARIOS if n not in RELATED_SCENARIOS]
    else:
        scenarios = args.scenarios
    wants_ancestry = bool(ANCESTRY_SCENARIOS.intersection(scenarios))
    if args.timeout is None:
        args.timeout = ANCESTRY_TIMEOUT if wants_ancestry else FLAT_TIMEOUT

    if not args.geno.with_suffix(".pgen").exists():
        parser.error(f"input not found: {args.geno}.pgen")
    if not args.stable_bin.exists():
        parser.error(
            f"pre-refactor CLI not found at {args.stable_bin}. Build it with: "
            "PYTHON=python3.11 bash tests/scripts/setup_stable_venv.sh origin/main"
        )
    if wants_ancestry:
        if not Path(f"{args.ref_panel}.bed").exists():
            parser.error(
                f"reference panel not found: {args.ref_panel}.bed. Fetch it with "
                "genotools-download, or pass --ref-panel."
            )
        if not args.ref_labels.exists():
            parser.error(f"reference labels not found: {args.ref_labels}")
    plink2 = _resolve_plink2()
    if plink2 is None:
        parser.error("plink2 not found (needed for genotype comparison)")

    _check_king_hint(args.stable_bin)

    args.out.mkdir(parents=True, exist_ok=True)
    # Copy the cohort so GWAS's {input}.pheno side-effect never touches source data.
    geno = args.out / "input"
    print(f"Copying cohort {args.geno} -> {geno}.*")
    _copy_pfileset(args.geno, geno)

    old_base = [str(args.stable_bin)]
    new_base = [sys.executable, "-m", "genotools"]

    results: list[tuple[str, bool, str]] = []
    for name in scenarios:
        tokens = SCENARIOS[name]
        scen_dir = args.out / name
        scen_dir.mkdir(parents=True, exist_ok=True)
        old_out = scen_dir / "old"
        new_out = scen_dir / "new"

        # Ancestry keeps no intermediates: ~10 GB per side of per-group
        # *_hwe_tmp.bed that the comparison never reads.
        full_output = name not in ANCESTRY_SCENARIOS
        cmd_kwargs = dict(
            ref_panel=args.ref_panel, ref_labels=args.ref_labels, full_output=full_output
        )

        print(f"\n=== scenario: {name}  ({' '.join(tokens)}) ===")
        old_res = _run(
            _build_cmd(old_base, geno, old_out, tokens, new=False, **cmd_kwargs),
            timeout=args.timeout,
        )
        new_res = _run(
            _build_cmd(new_base, geno, new_out, tokens, new=True, **cmd_kwargs),
            cwd=str(REPO_ROOT), timeout=args.timeout,
        )

        if name in RELATED_SCENARIOS:
            ok, msg = _parity_related(
                old_out, new_out, plink2, args.expect_duplicate_pairs
            )
        elif name in ANCESTRY_SCENARIOS:
            ok, msg = _parity_ancestry(
                old_out, new_out, plink2, args.release_json, args.keep
            )
        elif name in GWAS_SCENARIOS:
            ok, msg = _parity_gwas(old_out, new_out)
        else:
            ok, msg = _parity_qc(old_out, new_out, plink2)

        if not ok and (old_res.returncode != 0 or new_res.returncode != 0):
            msg += (
                f"\n    old rc={old_res.returncode} new rc={new_res.returncode}"
                f"\n    old stderr tail: {old_res.stderr.strip()[-300:]}"
                f"\n    new stderr tail: {new_res.stderr.strip()[-300:]}"
            )
        results.append((name, ok, msg))
        print(f"    {'PASS' if ok else 'FAIL'}: {msg}")

    print("\n" + "=" * 70)
    print("PARITY SUMMARY")
    print("=" * 70)
    for name, ok, msg in results:
        print(f"  [{'PASS' if ok else 'FAIL'}] {name}: {msg.splitlines()[0]}")

    failed = [n for n, ok, _ in results if not ok]
    if failed:
        print(f"\n{len(failed)} scenario(s) FAILED parity: {', '.join(failed)}")
        return 1
    print(f"\nAll {len(results)} scenario(s) passed parity.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
