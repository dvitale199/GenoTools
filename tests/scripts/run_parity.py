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
    "geno": ("--geno", "--geno"),
    "case_control": ("--case_control", "--case-control"),
    "haplotype": ("--haplotype", "--haplotype"),
    "hwe": ("--hwe", "--hwe"),
    "ld": ("--ld", "--ld"),
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
}
GWAS_SCENARIOS = {"gwas"}
DEFAULT_SCENARIOS = ["qc", "full-pipeline", "gwas"]


def _build_cmd(base: list[str], geno: Path, out: Path, tokens: list[str], new: bool) -> list[str]:
    cmd = list(base) + ["--pfile", str(geno), "--out", str(out)]
    for token in tokens:
        cmd.append(_STEP_FLAGS[token][1 if new else 0])
        cmd.extend(_STEP_ARGS.get(token, []))
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
    parser.add_argument("--timeout", type=int, default=3600, help="Per-CLI timeout in seconds")
    args = parser.parse_args()

    scenarios = list(SCENARIOS) if "all" in args.scenarios else args.scenarios

    if not args.geno.with_suffix(".pgen").exists():
        parser.error(f"input not found: {args.geno}.pgen")
    if not args.stable_bin.exists():
        parser.error(
            f"pre-refactor CLI not found at {args.stable_bin}. Build it with: "
            "PYTHON=python3.11 bash tests/scripts/setup_stable_venv.sh origin/main"
        )
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

        print(f"\n=== scenario: {name}  ({' '.join(tokens)}) ===")
        old_res = _run(_build_cmd(old_base, geno, old_out, tokens, new=False), timeout=args.timeout)
        new_res = _run(_build_cmd(new_base, geno, new_out, tokens, new=True), cwd=str(REPO_ROOT), timeout=args.timeout)

        if name in GWAS_SCENARIOS:
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
