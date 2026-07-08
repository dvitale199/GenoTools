"""Utilities for comparing genomic outputs in regression testing."""
from dataclasses import dataclass
from pathlib import Path
from typing import Optional
import pandas as pd
import json


@dataclass
class ComparisonResult:
    """Result of comparing two outputs."""
    equal: bool
    sample_diff: int
    variant_diff: int
    mismatched_samples: list[str]
    mismatched_variants: list[str]
    message: str


def compare_sample_ids(expected: Path, actual: Path) -> ComparisonResult:
    """
    Compare sample IDs between two psam/fam files.

    Args:
        expected: Path to expected psam/fam file
        actual: Path to actual psam/fam file

    Returns:
        ComparisonResult with details about differences
    """
    ext = expected.suffix
    if ext == ".psam":
        exp_df = pd.read_csv(expected, sep="\t", dtype=str)
        act_df = pd.read_csv(actual, sep="\t", dtype=str)
        exp_ids = set(exp_df["IID"])
        act_ids = set(act_df["IID"])
    else:  # .fam
        exp_df = pd.read_csv(expected, sep=r"\s+", header=None, dtype=str)
        act_df = pd.read_csv(actual, sep=r"\s+", header=None, dtype=str)
        exp_ids = set(exp_df[1])
        act_ids = set(act_df[1])

    missing = exp_ids - act_ids
    extra = act_ids - exp_ids

    return ComparisonResult(
        equal=len(missing) == 0 and len(extra) == 0,
        sample_diff=len(missing) + len(extra),
        variant_diff=0,
        mismatched_samples=sorted(list(missing | extra)),
        mismatched_variants=[],
        message=f"Missing: {len(missing)}, Extra: {len(extra)}"
    )


def compare_variant_ids(expected: Path, actual: Path) -> ComparisonResult:
    """
    Compare variant IDs between two pvar/bim files.

    Args:
        expected: Path to expected pvar/bim file
        actual: Path to actual pvar/bim file

    Returns:
        ComparisonResult with details about differences
    """
    ext = expected.suffix
    if ext == ".pvar":
        # PVAR files have header starting with #CHROM - we need to keep it
        # but skip any ## comment lines
        def read_pvar(path):
            # Skip lines starting with ## (VCF-style comments)
            with open(path) as f:
                skip_lines = 0
                for line in f:
                    if line.startswith("##"):
                        skip_lines += 1
                    else:
                        break
            df = pd.read_csv(path, sep="\t", dtype=str, skiprows=skip_lines)
            # Rename #CHROM to CHROM for easier access
            if "#CHROM" in df.columns:
                df = df.rename(columns={"#CHROM": "CHROM"})
            return df

        exp_df = read_pvar(expected)
        act_df = read_pvar(actual)
        exp_ids = set(exp_df["ID"])
        act_ids = set(act_df["ID"])
    else:  # .bim
        exp_df = pd.read_csv(expected, sep="\t", header=None, dtype=str)
        act_df = pd.read_csv(actual, sep="\t", header=None, dtype=str)
        exp_ids = set(exp_df[1])
        act_ids = set(act_df[1])

    missing = exp_ids - act_ids
    extra = act_ids - exp_ids

    return ComparisonResult(
        equal=len(missing) == 0 and len(extra) == 0,
        sample_diff=0,
        variant_diff=len(missing) + len(extra),
        mismatched_samples=[],
        mismatched_variants=sorted(list(missing | extra)),
        message=f"Missing: {len(missing)}, Extra: {len(extra)}"
    )


def compare_pfiles(expected_prefix: Path, actual_prefix: Path) -> ComparisonResult:
    """
    Compare two pfile sets (samples + variants).

    Args:
        expected_prefix: Path prefix for expected pfiles (without extension)
        actual_prefix: Path prefix for actual pfiles (without extension)

    Returns:
        ComparisonResult combining sample and variant differences
    """
    expected_prefix = Path(expected_prefix)
    actual_prefix = Path(actual_prefix)

    sample_result = compare_sample_ids(
        expected_prefix.with_suffix(".psam"),
        actual_prefix.with_suffix(".psam")
    )
    variant_result = compare_variant_ids(
        expected_prefix.with_suffix(".pvar"),
        actual_prefix.with_suffix(".pvar")
    )

    return ComparisonResult(
        equal=sample_result.equal and variant_result.equal,
        sample_diff=sample_result.sample_diff,
        variant_diff=variant_result.variant_diff,
        mismatched_samples=sample_result.mismatched_samples,
        mismatched_variants=variant_result.mismatched_variants,
        message=f"Samples: {sample_result.message}; Variants: {variant_result.message}"
    )


def compare_outlier_files(expected: Path, actual: Path) -> ComparisonResult:
    """
    Compare two outlier files (tab-separated with #FID header).

    Args:
        expected: Path to expected outliers file
        actual: Path to actual outliers file

    Returns:
        ComparisonResult with details about differences
    """
    exp_df = pd.read_csv(expected, sep="\t", dtype=str)
    act_df = pd.read_csv(actual, sep="\t", dtype=str)

    # Get sample IDs
    exp_ids = set(exp_df["IID"])
    act_ids = set(act_df["IID"])

    missing = exp_ids - act_ids
    extra = act_ids - exp_ids

    return ComparisonResult(
        equal=len(missing) == 0 and len(extra) == 0,
        sample_diff=len(missing) + len(extra),
        variant_diff=0,
        mismatched_samples=sorted(list(missing | extra)),
        mismatched_variants=[],
        message=f"Missing outliers: {len(missing)}, Extra outliers: {len(extra)}"
    )


def compare_metrics(
    expected: dict,
    actual: dict,
    numeric_tolerance: float = 0.0
) -> tuple[bool, str]:
    """
    Compare two metrics dictionaries.

    Args:
        expected: Expected metrics dictionary
        actual: Actual metrics dictionary
        numeric_tolerance: Tolerance for numeric comparisons

    Returns:
        Tuple of (equal, message)
    """
    differences = []

    for key in expected:
        if key not in actual:
            differences.append(f"Missing key: {key}")
            continue

        exp_val = expected[key]
        act_val = actual[key]

        if isinstance(exp_val, (int, float)) and isinstance(act_val, (int, float)):
            if abs(exp_val - act_val) > numeric_tolerance:
                differences.append(f"{key}: expected {exp_val}, got {act_val}")
        elif exp_val != act_val:
            differences.append(f"{key}: expected {exp_val}, got {act_val}")

    # Check for extra keys
    for key in actual:
        if key not in expected:
            differences.append(f"Extra key: {key}")

    if differences:
        return False, "; ".join(differences)
    return True, "Metrics match"


def compare_ancestry_predictions(
    expected: Path,
    actual: Path,
    tolerance: float = 0.0
) -> ComparisonResult:
    """
    Compare ancestry prediction outputs.

    Args:
        expected: Path to expected predictions file
        actual: Path to actual predictions file
        tolerance: Tolerance for probability comparisons (not used for label matching)

    Returns:
        ComparisonResult with details about prediction differences
    """
    exp_df = pd.read_csv(expected, sep="\t")
    act_df = pd.read_csv(actual, sep="\t")

    # Merge on sample ID
    merged = exp_df.merge(act_df, on="IID", suffixes=("_exp", "_act"))

    # Check predicted ancestry matches
    label_col_exp = "predicted_ancestry_exp" if "predicted_ancestry_exp" in merged.columns else "label_exp"
    label_col_act = "predicted_ancestry_act" if "predicted_ancestry_act" in merged.columns else "label_act"

    mismatched = merged[merged[label_col_exp] != merged[label_col_act]]

    return ComparisonResult(
        equal=len(mismatched) == 0,
        sample_diff=len(mismatched),
        variant_diff=0,
        mismatched_samples=list(mismatched["IID"]),
        mismatched_variants=[],
        message=f"{len(mismatched)} samples have different predictions"
    )


def load_golden_metrics(golden_dir: Path, step_name: str) -> Optional[dict]:
    """
    Load golden metrics for a specific step.

    Args:
        golden_dir: Path to golden files directory
        step_name: Name of the QC step

    Returns:
        Metrics dictionary or None if not found
    """
    metrics_path = golden_dir / step_name / "metrics.json"
    if not metrics_path.exists():
        return None

    with open(metrics_path) as f:
        return json.load(f)


# ---------------------------------------------------------------------------
# Genotype-content comparison
#
# compare_pfiles only diffs sample/variant ID *sets*. For true old-vs-new
# parity we must also confirm the surviving genotype matrix is identical --
# a regression that keeps the same IDs but corrupts genotypes or flips allele
# coding would otherwise pass. These helpers export pfiles to PLINK2 .traw
# (variant x sample additive-dosage matrix) and compare order-independently.
# ---------------------------------------------------------------------------

_TRAW_META_COLS = ["CHR", "SNP", "(C)M", "POS", "COUNTED", "ALT"]


def _compare_traw_files(expected: Path, actual: Path) -> ComparisonResult:
    """Compare two PLINK2 .traw genotype matrices, order-independently.

    Aligns on variant ID (SNP) and sample column, then compares allele coding
    (COUNTED/ALT) and per-cell genotype dosages. Matching missing values (NA in
    the same cell of both files) count as equal.

    Args:
        expected: Path to the expected .traw file.
        actual: Path to the actual .traw file.

    Returns:
        ComparisonResult. `mismatched_variants` lists variants that differ in
        coding or genotype; `variant_diff` aggregates ID-set and content diffs.
    """
    exp = pd.read_csv(expected, sep="\t", low_memory=False)
    act = pd.read_csv(actual, sep="\t", low_memory=False)

    exp_samples = [c for c in exp.columns if c not in _TRAW_META_COLS]
    act_samples = [c for c in act.columns if c not in _TRAW_META_COLS]

    exp = exp.set_index("SNP")
    act = act.set_index("SNP")

    exp_snps, act_snps = set(exp.index), set(act.index)
    missing_snps = exp_snps - act_snps
    extra_snps = act_snps - exp_snps
    common_snps = sorted(exp_snps & act_snps)

    exp_ids, act_ids = set(exp_samples), set(act_samples)
    missing_ids = exp_ids - act_ids
    extra_ids = act_ids - exp_ids
    common_ids = sorted(exp_ids & act_ids)

    mismatched_variants: set = set()
    genotype_mismatches = 0
    coding_mismatches = 0

    if common_snps:
        e = exp.loc[common_snps]
        a = act.loc[common_snps]

        # Allele-coding differences (COUNTED / ALT) per variant
        coding_diff = (
            (e["COUNTED"].astype(str).values != a["COUNTED"].astype(str).values)
            | (e["ALT"].astype(str).values != a["ALT"].astype(str).values)
        )
        coding_mismatches = int(coding_diff.sum())
        for i, snp in enumerate(common_snps):
            if coding_diff[i]:
                mismatched_variants.add(snp)

        # Per-cell genotype differences on shared samples
        if common_ids:
            eg = e[common_ids]
            ag = a[common_ids]
            neq = eg.values != ag.values
            both_nan = pd.isna(eg.values) & pd.isna(ag.values)
            cell_diff = neq & ~both_nan
            genotype_mismatches = int(cell_diff.sum())
            row_has_diff = cell_diff.any(axis=1)
            for i, snp in enumerate(common_snps):
                if row_has_diff[i]:
                    mismatched_variants.add(snp)

    sample_diff = len(missing_ids) + len(extra_ids)
    id_variant_diff = len(missing_snps) + len(extra_snps)

    equal = (
        sample_diff == 0
        and id_variant_diff == 0
        and genotype_mismatches == 0
        and coding_mismatches == 0
    )

    message = (
        f"variants: missing {len(missing_snps)}, extra {len(extra_snps)}; "
        f"samples: missing {len(missing_ids)}, extra {len(extra_ids)}; "
        f"coding mismatches {coding_mismatches}; "
        f"genotype cell mismatches {genotype_mismatches}"
    )

    return ComparisonResult(
        equal=equal,
        sample_diff=sample_diff,
        variant_diff=id_variant_diff + len(mismatched_variants),
        mismatched_samples=sorted(missing_ids | extra_ids),
        mismatched_variants=sorted(mismatched_variants) + sorted(missing_snps | extra_snps),
        message=message,
    )


def _resolve_plink2() -> Optional[str]:
    """Locate a plink2 executable (PATH, then the genotools cache)."""
    import shutil

    found = shutil.which("plink2")
    if found:
        return found
    candidate = Path.home() / ".genotools" / "misc" / "executables" / "plink2"
    return str(candidate) if candidate.exists() else None


def _export_traw(prefix: Path, plink2_exec: str, out_prefix: Path) -> Path:
    """Export a pfile set to a .traw genotype matrix via PLINK2."""
    import subprocess

    cmd = [
        plink2_exec,
        "--pfile", str(prefix),
        "--export", "A-transpose",
        "--out", str(out_prefix),
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    traw = out_prefix.with_suffix(".traw")
    if result.returncode != 0 or not traw.exists():
        raise RuntimeError(
            f"plink2 --export A-transpose failed for {prefix} "
            f"(returncode {result.returncode}):\n{result.stderr}"
        )
    return traw


def compare_genotypes(
    expected_prefix: Path,
    actual_prefix: Path,
    plink2_exec: Optional[str] = None,
) -> ComparisonResult:
    """Compare the genotype content of two pfile sets.

    Exports both to PLINK2 .traw and compares order-independently. Detects
    same-ID-but-different-genotype / flipped-coding regressions that
    compare_pfiles cannot see. Requires a PLINK2 executable.

    Args:
        expected_prefix: Path prefix for expected pfiles (without extension).
        actual_prefix: Path prefix for actual pfiles (without extension).
        plink2_exec: Path to plink2; auto-resolved from PATH/genotools cache
            if omitted.

    Returns:
        ComparisonResult over the genotype matrices.

    Raises:
        RuntimeError: If plink2 cannot be found or an export fails.
    """
    import tempfile

    expected_prefix = Path(expected_prefix)
    actual_prefix = Path(actual_prefix)

    plink2 = plink2_exec or _resolve_plink2()
    if plink2 is None:
        raise RuntimeError(
            "plink2 executable not found; pass plink2_exec=... explicitly"
        )

    with tempfile.TemporaryDirectory() as td:
        tdp = Path(td)
        exp_traw = _export_traw(expected_prefix, plink2, tdp / "expected")
        act_traw = _export_traw(actual_prefix, plink2, tdp / "actual")
        return _compare_traw_files(exp_traw, act_traw)
