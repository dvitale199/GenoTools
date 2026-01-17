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
