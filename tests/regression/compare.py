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


def _ext(prefix: Path, ext: str) -> Path:
    """Append an extension to a fileset prefix.

    Path.with_suffix() *replaces* the last dotted segment, so a legitimate
    prefix like "release.r12" or "cohort.v2/out" would lose part of its name.
    PLINK appends to whatever --out it is given, so this does too.
    """
    return prefix.parent / f"{prefix.name}{ext}"


def _read_psam(path: Path) -> pd.DataFrame:
    """Read a .psam, normalizing the leading '#' off the first column name.

    PLINK2 writes either "#FID<TAB>IID..." or, when there is no FID column,
    "#IID...". Stripping the '#' leaves an "IID" column either way.
    """
    df = pd.read_csv(path, sep="\t", dtype=str)
    if df.columns[0].startswith("#"):
        df = df.rename(columns={df.columns[0]: df.columns[0][1:]})
    return df


def _read_pvar(path: Path) -> pd.DataFrame:
    """Read a .pvar, skipping VCF-style ## meta lines and normalizing #CHROM."""
    with open(path) as fh:
        skip_lines = 0
        for line in fh:
            if line.startswith("##"):
                skip_lines += 1
            else:
                break
    df = pd.read_csv(path, sep="\t", dtype=str, skiprows=skip_lines)
    if "#CHROM" in df.columns:
        df = df.rename(columns={"#CHROM": "CHROM"})
    return df


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
        exp_df = _read_psam(expected)
        act_df = _read_psam(actual)
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
        exp_df = _read_pvar(expected)
        act_df = _read_pvar(actual)
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
        _ext(expected_prefix, ".psam"),
        _ext(actual_prefix, ".psam")
    )
    variant_result = compare_variant_ids(
        _ext(expected_prefix, ".pvar"),
        _ext(actual_prefix, ".pvar")
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
    """Locate a plink2 executable (PATH, then the genotools cache).

    Note this can differ from the plink2 GenoTools itself runs, which always
    comes from the genotools cache. The comparison operations used here are
    stable across versions, but flag semantics are not: --maj-ref writes
    different alleles under v2.00a5.10 than under v2.0.0-a.6.3, so a test that
    depends on a flag's effect should build its fixture explicitly rather than
    letting whichever plink2 is resolved decide what the fixture means.
    """
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
    traw = _ext(out_prefix, ".traw")
    if result.returncode != 0 or not traw.exists():
        raise RuntimeError(
            f"plink2 --export A-transpose failed for {prefix} "
            f"(returncode {result.returncode}):\n{result.stderr}"
        )
    return traw


# ---------------------------------------------------------------------------
# GWAS result comparison
#
# GWAS output is NOT a pfile: it is a PLINK2 --glm association table
# ({out}.PHENO1.glm.logistic.hybrid or .glm.linear). GWAS was broken until the
# last hardening PR, so old-vs-new GWAS parity had never been checked. These
# helpers align the two result tables on variant ID and compare per-variant
# p-values (and the genomic-inflation lambda derived from them) within a small
# tolerance.
# ---------------------------------------------------------------------------


def find_gwas_output(out_prefix: Path) -> Optional[Path]:
    """Locate a PLINK2 --glm output file for a run output prefix.

    Matches both logistic (``.glm.logistic.hybrid``) and linear
    (``.glm.linear``) results, regardless of phenotype name.

    Args:
        out_prefix: Run output prefix (without extension).

    Returns:
        Path to the glm results file, or None if not found.
    """
    out_prefix = Path(out_prefix)
    patterns = ["*.glm.logistic.hybrid", "*.glm.linear"]
    for pattern in patterns:
        matches = sorted(out_prefix.parent.glob(f"{out_prefix.name}{pattern}"))
        if matches:
            return matches[0]
    return None


def _lambda_gc(pvalues) -> float:
    """Compute the genomic inflation factor (lambda GC) from p-values.

    lambda_gc = median(chi2) / expected_median, with chi2 the 1-df chi-squared
    quantiles of (1 - p). Matches the approach used by the production inflation
    code so a comparison here is meaningful.
    """
    from scipy.stats import chi2

    p = pd.to_numeric(pd.Series(pvalues), errors="coerce").to_numpy()
    p = p[~pd.isna(p)]
    if len(p) == 0:
        return float("nan")
    chi2_stats = chi2.isf(p, 1)
    expected_median = chi2.ppf(0.5, 1)
    import numpy as np

    return float(np.nanmedian(chi2_stats) / expected_median)


def compare_gwas_results(
    expected: Path,
    actual: Path,
    p_tolerance: float = 1e-6,
    lambda_tolerance: float = 1e-3,
) -> ComparisonResult:
    """Compare two PLINK2 --glm association tables, order-independently.

    Aligns on variant ID (and TEST, when present), then compares per-variant
    p-values within ``p_tolerance`` and the derived lambda GC within
    ``lambda_tolerance``. ``mismatched_variants`` lists variants whose p-value
    differs beyond tolerance.

    Args:
        expected: Path to the expected .glm.* file.
        actual: Path to the actual .glm.* file.
        p_tolerance: Absolute tolerance for per-variant p-value comparison.
        lambda_tolerance: Absolute tolerance for lambda GC comparison.

    Returns:
        ComparisonResult over the association tables.
    """
    import numpy as np

    exp = pd.read_csv(expected, sep=r"\s+", dtype={"#CHROM": str}, low_memory=False)
    act = pd.read_csv(actual, sep=r"\s+", dtype={"#CHROM": str}, low_memory=False)

    # Restrict to the additive test if a TEST column is present.
    if "TEST" in exp.columns:
        exp = exp[exp["TEST"] == "ADD"]
    if "TEST" in act.columns:
        act = act[act["TEST"] == "ADD"]

    exp = exp.set_index("ID")
    act = act.set_index("ID")

    exp_ids, act_ids = set(exp.index), set(act.index)
    missing = exp_ids - act_ids
    extra = act_ids - exp_ids
    common = sorted(exp_ids & act_ids)

    mismatched_variants: list[str] = []
    if common:
        e = exp.loc[common]
        a = act.loc[common]
        e_p = pd.to_numeric(e["P"], errors="coerce").to_numpy()
        a_p = pd.to_numeric(a["P"], errors="coerce").to_numpy()
        both_nan = pd.isna(e_p) & pd.isna(a_p)
        # Treat matching NaNs as equal; compare finite p-values within tolerance.
        diff = ~both_nan & ~(np.abs(e_p - a_p) <= p_tolerance)
        mismatched_variants = [common[i] for i in np.where(diff)[0]]

    lambda_exp = _lambda_gc(exp["P"]) if "P" in exp.columns else float("nan")
    lambda_act = _lambda_gc(act["P"]) if "P" in act.columns else float("nan")
    lambda_ok = (
        (pd.isna(lambda_exp) and pd.isna(lambda_act))
        or abs(lambda_exp - lambda_act) <= lambda_tolerance
    )

    equal = (
        len(missing) == 0
        and len(extra) == 0
        and len(mismatched_variants) == 0
        and lambda_ok
    )

    message = (
        f"variants: missing {len(missing)}, extra {len(extra)}, "
        f"p-mismatches {len(mismatched_variants)}; "
        f"lambda expected {lambda_exp:.5f} vs actual {lambda_act:.5f} "
        f"(ok={lambda_ok})"
    )

    return ComparisonResult(
        equal=equal,
        sample_diff=0,
        variant_diff=len(missing) + len(extra) + len(mismatched_variants),
        mismatched_samples=[],
        mismatched_variants=sorted(missing | extra) + mismatched_variants,
        message=message,
    )


def compare_gwas(expected_prefix: Path, actual_prefix: Path, **kwargs) -> ComparisonResult:
    """Compare GWAS output for two run prefixes.

    Locates each run's .glm.* file and compares them via
    :func:`compare_gwas_results`.

    Args:
        expected_prefix: Expected run output prefix (without extension).
        actual_prefix: Actual run output prefix (without extension).
        **kwargs: Forwarded to compare_gwas_results (tolerances).

    Returns:
        ComparisonResult, non-equal with an explanatory message if either
        glm file is missing.
    """
    exp_glm = find_gwas_output(Path(expected_prefix))
    act_glm = find_gwas_output(Path(actual_prefix))

    if exp_glm is None or act_glm is None:
        return ComparisonResult(
            equal=False,
            sample_diff=0,
            variant_diff=0,
            mismatched_samples=[],
            mismatched_variants=[],
            message=(
                f"GWAS output missing: expected={exp_glm}, actual={act_glm}"
            ),
        )

    return compare_gwas_results(exp_glm, act_glm, **kwargs)


def _run_pgen_diff(
    expected_prefix: Path,
    actual_prefix: Path,
    plink2_exec: str,
    out_prefix: Path,
) -> tuple[int, set]:
    """Report genotype differences between two pfile sets via --pgen-diff.

    Unlike a .traw export this needs no full-matrix intermediate: PLINK2 walks
    both filesets and writes one row per differing (variant, sample) pair.

    Returns:
        (number of differing calls, set of variant IDs involved).

    Raises:
        RuntimeError: If PLINK2 fails.
    """
    import subprocess

    cmd = [
        plink2_exec,
        "--pfile", str(expected_prefix),
        # include-missing so a call present on one side and missing on the
        # other counts as a difference, matching the .traw comparison.
        "--pgen-diff", str(actual_prefix), "include-missing",
        "--out", str(out_prefix),
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    pdiff = _ext(out_prefix, ".pdiff")
    if result.returncode != 0:
        raise RuntimeError(
            f"plink2 --pgen-diff failed for {expected_prefix} vs "
            f"{actual_prefix} (returncode {result.returncode}):\n"
            f"{result.stderr}"
        )
    if not pdiff.exists():
        # PLINK2 writes a header-only .pdiff when there is nothing to report;
        # a missing file means it decided there was nothing to compare.
        return 0, set()

    ids: set = set()
    count = 0
    with open(pdiff) as fh:
        header = fh.readline()
        if not header.startswith("#ID"):
            raise RuntimeError(
                f"unexpected .pdiff header in {pdiff}: {header!r}"
            )
        for line in fh:
            if not line.strip():
                continue
            count += 1
            ids.add(line.split("\t", 1)[0])
    return count, ids


def compare_genotypes(
    expected_prefix: Path,
    actual_prefix: Path,
    plink2_exec: Optional[str] = None,
    method: str = "pgen-diff",
) -> ComparisonResult:
    """Compare the genotype content of two pfile sets.

    Detects same-ID-but-different-genotype and flipped-coding regressions that
    compare_pfiles cannot see. Requires a PLINK2 executable.

    Two methods, same verdict:

    ``pgen-diff`` (default)
        PLINK2 ``--pgen-diff`` walks both filesets in place and reports
        differing calls, plus a direct REF/ALT comparison read from the .pvar
        files. Needs no large intermediate, so it scales to real cohorts.
    ``traw``
        Exports both filesets to .traw text and diffs the matrices. Equivalent,
        but the export costs roughly (samples x variants x 2) bytes per side -
        tens of GB for a real cohort - so it is kept only for debugging, where
        the full matrix is what you want to look at.

    The .pvar allele comparison runs first because PLINK2 *refuses* to run
    --pgen-diff when the two filesets disagree on a variant's REF allele. On a
    coding mismatch the per-call diff is skipped and the coding difference is
    reported, which keeps a flipped fileset a clear failure rather than a
    PLINK error.

    Args:
        expected_prefix: Path prefix for expected pfiles (without extension).
        actual_prefix: Path prefix for actual pfiles (without extension).
        plink2_exec: Path to plink2; auto-resolved from PATH/genotools cache
            if omitted.
        method: "pgen-diff" (default) or "traw".

    Returns:
        ComparisonResult over the genotype content.

    Raises:
        RuntimeError: If plink2 cannot be found or a PLINK2 call fails.
        ValueError: If method is not recognized.
    """
    import tempfile

    expected_prefix = Path(expected_prefix)
    actual_prefix = Path(actual_prefix)

    if method not in ("pgen-diff", "traw"):
        raise ValueError(
            f"unknown method {method!r}; expected 'pgen-diff' or 'traw'"
        )

    plink2 = plink2_exec or _resolve_plink2()
    if plink2 is None:
        raise RuntimeError(
            "plink2 executable not found; pass plink2_exec=... explicitly"
        )

    if method == "traw":
        with tempfile.TemporaryDirectory() as td:
            tdp = Path(td)
            exp_traw = _export_traw(expected_prefix, plink2, tdp / "expected")
            act_traw = _export_traw(actual_prefix, plink2, tdp / "actual")
            return _compare_traw_files(exp_traw, act_traw)

    # ID sets: --pgen-diff only compares the intersection, so the sets have to
    # be checked separately or a dropped sample/variant would go unnoticed.
    exp_psam = _read_psam(_ext(expected_prefix, ".psam"))
    act_psam = _read_psam(_ext(actual_prefix, ".psam"))
    exp_ids, act_ids = set(exp_psam["IID"]), set(act_psam["IID"])
    missing_ids, extra_ids = exp_ids - act_ids, act_ids - exp_ids

    exp_pvar = _read_pvar(_ext(expected_prefix, ".pvar"))
    act_pvar = _read_pvar(_ext(actual_prefix, ".pvar"))
    exp_snps, act_snps = set(exp_pvar["ID"]), set(act_pvar["ID"])
    missing_snps, extra_snps = exp_snps - act_snps, act_snps - exp_snps

    # Allele coding on shared variants. --pgen-diff compares alleles, not
    # dosage codes, so a REF/ALT flip is invisible to it; compare directly.
    mismatched_variants: set = set()
    coding_mismatches = 0
    if {"REF", "ALT"} <= set(exp_pvar.columns) & set(act_pvar.columns):
        e = exp_pvar.set_index("ID")[["REF", "ALT"]]
        a = act_pvar.set_index("ID")[["REF", "ALT"]]
        shared = e.index.intersection(a.index)
        if len(shared):
            e, a = e.loc[shared], a.loc[shared]
            flipped = (e["REF"].values != a["REF"].values) | (
                e["ALT"].values != a["ALT"].values
            )
            coding_mismatches = int(flipped.sum())
            mismatched_variants.update(shared[flipped])

    if coding_mismatches:
        # PLINK2 refuses to run --pgen-diff when REF alleles conflict ("REF
        # allele on line N ... conflicts with loaded REF"), so the per-call
        # diff cannot run. The coding mismatch is already a failing verdict.
        genotype_mismatches = None
    else:
        with tempfile.TemporaryDirectory() as td:
            genotype_mismatches, diff_ids = _run_pgen_diff(
                expected_prefix, actual_prefix, plink2, Path(td) / "pdiff"
            )
        mismatched_variants.update(diff_ids)

    sample_diff = len(missing_ids) + len(extra_ids)
    id_variant_diff = len(missing_snps) + len(extra_snps)

    equal = (
        sample_diff == 0
        and id_variant_diff == 0
        and genotype_mismatches == 0
        and coding_mismatches == 0
    )
    calls = (
        "not compared (allele coding differs)"
        if genotype_mismatches is None
        else str(genotype_mismatches)
    )
    message = (
        f"variants: missing {len(missing_snps)}, extra {len(extra_snps)}; "
        f"samples: missing {len(missing_ids)}, extra {len(extra_ids)}; "
        f"coding mismatches {coding_mismatches}; "
        f"genotype cell mismatches {calls}"
    )

    return ComparisonResult(
        equal=equal,
        sample_diff=sample_diff,
        variant_diff=id_variant_diff + len(mismatched_variants),
        mismatched_samples=sorted(missing_ids | extra_ids),
        mismatched_variants=sorted(mismatched_variants)
        + sorted(missing_snps | extra_snps),
        message=message,
    )
