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

"""GWAS association analysis functions.

This module provides pure functions for genome-wide association analysis:
- Genomic inflation (lambda) calculation
- GWAS execution using PLINK2 --glm

Usage:
    from genotools.gwas.steps.association import run_gwas, calculate_inflation

    # Run GWAS with PCA covariates
    gwas_result = run_gwas(data, config, out_path, covar_path=eigenvec_path)

    # Calculate inflation from p-values
    inflation = calculate_inflation(pvalues, n_cases=500, n_controls=1000)
"""

import logging
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
from scipy.stats import ncx2

from genotools.core.exceptions import GWASError, ValidationError
from genotools.core.executors import get_plink2, run_command
from genotools.core.genotypes import GenotypeData
from genotools.core.logging import step_context
from genotools.gwas.config import GWASConfig, CovariateConfig
from genotools.gwas.results import GWASResult, InflationMetrics

logger = logging.getLogger(__name__)


def calculate_inflation(
    pvalues: np.ndarray,
    n_cases: Optional[int] = None,
    n_controls: Optional[int] = None,
) -> InflationMetrics:
    """Calculate genomic inflation factor (lambda GC).

    Computes the genomic inflation factor using the chi-squared
    approach: lambda = median(chi2) / expected_median.

    Args:
        pvalues: Array of p-values from association tests.
        n_cases: Number of cases (for sample-size normalization).
        n_controls: Number of controls (for sample-size normalization).

    Returns:
        InflationMetrics with lambda_gc and optionally lambda_1000
        if case/control counts are provided.

    Note:
        Lambda values close to 1.0 indicate well-controlled inflation.
        Values significantly > 1.0 suggest population stratification
        or other systematic bias.
    """
    # Filter out NaN p-values
    pvalues = pvalues[~np.isnan(pvalues)]

    if len(pvalues) == 0:
        return InflationMetrics(
            lambda_gc=np.nan,
            lambda_1000=np.nan,
            n_cases=n_cases,
            n_controls=n_controls,
        )

    # Convert p-values to chi-squared statistics
    # Using non-central chi-squared with df=1, nc=0
    chi2_stats = ncx2.ppf(1 - pvalues, 1, nc=0)
    expected_median = ncx2.ppf(0.5, 1, nc=0)

    # Calculate lambda GC
    lambda_gc = float(np.nanmedian(chi2_stats) / expected_median)

    # Calculate sample-size normalized lambda (lambda_1000)
    lambda_1000: Optional[float] = None
    if n_cases is not None and n_controls is not None and n_cases > 0 and n_controls > 0:
        # Normalize to 1000 cases + 1000 controls
        # Formula: lambda_1000 = 1 + (lambda - 1) * (1/n_cases + 1/n_controls) / (1/1000 + 1/1000)
        lambda_1000 = float(
            1 + (lambda_gc - 1) * (1/n_cases + 1/n_controls) / (1/1000 + 1/1000)
        )

    return InflationMetrics(
        lambda_gc=lambda_gc,
        lambda_1000=lambda_1000,
        n_cases=n_cases,
        n_controls=n_controls,
    )


def _extract_covariate_names(covar_path: Path) -> str:
    """Extract covariate column names from a covariate file.

    Reads the header of the covariate file and returns all column
    names except #FID/#IID (or just #IID) as a comma-separated string.

    Args:
        covar_path: Path to covariate file.

    Returns:
        Comma-separated string of covariate column names.
    """
    # Read just the header
    header_df = pd.read_csv(covar_path, sep=r"\s+", header=0, nrows=0)
    columns = header_df.columns.tolist()

    # Remove ID columns
    if "#FID" in columns:
        # Account for #FID and IID columns
        covar_names = columns[2:]
    else:
        # Just #IID column
        covar_names = columns[1:]

    return ",".join(covar_names)


def _create_phenotype_file(data: GenotypeData, pheno_name: str) -> Path:
    """Create phenotype file from psam.

    Extracts the phenotype column from the psam file and writes
    it as a separate phenotype file for PLINK2.

    Args:
        data: Input genotype data (must have .psam file).
        pheno_name: Name of phenotype column in psam.

    Returns:
        Path to created phenotype file.
    """
    psam_path = data.path.with_suffix(".psam")
    psam = pd.read_csv(psam_path, sep=r"\s+")

    pheno_path = data.path.parent / f"{data.path.name}.pheno"

    if "#FID" in psam.columns:
        psam[["#FID", "IID", pheno_name]].to_csv(pheno_path, index=False)
    else:
        psam[["#IID", pheno_name]].to_csv(pheno_path, index=False)

    return pheno_path


def _get_case_control_counts(data: GenotypeData) -> tuple[int, int]:
    """Get case and control counts from psam file.

    Args:
        data: Input genotype data.

    Returns:
        Tuple of (n_cases, n_controls).
    """
    psam_path = data.path.with_suffix(".psam")
    psam = pd.read_csv(psam_path, sep=r"\s+", dtype={"PHENO1": str})

    counts = psam["PHENO1"].value_counts()
    n_controls = counts.get("1", 0)
    n_cases = counts.get("2", 0)

    return int(n_cases), int(n_controls)


def run_gwas(
    data: GenotypeData,
    config: GWASConfig,
    out_path: Path,
    covar_path: Optional[Path] = None,
    covar_names: Optional[str] = None,
) -> GWASResult:
    """Run genome-wide association analysis.

    Performs association testing using PLINK2 --glm. Supports both
    binary (logistic) and quantitative (linear) phenotypes with
    optional covariates.

    Args:
        data: Input genotype data.
        config: GWAS configuration with phenotype and GLM options.
        out_path: Output path prefix (without extension).
            Will create {out_path}.PHENO1.glm.logistic.hybrid or
            {out_path}.PHENO1.glm.linear.
        covar_path: Optional path to covariate file. If provided,
            covariates will be included in the model.
        covar_names: Optional comma-separated list of covariate
            column names. If None and covar_path is provided,
            all non-ID columns will be used.

    Returns:
        GWASResult with output path and inflation metrics.

    Raises:
        ValidationError: If input data does not exist.
        GWASError: If GWAS fails.
        ExternalToolError: If PLINK2 fails.
    """
    step_name = "gwas"

    with step_context(step_name):
        logger.info(f"Running GWAS (phenotype={config.pheno_name})")

        # Validate input exists
        input_file = data.path.with_suffix(".pgen" if data.format == "pfile" else ".bed")
        if not input_file.exists():
            raise ValidationError(f"Input file does not exist: {input_file}")

        # Ensure output directory exists
        out_path.parent.mkdir(parents=True, exist_ok=True)

        # Create phenotype file
        pheno_path = _create_phenotype_file(data, config.pheno_name)
        logger.debug(f"Created phenotype file: {pheno_path}")

        # Build PLINK2 command. glm_options is a space-separated string of --glm
        # modifiers; each must be a separate argv token. The legacy code relied
        # on shell_do splitting the whole command on whitespace, but run_command
        # passes list elements verbatim, so an un-split string reaches PLINK2 as
        # one giant (invalid) --glm argument. allow-no-covars is itself a --glm
        # modifier, so it must stay in the --glm group when no covariates apply.
        plink2 = get_plink2()
        glm_args = config.glm_options.split()
        if covar_path is None:
            glm_args.append("allow-no-covars")

        cmd = [
            str(plink2),
            "--pfile" if data.format == "pfile" else "--bfile",
            str(data.path),
            "--glm", *glm_args,
            "--pheno-name", config.pheno_name,
            "--pheno", str(pheno_path),
        ]

        # Add covariates if provided
        if covar_path is not None:
            if not covar_path.exists():
                raise ValidationError(f"Covariate file does not exist: {covar_path}")

            # Get covariate names if not specified
            if covar_names is None:
                covar_names = _extract_covariate_names(covar_path)

            cmd.extend([
                "--covar", str(covar_path),
                "--covar-name", covar_names,
            ])

            if config.covariate_variance_standardize:
                cmd.append("--covar-variance-standardize")

        cmd.extend(["--out", str(out_path)])

        # Run GWAS
        result = run_command(cmd, tool_name="plink2")

        # Check for output files and determine phenotype type
        logistic_path = out_path.parent / f"{out_path.name}.{config.pheno_name}.glm.logistic.hybrid"
        linear_path = out_path.parent / f"{out_path.name}.{config.pheno_name}.glm.linear"

        if logistic_path.exists():
            return _process_logistic_results(
                output_path=logistic_path,
                data=data,
                config=config,
                log=result.stderr,
            )
        elif linear_path.exists():
            return _process_linear_results(
                output_path=linear_path,
                config=config,
                log=result.stderr,
            )
        else:
            logger.error(f"GWAS failed - no output file found at {out_path}")
            return GWASResult(
                output_path=None,
                phenotype_type="unknown",
                inflation=None,
                n_variants_tested=0,
                success=False,
                log=f"Check {out_path}.log for more information",
            )


def _process_logistic_results(
    output_path: Path,
    data: GenotypeData,
    config: GWASConfig,
    log: str,
) -> GWASResult:
    """Process logistic regression GWAS results.

    Args:
        output_path: Path to .glm.logistic.hybrid file.
        data: Input genotype data (for case/control counts).
        config: GWAS configuration.
        log: PLINK2 log output.

    Returns:
        GWASResult with inflation metrics.
    """
    logger.info("Processing logistic regression results")

    # Read GWAS results
    gwas_df = pd.read_csv(output_path, sep=r"\s+", dtype={"#CHROM": str})

    # Filter to additive test results
    gwas_df_add = gwas_df[gwas_df["TEST"] == "ADD"]
    n_variants = len(gwas_df_add)

    # Get case/control counts
    n_cases, n_controls = _get_case_control_counts(data)

    # Calculate inflation metrics
    pvalues = gwas_df_add["P"].values
    inflation = calculate_inflation(pvalues, n_cases=n_cases, n_controls=n_controls)

    # Calculate MAF-filtered inflation if requested
    lambda_gc_maf: Optional[float] = None
    lambda_1000_maf: Optional[float] = None

    if config.maf_lambdas:
        # Filter to common variants (MAF >= threshold)
        maf_mask = (
            (gwas_df_add["A1_FREQ"] >= config.maf_threshold) &
            (gwas_df_add["A1_FREQ"] <= 1 - config.maf_threshold) &
            (~gwas_df_add["P"].isna())
        )
        gwas_df_maf = gwas_df_add[maf_mask]

        if len(gwas_df_maf) > 0:
            maf_inflation = calculate_inflation(
                gwas_df_maf["P"].values,
                n_cases=n_cases,
                n_controls=n_controls,
            )
            lambda_gc_maf = maf_inflation.lambda_gc
            lambda_1000_maf = maf_inflation.lambda_1000

    # Update inflation with MAF metrics
    inflation = InflationMetrics(
        lambda_gc=inflation.lambda_gc,
        lambda_1000=inflation.lambda_1000,
        lambda_gc_maf=lambda_gc_maf,
        lambda_1000_maf=lambda_1000_maf,
        n_cases=n_cases,
        n_controls=n_controls,
    )

    lambda_1000_str = (
        f"{inflation.lambda_1000:.4f}"
        if inflation.lambda_1000 is not None
        else "N/A"
    )
    logger.info(
        f"GWAS complete: {n_variants} variants tested, "
        f"lambda={inflation.lambda_gc:.4f}, "
        f"lambda1000={lambda_1000_str}"
    )

    return GWASResult(
        output_path=output_path,
        phenotype_type="binary",
        inflation=inflation,
        n_variants_tested=n_variants,
        success=True,
        log=log,
    )


def _process_linear_results(
    output_path: Path,
    config: GWASConfig,
    log: str,
) -> GWASResult:
    """Process linear regression GWAS results.

    Args:
        output_path: Path to .glm.linear file.
        config: GWAS configuration.
        log: PLINK2 log output.

    Returns:
        GWASResult with inflation metrics.
    """
    logger.info("Processing linear regression results")

    # Read GWAS results
    gwas_df = pd.read_csv(output_path, sep=r"\s+", dtype={"#CHROM": str})

    # Filter to additive test results
    gwas_df_add = gwas_df[gwas_df["TEST"] == "ADD"]
    n_variants = len(gwas_df_add)

    # Calculate inflation (no normalization for quantitative traits)
    pvalues = gwas_df_add["P"].values
    inflation = calculate_inflation(pvalues)

    # For linear regression, we don't have case/control lambda_1000
    # Set metrics without sample-size normalization
    inflation = InflationMetrics(
        lambda_gc=inflation.lambda_gc,
        lambda_1000=None,  # Not applicable for quantitative traits
        n_cases=None,
        n_controls=None,
    )

    logger.info(f"GWAS complete: {n_variants} variants tested, lambda={inflation.lambda_gc:.4f}")

    return GWASResult(
        output_path=output_path,
        phenotype_type="quantitative",
        inflation=inflation,
        n_variants_tested=n_variants,
        success=True,
        log=log,
    )
