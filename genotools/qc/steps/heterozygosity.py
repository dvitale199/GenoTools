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

"""Heterozygosity filtering for samples.

Removes samples with extreme heterozygosity rates using PLINK2.
"""

from pathlib import Path
from typing import Any, Dict, Tuple

import pandas as pd

from genotools.core.exceptions import QCError, ValidationError
from genotools.core.executors import run_plink2
from genotools.core.genotypes import GenotypeData
from genotools.core.logging import RAW_LOG_HINT, get_logger, step_context
from genotools.qc.config import HetConfig
from genotools.qc.results import FilterResult

logger = get_logger(__name__)


def select_het_outliers(
    het: pd.DataFrame, config: HetConfig
) -> Tuple[pd.DataFrame, Dict[str, Any]]:
    """Pick the outlier rows from a PLINK ``--het`` table.

    Split out from :func:`filter_heterozygosity` so the bound arithmetic can be
    tested against a frame with known mean and sigma, without PLINK in the loop.

    Three modes, in precedence order:

    1. ``config.auto_detect`` - bounds at ``mean +/- auto_sd * sd`` of **F**.
       This is the CLI's ``--het sd [N]`` mode. Thresholding F (rather than the
       derived rate) means both CLI spellings of ``--het`` bound the same
       statistic.
    2. The legacy ``f_lower == f_upper == -1.0`` sentinel - bounds at 3 standard
       deviations of the derived **heterozygosity rate**. Preserved verbatim,
       including the ``HET`` column it adds to the output, for Python-API
       callers that already pass it. The CLI no longer produces this form.
    3. Otherwise - fixed bounds on F.

    Note that mean and sd are computed over every sample, outliers included, so
    a genuinely dispersed group widens its own bounds. That is deliberate: a
    real admixed subpopulation should not be clipped as contamination. The cost
    is that ``sd`` mode is a weak contamination screen for wide groups.

    Args:
        het: Parsed ``.het`` table, with F, OBS_CT and O(HOM) columns.
        config: Heterozygosity filtering configuration.

    Returns:
        Tuple of (outlier rows, metrics describing the mode and bounds).
    """
    if config.auto_detect:
        mean = het["F"].mean()
        std = het["F"].std()
        low = mean - (config.auto_sd * std)
        high = mean + (config.auto_sd * std)
        outliers = het[(het["F"] < low) | (het["F"] > high)]
        metrics = {
            "het_mode": "sd",
            "het_statistic": "F",
            "het_sd": config.auto_sd,
            "het_lower": low,
            "het_upper": high,
        }
    elif config.f_lower == -1.0 and config.f_upper == -1.0:
        het = het.assign(HET=(het["OBS_CT"] - het["O(HOM)"]) / het["OBS_CT"])
        mean = het["HET"].mean()
        std = het["HET"].std()
        low = mean - (3 * std)
        high = mean + (3 * std)
        outliers = het[(het["HET"] < low) | (het["HET"] > high)]
        metrics = {
            "het_mode": "sd",
            "het_statistic": "rate",
            "het_sd": 3.0,
            "het_lower": low,
            "het_upper": high,
        }
    else:
        outliers = het[
            (het["F"] <= config.f_lower) | (het["F"] >= config.f_upper)
        ]
        metrics = {
            "het_mode": "fixed",
            "het_statistic": "F",
            "het_sd": None,
            "het_lower": config.f_lower,
            "het_upper": config.f_upper,
        }
    return outliers, metrics


def filter_heterozygosity(
    data: GenotypeData,
    config: HetConfig,
    out_path: Path,
) -> FilterResult:
    """Filter samples by heterozygosity rate.

    Removes samples with extreme heterozygosity, which may indicate
    sample contamination or inbreeding. The filtering process:
    1. LD-prune variants (geno=0.01, maf=0.05, indep-pairwise 50 5 0.5)
    2. Calculate per-sample heterozygosity on pruned set
    3. Filter samples outside fixed F bounds, or outside bounds derived
       from the data - see :func:`select_het_outliers`

    Args:
        data: Input genotype data.
        config: Heterozygosity filtering configuration.
        out_path: Output path prefix (without extension).

    Returns:
        FilterResult with filtered output and metrics including:
            - samples_removed: Number of samples with extreme heterozygosity
            - pruned_samples_file: Path to .outliers file with removed sample IDs

    Raises:
        ValidationError: If input data does not exist.
        QCError: If heterozygosity calculation fails.
        ExternalToolError: If PLINK2 fails.
    """
    step_name = "het_prune"

    with step_context(step_name):
        if config.auto_detect:
            logger.info(
                f"Filtering samples by heterozygosity "
                f"(bounds derived per dataset at {config.auto_sd} sd of F)"
            )
        else:
            logger.info(
                f"Filtering samples by heterozygosity "
                f"(f_lower={config.f_lower}, f_upper={config.f_upper})"
            )

        # Validate input exists
        input_file = data.path.with_suffix(".pgen" if data.format == "pfile" else ".bed")
        if not input_file.exists():
            raise ValidationError(f"Input file does not exist: {input_file}")

        # Ensure output directory exists
        out_path.parent.mkdir(parents=True, exist_ok=True)

        # Create temp file names
        het_tmp = out_path.parent / f"{out_path.name}_het_tmp"
        het_tmp2 = out_path.parent / f"{out_path.name}_het_tmp2"
        het_tmp3 = out_path.parent / f"{out_path.name}_het_tmp3"
        outliers_out = out_path.parent / f"{out_path.name}.outliers"

        # Step 1: LD prune to get independent variants
        run_plink2(
            input_path=data.path,
            output_path=het_tmp,
            input_format=data.format,
            extra_args=[
                "--geno",
                "0.01",
                "--maf",
                "0.05",
                "--indep-pairwise",
                "50",
                "5",
                "0.5",
            ],
            make_pgen=False,
        )

        # Step 2: Extract pruned variants
        prune_in_file = het_tmp.with_suffix(".prune.in")
        if not prune_in_file.exists():
            raise QCError(
                f"LD pruning failed - {prune_in_file} not found. {RAW_LOG_HINT}"
            )

        run_plink2(
            input_path=data.path,
            output_path=het_tmp2,
            input_format=data.format,
            extra_args=["--extract", str(prune_in_file)],
        )

        # Step 3: Calculate heterozygosity
        run_plink2(
            input_path=het_tmp2,
            output_path=het_tmp3,
            input_format="pfile",
            extra_args=["--het"],
            make_pgen=False,
        )

        # Process het file
        het_file = het_tmp3.with_suffix(".het")
        if het_file.exists():
            het = pd.read_csv(het_file, sep=r"\s+")

            het_outliers, het_metrics = select_het_outliers(het, config)

            if het_metrics["het_mode"] == "sd":
                logger.info(
                    f"Derived bounds from the data: "
                    f"sd={het_metrics['het_sd']} on {het_metrics['het_statistic']}, "
                    f"bounds=[{het_metrics['het_lower']:.4f}, "
                    f"{het_metrics['het_upper']:.4f}]"
                )

            outlier_count = het_outliers.shape[0]
            het_outliers.to_csv(outliers_out, sep="\t", header=True, index=False)

            # Remove outliers from original data
            run_plink2(
                input_path=data.path,
                output_path=out_path,
                input_format=data.format,
                extra_args=["--remove", str(outliers_out)],
            )

            # Load output
            output = GenotypeData.from_path(out_path)

            # Cleanup temp files
            for temp_file in [
                prune_in_file,
                het_tmp.with_suffix(".prune.out"),
                het_tmp.with_suffix(".log"),
                het_tmp2.with_suffix(".pgen"),
                het_tmp2.with_suffix(".pvar"),
                het_tmp2.with_suffix(".psam"),
                het_tmp2.with_suffix(".log"),
                het_file,
                het_tmp3.with_suffix(".log"),
            ]:
                if temp_file.exists():
                    temp_file.unlink()

            samples_removed = data.sample_count - output.sample_count

            logger.info(
                f"Heterozygosity filtering complete: {samples_removed} samples removed"
            )

            return FilterResult(
                output=output,
                samples_removed=samples_removed,
                variants_removed=0,
                # Counts only. The JSON report renders every metric as a row
                # in a long (step, metric, pruned_count) table, so a
                # descriptive value like het_mode="sd" would land under a
                # column meaning "how many were pruned". The mode, statistic,
                # multiplier and derived bounds go to the log instead - see
                # the "Derived bounds from the data" line above.
                metrics={
                    "outlier_count": outlier_count,
                },
                log="",
                pruned_samples_file=outliers_out if outlier_count > 0 else None,
                step_name=step_name,
            )
        else:
            raise QCError(
                f"Heterozygosity calculation failed - {het_file} not found. "
                f"{RAW_LOG_HINT} "
                "Note: This can happen if there's only one sample in the genotype file."
            )
