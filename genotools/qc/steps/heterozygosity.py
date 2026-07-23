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

import pandas as pd

from genotools.core.exceptions import QCError, ValidationError
from genotools.core.executors import run_plink2
from genotools.core.genotypes import GenotypeData
from genotools.core.logging import get_logger, step_context
from genotools.qc.config import HetConfig
from genotools.qc.results import FilterResult

logger = get_logger(__name__)


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
    3. Filter samples outside F-statistic bounds or 3 std devs from mean

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
        logger.info(
            f"Filtering samples by heterozygosity "
            f"(f_lower={config.f_lower}, f_upper={config.f_upper}, "
            f"auto_detect={config.auto_detect})"
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
                f"LD pruning failed - {prune_in_file} not found. "
                f"Check {het_tmp}.log for details."
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

            # Determine if using auto-detect mode
            use_auto = config.auto_detect or (
                config.f_lower == -1.0 and config.f_upper == -1.0
            )

            if use_auto:
                # Calculate heterozygosity rate and use 3 std devs
                het["HET"] = (het["OBS_CT"] - het["O(HOM)"]) / het["OBS_CT"]
                het_mean = het["HET"].mean()
                het_std = het["HET"].std()
                het_low = het_mean - (3 * het_std)
                het_high = het_mean + (3 * het_std)
                het_outliers = het[(het["HET"] < het_low) | (het["HET"] > het_high)]
                logger.info(
                    f"Auto-detect mode: mean={het_mean:.4f}, std={het_std:.4f}, "
                    f"bounds=[{het_low:.4f}, {het_high:.4f}]"
                )
            else:
                # Use fixed F-statistic bounds
                het_outliers = het[
                    (het["F"] <= config.f_lower) | (het["F"] >= config.f_upper)
                ]

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
                f"Check {het_tmp}.log, {het_tmp2}.log, or {het_tmp3}.log for details. "
                "Note: This can happen if there's only one sample in the genotype file."
            )
