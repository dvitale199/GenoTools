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

"""Haplotype missingness filtering.

Removes variants with non-random patterns of missing data using
PLINK 1.9 --test-mishap.
"""

from pathlib import Path

import pandas as pd

from genotools.core.exceptions import QCError, ValidationError
from genotools.core.executors import run_plink, run_plink2
from genotools.core.genotypes import GenotypeData
from genotools.core.logging import RAW_LOG_HINT, get_logger, step_context
from genotools.qc.config import HaplotypeConfig
from genotools.qc.results import FilterResult

logger = get_logger(__name__)


def filter_haplotype(
    data: GenotypeData,
    config: HaplotypeConfig,
    out_path: Path,
) -> FilterResult:
    """Filter variants by haplotype missingness pattern.

    Removes variants with non-random patterns of missing data that
    suggest genotyping problems. Uses PLINK 1.9 --test-mishap.

    When significant associations are found, both the central variant
    and flanking variants are removed.

    Args:
        data: Input genotype data.
        config: Haplotype missingness configuration with p-value threshold.
        out_path: Output path prefix (without extension).

    Returns:
        FilterResult with filtered output and metrics including:
            - variants_removed: Number of variants filtered
            - haplotype_removed_count: Same as variants_removed (legacy metric name)

    Raises:
        ValidationError: If input data does not exist.
        QCError: If haplotype test fails.
        ExternalToolError: If PLINK fails.
    """
    step_name = "haplotype_prune"

    with step_context(step_name):
        # Determine effective p-value threshold
        p_threshold = config.p_threshold

        # Validate input exists
        input_file = data.path.with_suffix(".pgen" if data.format == "pfile" else ".bed")
        if not input_file.exists():
            raise ValidationError(f"Input file does not exist: {input_file}")

        # Ensure output directory exists
        out_path.parent.mkdir(parents=True, exist_ok=True)

        # Convert to bfile for PLINK 1.9 --test-mishap
        bfile_data = data.ensure_bfile(out_path.parent / f"{out_path.name}_bfile")

        # Get initial variant count
        initial_variant_count = data.variant_count

        # Adjust threshold for large sample sizes if configured
        sample_count = data.sample_count
        if config.adjust_for_sample_size and sample_count > 10000:
            p_threshold = 0.05 / sample_count
            logger.info(
                f"Adjusted p-threshold for large sample size: {p_threshold:.2e}"
            )

        logger.info(f"Filtering variants by haplotype missingness (p={p_threshold})")

        # Create temp file names
        hap_tmp = out_path.parent / f"{out_path.name}_hap_tmp"

        # Run PLINK 1.9 --test-mishap
        run_plink(
            input_path=bfile_data.path,
            output_path=hap_tmp,
            input_format="bfile",
            extra_args=["--maf", "0.05", "--test-mishap"],
        )

        mishap_file = hap_tmp.with_suffix(".missing.hap")
        if mishap_file.exists():
            # Read mishap file and identify variants to exclude
            mis_hap = pd.read_csv(mishap_file, sep=r"\s+")

            # Get reference SNPs (central variants) for significant associations
            ref_snps = list(mis_hap[mis_hap.P <= p_threshold].loc[:, "SNP"])

            # Get flanking SNPs for significant associations
            mis_hap_snps = list(
                mis_hap[mis_hap.P <= p_threshold].loc[:, "FLANKING"].str.split("|")
            )
            flanking_snps = [rsid for ls in mis_hap_snps for rsid in ls]

            # Combine reference and flanking SNPs
            all_snps = ref_snps + flanking_snps
            snp_df = pd.DataFrame({"snp": all_snps})
            exclude_file = hap_tmp.with_suffix(".exclude")
            snp_df["snp"].to_csv(exclude_file, sep="\t", header=False, index=False)

            # Exclude variants and convert back to pfile
            run_plink2(
                input_path=bfile_data.path,
                output_path=out_path,
                input_format="bfile",
                extra_args=["--exclude", str(exclude_file)],
            )

            # Load output
            output = GenotypeData.from_path(out_path)
            variants_removed = initial_variant_count - output.variant_count

            # Cleanup temp files
            for temp_file in [
                exclude_file,
                mishap_file,
                hap_tmp.with_suffix(".hh"),
                hap_tmp.with_suffix(".log"),
            ]:
                if temp_file.exists():
                    temp_file.unlink()

            # Cleanup bfiles if we created them
            if bfile_data.path != data.path:
                for ext in [".bed", ".bim", ".fam"]:
                    bfile = bfile_data.path.with_suffix(ext)
                    if bfile.exists():
                        bfile.unlink()

            logger.info(
                f"Haplotype missingness filtering complete: "
                f"{variants_removed} variants removed"
            )

            return FilterResult(
                output=output,
                samples_removed=0,
                variants_removed=variants_removed,
                metrics={
                    "haplotype_removed_count": variants_removed,
                },
                log="",
                pruned_samples_file=None,
                step_name=step_name,
            )
        else:
            # Cleanup bfiles if we created them
            if bfile_data.path != data.path:
                for ext in [".bed", ".bim", ".fam"]:
                    bfile = bfile_data.path.with_suffix(ext)
                    if bfile.exists():
                        bfile.unlink()

            raise QCError(
                f"Haplotype missingness test failed. {RAW_LOG_HINT}"
            )
