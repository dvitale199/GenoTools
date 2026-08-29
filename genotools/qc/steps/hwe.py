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

"""Hardy-Weinberg equilibrium filtering.

Removes variants violating Hardy-Weinberg equilibrium using PLINK 1.9.
"""

from pathlib import Path

from genotools.core.exceptions import QCError, ValidationError
from genotools.core.executors import run_plink, run_plink2
from genotools.core.genotypes import GenotypeData
from genotools.core.logging import RAW_LOG_HINT, get_logger, step_context
from genotools.qc.config import HWEConfig
from genotools.qc.results import FilterResult

logger = get_logger(__name__)


def filter_hwe(
    data: GenotypeData,
    config: HWEConfig,
    out_path: Path,
) -> FilterResult:
    """Filter variants by Hardy-Weinberg equilibrium.

    Removes variants that violate Hardy-Weinberg equilibrium, which may
    indicate genotyping errors. Uses PLINK 1.9 --hwe.

    Args:
        data: Input genotype data.
        config: HWE configuration with p-value threshold and control filter.
        out_path: Output path prefix (without extension).

    Returns:
        FilterResult with filtered output and metrics including:
            - variants_removed: Number of variants filtered
            - hwe_removed_count: Same as variants_removed (legacy metric name)

    Raises:
        ValidationError: If input data does not exist.
        QCError: If HWE test fails.
        ExternalToolError: If PLINK fails.
    """
    step_name = "hwe_prune"

    with step_context(step_name):
        logger.info(
            f"Filtering variants by HWE "
            f"(threshold={config.hwe_threshold}, filter_controls={config.filter_controls})"
        )

        # Validate input exists
        input_file = data.path.with_suffix(".pgen" if data.format == "pfile" else ".bed")
        if not input_file.exists():
            raise ValidationError(f"Input file does not exist: {input_file}")

        # Ensure output directory exists
        out_path.parent.mkdir(parents=True, exist_ok=True)

        # Convert to bfile for PLINK 1.9 --hwe
        bfile_data = data.ensure_bfile(out_path.parent / f"{out_path.name}_bfile")

        # Get initial variant count
        initial_variant_count = data.variant_count

        # Create temp file names
        hwe_tmp = out_path.parent / f"{out_path.name}_hwe_tmp"

        # Build PLINK 1.9 command
        extra_args = ["--hwe", str(config.hwe_threshold), "--write-snplist"]
        if config.filter_controls:
            extra_args.insert(0, "--filter-controls")

        # Run PLINK 1.9 --hwe
        run_plink(
            input_path=bfile_data.path,
            output_path=hwe_tmp,
            input_format="bfile",
            extra_args=extra_args,
        )

        snplist_file = hwe_tmp.with_suffix(".snplist")
        if snplist_file.exists():
            # Extract passing variants and convert back to pfile
            run_plink2(
                input_path=bfile_data.path,
                output_path=out_path,
                input_format="bfile",
                extra_args=["--extract", str(snplist_file)],
            )

            # Load output
            output = GenotypeData.from_path(out_path)
            variants_removed = initial_variant_count - output.variant_count

            # Cleanup temp files
            for temp_file in [
                snplist_file,
                hwe_tmp.with_suffix(".hh"),
                hwe_tmp.with_suffix(".log"),
            ]:
                if temp_file.exists():
                    temp_file.unlink()

            # Cleanup bfiles if we created them
            if bfile_data.path != data.path:
                for ext in [".bed", ".bim", ".fam"]:
                    bfile = bfile_data.path.with_suffix(ext)
                    if bfile.exists():
                        bfile.unlink()

            logger.info(f"HWE filtering complete: {variants_removed} variants removed")

            return FilterResult(
                output=output,
                samples_removed=0,
                variants_removed=variants_removed,
                metrics={
                    "hwe_removed_count": variants_removed,
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
                f"HWE test failed. {RAW_LOG_HINT}"
            )
