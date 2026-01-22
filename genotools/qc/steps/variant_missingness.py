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

"""Variant missingness filtering.

Removes variants with high missing genotype rates using PLINK2 --geno.
"""

import logging
from pathlib import Path

from genotools.core.exceptions import ValidationError
from genotools.core.executors import run_plink2
from genotools.core.genotypes import GenotypeData
from genotools.core.logging import step_context
from genotools.qc.config import GenoConfig
from genotools.qc.results import FilterResult

logger = logging.getLogger(__name__)


def filter_variant_missingness(
    data: GenotypeData,
    config: GenoConfig,
    out_path: Path,
) -> FilterResult:
    """Filter variants by missingness rate.

    Removes variants with missing genotype rate exceeding the threshold.
    Uses PLINK2 --geno flag.

    Args:
        data: Input genotype data.
        config: Variant missingness configuration with geno threshold.
        out_path: Output path prefix (without extension).

    Returns:
        FilterResult with filtered output and metrics including:
            - variants_removed: Number of variants filtered
            - geno_removed_count: Same as variants_removed (legacy metric name)

    Raises:
        ValidationError: If input data does not exist.
        ExternalToolError: If PLINK2 fails.
    """
    step_name = "geno_prune"

    with step_context(step_name):
        logger.info(f"Filtering variants by missingness (geno={config.geno})")

        # Validate input exists
        input_file = data.path.with_suffix(".pgen" if data.format == "pfile" else ".bed")
        if not input_file.exists():
            raise ValidationError(f"Input file does not exist: {input_file}")

        # Ensure output directory exists
        out_path.parent.mkdir(parents=True, exist_ok=True)

        # Get initial variant count
        initial_variant_count = data.variant_count

        # Run PLINK2 with --geno
        result = run_plink2(
            input_path=data.path,
            output_path=out_path,
            input_format=data.format,
            extra_args=["--geno", str(config.geno)],
        )

        # Load output to get counts
        output = GenotypeData.from_path(out_path)

        variants_removed = initial_variant_count - output.variant_count

        logger.info(
            f"Variant missingness filtering complete: {variants_removed} variants removed"
        )

        return FilterResult(
            output=output,
            samples_removed=0,
            variants_removed=variants_removed,
            metrics={
                "geno_removed_count": variants_removed,
            },
            log=result.stderr,
            pruned_samples_file=None,
            step_name=step_name,
        )
