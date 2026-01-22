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

"""Linkage disequilibrium pruning.

Removes variants in high LD using PLINK2 --indep-pairwise.
"""

import logging
from pathlib import Path

from genotools.core.exceptions import QCError, ValidationError
from genotools.core.executors import run_plink2
from genotools.core.genotypes import GenotypeData
from genotools.core.logging import step_context
from genotools.qc.config import LDConfig
from genotools.qc.results import FilterResult

logger = logging.getLogger(__name__)


def prune_ld(
    data: GenotypeData,
    config: LDConfig,
    out_path: Path,
) -> FilterResult:
    """Prune variants by linkage disequilibrium.

    Removes variants in high LD to create a set of approximately
    independent variants. Uses PLINK2 --indep-pairwise.

    Args:
        data: Input genotype data.
        config: LD pruning configuration with window, step, and r2 threshold.
        out_path: Output path prefix (without extension).

    Returns:
        FilterResult with filtered output and metrics including:
            - variants_removed: Number of variants pruned
            - ld_removed_count: Same as variants_removed (legacy metric name)

    Raises:
        ValidationError: If input data does not exist.
        QCError: If LD pruning fails.
        ExternalToolError: If PLINK2 fails.
    """
    step_name = "ld_prune"

    with step_context(step_name):
        logger.info(
            f"Pruning variants by LD "
            f"(window={config.window_size}, step={config.step_size}, "
            f"r2={config.r2_threshold})"
        )

        # Validate input exists
        input_file = data.path.with_suffix(".pgen" if data.format == "pfile" else ".bed")
        if not input_file.exists():
            raise ValidationError(f"Input file does not exist: {input_file}")

        # Ensure output directory exists
        out_path.parent.mkdir(parents=True, exist_ok=True)

        # Get initial variant count
        initial_variant_count = data.variant_count

        # Create temp file names
        ld_tmp = out_path.parent / f"{out_path.name}_ld_tmp"

        # Step 1: Calculate LD and get list of variants to keep
        run_plink2(
            input_path=data.path,
            output_path=ld_tmp,
            input_format=data.format,
            extra_args=[
                "--indep-pairwise",
                str(config.window_size),
                str(config.step_size),
                str(config.r2_threshold),
            ],
            make_pgen=False,
        )

        prune_in_file = ld_tmp.with_suffix(".prune.in")
        if prune_in_file.exists():
            # Step 2: Extract pruned variants
            run_plink2(
                input_path=data.path,
                output_path=out_path,
                input_format=data.format,
                extra_args=["--extract", str(prune_in_file)],
            )

            # Load output
            output = GenotypeData.from_path(out_path)
            variants_removed = initial_variant_count - output.variant_count

            # Cleanup temp files
            for temp_file in [
                prune_in_file,
                ld_tmp.with_suffix(".prune.out"),
                ld_tmp.with_suffix(".log"),
            ]:
                if temp_file.exists():
                    temp_file.unlink()

            logger.info(f"LD pruning complete: {variants_removed} variants removed")

            return FilterResult(
                output=output,
                samples_removed=0,
                variants_removed=variants_removed,
                metrics={
                    "ld_removed_count": variants_removed,
                },
                log="",
                pruned_samples_file=None,
                step_name=step_name,
            )
        else:
            raise QCError(
                f"LD pruning failed - {prune_in_file} not found. "
                f"Check {ld_tmp}.log or {out_path}.log for more information."
            )
