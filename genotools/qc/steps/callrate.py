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

"""Call rate filtering for samples.

Removes samples with high missing genotype rates using PLINK2 --mind.
"""

import logging
from pathlib import Path

import pandas as pd

from genotools.core.exceptions import ValidationError
from genotools.core.executors import run_plink2
from genotools.core.genotypes import GenotypeData
from genotools.core.logging import step_context
from genotools.qc.config import CallrateConfig
from genotools.qc.results import FilterResult

logger = logging.getLogger(__name__)


def filter_callrate(
    data: GenotypeData,
    config: CallrateConfig,
    out_path: Path,
) -> FilterResult:
    """Filter samples by call rate (missingness).

    Removes samples with missing genotype rate exceeding the threshold.
    Uses PLINK2 --mind flag.

    Args:
        data: Input genotype data.
        config: Callrate filtering configuration with mind threshold.
        out_path: Output path prefix (without extension).

    Returns:
        FilterResult with filtered output and metrics including:
            - samples_removed: Number of samples filtered
            - pruned_samples_file: Path to .outliers file with removed sample IDs

    Raises:
        ValidationError: If input data does not exist.
        ExternalToolError: If PLINK2 fails.
    """
    step_name = "callrate_prune"

    with step_context(step_name):
        logger.info(f"Filtering samples by callrate (mind={config.mind})")

        # Validate input exists
        input_file = data.path.with_suffix(".pgen" if data.format == "pfile" else ".bed")
        if not input_file.exists():
            raise ValidationError(f"Input file does not exist: {input_file}")

        # Ensure output directory exists
        out_path.parent.mkdir(parents=True, exist_ok=True)

        # Run PLINK2 with --mind
        result = run_plink2(
            input_path=data.path,
            output_path=out_path,
            input_format=data.format,
            extra_args=["--mind", str(config.mind)],
        )

        # Load output to get counts
        output = GenotypeData.from_path(out_path)

        # Process outlier file
        mindrem_file = out_path.with_suffix(".mindrem.id")
        outliers_file = out_path.parent / f"{out_path.name}.outliers"
        pruned_file: Path | None = None

        if mindrem_file.exists():
            outliers = pd.read_csv(mindrem_file, sep=r"\s+", dtype=str)
            # Ensure #FID header format (PLINK2 may output #FID or FID)
            if "FID" in outliers.columns and "#FID" not in outliers.columns:
                outliers = outliers.rename(columns={"FID": "#FID"})
            outliers.to_csv(outliers_file, sep="\t", header=True, index=False)
            outlier_count = len(outliers)
            pruned_file = outliers_file
            # Clean up mindrem.id file
            mindrem_file.unlink()
        else:
            outlier_count = 0

        samples_removed = data.sample_count - output.sample_count

        logger.info(f"Callrate filtering complete: {samples_removed} samples removed")

        return FilterResult(
            output=output,
            samples_removed=samples_removed,
            variants_removed=0,
            metrics={
                "outlier_count": outlier_count,
            },
            log=result.stderr,
            pruned_samples_file=pruned_file,
            step_name=step_name,
        )
