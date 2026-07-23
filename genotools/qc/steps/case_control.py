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

"""Case-control differential missingness filtering.

Removes variants with significantly different missing rates between
cases and controls using PLINK 1.9 --test-missing.
"""

from pathlib import Path

import pandas as pd

from genotools.core.exceptions import QCError, ValidationError
from genotools.core.executors import run_plink, run_plink2
from genotools.core.genotypes import GenotypeData
from genotools.core.logging import get_logger, step_context
from genotools.qc.config import CaseControlConfig
from genotools.qc.results import FilterResult

logger = get_logger(__name__)


def filter_case_control(
    data: GenotypeData,
    config: CaseControlConfig,
    out_path: Path,
) -> FilterResult:
    """Filter variants by case-control differential missingness.

    Removes variants where missing genotype rates differ significantly
    between cases and controls. Uses PLINK 1.9 --test-missing.

    Args:
        data: Input genotype data.
        config: Case-control missingness configuration with p-value threshold.
        out_path: Output path prefix (without extension).

    Returns:
        FilterResult with filtered output and metrics including:
            - variants_removed: Number of variants filtered
            - mis_removed_count: Same as variants_removed (legacy metric name)

    Raises:
        ValidationError: If input data does not exist.
        QCError: If data lacks both cases and controls, or test fails.
        ExternalToolError: If PLINK fails.
    """
    step_name = "case_control_missingness_prune"

    with step_context(step_name):
        logger.info(
            f"Filtering variants by case-control missingness (p={config.p_threshold})"
        )

        # Validate input exists
        input_file = data.path.with_suffix(".pgen" if data.format == "pfile" else ".bed")
        if not input_file.exists():
            raise ValidationError(f"Input file does not exist: {input_file}")

        # Ensure output directory exists
        out_path.parent.mkdir(parents=True, exist_ok=True)

        # Convert to bfile for PLINK 1.9 --test-missing
        bfile_data = data.ensure_bfile(out_path.parent / f"{out_path.name}_bfile")

        # Get initial variant count
        initial_variant_count = data.variant_count

        # Check if data contains both cases and controls
        # fam is guaranteed to exist after ensure_bfile()
        assert bfile_data.fam is not None, "bfile format must have .fam file"
        fam = pd.read_csv(
            bfile_data.fam, sep=r"\s+", header=None, usecols=[5], names=["case"]
        )
        if not all(x in fam["case"].unique() for x in [1, 2]):
            # Cleanup bfiles if we created them
            if bfile_data.path != data.path:
                for ext in [".bed", ".bim", ".fam"]:
                    bfile = bfile_data.path.with_suffix(ext)
                    if bfile.exists():
                        bfile.unlink()

            raise QCError(
                "Case-control missingness filtering requires both cases (2) and "
                "controls (1) in the phenotype column. Check your data."
            )

        # Create temp file names
        mis_tmp = out_path.parent / f"{out_path.name}_mis_tmp"

        # Run PLINK 1.9 --test-missing
        run_plink(
            input_path=bfile_data.path,
            output_path=mis_tmp,
            input_format="bfile",
            extra_args=["--test-missing"],
        )

        missing_file = mis_tmp.with_suffix(".missing")
        if missing_file.exists():
            # Read missing file and identify variants to exclude
            mis = pd.read_csv(missing_file, sep=r"\s+")
            exclude = mis[mis.P <= config.p_threshold].loc[:, "SNP"]
            exclude_file = mis_tmp.with_suffix(".exclude")
            exclude.to_csv(exclude_file, sep="\t", header=False, index=False)

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
                mis_tmp.with_suffix(".hh"),
                missing_file,
                mis_tmp.with_suffix(".log"),
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
                f"Case-control missingness filtering complete: "
                f"{variants_removed} variants removed"
            )

            return FilterResult(
                output=output,
                samples_removed=0,
                variants_removed=variants_removed,
                metrics={
                    "mis_removed_count": variants_removed,
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
                f"Case-control missingness test failed. "
                f"Check {mis_tmp}.log for more information."
            )
