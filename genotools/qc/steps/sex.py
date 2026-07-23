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

"""Sex check filtering for samples.

Removes samples with sex discrepancies using PLINK 1.9 --check-sex.
"""

from pathlib import Path

import pandas as pd

from genotools.core.exceptions import QCError, ValidationError
from genotools.core.executors import run_plink, run_plink2
from genotools.core.genotypes import GenotypeData
from genotools.core.logging import get_logger, step_context
from genotools.qc.config import SexConfig
from genotools.qc.results import FilterResult

logger = get_logger(__name__)


def filter_sex(
    data: GenotypeData,
    config: SexConfig,
    out_path: Path,
) -> FilterResult:
    """Filter samples by sex check discrepancies.

    Removes samples where reported sex doesn't match genetic sex based on
    X chromosome F-statistic. Uses two methods:
    1. Standard --check-sex on all variants
    2. --check-sex on X chromosome variants only (filtered for MAF, geno, HWE)

    Samples flagged by either method are removed.

    Args:
        data: Input genotype data.
        config: Sex check configuration with F-statistic thresholds.
        out_path: Output path prefix (without extension).

    Returns:
        FilterResult with filtered output and metrics including:
            - samples_removed: Number of samples with sex discrepancies
            - pruned_samples_file: Path to .outliers file with removed sample IDs

    Raises:
        ValidationError: If input data does not exist.
        QCError: If sex check fails.
        ExternalToolError: If PLINK fails.
    """
    step_name = "sex_prune"

    with step_context(step_name):
        logger.info(
            f"Filtering samples by sex check "
            f"(female_max_f={config.female_max_f}, male_min_f={config.male_min_f})"
        )

        # Validate input exists
        input_file = data.path.with_suffix(".pgen" if data.format == "pfile" else ".bed")
        if not input_file.exists():
            raise ValidationError(f"Input file does not exist: {input_file}")

        # Ensure output directory exists
        out_path.parent.mkdir(parents=True, exist_ok=True)

        # Convert to bfile for PLINK 1.9 --check-sex
        bfile_data = data.ensure_bfile(out_path.parent / f"{out_path.name}_bfile")

        # Create temp file names
        sex_tmp1 = out_path.parent / f"{out_path.name}_sex_tmp1"
        sex_tmp2 = out_path.parent / f"{out_path.name}_sex_tmp2"
        sex_fails = out_path.parent / f"{out_path.name}.outliers"

        # Method 1: Standard sex check
        run_plink(
            input_path=bfile_data.path,
            output_path=sex_tmp1,
            input_format="bfile",
            extra_args=[
                "--check-sex",
                str(config.female_max_f),
                str(config.male_min_f),
                "--maf",
                "0.05",
            ],
        )

        # Method 2: Sex check on X chromosome only with QC filters
        run_plink(
            input_path=bfile_data.path,
            output_path=sex_tmp2,
            input_format="bfile",
            extra_args=[
                "--chr",
                "23",
                "--from-bp",
                "2781479",
                "--to-bp",
                "155701383",
                "--maf",
                "0.05",
                "--geno",
                "0.05",
                "--hwe",
                "1E-5",
                "--check-sex",
                str(config.female_max_f),
                str(config.male_min_f),
            ],
        )

        # Process results
        sexcheck1 = sex_tmp1.with_suffix(".sexcheck")
        sexcheck2 = sex_tmp2.with_suffix(".sexcheck")

        if sexcheck1.exists() and sexcheck2.exists():
            # Grab fails from .sexcheck files
            sex1 = pd.read_csv(sexcheck1, sep=r"\s+")
            sex_fail1 = sex1[sex1.STATUS == "PROBLEM"]

            sex2 = pd.read_csv(sexcheck2, sep=r"\s+")
            sex_fail2 = sex2[sex2.STATUS == "PROBLEM"]

            # Combine and deduplicate
            sex_fail_df = pd.concat([sex_fail1, sex_fail2], ignore_index=True)
            sex_fail_ids = sex_fail_df.loc[:, ["FID", "IID"]].drop_duplicates(
                subset=["FID", "IID"]
            )
            sex_fail_ids = sex_fail_ids.rename(columns={"FID": "#FID"})
            sex_fail_count = sex_fail_ids.shape[0]
            sex_fail_ids.to_csv(sex_fails, sep="\t", header=True, index=False)

            # Remove sex fail samples from geno
            run_plink2(
                input_path=data.path,
                output_path=out_path,
                input_format=data.format,
                extra_args=["--remove", str(sex_fails)],
            )

            # Load output
            output = GenotypeData.from_path(out_path)

            # Cleanup temp files
            for temp_file in [
                sex_tmp1.with_suffix(".hh"),
                sexcheck1,
                sex_tmp2.with_suffix(".hh"),
                sexcheck2,
                sex_tmp1.with_suffix(".log"),
                sex_tmp2.with_suffix(".log"),
            ]:
                if temp_file.exists():
                    temp_file.unlink()

            # Cleanup bfiles if we created them
            if bfile_data.path != data.path:
                for ext in [".bed", ".bim", ".fam"]:
                    bfile = bfile_data.path.with_suffix(ext)
                    if bfile.exists():
                        bfile.unlink()

            samples_removed = data.sample_count - output.sample_count

            logger.info(f"Sex check filtering complete: {samples_removed} samples removed")

            return FilterResult(
                output=output,
                samples_removed=samples_removed,
                variants_removed=0,
                metrics={
                    "outlier_count": sex_fail_count,
                },
                log="",
                pruned_samples_file=sex_fails if sex_fail_count > 0 else None,
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
                f"Sex check failed. Check {sex_tmp1}.log and {sex_tmp2}.log "
                "for more information."
            )
