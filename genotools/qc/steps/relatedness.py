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

"""Relatedness filtering for samples.

Removes related and/or duplicated samples using PLINK2 KING-based methods.
"""

import logging
import platform
import shutil
from pathlib import Path

import numpy as np
import pandas as pd

from genotools.core.exceptions import QCError, ValidationError
from genotools.core.executors import run_command, run_plink2, get_king
from genotools.core.genotypes import GenotypeData
from genotools.core.logging import step_context
from genotools.qc.config import RelatedConfig
from genotools.qc.results import FilterResult

logger = logging.getLogger(__name__)


def filter_relatedness(
    data: GenotypeData,
    config: RelatedConfig,
    out_path: Path,
) -> FilterResult:
    """Filter related and/or duplicated samples.

    Uses PLINK2's KING-based relatedness methods to identify and remove
    related sample pairs. The process:
    1. Create KING kinship table with --make-king-table
    2. Identify samples to remove with --king-cutoff at related threshold
    3. Identify duplicates with --king-cutoff at duplicate threshold
    4. Remove selected samples based on config settings

    Kinship coefficient reference (Manichaikul et al. 2010):
        - >= 0.354: MZ twin/duplicate
        - 0.177-0.354: 1st degree (parent-child, full sibling)
        - 0.0884-0.177: 2nd degree (half-sibling, avuncular, grandparent)
        - 0.0442-0.0884: 3rd degree (first cousin)

    Args:
        data: Input genotype data.
        config: Relatedness filtering configuration.
        out_path: Output path prefix (without extension).

    Returns:
        FilterResult with filtered output and metrics including:
            - samples_removed: Number of related/duplicated samples removed
            - related_count: Number of related (non-duplicate) samples
            - duplicated_count: Number of duplicate samples
            - pruned_samples_file: Path to file with removed sample IDs
            - related_samples: Path to .related file with pair information

    Raises:
        ValidationError: If input data does not exist.
        QCError: If relatedness analysis fails.
        ExternalToolError: If PLINK2 fails.
    """
    step_name = "related_prune"

    with step_context(step_name):
        logger.info(
            f"Filtering samples by relatedness "
            f"(related_cutoff={config.related_cutoff}, "
            f"duplicated_cutoff={config.duplicated_cutoff})"
        )

        # Validate input exists
        input_file = data.path.with_suffix(".pgen" if data.format == "pfile" else ".bed")
        if not input_file.exists():
            raise ValidationError(f"Input file does not exist: {input_file}")

        # Ensure output directory exists
        out_path.parent.mkdir(parents=True, exist_ok=True)

        # Create temp file names
        king1 = out_path.parent / f"{out_path.name}_related_king"
        king2 = out_path.parent / f"{out_path.name}_duplicated_king"
        related_pairs = out_path.parent / f"{out_path.name}_pairs"
        related_out = related_pairs.with_suffix(".related")
        pruned_out: Path | None = out_path.parent / f"{out_path.name}.pruned"

        # Step 1: Create KING table with related pairs
        run_plink2(
            input_path=data.path,
            output_path=related_pairs,
            input_format=data.format,
            extra_args=[
                "--make-king-table",
                "--make-king",
                "triangle",
                "bin",
                "--king-table-filter",
                str(config.related_cutoff),
            ],
            make_pgen=False,
        )

        # Step 2: Get samples to remove at related threshold
        run_plink2(
            input_path=data.path,
            output_path=king1,
            input_format=data.format,
            extra_args=[
                "--king-cutoff",
                str(related_pairs),
                str(config.related_cutoff),
            ],
            make_pgen=False,
        )

        # Step 3: Get samples to remove at duplicate threshold
        run_plink2(
            input_path=data.path,
            output_path=king2,
            input_format=data.format,
            extra_args=[
                "--king-cutoff",
                str(related_pairs),
                str(config.duplicated_cutoff),
            ],
            make_pgen=False,
        )

        # Check output files exist
        kin0_file = related_pairs.with_suffix(".kin0")
        related_file = king1.with_name(f"{king1.name}.king.cutoff.out.id")
        duplicated_file = king2.with_name(f"{king2.name}.king.cutoff.out.id")

        if kin0_file.exists() and related_file.exists() and duplicated_file.exists():
            # Create .related file with relationship categories
            kinship = pd.read_csv(kin0_file, sep=r"\s+")
            kinship["REL"] = pd.cut(
                x=kinship["KINSHIP"],
                bins=[-np.inf, 0.0884, 0.177, 0.354, np.inf],
                labels=["unrel", "second_deg", "first_deg", "duplicate"],
            )
            kinship.to_csv(related_out, index=False)

            # Copy cutoff files
            shutil.copy(related_file, king1.with_suffix(".related"))
            related_count = sum(1 for _ in open(king1.with_suffix(".related"))) - 1

            shutil.copy(duplicated_file, king2.with_suffix(".duplicated"))
            duplicated_count = sum(1 for _ in open(king2.with_suffix(".duplicated"))) - 1

            # Adjust related count (duplicates are also counted as related)
            related_count = related_count - duplicated_count

            # Read duplicate samples
            duplicated = pd.read_csv(king2.with_suffix(".duplicated"), sep=r"\s+")

            # Determine which samples to remove based on config
            if config.prune_related and config.prune_duplicated:
                # Remove all related (includes duplicates)
                run_plink2(
                    input_path=data.path,
                    output_path=out_path,
                    input_format=data.format,
                    extra_args=["--remove", str(related_file)],
                )

                related = pd.read_csv(king1.with_suffix(".related"), sep=r"\s+")
                grm_pruned = pd.concat([related, duplicated], ignore_index=True)

            elif config.prune_duplicated and not config.prune_related:
                # Remove duplicates only
                run_plink2(
                    input_path=data.path,
                    output_path=out_path,
                    input_format=data.format,
                    extra_args=["--remove", str(duplicated_file)],
                )

                grm_pruned = duplicated
                related_count = 0

            else:
                # No pruning - just copy files
                for ext in ["pgen", "psam", "pvar"]:
                    src = data.path.with_suffix(f".{ext}")
                    dst = out_path.with_suffix(f".{ext}")
                    shutil.copy(src, dst)

                grm_pruned = pd.DataFrame()
                related_count = 0
                duplicated_count = 0

            # Format pruned samples file
            if not grm_pruned.empty:
                if "#FID" in grm_pruned.columns:
                    grm_pruned = grm_pruned.rename(columns={"#IID": "IID"})
                else:
                    grm_pruned["#FID"] = 0
                    if "#IID" in grm_pruned.columns:
                        grm_pruned = grm_pruned.rename(columns={"#IID": "IID"})
                grm_pruned = grm_pruned.drop_duplicates(subset=["#FID", "IID"], keep="last")
                grm_pruned.to_csv(pruned_out, sep="\t", header=True, index=False)
            else:
                pruned_out = None

            # Load output
            output = GenotypeData.from_path(out_path)

            # Cleanup temp files
            for temp_file in [
                king1.with_name(f"{king1.name}.king.cutoff.in.id"),
                related_file,
                king1.with_suffix(".related"),
                king1.with_suffix(".log"),
                king2.with_suffix(".duplicated"),
                king2.with_name(f"{king2.name}.king.cutoff.in.id"),
                duplicated_file,
                king2.with_suffix(".log"),
                related_pairs.with_suffix(".king.bin"),
                related_pairs.with_suffix(".king.id"),
                kin0_file,
                related_pairs.with_suffix(".log"),
            ]:
                if temp_file.exists():
                    temp_file.unlink()

            samples_removed = data.sample_count - output.sample_count

            logger.info(
                f"Relatedness filtering complete: {samples_removed} samples removed "
                f"({related_count} related, {duplicated_count} duplicated)"
            )

            return FilterResult(
                output=output,
                samples_removed=samples_removed,
                variants_removed=0,
                metrics={
                    "related_count": related_count,
                    "duplicated_count": duplicated_count,
                },
                log="",
                pruned_samples_file=pruned_out,
                step_name=step_name,
            )
        else:
            raise QCError(
                "Relatedness analysis failed. "
                "Check all_plink_logs.gtlog for more information."
            )


def verify_kinship(
    data: GenotypeData,
    out_path: Path,
) -> FilterResult:
    """Verify family relationships using KING.

    Identifies samples with discordant family IDs by running KING
    relatedness analysis. This does NOT remove samples, only reports
    discrepancies.

    Note: KING is only available on Linux.

    Args:
        data: Input genotype data.
        out_path: Output path prefix (without extension).

    Returns:
        FilterResult with metrics including:
            - same_fid_unrelated_count: Pairs with same FID but unrelated
            - diff_fid_related_count: Pairs with different FID but related
            - same_fid_unrelated: Path to file with same-FID unrelated pairs
            - diff_fid_related: Path to file with different-FID related pairs

    Raises:
        QCError: If not running on Linux or KING fails.
    """
    step_name = "confirming_kinship"

    with step_context(step_name):
        logger.info("Verifying kinship relationships using KING")

        # Check platform
        if platform.system() != "Linux":
            logger.warning("KING-based kinship verification only runs on Linux")
            return FilterResult(
                output=data,  # No changes made
                samples_removed=0,
                variants_removed=0,
                metrics={
                    "same_fid_unrelated_count": 0,
                    "diff_fid_related_count": 0,
                },
                log="Skipped: KING only available on Linux",
                pruned_samples_file=None,
                step_name=step_name,
            )

        # Validate input exists
        input_file = data.path.with_suffix(".pgen" if data.format == "pfile" else ".bed")
        if not input_file.exists():
            raise ValidationError(f"Input file does not exist: {input_file}")

        # Ensure output directory exists
        out_path.parent.mkdir(parents=True, exist_ok=True)

        # Convert to bfile if needed (KING requires bfiles)
        bfile_data = data.ensure_bfile(out_path.parent / f"{out_path.name}_bfile")

        # Output files
        temp = out_path.parent / f"{out_path.name}_temp"
        same_fid_unrelated = out_path.parent / f"{out_path.name}_same_fid.unrelated"
        diff_fid_related = out_path.parent / f"{out_path.name}_diff_fid.related"

        # Get KING executable
        king_exec = get_king()
        if king_exec is None:
            raise QCError("KING executable not found")

        # Run KING
        king_cmd = [
            str(king_exec),
            "-b",
            f"{bfile_data.path}.bed",
            "--related",
            "--build",
            "--degree",
            "3",
            "--prefix",
            str(temp),
        ]
        run_command(king_cmd, tool_name="king", check=True)

        # Process results
        third_deg_cutoff = 0.0442
        num_same_fid_unrelated = 0
        num_diff_fid_related = 0

        kin0_file = temp.with_suffix(".kin0")
        kin_file = temp.with_suffix(".kin")

        if kin0_file.exists() or kin_file.exists():
            # Different FID but related
            if kin0_file.exists():
                kin0_file.rename(diff_fid_related)
                num_diff_fid_related = sum(1 for _ in open(diff_fid_related)) - 1

            # Same FID but unrelated
            if kin_file.exists():
                kin = pd.read_csv(kin_file, sep=r"\s+")
                unrelated = kin[kin["Kinship"] <= third_deg_cutoff]
                unrelated.to_csv(same_fid_unrelated, sep="\t", header=True, index=False)
                num_same_fid_unrelated = len(unrelated)

            # Cleanup temp files
            for temp_file in [
                kin_file,
                temp.with_suffix(".seg"),
                temp.parent / f"{temp.name}X.kin",
                temp.parent / f"{temp.name}X.kin0",
                temp.with_suffix(".segments.gz"),
                temp.parent / f"{temp.name}allsegs.txt",
                temp.parent / f"{temp.name}build.log",
                temp.parent / f"{temp.name}updateids.txt",
                temp.parent / f"{temp.name}updateparents.txt",
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
                f"Kinship verification complete: "
                f"{num_same_fid_unrelated} same-FID unrelated pairs, "
                f"{num_diff_fid_related} different-FID related pairs"
            )

            return FilterResult(
                output=data,  # No samples removed
                samples_removed=0,
                variants_removed=0,
                metrics={
                    "same_fid_unrelated_count": num_same_fid_unrelated,
                    "diff_fid_related_count": num_diff_fid_related,
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

            raise QCError("KING relatedness analysis failed")
