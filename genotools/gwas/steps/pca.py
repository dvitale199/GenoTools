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

"""PCA computation functions.

This module provides pure functions for PCA-related operations:
- Variant pruning for PCA (filtering + LD pruning)
- PCA computation using PLINK2

Usage:
    from genotools.gwas.steps.pca import run_pca_pruning, run_pca

    # Just pruning
    pruning_result = run_pca_pruning(data, config.pruning, out_path, config.build)

    # Full PCA (includes pruning)
    pca_result = run_pca(data, config, out_path)
"""

from pathlib import Path
from typing import Optional

from genotools.core.exceptions import ValidationError
from genotools.core.executors import get_plink2, run_command, run_plink2
from genotools.core.genotypes import GenotypeData
from genotools.core.logging import RAW_LOG_HINT, get_logger, step_context
from genotools.gwas.config import PCAConfig, PCAPruningConfig, get_exclusion_regions
from genotools.gwas.results import PCAResult, PruningResult

logger = get_logger(__name__)


def _write_exclusion_file(out_path: Path, build: str) -> Path:
    """Write exclusion regions file for the specified build.

    Args:
        out_path: Base output path (file will be named {out_path}_{build}.exclusion).
        build: Genome build ('hg19' or 'hg38').

    Returns:
        Path to the written exclusion file.
    """
    exclusion_content = get_exclusion_regions(build)
    exclusion_file = out_path.parent / f"{out_path.name}_{build}.exclusion"
    with open(exclusion_file, "w") as f:
        f.write(exclusion_content)
    return exclusion_file


def _cleanup_temp_files(out_path: Path, extensions: list[str]) -> None:
    """Clean up temporary files by extension.

    Args:
        out_path: Base path prefix.
        extensions: List of extensions to remove (including the dot).
    """
    for ext in extensions:
        temp_file = out_path.parent / f"{out_path.name}{ext}"
        if temp_file.exists():
            try:
                temp_file.unlink()
            except OSError:
                pass


def run_pca_pruning(
    data: GenotypeData,
    config: PCAPruningConfig,
    out_path: Path,
    build: str = "hg38",
) -> PruningResult:
    """Prune variants for PCA analysis.

    Applies quality filtering and LD pruning to prepare a clean
    set of independent variants for PCA computation. Excludes
    problematic genomic regions (MHC, centromeres).

    Args:
        data: Input genotype data.
        config: Pruning configuration with MAF, missingness, HWE,
            and LD pruning parameters.
        out_path: Output path prefix (without extension).
        build: Genome build for exclusion regions ('hg19' or 'hg38').

    Returns:
        PruningResult with pruned output path and variant counts.

    Raises:
        ValidationError: If input data does not exist.
        ExternalToolError: If PLINK2 fails.
    """
    step_name = "pca_pruning"

    with step_context(step_name):
        logger.info(
            f"Pruning variants for PCA "
            f"(maf={config.maf}, geno={config.geno}, hwe={config.hwe})"
        )

        # Validate input exists
        input_file = data.pgen if data.format == "pfile" else data.bed
        if not input_file.exists():
            raise ValidationError(f"Input file does not exist: {input_file}")

        # Ensure output directory exists
        out_path.parent.mkdir(parents=True, exist_ok=True)

        # Get initial variant count
        initial_variant_count = data.variant_count

        # Create temp file names
        tmp_path = out_path.parent / f"{out_path.name}_tmp"
        pruned_path = out_path.parent / f"{out_path.name}_pruned"

        # Write exclusion file for problematic regions. From here on every
        # intermediate is removed however PLINK exits: a failed command raises
        # out of this function, and these files used to be left beside the
        # run's real outputs (#329).
        exclusion_file = _write_exclusion_file(out_path, build)
        try:
            return _prune(data, config, out_path, tmp_path, pruned_path, exclusion_file,
                          initial_variant_count)
        finally:
            _cleanup_temp_files(tmp_path, [".pgen", ".pvar", ".psam", ".log"])
            _cleanup_temp_files(pruned_path, [".prune.in", ".prune.out", ".log"])
            exclusion_file.unlink(missing_ok=True)


def _prune(
    data: GenotypeData,
    config: PCAPruningConfig,
    out_path: Path,
    tmp_path: Path,
    pruned_path: Path,
    exclusion_file: Path,
    initial_variant_count: int,
) -> PruningResult:
    """Filter, LD-prune and extract; the caller removes the intermediates."""
    plink2 = get_plink2()

    # Step 1: Filter variants by MAF, missingness, HWE, and exclude regions
    filter_cmd = [
        str(plink2),
        "--pfile" if data.format == "pfile" else "--bfile",
        str(data.path),
        "--maf", str(config.maf),
        "--geno", str(config.geno),
        "--hwe", str(config.hwe),
        "--autosome",
        "--exclude", "range", str(exclusion_file),
        "--make-pgen", "psam-cols=fid,parents,sex,pheno1,phenos",
        "--out", str(tmp_path),
    ]
    run_command(filter_cmd, tool_name="plink2")

    # Step 2: LD pruning with --indep-pairwise
    window, step, r2 = config.indep_pairwise_args
    prune_cmd = [
        str(plink2),
        "--pfile", str(tmp_path),
        "--indep-pairwise", str(window), str(step), str(r2),
        "--autosome",
        "--out", str(pruned_path),
    ]
    run_command(prune_cmd, tool_name="plink2")

    prune_in_file = Path(f"{pruned_path}.prune.in")
    if not prune_in_file.exists():
        # PLINK exited 0 but wrote no prune list. Too few samples cannot get
        # here: --indep-pairwise exits 13 below 50, and the runner skips such
        # a dataset before the step runs.
        logger.warning(f"PCA pruning failed - {prune_in_file} not found.")
        return PruningResult(
            output_path=out_path,
            variants_before=initial_variant_count,
            variants_after=0,
            exclusion_file=None,
            success=False,
            log=RAW_LOG_HINT,
        )

    # Step 3: Extract pruned variants
    extract_cmd = [
        str(plink2),
        "--pfile", str(tmp_path),
        "--extract", str(prune_in_file),
        "--make-pgen", "psam-cols=fid,parents,sex,pheno1,phenos",
        "--out", str(out_path),
    ]
    run_command(extract_cmd, tool_name="plink2")

    variants_after = GenotypeData.from_path(out_path).variant_count
    logger.info(
        f"PCA pruning complete: {initial_variant_count} -> {variants_after} variants "
        f"({initial_variant_count - variants_after} removed)"
    )
    return PruningResult(
        output_path=out_path,
        variants_before=initial_variant_count,
        variants_after=variants_after,
        exclusion_file=None,  # Removed by the caller
        success=True,
        log="",
    )


def run_pca(
    data: GenotypeData,
    config: PCAConfig,
    out_path: Path,
) -> PCAResult:
    """Compute principal components.

    Runs the full PCA workflow:
    1. Prune variants (filter + LD prune)
    2. Compute PCs using PLINK2 --pca

    Args:
        data: Input genotype data.
        config: PCA configuration with number of PCs and pruning settings.
        out_path: Output path prefix (without extension).
            Will create {out_path}.eigenvec and {out_path}.eigenval.

    Returns:
        PCAResult with paths to eigenvector/eigenvalue files.

    Raises:
        ValidationError: If input data does not exist.
        ExternalToolError: If PLINK2 fails.
    """
    step_name = "pca"

    with step_context(step_name):
        logger.info(f"Computing {config.n_pcs} principal components")

        # Validate input exists
        input_file = data.pgen if data.format == "pfile" else data.bed
        if not input_file.exists():
            raise ValidationError(f"Input file does not exist: {input_file}")

        # Ensure output directory exists
        out_path.parent.mkdir(parents=True, exist_ok=True)

        # Step 1: Prune variants
        pruned_path = out_path.parent / f"{out_path.name}_pca_pruned"
        pruning_result = run_pca_pruning(
            data=data,
            config=config.pruning,
            out_path=pruned_path,
            build=config.build,
        )

        if not pruning_result.success:
            logger.error("PCA failed: variant pruning did not complete successfully")
            return PCAResult(
                eigenvec_path=None,
                eigenval_path=None,
                n_pcs=config.n_pcs,
                n_samples=0,
                pruning_result=pruning_result,
                success=False,
                log=pruning_result.log,
            )

        # Step 2: Compute PCA
        plink2 = get_plink2()
        pca_cmd = [
            str(plink2),
            "--pfile", str(pruning_result.output_path),
            "--pca", str(config.n_pcs),
            "--out", str(out_path),
        ]
        try:
            result = run_command(pca_cmd, tool_name="plink2")
        finally:
            # The pruned pfiles are only PCA's input; drop them however it ends.
            _cleanup_temp_files(pruned_path, [".pgen", ".pvar", ".psam", ".log"])

        eigenvec_path = Path(f"{out_path}.eigenvec")
        eigenval_path = Path(f"{out_path}.eigenval")

        if eigenvec_path.exists():
            # Count samples in eigenvec file
            with open(eigenvec_path) as f:
                n_samples = sum(1 for _ in f) - 1  # Subtract header

            logger.info(
                f"PCA complete: {config.n_pcs} PCs computed for {n_samples} samples"
            )

            return PCAResult(
                eigenvec_path=eigenvec_path,
                eigenval_path=eigenval_path if eigenval_path.exists() else None,
                n_pcs=config.n_pcs,
                n_samples=n_samples,
                pruning_result=pruning_result,
                success=True,
                log=result.stderr,
            )
        else:
            logger.error(f"PCA failed - {eigenvec_path} not found")

            return PCAResult(
                eigenvec_path=None,
                eigenval_path=None,
                n_pcs=config.n_pcs,
                n_samples=0,
                pruning_result=pruning_result,
                success=False,
                log=RAW_LOG_HINT,
            )
