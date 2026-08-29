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

"""Association analysis pipeline orchestration.

This module provides the AssocPipeline class that orchestrates
PCA and GWAS analysis, handling the flow of data between steps.

Usage:
    from genotools.gwas.pipeline import AssocPipeline
    from genotools.gwas.config import AssocConfig

    # Create pipeline with config
    config = AssocConfig(run_pca=True, run_gwas=True)
    pipeline = AssocPipeline(config)

    # Run pipeline
    result = pipeline.run(data, out_path)

    # Access results
    print(result.inflation_metrics)
    print(result.to_dict())  # Legacy format
"""

import warnings
from pathlib import Path
from typing import Optional

from genotools.core.exceptions import GWASError, ValidationError
from genotools.core.genotypes import GenotypeData
from genotools.core.logging import get_logger, step_context
from genotools.gwas.config import AssocConfig, CovariateConfig
from genotools.gwas.results import AssocResult, PCAResult, GWASResult
from genotools.gwas.steps.pca import run_pca
from genotools.gwas.steps.association import run_gwas

logger = get_logger(__name__)


class AssocPipeline:
    """Orchestrates PCA and GWAS analysis.

    This class manages the flow of association analysis:
    1. PCA computation (optional)
    2. GWAS execution (optional)
    3. Covariate handling (PCA results or external file)

    The pipeline ensures proper data flow between steps and
    handles the complex logic around covariate sources.

    Attributes:
        config: AssocConfig with PCA and GWAS settings.
    """

    def __init__(self, config: Optional[AssocConfig] = None):
        """Initialize the association pipeline.

        Args:
            config: Association configuration. Uses defaults if not provided.
        """
        self.config = config or AssocConfig()

    def run(
        self,
        data: GenotypeData,
        out_path: Path,
    ) -> AssocResult:
        """Run the association pipeline.

        Executes PCA and/or GWAS based on configuration.
        Handles covariate logic:
        - If PCA runs and external covariates provided: warn and use external
        - If PCA runs and no external covariates: use PCA as covariates
        - If no PCA and external covariates: use external
        - If no PCA and no external covariates: run GWAS without covariates

        Args:
            data: Input genotype data.
            out_path: Output path prefix (without extension).

        Returns:
            AssocResult with PCA and/or GWAS results.

        Raises:
            ValidationError: If input data does not exist.
            GWASError: If a required step fails.
        """
        with step_context("association_pipeline"):
            logger.info("Starting association analysis pipeline")

            # Validate input exists
            input_file = data.path.with_suffix(".pgen" if data.format == "pfile" else ".bed")
            if not input_file.exists():
                raise ValidationError(f"Input file does not exist: {input_file}")

            # Ensure output directory exists
            out_path.parent.mkdir(parents=True, exist_ok=True)

            pca_result: Optional[PCAResult] = None
            gwas_result: Optional[GWASResult] = None

            # Step 1: Run PCA if configured
            if self.config.run_pca:
                pca_result = run_pca(
                    data=data,
                    config=self.config.pca,
                    out_path=out_path,
                )

                if not pca_result.success:
                    logger.warning("PCA failed, continuing without PCA results")

            # Step 2: Run GWAS if configured
            if self.config.run_gwas:
                # Determine covariate source
                covar_path, covar_names = self._resolve_covariates(pca_result, out_path)

                gwas_result = run_gwas(
                    data=data,
                    config=self.config.gwas,
                    out_path=out_path,
                    covar_path=covar_path,
                    covar_names=covar_names,
                )

            # Determine overall success
            success = True
            if self.config.run_pca and self.config.run_gwas:
                success = (
                    (pca_result is not None and pca_result.success) and
                    (gwas_result is not None and gwas_result.success)
                )
            elif self.config.run_pca:
                success = pca_result is not None and pca_result.success
            elif self.config.run_gwas:
                success = gwas_result is not None and gwas_result.success

            logger.info(f"Association analysis complete (success={success})")

            return AssocResult(
                pca_result=pca_result,
                gwas_result=gwas_result,
                success=success,
            )

    def _resolve_covariates(
        self,
        pca_result: Optional[PCAResult],
        out_path: Path,
    ) -> tuple[Optional[Path], Optional[str]]:
        """Resolve which covariates to use for GWAS.

        Implements the covariate selection logic:
        1. External covariates take precedence (with warning if PCA also ran)
        2. PCA eigenvectors used if available and no external covariates
        3. No covariates if neither available

        Args:
            pca_result: PCA result, if PCA was run.
            out_path: Output path for locating eigenvector file.

        Returns:
            Tuple of (covar_path, covar_names) or (None, None).
        """
        covar_config = self.config.covariates
        eigenvec_path = out_path.with_suffix(".eigenvec")

        # Check if PCA produced eigenvectors
        pca_available = (
            pca_result is not None and
            pca_result.success and
            eigenvec_path.exists()
        )

        # Check if external covariates are configured
        has_external = covar_config.has_external_covariates()

        if has_external:
            if pca_available:
                warnings.warn(
                    "PCA ran and external covariates provided! "
                    "Defaulting to external covariates. "
                    "Recommend merging PCA eigenvectors with external covariates and rerunning."
                )
            logger.info(f"Using external covariates: {covar_config.covar_path}")
            return Path(covar_config.covar_path), covar_config.covar_names  # type: ignore

        elif pca_available and covar_config.use_pca_as_covariates:
            logger.info("Using PCA eigenvectors as covariates")
            return eigenvec_path, None

        else:
            logger.info("Running GWAS without covariates")
            return None, None


def run_association(
    data: GenotypeData,
    out_path: Path,
    config: Optional[AssocConfig] = None,
) -> AssocResult:
    """Convenience function to run association analysis.

    This is a simple wrapper around AssocPipeline for cases
    where you don't need to reuse the pipeline object.

    Args:
        data: Input genotype data.
        out_path: Output path prefix.
        config: Association configuration. Uses defaults if not provided.

    Returns:
        AssocResult with PCA and/or GWAS results.

    Example:
        from genotools.gwas import run_association
        from genotools.core import GenotypeData

        data = GenotypeData.from_path("/path/to/genotypes")
        result = run_association(data, Path("/path/to/output"))
        print(result.inflation_metrics)
    """
    pipeline = AssocPipeline(config)
    return pipeline.run(data, out_path)
