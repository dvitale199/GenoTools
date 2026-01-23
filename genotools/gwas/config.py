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

"""GWAS/Association configuration dataclasses.

This module provides configuration classes for PCA and GWAS analysis.
All configs are frozen dataclasses with validation.

Usage:
    from genotools.gwas.config import PCAConfig, GWASConfig, AssocConfig

    # Use defaults
    pca_config = PCAConfig()

    # Custom configuration
    pca_config = PCAConfig(n_pcs=20, maf=0.05)
    gwas_config = GWASConfig(maf_lambdas=True)

    # Combined config
    assoc_config = AssocConfig(pca=pca_config, gwas=gwas_config)
"""

from dataclasses import dataclass, field
from typing import Literal, Optional, Tuple

from genotools.core.config import ThresholdConfig, BaseConfig


# Build-specific exclusion regions for PCA
# These regions are excluded to avoid inflation from population structure
# Includes: MHC region, centromeres, and other high-LD regions

HG19_EXCLUSION_REGIONS: str = """5 44000000 51500000 r1
6 25000000 33500000 r2
8 8000000 12000000 r3
11 45000000 57000000 r4
"""

HG38_EXCLUSION_REGIONS: str = """1 47534328 51534328 r1
2 133742429 137242430 r2
2 182135273 189135274 r3
3 47458510 49962567 r4
3 83450849 86950850 r5
5 98664296 101164296 r6
5 129664307 132664308 r7
5 136164311 139164311 r8
6 24999772 35032223 r9
6 139678863 142178863 r10
8 7142478 13142491 r11
8 110987771 113987771 r12
11 87789108 90766832 r13
12 109062195 111562196 r14
20 33412194 35912078 r15
"""


def get_exclusion_regions(build: str) -> str:
    """Get exclusion regions for a genome build.

    Args:
        build: Genome build ('hg19' or 'hg38').

    Returns:
        String with exclusion regions in BED-like format.

    Raises:
        ValueError: If build is not recognized.
    """
    if build == "hg19":
        return HG19_EXCLUSION_REGIONS
    elif build == "hg38":
        return HG38_EXCLUSION_REGIONS
    else:
        raise ValueError(f"Invalid build: {build}. Must be 'hg19' or 'hg38'.")


@dataclass(frozen=True)
class PCAPruningConfig(ThresholdConfig):
    """Configuration for variant pruning before PCA.

    Controls quality filtering and LD pruning applied before
    running PCA to ensure only high-quality, independent variants
    are used for ancestry/population structure inference.

    Attributes:
        maf: Minimum allele frequency threshold. Variants with
            MAF < this are excluded. Default is 0.01 (1%).
        geno: Maximum variant missingness. Variants with missing
            rate > this are excluded. Default is 0.01 (1%).
        hwe: Hardy-Weinberg equilibrium p-value threshold. Variants
            with HWE p < this are excluded. Default is 5e-6.
        ld_window_size: LD pruning window size in kb. Default is 1000.
        ld_step_size: LD pruning step size in variants. Default is 10.
        ld_r2_threshold: LD pruning r² threshold. Pairs with r² above
            this have one variant removed. Default is 0.02.
    """

    maf: float = 0.01
    geno: float = 0.01
    hwe: float = 5e-6
    ld_window_size: int = 1000
    ld_step_size: int = 10
    ld_r2_threshold: float = 0.02

    def __post_init__(self) -> None:
        self._validate_proportion(self.maf, "maf")
        self._validate_proportion(self.geno, "geno")
        self._validate_proportion(self.hwe, "hwe")
        self._validate_positive(float(self.ld_window_size), "ld_window_size")
        self._validate_positive(float(self.ld_step_size), "ld_step_size")
        self._validate_proportion(self.ld_r2_threshold, "ld_r2_threshold")

    @property
    def indep_pairwise_args(self) -> Tuple[int, int, float]:
        """Get PLINK --indep-pairwise arguments as tuple."""
        return (self.ld_window_size, self.ld_step_size, self.ld_r2_threshold)


@dataclass(frozen=True)
class PCAConfig(ThresholdConfig):
    """Configuration for PCA computation.

    Controls the number of principal components to compute and
    the variant filtering/pruning parameters.

    Attributes:
        n_pcs: Number of principal components to calculate.
            Default is 10.
        build: Genome build for exclusion regions ('hg19' or 'hg38').
            Default is 'hg38'.
        pruning: Variant pruning configuration. Uses default
            PCAPruningConfig if not specified.
    """

    n_pcs: int = 10
    build: Literal["hg19", "hg38"] = "hg38"
    pruning: PCAPruningConfig = field(default_factory=PCAPruningConfig)

    def __post_init__(self) -> None:
        self._validate_positive(float(self.n_pcs), "n_pcs")
        if self.build not in ("hg19", "hg38"):
            raise ValueError(f"build must be 'hg19' or 'hg38', got '{self.build}'")

    def get_exclusion_regions(self) -> str:
        """Get exclusion regions for this config's build."""
        return get_exclusion_regions(self.build)


# PLINK2 GLM options for GWAS
# These options control output columns and regression settings
PLINK2_GLM_OPTIONS: str = (
    "hide-covar "
    "firth-fallback "
    "no-x-sex "
    "cols=+a1freq,"
    "+a1freqcc,"
    "+a1count,"
    "+totallele,"
    "+a1countcc,"
    "+totallelecc,"
    "+err"
)


@dataclass(frozen=True)
class GWASConfig(ThresholdConfig):
    """Configuration for GWAS analysis.

    Controls association testing parameters including phenotype
    handling, covariate usage, and inflation calculation.

    Attributes:
        pheno_name: Name of phenotype column in psam file.
            Default is 'PHENO1'.
        use_covariates: Whether to include covariates in association
            testing. Default is True.
        covariate_variance_standardize: Whether to standardize
            covariates to zero mean and unit variance. Default is True.
        maf_lambdas: If True, calculate additional inflation metrics
            restricted to variants with MAF >= 0.01. Default is False.
        maf_threshold: MAF threshold for maf_lambdas filtering.
            Default is 0.01.
        glm_options: PLINK2 --glm options string.
            Default uses standard options with Firth fallback.
    """

    pheno_name: str = "PHENO1"
    use_covariates: bool = True
    covariate_variance_standardize: bool = True
    maf_lambdas: bool = False
    maf_threshold: float = 0.01
    glm_options: str = PLINK2_GLM_OPTIONS

    def __post_init__(self) -> None:
        self._validate_proportion(self.maf_threshold, "maf_threshold")


@dataclass(frozen=True)
class CovariateConfig(BaseConfig):
    """Configuration for covariates in GWAS.

    Specifies external covariate file and which columns to use.
    If not provided, PCA results can be used as covariates.

    Attributes:
        covar_path: Path to covariate file (tab-separated with
            #FID and IID columns). None means no external covariates.
        covar_names: Comma-separated list of covariate column names
            to use. None means use all columns except #FID and IID.
        use_pca_as_covariates: If True and covar_path is None, use
            PCA eigenvectors as covariates. Default is True.
    """

    covar_path: Optional[str] = None
    covar_names: Optional[str] = None
    use_pca_as_covariates: bool = True

    def has_external_covariates(self) -> bool:
        """Check if external covariates are configured."""
        return self.covar_path is not None


@dataclass(frozen=True)
class AssocConfig(BaseConfig):
    """Combined configuration for association analysis.

    Combines PCA and GWAS configurations with flags controlling
    which analyses to run.

    Attributes:
        pca: PCA configuration. Uses defaults if not specified.
        gwas: GWAS configuration. Uses defaults if not specified.
        covariates: Covariate configuration. Uses defaults if not specified.
        run_pca: Whether to run PCA. Default is True.
        run_gwas: Whether to run GWAS. Default is True.
    """

    pca: PCAConfig = field(default_factory=PCAConfig)
    gwas: GWASConfig = field(default_factory=GWASConfig)
    covariates: CovariateConfig = field(default_factory=CovariateConfig)
    run_pca: bool = True
    run_gwas: bool = True

    def __post_init__(self) -> None:
        if not self.run_pca and not self.run_gwas:
            raise ValueError(
                "At least one of run_pca or run_gwas must be True"
            )
