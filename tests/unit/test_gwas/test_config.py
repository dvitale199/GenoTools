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

"""Tests for GWAS configuration classes."""

import pytest

from genotools.gwas.config import (
    PCAPruningConfig,
    PCAConfig,
    GWASConfig,
    CovariateConfig,
    AssocConfig,
    get_exclusion_regions,
    HG19_EXCLUSION_REGIONS,
    HG38_EXCLUSION_REGIONS,
)


class TestPCAPruningConfig:
    """Tests for PCAPruningConfig."""

    def test_default_values(self) -> None:
        """Default values are set correctly."""
        config = PCAPruningConfig()
        assert config.maf == 0.01
        assert config.geno == 0.01
        assert config.hwe == 5e-6
        assert config.ld_window_size == 1000
        assert config.ld_step_size == 10
        assert config.ld_r2_threshold == 0.02

    def test_custom_values(self) -> None:
        """Custom values are accepted."""
        config = PCAPruningConfig(
            maf=0.05,
            geno=0.02,
            hwe=1e-5,
            ld_window_size=500,
            ld_step_size=5,
            ld_r2_threshold=0.1,
        )
        assert config.maf == 0.05
        assert config.geno == 0.02
        assert config.hwe == 1e-5
        assert config.ld_window_size == 500
        assert config.ld_step_size == 5
        assert config.ld_r2_threshold == 0.1

    def test_indep_pairwise_args(self) -> None:
        """indep_pairwise_args property returns correct tuple."""
        config = PCAPruningConfig(ld_window_size=500, ld_step_size=5, ld_r2_threshold=0.1)
        assert config.indep_pairwise_args == (500, 5, 0.1)

    def test_invalid_maf_negative(self) -> None:
        """Negative maf raises ValueError."""
        with pytest.raises(ValueError, match="maf"):
            PCAPruningConfig(maf=-0.01)

    def test_invalid_maf_too_large(self) -> None:
        """maf > 1 raises ValueError."""
        with pytest.raises(ValueError, match="maf"):
            PCAPruningConfig(maf=1.5)

    def test_invalid_geno_negative(self) -> None:
        """Negative geno raises ValueError."""
        with pytest.raises(ValueError, match="geno"):
            PCAPruningConfig(geno=-0.01)

    def test_invalid_hwe_negative(self) -> None:
        """Negative hwe raises ValueError."""
        with pytest.raises(ValueError, match="hwe"):
            PCAPruningConfig(hwe=-1e-6)

    def test_invalid_ld_window_size_negative(self) -> None:
        """Negative ld_window_size raises ValueError."""
        with pytest.raises(ValueError, match="ld_window_size"):
            PCAPruningConfig(ld_window_size=-100)

    def test_invalid_ld_step_size_zero(self) -> None:
        """Zero ld_step_size raises ValueError."""
        with pytest.raises(ValueError, match="ld_step_size"):
            PCAPruningConfig(ld_step_size=0)

    def test_invalid_ld_r2_threshold_too_large(self) -> None:
        """ld_r2_threshold > 1 raises ValueError."""
        with pytest.raises(ValueError, match="ld_r2_threshold"):
            PCAPruningConfig(ld_r2_threshold=1.5)

    def test_frozen(self) -> None:
        """Config is immutable."""
        config = PCAPruningConfig()
        with pytest.raises(AttributeError):
            config.maf = 0.05  # type: ignore[misc]


class TestPCAConfig:
    """Tests for PCAConfig."""

    def test_default_values(self) -> None:
        """Default values are set correctly."""
        config = PCAConfig()
        assert config.n_pcs == 10
        assert config.build == "hg38"
        assert isinstance(config.pruning, PCAPruningConfig)

    def test_custom_n_pcs(self) -> None:
        """Custom n_pcs is accepted."""
        config = PCAConfig(n_pcs=20)
        assert config.n_pcs == 20

    def test_custom_build(self) -> None:
        """Custom build is accepted."""
        config = PCAConfig(build="hg19")
        assert config.build == "hg19"

    def test_custom_pruning(self) -> None:
        """Custom pruning config is accepted."""
        pruning = PCAPruningConfig(maf=0.05)
        config = PCAConfig(pruning=pruning)
        assert config.pruning.maf == 0.05

    def test_invalid_n_pcs_zero(self) -> None:
        """Zero n_pcs raises ValueError."""
        with pytest.raises(ValueError, match="n_pcs"):
            PCAConfig(n_pcs=0)

    def test_invalid_n_pcs_negative(self) -> None:
        """Negative n_pcs raises ValueError."""
        with pytest.raises(ValueError, match="n_pcs"):
            PCAConfig(n_pcs=-5)

    def test_invalid_build(self) -> None:
        """Invalid build raises ValueError."""
        with pytest.raises(ValueError, match="build"):
            PCAConfig(build="hg37")  # type: ignore[arg-type]

    def test_get_exclusion_regions_hg38(self) -> None:
        """get_exclusion_regions returns correct regions for hg38."""
        config = PCAConfig(build="hg38")
        regions = config.get_exclusion_regions()
        assert "6 24999772 35032223" in regions  # MHC region

    def test_get_exclusion_regions_hg19(self) -> None:
        """get_exclusion_regions returns correct regions for hg19."""
        config = PCAConfig(build="hg19")
        regions = config.get_exclusion_regions()
        assert "6 25000000 33500000" in regions  # MHC region

    def test_frozen(self) -> None:
        """Config is immutable."""
        config = PCAConfig()
        with pytest.raises(AttributeError):
            config.n_pcs = 20  # type: ignore[misc]


class TestGWASConfig:
    """Tests for GWASConfig."""

    def test_default_values(self) -> None:
        """Default values are set correctly."""
        config = GWASConfig()
        assert config.pheno_name == "PHENO1"
        assert config.use_covariates is True
        assert config.covariate_variance_standardize is True
        assert config.maf_lambdas is False
        assert config.maf_threshold == 0.01
        assert "firth-fallback" in config.glm_options

    def test_custom_pheno_name(self) -> None:
        """Custom pheno_name is accepted."""
        config = GWASConfig(pheno_name="MY_PHENO")
        assert config.pheno_name == "MY_PHENO"

    def test_enable_maf_lambdas(self) -> None:
        """maf_lambdas can be enabled."""
        config = GWASConfig(maf_lambdas=True)
        assert config.maf_lambdas is True

    def test_custom_maf_threshold(self) -> None:
        """Custom maf_threshold is accepted."""
        config = GWASConfig(maf_threshold=0.05)
        assert config.maf_threshold == 0.05

    def test_invalid_maf_threshold_negative(self) -> None:
        """Negative maf_threshold raises ValueError."""
        with pytest.raises(ValueError, match="maf_threshold"):
            GWASConfig(maf_threshold=-0.01)

    def test_invalid_maf_threshold_too_large(self) -> None:
        """maf_threshold > 1 raises ValueError."""
        with pytest.raises(ValueError, match="maf_threshold"):
            GWASConfig(maf_threshold=1.5)

    def test_frozen(self) -> None:
        """Config is immutable."""
        config = GWASConfig()
        with pytest.raises(AttributeError):
            config.pheno_name = "OTHER"  # type: ignore[misc]


class TestCovariateConfig:
    """Tests for CovariateConfig."""

    def test_default_values(self) -> None:
        """Default values are set correctly."""
        config = CovariateConfig()
        assert config.covar_path is None
        assert config.covar_names is None
        assert config.use_pca_as_covariates is True

    def test_with_covar_path(self) -> None:
        """Covariate path can be set."""
        config = CovariateConfig(covar_path="/path/to/covars.txt")
        assert config.covar_path == "/path/to/covars.txt"

    def test_with_covar_names(self) -> None:
        """Covariate names can be set."""
        config = CovariateConfig(
            covar_path="/path/to/covars.txt",
            covar_names="PC1,PC2,AGE",
        )
        assert config.covar_names == "PC1,PC2,AGE"

    def test_has_external_covariates_true(self) -> None:
        """has_external_covariates returns True when path is set."""
        config = CovariateConfig(covar_path="/path/to/covars.txt")
        assert config.has_external_covariates() is True

    def test_has_external_covariates_false(self) -> None:
        """has_external_covariates returns False when path is None."""
        config = CovariateConfig()
        assert config.has_external_covariates() is False

    def test_disable_pca_as_covariates(self) -> None:
        """use_pca_as_covariates can be disabled."""
        config = CovariateConfig(use_pca_as_covariates=False)
        assert config.use_pca_as_covariates is False


class TestAssocConfig:
    """Tests for AssocConfig."""

    def test_default_values(self) -> None:
        """Default values are set correctly."""
        config = AssocConfig()
        assert config.run_pca is True
        assert config.run_gwas is True
        assert isinstance(config.pca, PCAConfig)
        assert isinstance(config.gwas, GWASConfig)
        assert isinstance(config.covariates, CovariateConfig)

    def test_run_pca_only(self) -> None:
        """Can run PCA only."""
        config = AssocConfig(run_pca=True, run_gwas=False)
        assert config.run_pca is True
        assert config.run_gwas is False

    def test_run_gwas_only(self) -> None:
        """Can run GWAS only."""
        config = AssocConfig(run_pca=False, run_gwas=True)
        assert config.run_pca is False
        assert config.run_gwas is True

    def test_invalid_neither_enabled(self) -> None:
        """Must enable at least one of run_pca or run_gwas."""
        with pytest.raises(ValueError, match="At least one"):
            AssocConfig(run_pca=False, run_gwas=False)

    def test_custom_pca_config(self) -> None:
        """Custom PCA config is accepted."""
        pca = PCAConfig(n_pcs=20)
        config = AssocConfig(pca=pca)
        assert config.pca.n_pcs == 20

    def test_custom_gwas_config(self) -> None:
        """Custom GWAS config is accepted."""
        gwas = GWASConfig(maf_lambdas=True)
        config = AssocConfig(gwas=gwas)
        assert config.gwas.maf_lambdas is True


class TestExclusionRegions:
    """Tests for exclusion region functions and constants."""

    def test_hg19_exclusion_regions(self) -> None:
        """HG19 exclusion regions contain expected data."""
        assert "5 44000000 51500000" in HG19_EXCLUSION_REGIONS
        assert "6 25000000 33500000" in HG19_EXCLUSION_REGIONS  # MHC

    def test_hg38_exclusion_regions(self) -> None:
        """HG38 exclusion regions contain expected data."""
        assert "6 24999772 35032223" in HG38_EXCLUSION_REGIONS  # MHC
        assert "8 7142478 13142491" in HG38_EXCLUSION_REGIONS

    def test_get_exclusion_regions_hg19(self) -> None:
        """get_exclusion_regions returns hg19 regions."""
        regions = get_exclusion_regions("hg19")
        assert regions == HG19_EXCLUSION_REGIONS

    def test_get_exclusion_regions_hg38(self) -> None:
        """get_exclusion_regions returns hg38 regions."""
        regions = get_exclusion_regions("hg38")
        assert regions == HG38_EXCLUSION_REGIONS

    def test_get_exclusion_regions_invalid_build(self) -> None:
        """get_exclusion_regions raises ValueError for invalid build."""
        with pytest.raises(ValueError, match="Invalid build"):
            get_exclusion_regions("hg37")


class TestConfigImmutability:
    """Tests that all configs are immutable (frozen)."""

    def test_all_configs_frozen(self) -> None:
        """All config classes are frozen dataclasses."""
        configs = [
            PCAPruningConfig(),
            PCAConfig(),
            GWASConfig(),
            CovariateConfig(),
            AssocConfig(),
        ]

        for config in configs:
            # Try to modify an attribute - should raise
            attr = list(config.__dataclass_fields__.keys())[0]
            with pytest.raises(AttributeError):
                setattr(config, attr, "new_value")
