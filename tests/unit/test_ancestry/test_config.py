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

"""Tests for ancestry configuration classes."""

import pytest

from genotools.ancestry.config import (
    AdmixedConfig,
    AncestryConfig,
    ClassifierConfig,
    GridSearchConfig,
    InferenceConfig,
    InferenceMode,
    PCAConfig,
    TrainingConfig,
    UMAPConfig,
    DEFAULT_ANCESTRY_LABELS,
)


class TestPCAConfig:
    """Tests for PCAConfig."""

    def test_default_values(self) -> None:
        """Default values are set correctly."""
        config = PCAConfig()
        assert config.n_components == 50
        assert config.maf == 0.01
        assert config.geno == 0.1
        assert config.hwe == 5e-6
        assert config.random_state == 123

    def test_custom_values(self) -> None:
        """Custom values are accepted."""
        config = PCAConfig(n_components=20, maf=0.05)
        assert config.n_components == 20
        assert config.maf == 0.05

    def test_invalid_n_components_too_small(self) -> None:
        """n_components < 2 raises ValueError."""
        with pytest.raises(ValueError, match="n_components must be >= 2"):
            PCAConfig(n_components=1)

    def test_invalid_n_components_too_large(self) -> None:
        """n_components > 1000 raises ValueError."""
        with pytest.raises(ValueError, match="n_components must be <= 1000"):
            PCAConfig(n_components=1500)

    def test_invalid_maf_negative(self) -> None:
        """Negative maf raises ValueError."""
        with pytest.raises(ValueError, match="maf"):
            PCAConfig(maf=-0.1)

    def test_invalid_maf_too_large(self) -> None:
        """maf > 1 raises ValueError."""
        with pytest.raises(ValueError, match="maf"):
            PCAConfig(maf=1.5)

    def test_invalid_hwe_zero(self) -> None:
        """hwe = 0 raises ValueError."""
        with pytest.raises(ValueError, match="hwe must be in"):
            PCAConfig(hwe=0)

    def test_invalid_hwe_negative(self) -> None:
        """Negative hwe raises ValueError."""
        with pytest.raises(ValueError, match="hwe must be in"):
            PCAConfig(hwe=-1e-6)

    def test_frozen(self) -> None:
        """Config is immutable."""
        config = PCAConfig()
        with pytest.raises(AttributeError):
            config.n_components = 30  # type: ignore[misc]


class TestUMAPConfig:
    """Tests for UMAPConfig."""

    def test_default_values(self) -> None:
        """Default values are set correctly."""
        config = UMAPConfig()
        assert config.n_neighbors == 15
        assert config.n_components == 2
        assert config.min_dist == 0.1
        assert config.metric == "euclidean"
        assert config.a is None
        assert config.b is None
        assert config.random_state == 123

    def test_custom_values(self) -> None:
        """Custom values are accepted."""
        config = UMAPConfig(n_neighbors=10, n_components=3, a=1.0, b=0.5)
        assert config.n_neighbors == 10
        assert config.n_components == 3
        assert config.a == 1.0
        assert config.b == 0.5

    def test_invalid_n_neighbors_too_small(self) -> None:
        """n_neighbors < 2 raises ValueError."""
        with pytest.raises(ValueError, match="n_neighbors must be >= 2"):
            UMAPConfig(n_neighbors=1)

    def test_invalid_n_components_too_small(self) -> None:
        """n_components < 1 raises ValueError."""
        with pytest.raises(ValueError, match="n_components must be >= 1"):
            UMAPConfig(n_components=0)

    def test_invalid_a_not_positive(self) -> None:
        """a <= 0 raises ValueError."""
        with pytest.raises(ValueError, match="a must be positive"):
            UMAPConfig(a=0)

    def test_invalid_b_not_positive(self) -> None:
        """b <= 0 raises ValueError."""
        with pytest.raises(ValueError, match="b must be positive"):
            UMAPConfig(b=-0.5)


class TestClassifierConfig:
    """Tests for ClassifierConfig."""

    def test_default_values(self) -> None:
        """Default values are set correctly."""
        config = ClassifierConfig()
        assert config.n_estimators == 100
        assert config.max_depth == 6
        assert config.learning_rate == 0.1
        assert config.booster == "gblinear"
        assert config.reg_lambda == 1.0
        assert config.random_state == 123

    def test_custom_values(self) -> None:
        """Custom values are accepted."""
        config = ClassifierConfig(n_estimators=200, booster="gbtree")
        assert config.n_estimators == 200
        assert config.booster == "gbtree"

    def test_invalid_n_estimators(self) -> None:
        """n_estimators < 1 raises ValueError."""
        with pytest.raises(ValueError, match="n_estimators must be >= 1"):
            ClassifierConfig(n_estimators=0)

    def test_invalid_max_depth(self) -> None:
        """max_depth < 1 raises ValueError."""
        with pytest.raises(ValueError, match="max_depth must be >= 1"):
            ClassifierConfig(max_depth=0)

    def test_invalid_learning_rate(self) -> None:
        """Negative learning_rate raises ValueError."""
        with pytest.raises(ValueError, match="learning_rate must be positive"):
            ClassifierConfig(learning_rate=-0.1)

    def test_invalid_booster(self) -> None:
        """Invalid booster raises ValueError."""
        with pytest.raises(ValueError, match="booster must be"):
            ClassifierConfig(booster="invalid")

    def test_valid_boosters(self) -> None:
        """All valid boosters are accepted."""
        for booster in ("gbtree", "gblinear", "dart"):
            config = ClassifierConfig(booster=booster)
            assert config.booster == booster


class TestGridSearchConfig:
    """Tests for GridSearchConfig."""

    def test_default_values(self) -> None:
        """Default values are set correctly."""
        config = GridSearchConfig()
        assert config.umap_n_neighbors == (5, 20)
        assert config.umap_n_components == (15, 25)
        assert config.umap_a == (0.75, 1.0, 1.5)
        assert config.umap_b == (0.25, 0.5, 0.75)
        assert config.cv_folds == 5
        assert config.scoring == "balanced_accuracy"

    def test_invalid_cv_folds(self) -> None:
        """cv_folds < 2 raises ValueError."""
        with pytest.raises(ValueError, match="cv_folds must be >= 2"):
            GridSearchConfig(cv_folds=1)

    def test_empty_umap_n_neighbors(self) -> None:
        """Empty umap_n_neighbors raises ValueError."""
        with pytest.raises(ValueError, match="umap_n_neighbors cannot be empty"):
            GridSearchConfig(umap_n_neighbors=())


class TestAdmixedConfig:
    """Tests for AdmixedConfig."""

    def test_default_values(self) -> None:
        """Default values are set correctly."""
        config = AdmixedConfig()
        assert config.n_clusters_cas == 2
        assert config.birch_threshold is None

    def test_invalid_n_clusters(self) -> None:
        """n_clusters_cas < 1 raises ValueError."""
        with pytest.raises(ValueError, match="n_clusters_cas must be >= 1"):
            AdmixedConfig(n_clusters_cas=0)


class TestInferenceConfig:
    """Tests for InferenceConfig."""

    def test_default_values(self) -> None:
        """Default values are set correctly."""
        config = InferenceConfig()
        assert config.mode == InferenceMode.LOCAL
        assert config.container_image == "mkoretsky1/genotools_ancestry:python3.11"
        assert config.cloud_project is None
        assert config.cloud_region == "us-central1"

    def test_all_modes(self) -> None:
        """All inference modes are valid."""
        for mode in InferenceMode:
            if mode == InferenceMode.CLOUD:
                config = InferenceConfig(mode=mode, cloud_project="test-project")
            else:
                config = InferenceConfig(mode=mode)
            assert config.mode == mode

    def test_cloud_mode_requires_project(self) -> None:
        """Cloud mode without project raises ValueError."""
        with pytest.raises(ValueError, match="cloud_project is required"):
            InferenceConfig(mode=InferenceMode.CLOUD)

    def test_cloud_mode_with_project(self) -> None:
        """Cloud mode with project is valid."""
        config = InferenceConfig(
            mode=InferenceMode.CLOUD, cloud_project="my-project"
        )
        assert config.cloud_project == "my-project"


class TestTrainingConfig:
    """Tests for TrainingConfig."""

    def test_default_values(self) -> None:
        """Default values are set correctly."""
        config = TrainingConfig()
        assert config.test_size == 0.2
        assert config.random_state == 123
        assert config.n_jobs is None
        assert config.gb_per_worker == 3.0

    def test_invalid_test_size_zero(self) -> None:
        """test_size = 0 raises ValueError."""
        with pytest.raises(ValueError, match="test_size must be in"):
            TrainingConfig(test_size=0)

    def test_invalid_test_size_one(self) -> None:
        """test_size = 1 raises ValueError."""
        with pytest.raises(ValueError, match="test_size must be in"):
            TrainingConfig(test_size=1.0)

    def test_invalid_test_size_negative(self) -> None:
        """Negative test_size raises ValueError."""
        with pytest.raises(ValueError, match="test_size must be in"):
            TrainingConfig(test_size=-0.1)


class TestAncestryConfig:
    """Tests for AncestryConfig."""

    def test_default_values(self) -> None:
        """Default values are set correctly."""
        config = AncestryConfig()
        assert isinstance(config.pca, PCAConfig)
        assert isinstance(config.umap, UMAPConfig)
        assert isinstance(config.classifier, ClassifierConfig)
        assert isinstance(config.grid_search, GridSearchConfig)
        assert isinstance(config.admixed, AdmixedConfig)
        assert isinstance(config.inference, InferenceConfig)
        assert isinstance(config.training, TrainingConfig)
        assert config.labels == DEFAULT_ANCESTRY_LABELS
        assert config.min_samples_per_ancestry == 0

    def test_nested_config(self) -> None:
        """Nested configs can be customized."""
        config = AncestryConfig(
            pca=PCAConfig(n_components=20),
            umap=UMAPConfig(n_neighbors=10),
        )
        assert config.pca.n_components == 20
        assert config.umap.n_neighbors == 10

    def test_empty_labels(self) -> None:
        """Empty labels raises ValueError."""
        with pytest.raises(ValueError, match="labels cannot be empty"):
            AncestryConfig(labels=())

    def test_invalid_min_samples(self) -> None:
        """Negative min_samples_per_ancestry raises ValueError."""
        with pytest.raises(ValueError, match="min_samples_per_ancestry must be >= 0"):
            AncestryConfig(min_samples_per_ancestry=-1)

    def test_default_ancestry_labels(self) -> None:
        """Default ancestry labels are correct."""
        expected = (
            "AFR", "SAS", "EAS", "EUR", "AMR",
            "AJ", "CAS", "MDE", "FIN", "AAC",
        )
        assert DEFAULT_ANCESTRY_LABELS == expected


class TestInferenceMode:
    """Tests for InferenceMode enum."""

    def test_all_modes_have_values(self) -> None:
        """All modes have string values."""
        assert InferenceMode.LOCAL.value == "local"
        assert InferenceMode.CONTAINER.value == "container"
        assert InferenceMode.SINGULARITY.value == "singularity"
        assert InferenceMode.CLOUD.value == "cloud"


class TestConfigImmutability:
    """Tests that all configs are immutable (frozen)."""

    def test_all_configs_frozen(self) -> None:
        """All config classes are frozen dataclasses."""
        configs = [
            PCAConfig(),
            UMAPConfig(),
            ClassifierConfig(),
            GridSearchConfig(),
            AdmixedConfig(),
            InferenceConfig(),
            TrainingConfig(),
            AncestryConfig(),
        ]

        for config in configs:
            # Get first attribute
            attr = list(config.__dataclass_fields__.keys())[0]
            with pytest.raises(AttributeError):
                setattr(config, attr, "new_value")
