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

"""Ancestry prediction configuration dataclasses.

This module provides typed configuration for the ancestry prediction pipeline,
including PCA, UMAP dimensionality reduction, and XGBoost classification.

CRITICAL: umap-learn must be pinned to version 0.5.3 for model compatibility.
Old pickled models may fail with newer UMAP versions.
See requirements.txt: umap-learn==0.5.3

Example:
    >>> from genotools.ancestry.config import AncestryConfig
    >>> config = AncestryConfig()
    >>> config.pca.n_components
    50
    >>> # Custom configuration
    >>> custom_config = AncestryConfig(
    ...     pca=PCAConfig(n_components=20),
    ...     classifier=ClassifierConfig(n_estimators=200)
    ... )
"""

from dataclasses import dataclass, field
from enum import Enum
from typing import Optional, Tuple

from genotools.core.config import BaseConfig, ThresholdConfig


class InferenceMode(Enum):
    """Ancestry inference execution modes.

    Determines where the prediction computation runs:
    - LOCAL: In-process Python execution (default)
    - CONTAINER: Docker container with bundled model
    - SINGULARITY: Singularity container for HPC environments
    - CLOUD: Google Cloud AI Platform prediction
    """

    LOCAL = "local"
    CONTAINER = "container"
    SINGULARITY = "singularity"
    CLOUD = "cloud"


# Default ancestry labels supported by GenoTools
# These are the 10 labels from the current implementation
DEFAULT_ANCESTRY_LABELS: Tuple[str, ...] = (
    "AFR",  # African
    "SAS",  # South Asian
    "EAS",  # East Asian
    "EUR",  # European
    "AMR",  # American/Latino (admixed)
    "AJ",   # Ashkenazi Jewish
    "CAS",  # Central Asian
    "MDE",  # Middle Eastern
    "FIN",  # Finnish
    "AAC",  # African American (admixed)
)


@dataclass(frozen=True)
class PCAConfig(ThresholdConfig):
    """Configuration for PCA dimensionality reduction.

    PCA is the first step in the ancestry prediction pipeline. It reduces
    the high-dimensional genotype data to a smaller set of principal
    components that capture population structure.

    The default of 50 components provides sufficient resolution for
    distinguishing 10 ancestry groups while maintaining computational
    efficiency.

    Attributes:
        n_components: Number of principal components to compute.
            More components capture more variance but increase computation
            time. Default is 50 (matching current implementation).
        maf: Minor allele frequency threshold for variant filtering before PCA.
            Variants with MAF < threshold are excluded. Default is 0.01.
        geno: Variant missingness threshold. Variants with missing rate
            > threshold are excluded before PCA. Default is 0.1.
        hwe: Hardy-Weinberg equilibrium p-value threshold. Variants with
            HWE p < threshold are excluded. Default is 5e-6.
        random_state: Random seed for reproducibility. Default is 123.
    """

    n_components: int = 50
    maf: float = 0.01
    geno: float = 0.1
    hwe: float = 5e-6
    random_state: int = 123

    def __post_init__(self) -> None:
        """Validate PCA configuration."""
        if self.n_components < 2:
            raise ValueError(
                f"n_components must be >= 2, got {self.n_components}"
            )
        if self.n_components > 1000:
            raise ValueError(
                f"n_components must be <= 1000, got {self.n_components}"
            )
        self._validate_proportion(self.maf, "maf")
        self._validate_proportion(self.geno, "geno")
        # HWE can be very small, so allow scientific notation range
        if self.hwe <= 0 or self.hwe > 1:
            raise ValueError(
                f"hwe must be in (0, 1], got {self.hwe}"
            )


@dataclass(frozen=True)
class UMAPConfig(ThresholdConfig):
    """Configuration for UMAP dimensionality reduction.

    UMAP (Uniform Manifold Approximation and Projection) is used to further
    reduce the PCA components while preserving local structure. This makes
    the ancestry clusters more separable for classification.

    CRITICAL: Requires umap-learn==0.5.3 for model compatibility. Models
    trained with one UMAP version may not load correctly with another.

    The parameters below are the GridSearchCV hyperparameter ranges from
    the current implementation. Final values are determined during training.

    Attributes:
        n_neighbors: Number of neighbors for local structure. Higher values
            capture more global structure. Default is 15 (UMAP default).
        n_components: Output dimensionality. Default is 2 for visualization,
            but grid search may select higher values.
        min_dist: Minimum distance between points in embedding. Lower values
            allow tighter clusters. Default is 0.1.
        metric: Distance metric. Default is "euclidean".
        a: UMAP curve parameter controlling how tight clusters form.
            Default is None (UMAP calculates from spread).
        b: UMAP curve parameter controlling cluster density.
            Default is None (UMAP calculates from min_dist).
        random_state: Random seed for reproducibility. Default is 123.
    """

    n_neighbors: int = 15
    n_components: int = 2
    min_dist: float = 0.1
    metric: str = "euclidean"
    a: Optional[float] = None
    b: Optional[float] = None
    random_state: int = 123

    def __post_init__(self) -> None:
        """Validate UMAP configuration."""
        if self.n_neighbors < 2:
            raise ValueError(
                f"n_neighbors must be >= 2, got {self.n_neighbors}"
            )
        if self.n_components < 1:
            raise ValueError(
                f"n_components must be >= 1, got {self.n_components}"
            )
        self._validate_proportion(self.min_dist, "min_dist", 0.0, 1.0)
        if self.a is not None and self.a <= 0:
            raise ValueError(f"a must be positive, got {self.a}")
        if self.b is not None and self.b <= 0:
            raise ValueError(f"b must be positive, got {self.b}")


@dataclass(frozen=True)
class ClassifierConfig(ThresholdConfig):
    """Configuration for XGBoost classifier.

    The classifier predicts ancestry labels from the UMAP embedding.
    Uses XGBoost with linear booster (gblinear) for efficiency and
    interpretability.

    Hyperparameter ranges from current implementation's GridSearchCV:
    - lambda: [10^-3, 10^-2, 10^-1, 1, 10, 100]
    - n_neighbors (UMAP): [5, 20]
    - n_components (UMAP): [15, 25]
    - a (UMAP): [0.75, 1.0, 1.5]
    - b (UMAP): [0.25, 0.5, 0.75]

    Attributes:
        n_estimators: Number of boosting rounds. Default is 100.
        max_depth: Maximum tree depth. Not used with gblinear booster.
            Default is 6 for tree boosters.
        learning_rate: Step size shrinkage. Default is 0.1.
        booster: Boosting algorithm. Default is "gblinear" for linear
            booster matching current implementation.
        reg_lambda: L2 regularization term. Default is 1.0.
        random_state: Random seed for reproducibility. Default is 123.
    """

    n_estimators: int = 100
    max_depth: int = 6
    learning_rate: float = 0.1
    booster: str = "gblinear"
    reg_lambda: float = 1.0
    random_state: int = 123

    def __post_init__(self) -> None:
        """Validate classifier configuration."""
        if self.n_estimators < 1:
            raise ValueError(
                f"n_estimators must be >= 1, got {self.n_estimators}"
            )
        if self.max_depth < 1:
            raise ValueError(
                f"max_depth must be >= 1, got {self.max_depth}"
            )
        self._validate_positive(self.learning_rate, "learning_rate")
        if self.booster not in ("gbtree", "gblinear", "dart"):
            raise ValueError(
                f"booster must be 'gbtree', 'gblinear', or 'dart', "
                f"got '{self.booster}'"
            )


@dataclass(frozen=True)
class GridSearchConfig(BaseConfig):
    """Configuration for hyperparameter grid search.

    Controls the grid search parameters used during model training.
    These are the exact parameter ranges from the current implementation.

    Attributes:
        umap_n_neighbors: Grid values for UMAP n_neighbors. Default is (5, 20).
        umap_n_components: Grid values for UMAP n_components. Default is (15, 25).
        umap_a: Grid values for UMAP a parameter. Default is (0.75, 1.0, 1.5).
        umap_b: Grid values for UMAP b parameter. Default is (0.25, 0.5, 0.75).
        xgb_lambda: Grid values for XGBoost L2 regularization.
            Default is powers of 10 from 10^-3 to 10^2.
        cv_folds: Number of cross-validation folds. Default is 5.
        scoring: Scoring metric for grid search. Default is "balanced_accuracy".
    """

    umap_n_neighbors: Tuple[int, ...] = (5, 20)
    umap_n_components: Tuple[int, ...] = (15, 25)
    umap_a: Tuple[float, ...] = (0.75, 1.0, 1.5)
    umap_b: Tuple[float, ...] = (0.25, 0.5, 0.75)
    xgb_lambda: Tuple[float, ...] = (0.001, 0.01, 0.1, 1.0, 10.0, 100.0)
    cv_folds: int = 5
    scoring: str = "balanced_accuracy"

    def __post_init__(self) -> None:
        """Validate grid search configuration."""
        if self.cv_folds < 2:
            raise ValueError(
                f"cv_folds must be >= 2, got {self.cv_folds}"
            )
        if not self.umap_n_neighbors:
            raise ValueError("umap_n_neighbors cannot be empty")
        if not self.umap_n_components:
            raise ValueError("umap_n_components cannot be empty")


@dataclass(frozen=True)
class AdmixedConfig(ThresholdConfig):
    """Configuration for admixed sample detection.

    Samples that don't clearly belong to a single ancestry group are
    classified as "Complex Admixture History" (CAH). This is determined
    by measuring the distance to each ancestry's centroid in PC space.

    Samples closest to the global centroid (ALL) rather than any specific
    ancestry centroid are classified as CAH.

    Attributes:
        n_clusters_cas: Number of Birch clusters for CAS (Central Asian)
            splitting. CAS shows bimodal distribution, so it's split into
            CAS and CAS2 subclusters. Default is 2.
        birch_threshold: Birch clustering threshold. Default is None
            (use sklearn default).
    """

    n_clusters_cas: int = 2
    birch_threshold: Optional[float] = None

    def __post_init__(self) -> None:
        """Validate admixed detection configuration."""
        if self.n_clusters_cas < 1:
            raise ValueError(
                f"n_clusters_cas must be >= 1, got {self.n_clusters_cas}"
            )


@dataclass(frozen=True)
class InferenceConfig(BaseConfig):
    """Configuration for ancestry inference execution.

    Controls where and how the prediction computation runs. Supports
    local execution, Docker containers, Singularity containers (for HPC),
    and Google Cloud AI Platform.

    Attributes:
        mode: Execution mode (local, container, singularity, cloud).
            Default is LOCAL.
        container_image: Docker/Singularity image name for container modes.
            Default is the official GenoTools ancestry image.
        cloud_project: Google Cloud project ID for cloud mode.
        cloud_region: Google Cloud region for cloud mode.
            Default is "us-central1".
        cloud_endpoint: Cloud AI Platform endpoint ID for cloud mode.
    """

    mode: InferenceMode = InferenceMode.LOCAL
    container_image: str = "mkoretsky1/genotools_ancestry:python3.11"
    cloud_project: Optional[str] = None
    cloud_region: str = "us-central1"
    cloud_endpoint: Optional[str] = None

    def __post_init__(self) -> None:
        """Validate inference configuration."""
        if self.mode == InferenceMode.CLOUD:
            if not self.cloud_project:
                raise ValueError(
                    "cloud_project is required when mode is CLOUD"
                )


@dataclass(frozen=True)
class TrainingConfig(BaseConfig):
    """Configuration for model training.

    Controls the training process including train/test split and
    resource allocation.

    Attributes:
        test_size: Proportion of reference samples for testing.
            Default is 0.2 (20% test, 80% train).
        random_state: Random seed for train/test split. Default is 123.
        n_jobs: Number of parallel jobs for grid search. None means
            auto-detect based on available resources. Default is None.
        gb_per_worker: Memory per worker for job calculation.
            Default is 3 GB.
    """

    test_size: float = 0.2
    random_state: int = 123
    n_jobs: Optional[int] = None
    gb_per_worker: float = 3.0

    def __post_init__(self) -> None:
        """Validate training configuration."""
        if self.test_size <= 0 or self.test_size >= 1:
            raise ValueError(
                f"test_size must be in (0, 1), got {self.test_size}"
            )


@dataclass(frozen=True)
class AncestryConfig(BaseConfig):
    """Full ancestry prediction configuration.

    This is the top-level configuration class that combines all
    sub-configurations for the ancestry prediction pipeline.

    Example:
        >>> config = AncestryConfig()
        >>> config.pca.n_components
        50
        >>> # With custom settings
        >>> config = AncestryConfig(
        ...     pca=PCAConfig(n_components=20),
        ...     inference=InferenceConfig(mode=InferenceMode.CONTAINER)
        ... )

    Attributes:
        pca: PCA dimensionality reduction configuration.
        umap: UMAP configuration (defaults used if grid searching).
        classifier: XGBoost classifier configuration.
        grid_search: Hyperparameter grid search configuration.
        admixed: Admixed sample detection configuration.
        inference: Inference execution configuration.
        training: Model training configuration.
        labels: Supported ancestry labels. Default includes all 10 labels.
        min_samples_per_ancestry: Minimum samples required per ancestry
            to include in output. Groups below this are pruned.
            Default is 0 (include all).
    """

    pca: PCAConfig = field(default_factory=PCAConfig)
    umap: UMAPConfig = field(default_factory=UMAPConfig)
    classifier: ClassifierConfig = field(default_factory=ClassifierConfig)
    grid_search: GridSearchConfig = field(default_factory=GridSearchConfig)
    admixed: AdmixedConfig = field(default_factory=AdmixedConfig)
    inference: InferenceConfig = field(default_factory=InferenceConfig)
    training: TrainingConfig = field(default_factory=TrainingConfig)
    labels: Tuple[str, ...] = DEFAULT_ANCESTRY_LABELS
    min_samples_per_ancestry: int = 0

    def __post_init__(self) -> None:
        """Validate ancestry configuration."""
        if not self.labels:
            raise ValueError("labels cannot be empty")
        if self.min_samples_per_ancestry < 0:
            raise ValueError(
                f"min_samples_per_ancestry must be >= 0, "
                f"got {self.min_samples_per_ancestry}"
            )
