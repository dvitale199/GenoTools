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

"""Ancestry prediction module for GenoTools.

This module provides a clean ML-pattern interface for ancestry prediction
using PCA, UMAP, and XGBoost classification.

Basic Usage:
    >>> from genotools.ancestry import AncestryModel, ReferencePanel, AncestryConfig
    >>>
    >>> # Load reference panel
    >>> ref = ReferencePanel.load(
    ...     geno_path=Path("/path/to/ref_panel"),
    ...     labels_path=Path("/path/to/labels.txt")
    ... )
    >>>
    >>> # Create and train model
    >>> model = AncestryModel()
    >>> model.fit(ref)
    >>>
    >>> # Predict ancestry
    >>> predictions = model.predict(test_data)
    >>> print(predictions.ancestry_counts)
    {'EUR': 100, 'AFR': 50, ...}
    >>>
    >>> # Save model for later use
    >>> model.save(Path("model.pkl"))

Configuration:
    >>> from genotools.ancestry import AncestryConfig, PCAConfig
    >>>
    >>> # Custom configuration
    >>> config = AncestryConfig(
    ...     pca=PCAConfig(n_components=30),
    ...     min_samples_per_ancestry=50
    ... )
    >>> model = AncestryModel(config=config)

Inference Modes:
    The module supports multiple inference modes for different environments:
    - LOCAL: In-process Python execution (default)
    - CONTAINER: Docker container with bundled model
    - SINGULARITY: Singularity container for HPC environments
    - CLOUD: Google Cloud AI Platform

    >>> from genotools.ancestry import InferenceMode, InferenceConfig
    >>> config = AncestryConfig(
    ...     inference=InferenceConfig(mode=InferenceMode.CONTAINER)
    ... )

UMAP Version Note:
    CRITICAL: This module requires umap-learn==0.5.3 for compatibility with
    existing trained models. Models trained with different UMAP versions may
    fail to load or produce incorrect predictions.
"""

# Configuration classes
from genotools.ancestry.config import (
    AncestryConfig,
    PCAConfig,
    UMAPConfig,
    ClassifierConfig,
    GridSearchConfig,
    AdmixedConfig,
    InferenceConfig,
    TrainingConfig,
    InferenceMode,
    DEFAULT_ANCESTRY_LABELS,
)

# Result classes
from genotools.ancestry.results import (
    AncestryPredictions,
    AncestryResult,
    TrainingMetrics,
    PCAResult,
    UMAPResult,
    SplitResult,
)

# Reference panel management
from genotools.ancestry.reference import (
    ReferencePanel,
    get_default_model_path,
    validate_model_files,
)

# Main model class
from genotools.ancestry.model import (
    AncestryModel,
    load_trained_pipeline,
)

# Reducers (for advanced usage)
from genotools.ancestry.reducers import (
    PCAReducer,
    UMAPReducer,
    flashpca_scale,
    run_pca,
    run_umap,
)


__all__ = [
    # Main classes
    "AncestryModel",
    "ReferencePanel",
    # Configuration
    "AncestryConfig",
    "PCAConfig",
    "UMAPConfig",
    "ClassifierConfig",
    "GridSearchConfig",
    "AdmixedConfig",
    "InferenceConfig",
    "TrainingConfig",
    "InferenceMode",
    "DEFAULT_ANCESTRY_LABELS",
    # Results
    "AncestryPredictions",
    "AncestryResult",
    "TrainingMetrics",
    "PCAResult",
    "UMAPResult",
    "SplitResult",
    # Utilities
    "get_default_model_path",
    "validate_model_files",
    "load_trained_pipeline",
    # Reducers
    "PCAReducer",
    "UMAPReducer",
    "flashpca_scale",
    "run_pca",
    "run_umap",
]
