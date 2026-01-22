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

"""UMAP dimensionality reduction for ancestry prediction.

This module provides UMAP (Uniform Manifold Approximation and Projection)
transformation for the ancestry prediction pipeline. UMAP is used to
reduce PCA components to a lower-dimensional space that preserves local
structure, making ancestry clusters more separable.

CRITICAL: Requires umap-learn==0.5.3 for model compatibility.
Models trained with one UMAP version may not load correctly with another.

Example:
    >>> from genotools.ancestry.reducers.umap_reducer import UMAPReducer
    >>> reducer = UMAPReducer(n_neighbors=15, n_components=2)
    >>> reducer.fit(X_train_pca)
    >>> X_train_umap = reducer.transform(X_train_pca)
    >>> X_new_umap = reducer.transform(X_new_pca)
"""

from dataclasses import dataclass, field
from typing import Any, Dict, Optional

import numpy as np
import pandas as pd
from umap import UMAP

from genotools.ancestry.config import UMAPConfig
from genotools.core.exceptions import AncestryError


@dataclass
class UMAPReducer:
    """UMAP dimensionality reducer for ancestry prediction.

    This class wraps the UMAP library with a fit/transform interface
    for use in the ancestry pipeline. It can be used standalone or
    as part of a sklearn Pipeline with XGBoost classifier.

    CRITICAL: Requires umap-learn==0.5.3 for compatibility with
    existing trained models.

    Attributes:
        config: UMAP configuration parameters.
        umap: Fitted UMAP object.
    """

    config: UMAPConfig = field(default_factory=UMAPConfig)
    umap: Optional[UMAP] = field(default=None, repr=False)
    _is_fitted: bool = field(default=False, repr=False)

    @property
    def is_fitted(self) -> bool:
        """Whether the reducer has been fitted."""
        return self._is_fitted

    @property
    def n_components(self) -> int:
        """Output dimensionality."""
        return self.config.n_components

    def _create_umap(self) -> UMAP:
        """Create UMAP instance from config.

        Returns:
            Configured UMAP instance.
        """
        kwargs: Dict[str, Any] = {
            "n_neighbors": self.config.n_neighbors,
            "n_components": self.config.n_components,
            "min_dist": self.config.min_dist,
            "metric": self.config.metric,
            "random_state": self.config.random_state,
        }

        # Only add a/b parameters if explicitly set
        if self.config.a is not None:
            kwargs["a"] = self.config.a
        if self.config.b is not None:
            kwargs["b"] = self.config.b

        return UMAP(**kwargs)

    def fit(self, X: np.ndarray) -> "UMAPReducer":  # type: ignore[type-arg]
        """Fit UMAP on input data.

        Args:
            X: Input data (samples x features), typically PCA components.

        Returns:
            Self for method chaining.

        Raises:
            AncestryError: If input data is too small.
        """
        if X.shape[0] < self.config.n_neighbors:
            raise AncestryError(
                f"Input has {X.shape[0]} samples but UMAP requires at least "
                f"{self.config.n_neighbors} samples (n_neighbors). "
                f"Reduce n_neighbors or use more samples."
            )

        self.umap = self._create_umap()
        self.umap.fit(X)
        self._is_fitted = True
        return self

    def transform(
        self,
        X: np.ndarray,  # type: ignore[type-arg]
        return_dataframe: bool = False,
    ) -> np.ndarray | pd.DataFrame:  # type: ignore[type-arg]
        """Transform data using fitted UMAP.

        Args:
            X: Input data (samples x features) to transform.
            return_dataframe: If True, return DataFrame with UMAP column names.

        Returns:
            UMAP-transformed data (samples x n_components).

        Raises:
            AncestryError: If UMAP not fitted.
        """
        if not self._is_fitted or self.umap is None:
            raise AncestryError("UMAPReducer must be fitted before transform")

        X_umap = self.umap.transform(X)

        if return_dataframe:
            columns = [f"UMAP{i+1}" for i in range(self.config.n_components)]
            return pd.DataFrame(X_umap, columns=columns)
        return X_umap

    def fit_transform(
        self,
        X: np.ndarray,  # type: ignore[type-arg]
        return_dataframe: bool = False,
    ) -> np.ndarray | pd.DataFrame:  # type: ignore[type-arg]
        """Fit UMAP and transform data in one step.

        Args:
            X: Input data (samples x features).
            return_dataframe: If True, return DataFrame with UMAP column names.

        Returns:
            UMAP-transformed data (samples x n_components).
        """
        self.fit(X)
        return self.transform(X, return_dataframe=return_dataframe)

    @classmethod
    def from_params(
        cls,
        params: Dict[str, Any],
        random_state: int = 123,
    ) -> "UMAPReducer":
        """Create UMAPReducer from parameter dictionary.

        This is useful for creating a reducer with parameters from
        a grid search result (e.g., best_params_ from GridSearchCV).

        Args:
            params: Dictionary with UMAP parameters. Keys should be
                prefixed with "umap__" (sklearn Pipeline convention)
                or be plain parameter names.
            random_state: Random seed for reproducibility.

        Returns:
            Configured UMAPReducer.
        """
        # Handle sklearn Pipeline parameter naming convention
        clean_params: Dict[str, Any] = {}
        for key, value in params.items():
            if key.startswith("umap__"):
                clean_key = key[6:]  # Remove "umap__" prefix
                clean_params[clean_key] = value
            elif key in ("n_neighbors", "n_components", "min_dist", "metric", "a", "b"):
                clean_params[key] = value

        # Build config
        config = UMAPConfig(
            n_neighbors=clean_params.get("n_neighbors", 15),
            n_components=clean_params.get("n_components", 2),
            min_dist=clean_params.get("min_dist", 0.1),
            metric=clean_params.get("metric", "euclidean"),
            a=clean_params.get("a"),
            b=clean_params.get("b"),
            random_state=random_state,
        )

        return cls(config=config)


def run_umap(
    ref_pca: pd.DataFrame,
    new_pca: pd.DataFrame,
    new_labels: pd.Series,  # type: ignore[type-arg]
    params: Optional[Dict[str, Any]] = None,
    config: Optional[UMAPConfig] = None,
) -> Dict[str, Any]:
    """Run UMAP on reference PCA and project new samples.

    This is a convenience function that matches the interface of the
    original umap_transform_with_fitted method.

    Args:
        ref_pca: DataFrame with reference sample PCA components.
            Should have columns: FID, IID, PC1...PCn, label
        new_pca: DataFrame with new sample PCA components.
            Should have columns: FID, IID, PC1...PCn
        new_labels: Series with predicted labels for new samples.
        params: Optional dictionary with UMAP parameters from grid search.
            If provided, overrides config.
        config: Optional UMAP configuration. Default is UMAPConfig().

    Returns:
        Dictionary with:
            - reducer: Fitted UMAPReducer
            - ref_umap: DataFrame with reference UMAP + labels + dataset='ref'
            - new_umap: DataFrame with new samples UMAP + labels + dataset='predicted'
            - total_umap: DataFrame with combined ref + new UMAP
    """
    # Get reference data without IDs and labels
    ref_labels = ref_pca["label"].values
    ref_ids = ref_pca[["FID", "IID"]]
    X_ref = ref_pca.drop(columns=["FID", "IID", "label"]).values

    # Get new sample data without IDs
    new_ids = new_pca[["FID", "IID"]]
    X_new = new_pca.drop(columns=["FID", "IID", "label"], errors="ignore").values

    # Create reducer
    if params is not None:
        reducer = UMAPReducer.from_params(params)
    elif config is not None:
        reducer = UMAPReducer(config=config)
    else:
        reducer = UMAPReducer()

    # Fit on reference data
    reducer.fit(X_ref)

    # Transform both datasets
    ref_umap_arr = reducer.transform(X_ref, return_dataframe=True)
    new_umap_arr = reducer.transform(X_new, return_dataframe=True)

    assert isinstance(ref_umap_arr, pd.DataFrame)
    assert isinstance(new_umap_arr, pd.DataFrame)

    # Build output DataFrames
    ref_umap = ref_umap_arr.reset_index(drop=True)
    ref_umap["label"] = ref_labels
    ref_umap["dataset"] = "ref"

    new_umap = new_umap_arr.reset_index(drop=True)
    new_umap["label"] = new_labels.values if hasattr(new_labels, "values") else new_labels
    new_umap["dataset"] = "predicted"

    # Combine
    total_umap = pd.concat([ref_umap, new_umap], ignore_index=True)

    return {
        "reducer": reducer,
        "ref_umap": ref_umap,
        "new_umap": new_umap,
        "total_umap": total_umap,
    }
