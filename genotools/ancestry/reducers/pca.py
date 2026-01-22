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

"""PCA dimensionality reduction for ancestry prediction.

This module provides PCA computation and transformation for the ancestry
prediction pipeline. It uses flashPCA-style scaling for compatibility with
the existing implementation and reference panels.

The flashPCA scaling computes standard deviation as sqrt(p/2 * (1 - p/2))
where p is the mean allele frequency. This differs from standard scaling
and is important for reproducibility with existing trained models.

Reference: https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0093766

Example:
    >>> from genotools.ancestry.reducers.pca import PCAReducer
    >>> pca = PCAReducer(n_components=50)
    >>> pca.fit(X_train)
    >>> X_train_pca = pca.transform(X_train)
    >>> X_new_pca = pca.transform(X_new)
"""

from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, Tuple, Any

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.impute import SimpleImputer

from genotools.ancestry.config import PCAConfig
from genotools.core.exceptions import AncestryError


def flashpca_scale(
    data: np.ndarray,  # type: ignore[type-arg]
    mean: Optional[np.ndarray] = None,  # type: ignore[type-arg]
    sd: Optional[np.ndarray] = None,  # type: ignore[type-arg]
    compute_stats: bool = False,
    eps: float = 1e-12,
) -> Tuple[
    np.ndarray,  # type: ignore[type-arg]
    np.ndarray,  # type: ignore[type-arg]
    np.ndarray,  # type: ignore[type-arg]
    np.ndarray,  # type: ignore[type-arg]
]:
    """Apply flashPCA-style scaling to genotype data.

    FlashPCA uses a specific scaling formula where the standard deviation
    is computed as sqrt(mean/2 * (1 - mean/2)), which is appropriate for
    genotype data where values represent allele counts (0, 1, 2).

    Reference: Abraham & Inouye, 2014, PLOS ONE
    https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0093766

    Args:
        data: Genotype matrix (samples x variants) with values 0, 1, 2.
        mean: Pre-computed means. If None and compute_stats=True, computed.
        sd: Pre-computed standard deviations. If None and compute_stats=True,
            computed using flashPCA formula.
        compute_stats: If True, compute mean and sd from data.
        eps: Small value to avoid division by zero. Default is 1e-12.

    Returns:
        Tuple of:
            - scaled_data: Scaled genotype matrix
            - mean: Mean values per variant
            - sd: Standard deviation values per variant
            - keep_mask: Boolean mask of variants to keep (sd > eps)

    Raises:
        ValueError: If mean/sd not provided and compute_stats is False.
    """
    if compute_stats:
        mean = data.mean(axis=0)
        # flashPCA-style standard deviation
        sd = np.sqrt((mean / 2) * (1 - mean / 2))
    elif mean is None or sd is None:
        raise ValueError(
            "mean and sd must be provided if compute_stats is False"
        )

    # Create mask for non-zero SD variants
    keep_mask = sd > eps

    # Filter data and stats to kept variants
    data_filtered = data[:, keep_mask]
    mean_filtered = mean[keep_mask]
    sd_filtered = sd[keep_mask]

    # Scale data
    scaled_data = (data_filtered - mean_filtered) / sd_filtered

    return scaled_data, mean, sd, keep_mask


@dataclass
class PCAReducer:
    """PCA dimensionality reducer for ancestry prediction.

    This class wraps sklearn's PCA with flashPCA-style scaling and
    provides a fit/transform interface for use in the ancestry pipeline.

    Attributes:
        config: PCA configuration parameters.
        pca: Fitted sklearn PCA object.
        mean: Training data mean for scaling.
        sd: Training data standard deviation for scaling.
        keep_mask: Boolean mask of variants kept after SD filtering.
        imputer: SimpleImputer for handling missing values.
    """

    config: PCAConfig = field(default_factory=PCAConfig)
    pca: Optional[PCA] = field(default=None, repr=False)
    mean: Optional[np.ndarray] = field(default=None, repr=False)  # type: ignore[type-arg]
    sd: Optional[np.ndarray] = field(default=None, repr=False)  # type: ignore[type-arg]
    keep_mask: Optional[np.ndarray] = field(default=None, repr=False)  # type: ignore[type-arg]
    imputer: Optional[SimpleImputer] = field(default=None, repr=False)
    _is_fitted: bool = field(default=False, repr=False)

    @property
    def n_components(self) -> int:
        """Number of PCA components."""
        return self.config.n_components

    @property
    def is_fitted(self) -> bool:
        """Whether the reducer has been fitted."""
        return self._is_fitted

    @property
    def n_variants_used(self) -> int:
        """Number of variants used after SD filtering."""
        if self.keep_mask is None:
            return 0
        return int(self.keep_mask.sum())

    @property
    def n_variants_dropped(self) -> int:
        """Number of variants dropped due to zero SD."""
        if self.keep_mask is None:
            return 0
        return int((~self.keep_mask).sum())

    @property
    def eigenvalues(self) -> Optional[np.ndarray]:  # type: ignore[type-arg]
        """Eigenvalues (explained variance) for each component."""
        if self.pca is None:
            return None
        return self.pca.explained_variance_

    @property
    def explained_variance_ratio(self) -> Optional[np.ndarray]:  # type: ignore[type-arg]
        """Proportion of variance explained by each component."""
        if self.pca is None:
            return None
        return self.pca.explained_variance_ratio_

    def _get_pc_names(self) -> list[str]:
        """Get column names for PC output."""
        return [f"PC{i+1}" for i in range(self.config.n_components)]

    def fit(self, X: np.ndarray) -> "PCAReducer":  # type: ignore[type-arg]
        """Fit PCA on training data.

        Args:
            X: Genotype matrix (samples x variants) with values 0, 1, 2.
                Missing values (NaN) are imputed with mean.

        Returns:
            Self for method chaining.

        Raises:
            AncestryError: If insufficient samples or variants for PCA.
        """
        # Validate input dimensions
        if X.shape[0] < 50:
            raise AncestryError(
                f"Training data has only {X.shape[0]} samples, "
                f"which is insufficient for PCA. Need at least 50."
            )
        if X.shape[1] < 50:
            raise AncestryError(
                f"Training data has only {X.shape[1]} SNPs, "
                f"which is insufficient for PCA. Need at least 50."
            )

        # Impute missing values
        self.imputer = SimpleImputer(missing_values=np.nan, strategy="mean")
        X_imputed = self.imputer.fit_transform(X)

        # Apply flashPCA scaling
        X_scaled, self.mean, self.sd, self.keep_mask = flashpca_scale(
            X_imputed, compute_stats=True
        )

        # Check we still have enough variants after filtering
        if X_scaled.shape[1] < 50:
            raise AncestryError(
                f"After removing {self.n_variants_dropped} zero-SD SNPs, "
                f"only {self.n_variants_used} SNPs remain. "
                f"Insufficient for PCA. Consider relaxing filters."
            )

        # Fit PCA
        self.pca = PCA(
            n_components=self.config.n_components,
            svd_solver="full",
            random_state=self.config.random_state,
        )
        self.pca.fit(X_scaled)

        self._is_fitted = True
        return self

    def transform(
        self,
        X: np.ndarray,  # type: ignore[type-arg]
        return_dataframe: bool = False,
    ) -> np.ndarray | pd.DataFrame:  # type: ignore[type-arg]
        """Transform data using fitted PCA.

        Args:
            X: Genotype matrix (samples x variants) to transform.
            return_dataframe: If True, return DataFrame with PC column names.

        Returns:
            PCA-transformed data (samples x n_components).

        Raises:
            AncestryError: If PCA not fitted.
        """
        if not self._is_fitted:
            raise AncestryError("PCAReducer must be fitted before transform")

        if self.imputer is None or self.pca is None:
            raise AncestryError("PCAReducer is in invalid state")

        if self.mean is None or self.sd is None or self.keep_mask is None:
            raise AncestryError("PCAReducer scaling parameters not set")

        # Impute missing values
        X_imputed = self.imputer.transform(X)

        # Apply flashPCA scaling using fitted parameters
        X_scaled, _, _, _ = flashpca_scale(
            X_imputed,
            mean=self.mean,
            sd=self.sd,
            compute_stats=False,
        )

        # Apply keep_mask (filter to same variants as training)
        X_filtered = X_scaled

        # Transform with PCA
        X_pca = self.pca.transform(X_filtered)

        if return_dataframe:
            return pd.DataFrame(X_pca, columns=self._get_pc_names())
        return X_pca

    def fit_transform(
        self,
        X: np.ndarray,  # type: ignore[type-arg]
        return_dataframe: bool = False,
    ) -> np.ndarray | pd.DataFrame:  # type: ignore[type-arg]
        """Fit PCA and transform data in one step.

        Args:
            X: Genotype matrix (samples x variants).
            return_dataframe: If True, return DataFrame with PC column names.

        Returns:
            PCA-transformed data (samples x n_components).
        """
        self.fit(X)
        return self.transform(X, return_dataframe=return_dataframe)

    def save_eigenvalues(self, path: Path) -> Path:
        """Save eigenvalues and explained variance to file.

        Args:
            path: Output file path.

        Returns:
            Path to saved file.

        Raises:
            AncestryError: If PCA not fitted.
        """
        if not self._is_fitted or self.pca is None:
            raise AncestryError("PCAReducer must be fitted before saving")

        ev_df = pd.DataFrame({
            "PC": self._get_pc_names(),
            "eigenvalue": self.pca.explained_variance_,
            "explained_variance_ratio": self.pca.explained_variance_ratio_,
        })
        ev_df.to_csv(path, sep="\t", index=False)
        return path


def run_pca(
    X_train: np.ndarray,  # type: ignore[type-arg]
    X_test: np.ndarray,  # type: ignore[type-arg]
    X_new: np.ndarray,  # type: ignore[type-arg]
    train_ids: pd.DataFrame,
    test_ids: pd.DataFrame,
    new_ids: pd.DataFrame,
    train_labels: np.ndarray,  # type: ignore[type-arg]
    test_labels: np.ndarray,  # type: ignore[type-arg]
    config: Optional[PCAConfig] = None,
    out_path: Optional[Path] = None,
) -> dict[str, Any]:
    """Run PCA on training data and project test/new samples.

    This is a convenience function that wraps PCAReducer to match the
    interface of the original calculate_pcs method.

    Args:
        X_train: Training genotype matrix (samples x variants).
        X_test: Test genotype matrix.
        X_new: New samples genotype matrix.
        train_ids: DataFrame with FID, IID for training samples.
        test_ids: DataFrame with FID, IID for test samples.
        new_ids: DataFrame with FID, IID for new samples.
        train_labels: Ancestry labels for training samples.
        test_labels: Ancestry labels for test samples.
        config: PCA configuration. Default is PCAConfig().
        out_path: Optional path prefix for saving outputs.

    Returns:
        Dictionary with:
            - reducer: Fitted PCAReducer
            - X_train: Transformed training data
            - X_test: Transformed test data
            - train_pca: DataFrame with train IDs, PCs, and labels
            - test_pca: DataFrame with test IDs, PCs, and labels
            - ref_pca: DataFrame with combined train+test (full reference)
            - projected_pca: DataFrame with new samples IDs and PCs
            - eigenvalues: Eigenvalues array
            - explained_variance_ratio: Explained variance ratio array
    """
    if config is None:
        config = PCAConfig()

    # Create and fit reducer
    reducer = PCAReducer(config=config)
    reducer.fit(X_train)

    # Transform all datasets
    X_train_pca = reducer.transform(X_train, return_dataframe=True)
    X_test_pca = reducer.transform(X_test, return_dataframe=True)
    X_new_pca = reducer.transform(X_new, return_dataframe=True)

    # Build output DataFrames with IDs and labels
    train_ids_reset = train_ids.reset_index(drop=True)
    test_ids_reset = test_ids.reset_index(drop=True)
    new_ids_reset = new_ids.reset_index(drop=True)

    # Ensure PCA DataFrames are properly typed
    assert isinstance(X_train_pca, pd.DataFrame)
    assert isinstance(X_test_pca, pd.DataFrame)
    assert isinstance(X_new_pca, pd.DataFrame)

    train_pca = pd.concat([train_ids_reset, X_train_pca], axis=1)
    train_pca["label"] = train_labels

    test_pca = pd.concat([test_ids_reset, X_test_pca], axis=1)
    test_pca["label"] = test_labels

    ref_pca = pd.concat([train_pca, test_pca], ignore_index=True)

    projected_pca = pd.concat([new_ids_reset, X_new_pca], axis=1)
    projected_pca["label"] = "new"

    # Save outputs if path provided
    if out_path:
        train_pca.to_csv(
            f"{out_path}_labeled_train_pca.txt", sep="\t", index=False
        )
        ref_pca.to_csv(f"{out_path}_labeled_ref_pca.txt", sep="\t", index=False)
        projected_pca.to_csv(
            f"{out_path}_projected_new_pca.txt", sep="\t", index=False
        )
        reducer.save_eigenvalues(Path(f"{out_path}_pca_eigenvalues.txt"))

    return {
        "reducer": reducer,
        "X_train": X_train_pca.values,
        "X_test": X_test_pca.values,
        "train_pca": train_pca,
        "test_pca": test_pca,
        "ref_pca": ref_pca,
        "projected_pca": projected_pca,
        "eigenvalues": reducer.eigenvalues,
        "explained_variance_ratio": reducer.explained_variance_ratio,
    }
