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

"""Tests for ancestry reducers (PCA, UMAP)."""

from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pytest

from genotools.ancestry.config import PCAConfig, UMAPConfig
from genotools.ancestry.reducers.pca import PCAReducer, flashpca_scale
from genotools.ancestry.reducers.umap_reducer import UMAPReducer
from genotools.core.exceptions import AncestryError


class TestFlashPCAScale:
    """Tests for flashpca_scale function."""

    def test_compute_stats(self) -> None:
        """flashpca_scale computes mean and sd correctly."""
        # Create genotype-like data (0, 1, 2 values)
        np.random.seed(42)
        data = np.random.choice([0, 1, 2], size=(100, 50), p=[0.25, 0.5, 0.25])
        data = data.astype(float)

        scaled, mean, sd, keep_mask = flashpca_scale(data, compute_stats=True)

        # Mean should be computed
        assert mean is not None
        assert len(mean) == 50

        # SD should use flashPCA formula: sqrt(mean/2 * (1 - mean/2))
        expected_sd = np.sqrt((mean / 2) * (1 - mean / 2))
        np.testing.assert_array_almost_equal(sd, expected_sd)

        # Most SNPs should be kept
        assert keep_mask.sum() > 40

    def test_apply_existing_stats(self) -> None:
        """flashpca_scale applies provided mean and sd."""
        np.random.seed(42)
        data = np.random.choice([0, 1, 2], size=(50, 20), p=[0.25, 0.5, 0.25])
        data = data.astype(float)

        # First compute stats
        _, mean, sd, keep_mask = flashpca_scale(data, compute_stats=True)

        # Apply to new data
        new_data = np.random.choice([0, 1, 2], size=(30, 20), p=[0.25, 0.5, 0.25])
        new_data = new_data.astype(float)

        scaled, _, _, _ = flashpca_scale(
            new_data, mean=mean, sd=sd, compute_stats=False
        )

        # Output should be filtered to kept SNPs
        assert scaled.shape[1] == keep_mask.sum()

    def test_missing_stats_raises(self) -> None:
        """flashpca_scale raises without stats when not computing."""
        data = np.random.randn(10, 5)

        with pytest.raises(ValueError, match="mean and sd must be provided"):
            flashpca_scale(data, compute_stats=False)


class TestPCAReducer:
    """Tests for PCAReducer class."""

    @pytest.fixture
    def sample_data(self) -> np.ndarray:  # type: ignore[type-arg]
        """Create sample genotype data."""
        np.random.seed(42)
        # 100 samples, 200 SNPs
        data = np.random.choice([0, 1, 2], size=(100, 200), p=[0.25, 0.5, 0.25])
        return data.astype(float)

    def test_default_config(self) -> None:
        """PCAReducer uses default config."""
        reducer = PCAReducer()
        assert reducer.n_components == 50
        assert not reducer.is_fitted

    def test_custom_config(self) -> None:
        """PCAReducer accepts custom config."""
        config = PCAConfig(n_components=20)
        reducer = PCAReducer(config=config)
        assert reducer.n_components == 20

    def test_fit(self, sample_data: np.ndarray) -> None:  # type: ignore[type-arg]
        """PCAReducer fits on data."""
        reducer = PCAReducer(config=PCAConfig(n_components=10))
        reducer.fit(sample_data)

        assert reducer.is_fitted
        assert reducer.n_variants_used > 0
        assert reducer.pca is not None

    def test_transform_before_fit_raises(
        self, sample_data: np.ndarray  # type: ignore[type-arg]
    ) -> None:
        """Transform before fit raises error."""
        reducer = PCAReducer()

        with pytest.raises(AncestryError, match="must be fitted"):
            reducer.transform(sample_data)

    def test_transform(self, sample_data: np.ndarray) -> None:  # type: ignore[type-arg]
        """PCAReducer transforms data after fit."""
        config = PCAConfig(n_components=10)
        reducer = PCAReducer(config=config)
        reducer.fit(sample_data)

        transformed = reducer.transform(sample_data)
        assert isinstance(transformed, np.ndarray)
        assert transformed.shape == (100, 10)

    def test_transform_returns_dataframe(
        self, sample_data: np.ndarray  # type: ignore[type-arg]
    ) -> None:
        """Transform can return DataFrame."""
        config = PCAConfig(n_components=10)
        reducer = PCAReducer(config=config)
        reducer.fit(sample_data)

        result = reducer.transform(sample_data, return_dataframe=True)
        assert isinstance(result, pd.DataFrame)
        assert list(result.columns) == [f"PC{i+1}" for i in range(10)]

    def test_fit_transform(self, sample_data: np.ndarray) -> None:  # type: ignore[type-arg]
        """fit_transform works in one step."""
        config = PCAConfig(n_components=10)
        reducer = PCAReducer(config=config)

        result = reducer.fit_transform(sample_data)
        assert reducer.is_fitted
        assert result.shape == (100, 10)

    def test_eigenvalues(self, sample_data: np.ndarray) -> None:  # type: ignore[type-arg]
        """Eigenvalues are available after fit."""
        config = PCAConfig(n_components=10)
        reducer = PCAReducer(config=config)
        reducer.fit(sample_data)

        assert reducer.eigenvalues is not None
        assert len(reducer.eigenvalues) == 10
        # Eigenvalues should be in descending order
        assert all(
            reducer.eigenvalues[i] >= reducer.eigenvalues[i + 1]
            for i in range(len(reducer.eigenvalues) - 1)
        )

    def test_explained_variance_ratio(
        self, sample_data: np.ndarray  # type: ignore[type-arg]
    ) -> None:
        """Explained variance ratio is available after fit."""
        config = PCAConfig(n_components=10)
        reducer = PCAReducer(config=config)
        reducer.fit(sample_data)

        assert reducer.explained_variance_ratio is not None
        assert len(reducer.explained_variance_ratio) == 10
        # Sum should be <= 1
        assert reducer.explained_variance_ratio.sum() <= 1.0

    def test_save_eigenvalues(
        self, sample_data: np.ndarray, tmp_path: Path  # type: ignore[type-arg]
    ) -> None:
        """Eigenvalues can be saved to file."""
        config = PCAConfig(n_components=10)
        reducer = PCAReducer(config=config)
        reducer.fit(sample_data)

        output = tmp_path / "eigenvalues.txt"
        reducer.save_eigenvalues(output)

        assert output.exists()
        df = pd.read_csv(output, sep="\t")
        assert len(df) == 10
        assert "PC" in df.columns
        assert "eigenvalue" in df.columns

    def test_insufficient_samples_raises(self) -> None:
        """Fit with too few samples raises error."""
        # Only 20 samples
        data = np.random.choice([0, 1, 2], size=(20, 100))
        data = data.astype(float)

        reducer = PCAReducer()
        with pytest.raises(AncestryError, match="insufficient for PCA"):
            reducer.fit(data)

    def test_insufficient_snps_raises(self) -> None:
        """Fit with too few SNPs raises error."""
        # Only 20 SNPs
        data = np.random.choice([0, 1, 2], size=(100, 20))
        data = data.astype(float)

        reducer = PCAReducer()
        with pytest.raises(AncestryError, match="insufficient for PCA"):
            reducer.fit(data)

    def test_handles_missing_values(
        self, sample_data: np.ndarray  # type: ignore[type-arg]
    ) -> None:
        """PCAReducer handles NaN values via imputation."""
        # Add some missing values
        data_with_nan = sample_data.copy()
        data_with_nan[0, 0] = np.nan
        data_with_nan[10, 5] = np.nan

        config = PCAConfig(n_components=10)
        reducer = PCAReducer(config=config)
        reducer.fit(data_with_nan)

        # Should fit without error
        assert reducer.is_fitted


class TestUMAPReducer:
    """Tests for UMAPReducer class."""

    @pytest.fixture
    def sample_pca_data(self) -> np.ndarray:  # type: ignore[type-arg]
        """Create sample PCA-like data."""
        np.random.seed(42)
        # 100 samples, 20 PCs
        return np.random.randn(100, 20)

    def test_default_config(self) -> None:
        """UMAPReducer uses default config."""
        reducer = UMAPReducer()
        assert reducer.n_components == 2
        assert not reducer.is_fitted

    def test_custom_config(self) -> None:
        """UMAPReducer accepts custom config."""
        config = UMAPConfig(n_neighbors=10, n_components=3)
        reducer = UMAPReducer(config=config)
        assert reducer.n_components == 3

    def test_fit(self, sample_pca_data: np.ndarray) -> None:  # type: ignore[type-arg]
        """UMAPReducer fits on data."""
        reducer = UMAPReducer()
        reducer.fit(sample_pca_data)

        assert reducer.is_fitted
        assert reducer.umap is not None

    def test_transform(self, sample_pca_data: np.ndarray) -> None:  # type: ignore[type-arg]
        """UMAPReducer transforms data."""
        reducer = UMAPReducer()
        reducer.fit(sample_pca_data)

        transformed = reducer.transform(sample_pca_data)
        assert isinstance(transformed, np.ndarray)
        assert transformed.shape == (100, 2)

    def test_transform_returns_dataframe(
        self, sample_pca_data: np.ndarray  # type: ignore[type-arg]
    ) -> None:
        """Transform can return DataFrame."""
        reducer = UMAPReducer()
        reducer.fit(sample_pca_data)

        result = reducer.transform(sample_pca_data, return_dataframe=True)
        assert isinstance(result, pd.DataFrame)
        assert list(result.columns) == ["UMAP1", "UMAP2"]

    def test_fit_transform(self, sample_pca_data: np.ndarray) -> None:  # type: ignore[type-arg]
        """fit_transform works in one step."""
        reducer = UMAPReducer()
        result = reducer.fit_transform(sample_pca_data)

        assert reducer.is_fitted
        assert result.shape == (100, 2)

    def test_transform_before_fit_raises(
        self, sample_pca_data: np.ndarray  # type: ignore[type-arg]
    ) -> None:
        """Transform before fit raises error."""
        reducer = UMAPReducer()

        with pytest.raises(AncestryError, match="must be fitted"):
            reducer.transform(sample_pca_data)

    def test_from_params(self) -> None:
        """UMAPReducer can be created from parameter dict."""
        params = {
            "umap__n_neighbors": 20,
            "umap__n_components": 3,
            "umap__a": 1.5,
            "umap__b": 0.5,
        }
        reducer = UMAPReducer.from_params(params)

        assert reducer.config.n_neighbors == 20
        assert reducer.config.n_components == 3
        assert reducer.config.a == 1.5
        assert reducer.config.b == 0.5

    def test_from_params_plain_keys(self) -> None:
        """from_params works with plain key names."""
        params = {
            "n_neighbors": 20,
            "n_components": 3,
        }
        reducer = UMAPReducer.from_params(params)

        assert reducer.config.n_neighbors == 20
        assert reducer.config.n_components == 3

    def test_insufficient_samples_raises(self) -> None:
        """Fit with too few samples raises error."""
        # Only 5 samples, but n_neighbors defaults to 15
        data = np.random.randn(5, 20)

        reducer = UMAPReducer()
        with pytest.raises(AncestryError, match="requires at least"):
            reducer.fit(data)


class _RecordingReducer:
    """Stands in for UMAPReducer, recording the feature width it is handed."""

    widths: list = []

    def fit(self, X: np.ndarray) -> "_RecordingReducer":
        self.widths.append(("fit", X.shape[1]))
        return self

    def transform(self, X: np.ndarray, return_dataframe: bool = False):
        self.widths.append(("transform", X.shape[1]))
        frame = pd.DataFrame(np.asarray(X)[:, :2], columns=["UMAP1", "UMAP2"])
        return frame if return_dataframe else frame.values

    @classmethod
    def from_params(cls, params: Any) -> "_RecordingReducer":
        return cls()


def _pca_frame(n: int = 6, label: str = "EUR") -> pd.DataFrame:
    frame = pd.DataFrame({
        "FID": [f"S{i}" for i in range(n)],
        "IID": [f"S{i}" for i in range(n)],
        "PC1": np.linspace(0, 1, n),
        "PC2": np.linspace(1, 0, n),
        "PC3": np.linspace(0, 2, n),
        "label": [label] * n,
    })
    return frame


def test_run_umap_embeds_the_pcs_and_nothing_else(monkeypatch):
    """An extra column in the projected frame must not enter the embedding.

    `run_umap` used to select features by dropping FID/IID/label, so anything
    else written into `_projected_new_pca.txt` became a feature and moved every
    point -- silently, since UMAP accepts any width. Diagnostics write per-sample
    distances now, which is exactly the column that would have leaked in.
    """
    from genotools.ancestry.reducers import umap_reducer

    _RecordingReducer.widths = []
    monkeypatch.setattr(umap_reducer, "UMAPReducer", _RecordingReducer)

    new_pca = _pca_frame()
    new_pca["dist_to_all"] = 42.0
    new_pca["label_pre_admixture"] = "AFR"

    result = umap_reducer.run_umap(
        ref_pca=_pca_frame(),
        new_pca=new_pca,
        new_labels=new_pca["label"],
    )

    assert {width for _, width in _RecordingReducer.widths} == {3}
    assert list(result["new_umap"].columns) == ["UMAP1", "UMAP2", "label", "dataset"]


def test_run_umap_falls_back_when_components_are_not_named_pc(monkeypatch):
    """A frame naming its components something else still embeds them all."""
    from genotools.ancestry.reducers import umap_reducer

    _RecordingReducer.widths = []
    monkeypatch.setattr(umap_reducer, "UMAPReducer", _RecordingReducer)

    frame = pd.DataFrame({
        "FID": ["a", "b", "c"], "IID": ["a", "b", "c"],
        "C1": [0.0, 1.0, 2.0], "C2": [2.0, 1.0, 0.0], "label": ["EUR"] * 3,
    })
    umap_reducer.run_umap(ref_pca=frame, new_pca=frame.copy(), new_labels=frame["label"])
    assert {width for _, width in _RecordingReducer.widths} == {2}
