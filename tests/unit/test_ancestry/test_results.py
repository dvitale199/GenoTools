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

"""Tests for ancestry result classes."""

from pathlib import Path
from typing import Any, Dict

import numpy as np
import pandas as pd
import pytest

from genotools.ancestry.results import (
    AncestryPredictions,
    AncestryResult,
    PCAResult,
    SplitResult,
    TrainingMetrics,
    UMAPResult,
)


class TestAncestryPredictions:
    """Tests for AncestryPredictions class."""

    @pytest.fixture
    def sample_predictions(self) -> pd.DataFrame:
        """Create sample predictions DataFrame."""
        return pd.DataFrame({
            "FID": ["FAM1", "FAM1", "FAM2", "FAM2", "FAM3"],
            "IID": ["SAMP1", "SAMP2", "SAMP3", "SAMP4", "SAMP5"],
            "predicted_ancestry": ["EUR", "EUR", "AFR", "AFR", "EAS"],
        })

    def test_create_predictions(self, sample_predictions: pd.DataFrame) -> None:
        """Predictions can be created from DataFrame."""
        preds = AncestryPredictions(predictions=sample_predictions)
        assert preds.sample_count == 5
        assert len(preds.unique_ancestries) == 3

    def test_missing_columns_raises(self) -> None:
        """Missing required columns raises ValueError."""
        df = pd.DataFrame({
            "FID": ["FAM1"],
            "IID": ["SAMP1"],
            # Missing predicted_ancestry
        })
        with pytest.raises(ValueError, match="missing required columns"):
            AncestryPredictions(predictions=df)

    def test_get_ancestry(self, sample_predictions: pd.DataFrame) -> None:
        """get_ancestry returns correct ancestry for sample."""
        preds = AncestryPredictions(predictions=sample_predictions)
        assert preds.get_ancestry("SAMP1") == "EUR"
        assert preds.get_ancestry("SAMP3") == "AFR"
        assert preds.get_ancestry("SAMP5") == "EAS"

    def test_get_ancestry_not_found(
        self, sample_predictions: pd.DataFrame
    ) -> None:
        """get_ancestry returns None for unknown sample."""
        preds = AncestryPredictions(predictions=sample_predictions)
        assert preds.get_ancestry("UNKNOWN") is None

    def test_filter_by_ancestry(self, sample_predictions: pd.DataFrame) -> None:
        """filter_by_ancestry returns correct sample IDs."""
        preds = AncestryPredictions(predictions=sample_predictions)
        eur_samples = preds.filter_by_ancestry("EUR")
        assert set(eur_samples) == {"SAMP1", "SAMP2"}

        afr_samples = preds.filter_by_ancestry("AFR")
        assert set(afr_samples) == {"SAMP3", "SAMP4"}

    def test_filter_by_ancestry_empty(
        self, sample_predictions: pd.DataFrame
    ) -> None:
        """filter_by_ancestry returns empty list for missing ancestry."""
        preds = AncestryPredictions(predictions=sample_predictions)
        assert preds.filter_by_ancestry("SAS") == []

    def test_ancestry_counts(self, sample_predictions: pd.DataFrame) -> None:
        """ancestry_counts returns correct counts."""
        preds = AncestryPredictions(predictions=sample_predictions)
        counts = preds.ancestry_counts
        assert counts["EUR"] == 2
        assert counts["AFR"] == 2
        assert counts["EAS"] == 1

    def test_to_dataframe(self, sample_predictions: pd.DataFrame) -> None:
        """to_dataframe returns copy of predictions."""
        preds = AncestryPredictions(predictions=sample_predictions)
        df = preds.to_dataframe()

        # Should be a copy, not the same object
        assert df is not preds.predictions
        assert len(df) == 5

    def test_to_dict_format(self, sample_predictions: pd.DataFrame) -> None:
        """to_dict returns legacy format."""
        preds = AncestryPredictions(predictions=sample_predictions)
        result = preds.to_dict()

        assert "data" in result
        assert "metrics" in result
        assert "output" in result
        assert "ids" in result["data"]
        assert "predicted_counts" in result["metrics"]

    def test_save(
        self, sample_predictions: pd.DataFrame, tmp_path: Path
    ) -> None:
        """save writes predictions to file."""
        preds = AncestryPredictions(predictions=sample_predictions)
        output_path = tmp_path / "predictions.txt"
        preds.save(output_path)

        assert output_path.exists()
        # Read back and verify
        df = pd.read_csv(output_path, sep="\t")
        assert "FID" in df.columns
        assert "IID" in df.columns
        assert "label" in df.columns  # renamed from predicted_ancestry
        assert len(df) == 5


class TestTrainingMetrics:
    """Tests for TrainingMetrics class."""

    def test_create_metrics(self) -> None:
        """TrainingMetrics can be created."""
        metrics = TrainingMetrics(
            train_accuracy=0.95,
            test_accuracy=0.90,
            train_accuracy_ci=(0.93, 0.97),
            test_accuracy_ci=(0.85, 0.95),
            confusion_matrix=np.array([[10, 2], [1, 15]]),
            best_params={"umap__n_neighbors": 15},
            label_encoder_classes=["EUR", "AFR"],
        )
        assert metrics.train_accuracy == 0.95
        assert metrics.test_accuracy == 0.90

    def test_to_dict(self) -> None:
        """to_dict returns all metrics."""
        metrics = TrainingMetrics(
            train_accuracy=0.95,
            test_accuracy=0.90,
            train_accuracy_ci=(0.93, 0.97),
            test_accuracy_ci=(0.85, 0.95),
            confusion_matrix=np.array([[10, 2], [1, 15]]),
            best_params={"umap__n_neighbors": 15},
            label_encoder_classes=["EUR", "AFR"],
        )
        result = metrics.to_dict()

        assert result["train_accuracy"] == 0.95
        assert result["test_accuracy"] == 0.90
        assert "confusion_matrix" in result
        assert result["best_params"] == {"umap__n_neighbors": 15}


class TestPCAResult:
    """Tests for PCAResult class."""

    @pytest.fixture
    def sample_pca_result(self) -> PCAResult:
        """Create sample PCA result."""
        n_samples = 10
        n_components = 5

        return PCAResult(
            train_pca=pd.DataFrame({
                "FID": [f"F{i}" for i in range(n_samples)],
                "IID": [f"S{i}" for i in range(n_samples)],
                **{f"PC{i+1}": np.random.randn(n_samples) for i in range(n_components)},
                "label": ["EUR"] * (n_samples // 2) + ["AFR"] * (n_samples // 2),
            }),
            test_pca=pd.DataFrame({
                "FID": [f"F{i}" for i in range(n_samples)],
                "IID": [f"S{i}" for i in range(n_samples)],
                **{f"PC{i+1}": np.random.randn(n_samples) for i in range(n_components)},
                "label": ["EUR"] * (n_samples // 2) + ["AFR"] * (n_samples // 2),
            }),
            ref_pca=pd.DataFrame(),  # Combined train + test
            projected_pca=pd.DataFrame(),  # New samples
            eigenvalues=np.array([10.0, 5.0, 2.0, 1.0, 0.5]),
            explained_variance_ratio=np.array([0.5, 0.25, 0.1, 0.05, 0.025]),
            n_components=n_components,
            n_snps_used=10000,
            train_mean=np.zeros(10000),
            train_sd=np.ones(10000),
        )

    def test_save_eigenvalues(
        self, sample_pca_result: PCAResult, tmp_path: Path
    ) -> None:
        """save_eigenvalues writes to file."""
        output_path = tmp_path / "eigenvalues.txt"
        sample_pca_result.save_eigenvalues(output_path)

        assert output_path.exists()
        df = pd.read_csv(output_path, sep="\t")
        assert "PC" in df.columns
        assert "eigenvalue" in df.columns
        assert "explained_variance_ratio" in df.columns
        assert len(df) == 5


class TestUMAPResult:
    """Tests for UMAPResult class."""

    def test_create_result(self) -> None:
        """UMAPResult can be created."""
        result = UMAPResult(
            ref_umap=pd.DataFrame({
                "UMAP1": [0.1, 0.2],
                "UMAP2": [0.3, 0.4],
                "label": ["EUR", "AFR"],
            }),
            new_umap=pd.DataFrame({
                "UMAP1": [0.15],
                "UMAP2": [0.35],
                "label": ["EUR"],
            }),
            total_umap=pd.DataFrame(),
            params={"umap__n_neighbors": 15},
        )
        assert len(result.ref_umap) == 2
        assert len(result.new_umap) == 1


class TestSplitResult:
    """Tests for SplitResult class."""

    def test_create_result(self) -> None:
        """SplitResult can be created."""
        result = SplitResult(
            ancestry_paths={
                "EUR": Path("/path/to/eur"),
                "AFR": Path("/path/to/afr"),
            },
            ancestry_labels=["EUR", "AFR"],
            pruned_samples=pd.DataFrame({
                "FID": ["F1"],
                "IID": ["S1"],
                "step": ["insufficient_ancestry_sample_n"],
                "label": ["SAS"],
            }),
        )
        assert result.n_ancestries == 2
        assert result.n_pruned == 1


class TestAncestryResult:
    """Tests for AncestryResult class."""

    @pytest.fixture
    def sample_predictions(self) -> AncestryPredictions:
        """Create sample predictions."""
        return AncestryPredictions(
            predictions=pd.DataFrame({
                "FID": ["FAM1", "FAM2"],
                "IID": ["SAMP1", "SAMP2"],
                "predicted_ancestry": ["EUR", "AFR"],
            })
        )

    def test_create_minimal_result(
        self, sample_predictions: AncestryPredictions
    ) -> None:
        """AncestryResult can be created with just predictions."""
        result = AncestryResult(predictions=sample_predictions)
        assert result.predictions.sample_count == 2
        assert result.training_metrics is None
        assert result.pca_result is None

    def test_to_dict_format(
        self, sample_predictions: AncestryPredictions
    ) -> None:
        """to_dict returns legacy format."""
        result = AncestryResult(predictions=sample_predictions)
        legacy = result.to_dict()

        assert legacy["step"] == "predict_ancestry"
        assert "data" in legacy
        assert "metrics" in legacy
        assert "output" in legacy
