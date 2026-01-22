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

"""Ancestry prediction result dataclasses.

This module provides typed result classes for ancestry prediction outputs,
including predictions, model training metrics, and pipeline results.

Example:
    >>> from genotools.ancestry.results import AncestryPredictions
    >>> # Access predictions
    >>> predictions.get_ancestry("SAMPLE001")
    'EUR'
    >>> predictions.filter_by_ancestry("EUR")
    ['SAMPLE001', 'SAMPLE002', ...]
    >>> predictions.ancestry_counts
    {'EUR': 100, 'AFR': 50, ...}
"""

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd


@dataclass
class AncestryPredictions:
    """Ancestry prediction results for a set of samples.

    This class holds the predicted ancestry labels and associated metadata
    for all samples processed through the ancestry prediction pipeline.

    The predictions DataFrame contains:
    - FID: Family ID
    - IID: Individual ID
    - predicted_ancestry: Predicted ancestry label (e.g., "EUR", "AFR")
    - Additional columns for per-class probabilities (if available)

    Attributes:
        predictions: DataFrame with sample IDs and predicted ancestries.
            Must contain 'FID', 'IID', and 'predicted_ancestry' columns.
        probabilities: Optional DataFrame with per-class prediction
            probabilities. Columns are ancestry labels, rows are samples.
        output_path: Path where predictions were saved, if any.
    """

    predictions: pd.DataFrame
    probabilities: Optional[pd.DataFrame] = None
    output_path: Optional[Path] = None

    def __post_init__(self) -> None:
        """Validate predictions DataFrame."""
        required_cols = {"FID", "IID", "predicted_ancestry"}
        missing = required_cols - set(self.predictions.columns)
        if missing:
            raise ValueError(
                f"predictions DataFrame missing required columns: {missing}"
            )

    def get_ancestry(self, sample_id: str) -> Optional[str]:
        """Get predicted ancestry for a specific sample.

        Args:
            sample_id: Individual ID (IID) to look up.

        Returns:
            Predicted ancestry label, or None if sample not found.
        """
        mask = self.predictions["IID"] == sample_id
        if not mask.any():
            return None
        return str(self.predictions.loc[mask, "predicted_ancestry"].iloc[0])

    def filter_by_ancestry(self, ancestry: str) -> List[str]:
        """Get all sample IDs for a given ancestry.

        Args:
            ancestry: Ancestry label to filter by (e.g., "EUR").

        Returns:
            List of individual IDs (IID) with the specified ancestry.
        """
        mask = self.predictions["predicted_ancestry"] == ancestry
        return list(self.predictions.loc[mask, "IID"])

    @property
    def ancestry_counts(self) -> Dict[str, int]:
        """Count of samples per predicted ancestry.

        Returns:
            Dictionary mapping ancestry labels to sample counts.
        """
        counts = self.predictions["predicted_ancestry"].value_counts()
        return dict(counts)

    @property
    def sample_count(self) -> int:
        """Total number of samples with predictions."""
        return len(self.predictions)

    @property
    def unique_ancestries(self) -> List[str]:
        """List of unique ancestry labels in predictions."""
        return list(self.predictions["predicted_ancestry"].unique())

    def to_dataframe(self) -> pd.DataFrame:
        """Return predictions as DataFrame.

        Returns:
            Copy of the predictions DataFrame.
        """
        return self.predictions.copy()

    def to_dict(self) -> Dict[str, Any]:
        """Convert to legacy dictionary format for backward compatibility.

        Returns the format expected by existing pipeline code:
        {
            'data': {'ids': DataFrame, 'y_pred': array, ...},
            'metrics': {'predicted_counts': Series, ...},
            'output': {'labels_outpath': str}
        }

        Returns:
            Dictionary in legacy format.
        """
        return {
            "data": {
                "ids": self.predictions[["FID", "IID", "predicted_ancestry"]],
                "y_pred": self.predictions["predicted_ancestry"].values,
            },
            "metrics": {
                "predicted_counts": self.predictions[
                    "predicted_ancestry"
                ].value_counts()
            },
            "output": {
                "labels_outpath": (
                    str(self.output_path) if self.output_path else None
                )
            },
        }

    def save(self, path: Path) -> Path:
        """Save predictions to a tab-separated file.

        Args:
            path: Output file path.

        Returns:
            Path to the saved file.
        """
        self.predictions[["FID", "IID", "predicted_ancestry"]].rename(
            columns={"predicted_ancestry": "label"}
        ).to_csv(path, sep="\t", index=False)
        return path


@dataclass
class TrainingMetrics:
    """Metrics from model training.

    Captures accuracy, confusion matrix, and hyperparameters from
    the training process for model evaluation and diagnostics.

    Attributes:
        train_accuracy: Cross-validated training accuracy.
        test_accuracy: Accuracy on held-out test set.
        train_accuracy_ci: 95% confidence interval for train accuracy.
        test_accuracy_ci: 95% confidence interval for test accuracy.
        confusion_matrix: Confusion matrix as numpy array.
        best_params: Best hyperparameters from grid search.
        label_encoder_classes: Ordered list of class labels.
    """

    train_accuracy: float
    test_accuracy: float
    train_accuracy_ci: Tuple[float, float]
    test_accuracy_ci: Tuple[float, float]
    confusion_matrix: np.ndarray  # type: ignore[type-arg]
    best_params: Dict[str, Any]
    label_encoder_classes: List[str] = field(default_factory=list)

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for serialization.

        Returns:
            Dictionary with all metrics.
        """
        return {
            "train_accuracy": self.train_accuracy,
            "test_accuracy": self.test_accuracy,
            "train_accuracy_ci": self.train_accuracy_ci,
            "test_accuracy_ci": self.test_accuracy_ci,
            "confusion_matrix": self.confusion_matrix.tolist(),
            "best_params": self.best_params,
            "label_encoder_classes": self.label_encoder_classes,
        }


@dataclass
class PCAResult:
    """Results from PCA transformation.

    Contains the transformed data and metadata from PCA computation,
    including eigenvalues for variance analysis.

    Attributes:
        train_pca: PCA-transformed training data with labels.
        test_pca: PCA-transformed test data with labels.
        ref_pca: Full reference panel PCA (train + test combined).
        projected_pca: PCA-projected new samples.
        eigenvalues: Eigenvalues for each principal component.
        explained_variance_ratio: Proportion of variance explained.
        n_components: Number of PCA components.
        n_snps_used: Number of SNPs used in PCA (after filtering).
        train_mean: Training data mean (for projection).
        train_sd: Training data standard deviation (for projection).
    """

    train_pca: pd.DataFrame
    test_pca: pd.DataFrame
    ref_pca: pd.DataFrame
    projected_pca: pd.DataFrame
    eigenvalues: np.ndarray  # type: ignore[type-arg]
    explained_variance_ratio: np.ndarray  # type: ignore[type-arg]
    n_components: int
    n_snps_used: int
    train_mean: np.ndarray  # type: ignore[type-arg]
    train_sd: np.ndarray  # type: ignore[type-arg]

    def save_eigenvalues(self, path: Path) -> Path:
        """Save eigenvalues to a file.

        Args:
            path: Output file path.

        Returns:
            Path to the saved file.
        """
        col_names = [f"PC{i+1}" for i in range(self.n_components)]
        ev_df = pd.DataFrame({
            "PC": col_names,
            "eigenvalue": self.eigenvalues,
            "explained_variance_ratio": self.explained_variance_ratio,
        })
        ev_df.to_csv(path, sep="\t", index=False)
        return path


@dataclass
class UMAPResult:
    """Results from UMAP transformation.

    Contains the UMAP embeddings for visualization and the fitted
    UMAP model parameters.

    Attributes:
        ref_umap: UMAP embedding of reference samples with labels.
        new_umap: UMAP embedding of new samples with predicted labels.
        total_umap: Combined UMAP embedding (ref + new) with dataset labels.
        params: UMAP parameters used (may be from grid search).
    """

    ref_umap: pd.DataFrame
    new_umap: pd.DataFrame
    total_umap: pd.DataFrame
    params: Dict[str, Any]


@dataclass
class SplitResult:
    """Results from splitting cohort by ancestry.

    Contains the paths to ancestry-specific genotype files and
    information about samples that were pruned due to small group sizes.

    Attributes:
        ancestry_paths: Dictionary mapping ancestry labels to output paths.
        ancestry_labels: List of ancestries that met minimum sample threshold.
        pruned_samples: DataFrame of samples pruned due to small ancestry
            group size. Contains FID, IID, step, label columns.
    """

    ancestry_paths: Dict[str, Path]
    ancestry_labels: List[str]
    pruned_samples: pd.DataFrame

    @property
    def n_ancestries(self) -> int:
        """Number of ancestry groups in output."""
        return len(self.ancestry_labels)

    @property
    def n_pruned(self) -> int:
        """Number of samples pruned due to small group size."""
        return len(self.pruned_samples)


@dataclass
class AncestryResult:
    """Complete result of ancestry prediction pipeline.

    This is the top-level result class that combines all outputs from
    running ancestry prediction, including predictions, training metrics,
    PCA/UMAP results, and split cohort information.

    Attributes:
        predictions: Ancestry predictions for all samples.
        training_metrics: Model training metrics (if model was trained).
        pca_result: PCA transformation results.
        umap_result: UMAP transformation results (if computed).
        split_result: Cohort splitting results.
        model_path: Path to saved model file, if any.
        common_snps_path: Path to common SNPs file used for prediction.
    """

    predictions: AncestryPredictions
    training_metrics: Optional[TrainingMetrics] = None
    pca_result: Optional[PCAResult] = None
    umap_result: Optional[UMAPResult] = None
    split_result: Optional[SplitResult] = None
    model_path: Optional[Path] = None
    common_snps_path: Optional[Path] = None

    def to_dict(self) -> Dict[str, Any]:
        """Convert to legacy dictionary format for backward compatibility.

        Returns the format expected by existing pipeline code:
        {
            'step': 'predict_ancestry',
            'data': {...},
            'metrics': {...},
            'output': {...}
        }

        Returns:
            Dictionary in legacy format.
        """
        data_dict: Dict[str, Any] = {
            "predict_data": self.predictions.to_dict()["data"],
        }

        if self.training_metrics:
            data_dict["confusion_matrix"] = self.training_metrics.confusion_matrix
            data_dict["label_encoder"] = self.training_metrics.label_encoder_classes

        if self.pca_result:
            data_dict["train_pcs"] = self.pca_result.train_pca
            data_dict["ref_pcs"] = self.pca_result.ref_pca
            data_dict["projected_pcs"] = self.pca_result.projected_pca

        if self.umap_result:
            data_dict["total_umap"] = self.umap_result.total_umap
            data_dict["ref_umap"] = self.umap_result.ref_umap
            data_dict["new_samples_umap"] = self.umap_result.new_umap

        if self.split_result:
            data_dict["labels_list"] = self.split_result.ancestry_labels
            data_dict["pruned_samples"] = self.split_result.pruned_samples

        metrics_dict: Dict[str, Any] = {
            "predicted_counts": self.predictions.ancestry_counts,
        }

        if self.training_metrics:
            metrics_dict["test_accuracy"] = self.training_metrics.test_accuracy

        outfiles_dict: Dict[str, Any] = {
            "predicted_labels": self.predictions.to_dict()["output"],
        }

        if self.split_result:
            outfiles_dict["split_paths"] = list(
                self.split_result.ancestry_paths.values()
            )

        return {
            "step": "predict_ancestry",
            "data": data_dict,
            "metrics": metrics_dict,
            "output": outfiles_dict,
        }
