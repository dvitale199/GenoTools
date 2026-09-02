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

"""Ancestry prediction model with fit/predict interface.

This module provides the AncestryModel class which implements a scikit-learn
style fit/predict interface for ancestry prediction. The model uses:

1. PCA for initial dimensionality reduction
2. UMAP for manifold learning
3. XGBoost for classification

Example:
    >>> from genotools.ancestry.model import AncestryModel
    >>>
    >>> # Create and train model on reference genotypes + labels
    >>> model = AncestryModel()
    >>> model.fit(raw_ref_data, labels)
    >>>
    >>> # Predict ancestry
    >>> predictions = model.predict(test_data)
    >>> predictions.ancestry_counts
    {'EUR': 100, 'AFR': 50, ...}
    >>>
    >>> # Save and load model
    >>> model.save(Path("model.pkl"))
    >>> loaded = AncestryModel.load(Path("model.pkl"))
"""

import json
import os
import pickle
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import psutil
from sklearn import preprocessing
from sklearn.cluster import Birch
from sklearn.model_selection import GridSearchCV, StratifiedKFold, train_test_split
from sklearn.pipeline import Pipeline
from umap import UMAP
from xgboost import XGBClassifier

from genotools.ancestry.config import (
    AncestryConfig,
    InferenceMode,
    PCAConfig,
)
from genotools.ancestry.reducers.pca import PCAReducer
from genotools.ancestry.reducers.umap_reducer import UMAPReducer
from genotools.ancestry.results import (
    AncestryPredictions,
    AncestryResult,
    PCAResult,
    SplitResult,
    TrainingMetrics,
    UMAPResult,
)
from genotools.core.exceptions import AncestryError
from genotools.core.genotypes import GenotypeData
from genotools.core.logging import get_logger
from genotools.core.provenance import (
    package_versions,
    requirements_lines,
    version_drift,
)


logger = get_logger(__name__)


def _warn_on_version_drift(
    model: "AncestryModel", path: Optional[Path] = None
) -> None:
    """Say so when a model was fitted under different libraries than are loaded.

    This is the whole point of recording versions. A model fitted under
    ``umap-learn==0.5.3`` and loaded under a newer umap unpickles *successfully*
    and produces a different embedding, so the run finishes clean while making
    different ancestry calls. Nothing here blocks the load — the drift is often
    harmless, and a hard failure would strand every model the moment a floating
    dependency moved — but it stops being invisible.

    A model with no recorded versions gets its own warning rather than silence:
    "cannot tell" and "no drift" are different answers.
    """
    recorded = getattr(model, "versions", None)
    if not recorded:
        logger.warning(
            "Model provenance unknown: this model records no library versions, "
            "so GenoTools cannot tell whether it was fitted under different "
            "ones than are installed here. Ancestry calls can differ silently "
            "across library versions; retrain to record provenance."
        )
        return

    drift = version_drift(recorded)
    if not drift:
        return

    changes = "; ".join(f"{name} {was} -> {now}" for name, was, now in drift)
    how = "Reinstall the recorded versions, or retrain, to reproduce them."
    if path is not None and (Path(path) / "requirements.txt").exists():
        how = (
            f"To reproduce the original calls: "
            f"pip install -r {Path(path) / 'requirements.txt'} -- or retrain."
        )
    logger.warning(
        f"Model version drift: {changes}. This model was fitted under the "
        f"recorded versions, and the embedding can differ under different "
        f"ones, so ancestry calls may not match what this model was validated "
        f"on. {how}"
    )


@dataclass
class AncestryModel:
    """Ancestry prediction model with fit/predict interface.

    This class provides a scikit-learn style interface for ancestry prediction.
    The model combines PCA, UMAP, and XGBoost into a prediction pipeline.

    The typical workflow is:
    1. Create model: AncestryModel(config)
    2. Fit on reference genotypes + labels: model.fit(raw_ref_data, labels)
    3. Predict on new samples: model.predict(new_data)
    4. Optionally save: model.save(path)

    Attributes:
        config: Full ancestry configuration.
        pca_reducer: Fitted PCA reducer (set after fit).
        pipeline: Fitted sklearn Pipeline with UMAP + XGBoost (set after fit).
        label_encoder: Fitted LabelEncoder (set after fit).
        best_params: Best hyperparameters from grid search (set after fit).
        training_metrics: Metrics from training (set after fit).
        versions: Library versions present when the model was fitted (set
            after fit); None for a model saved before they were recorded.
    """

    config: AncestryConfig = field(default_factory=AncestryConfig)
    pca_reducer: Optional[PCAReducer] = field(default=None, repr=False)
    pipeline: Optional[Pipeline] = field(default=None, repr=False)
    label_encoder: Optional[preprocessing.LabelEncoder] = field(
        default=None, repr=False
    )
    best_params: Optional[Dict[str, Any]] = field(default=None, repr=False)
    training_metrics: Optional[TrainingMetrics] = field(default=None, repr=False)
    _is_fitted: bool = field(default=False, repr=False)
    _train_pca: Optional[pd.DataFrame] = field(default=None, repr=False)
    common_snps: Optional[List[str]] = field(default=None, repr=False)
    # Library versions behind the fit, captured by fit() and checked by load().
    # None on a model written before they were recorded, which is why load()
    # distinguishes "no drift" from "cannot tell".
    versions: Optional[Dict[str, str]] = field(default=None, repr=False)

    @property
    def is_fitted(self) -> bool:
        """Whether the model has been fitted."""
        return self._is_fitted

    def _prepare_training_data(
        self,
        raw_data: pd.DataFrame,
        labels: pd.Series,  # type: ignore[type-arg]
    ) -> Dict[str, Any]:
        """Prepare data for training with train/test split.

        Args:
            raw_data: DataFrame with FID, IID, and SNP columns.
            labels: Series mapping IID to ancestry label.

        Returns:
            Dictionary with train/test split data.
        """
        # Separate IDs and SNPs
        ids = raw_data[["FID", "IID"]]
        snps = raw_data.drop(columns=["FID", "IID"], errors="ignore")

        # Get labels for samples in order
        y = ids["IID"].map(labels)

        # Encode labels
        self.label_encoder = preprocessing.LabelEncoder()
        y_encoded = self.label_encoder.fit_transform(y)

        # Train/test split
        X_train, X_test, y_train, y_test, idx_train, idx_test = train_test_split(
            snps.values,
            y_encoded,
            np.arange(len(snps)),
            test_size=self.config.training.test_size,
            random_state=self.config.training.random_state,
            stratify=y_encoded,
        )

        return {
            "X_train": X_train,
            "X_test": X_test,
            "y_train": y_train,
            "y_test": y_test,
            "train_ids": ids.iloc[idx_train].reset_index(drop=True),
            "test_ids": ids.iloc[idx_test].reset_index(drop=True),
            "train_labels": self.label_encoder.inverse_transform(y_train),
            "test_labels": self.label_encoder.inverse_transform(y_test),
            "label_encoder": self.label_encoder,
        }

    def _run_pca(
        self,
        X_train: np.ndarray,  # type: ignore[type-arg]
        X_test: np.ndarray,  # type: ignore[type-arg]
        train_ids: pd.DataFrame,
        test_ids: pd.DataFrame,
        train_labels: np.ndarray,  # type: ignore[type-arg]
        test_labels: np.ndarray,  # type: ignore[type-arg]
    ) -> Dict[str, Any]:
        """Run PCA on training data and transform test data.

        Args:
            X_train: Training genotype matrix.
            X_test: Test genotype matrix.
            train_ids: DataFrame with FID, IID for training samples.
            test_ids: DataFrame with FID, IID for test samples.
            train_labels: Labels for training samples.
            test_labels: Labels for test samples.

        Returns:
            Dictionary with PCA results.
        """
        # Create and fit PCA reducer
        self.pca_reducer = PCAReducer(config=self.config.pca)
        self.pca_reducer.fit(X_train)

        # Transform data
        X_train_pca = self.pca_reducer.transform(X_train, return_dataframe=True)
        X_test_pca = self.pca_reducer.transform(X_test, return_dataframe=True)

        assert isinstance(X_train_pca, pd.DataFrame)
        assert isinstance(X_test_pca, pd.DataFrame)

        # Build labeled DataFrames
        train_pca = pd.concat(
            [train_ids.reset_index(drop=True), X_train_pca], axis=1
        )
        train_pca["label"] = train_labels

        test_pca = pd.concat([test_ids.reset_index(drop=True), X_test_pca], axis=1)
        test_pca["label"] = test_labels

        ref_pca = pd.concat([train_pca, test_pca], ignore_index=True)

        self._train_pca = train_pca

        return {
            "X_train": X_train_pca.values,
            "X_test": X_test_pca.values,
            "train_pca": train_pca,
            "test_pca": test_pca,
            "ref_pca": ref_pca,
        }

    def _train_classifier(
        self,
        X_train: np.ndarray,  # type: ignore[type-arg]
        X_test: np.ndarray,  # type: ignore[type-arg]
        y_train: np.ndarray,  # type: ignore[type-arg]
        y_test: np.ndarray,  # type: ignore[type-arg]
    ) -> Dict[str, Any]:
        """Train UMAP + XGBoost classifier with grid search.

        Args:
            X_train: Training PCA components.
            X_test: Test PCA components.
            y_train: Training labels (encoded).
            y_test: Test labels (encoded).

        Returns:
            Dictionary with trained classifier and metrics.
        """
        # Build parameter grid from config
        param_grid = {
            "umap__n_neighbors": list(self.config.grid_search.umap_n_neighbors),
            "umap__n_components": list(self.config.grid_search.umap_n_components),
            "umap__a": list(self.config.grid_search.umap_a),
            "umap__b": list(self.config.grid_search.umap_b),
            "xgb__lambda": list(self.config.grid_search.xgb_lambda),
        }

        # Calculate optimal number of parallel jobs
        n_jobs = self.config.training.n_jobs
        if n_jobs is None:
            available_ram_gb = psutil.virtual_memory().total / (1024**3)
            gb_per_worker = self.config.training.gb_per_worker
            max_workers_by_ram = max(1, int((available_ram_gb * 0.8) / gb_per_worker))
            cpu_count = os.cpu_count() or 1
            n_jobs = min(cpu_count, max_workers_by_ram)

        logger.info(
            f"Using {n_jobs} parallel workers for GridSearchCV "
            f"(RAM: {psutil.virtual_memory().total / (1024**3):.1f}GB)"
        )

        # Build pipeline
        umap = UMAP(random_state=self.config.umap.random_state)
        xgb = XGBClassifier(
            booster=self.config.classifier.booster,
            random_state=self.config.classifier.random_state,
        )
        pipeline = Pipeline([("umap", umap), ("xgb", xgb)])

        # Grid search with cross-validation
        cv = StratifiedKFold(
            n_splits=self.config.grid_search.cv_folds,
            shuffle=True,
            random_state=self.config.training.random_state,
        )
        grid_search = GridSearchCV(
            pipeline,
            param_grid,
            cv=cv,
            scoring=self.config.grid_search.scoring,
            n_jobs=n_jobs,
            verbose=1,
        )
        grid_search.fit(X_train, y_train)

        # Get results
        results_df = pd.DataFrame(grid_search.cv_results_)
        top_results = results_df[results_df["rank_test_score"] == 1]

        train_acc = grid_search.best_score_
        train_ci_interval = 1.96 * float(top_results["std_test_score"].iloc[0])

        self.pipeline = grid_search.best_estimator_
        self.best_params = grid_search.best_params_

        # Evaluate on test set
        test_acc = float(self.pipeline.score(X_test, y_test))
        test_margin = 1.96 * np.sqrt((test_acc * (1 - test_acc)) / len(y_test))

        # Confusion matrix
        y_pred = self.pipeline.predict(X_test)
        from sklearn import metrics

        confusion_matrix = metrics.confusion_matrix(y_test, y_pred)

        logger.info(f"Training Balanced Accuracy: {train_acc:.4f}")
        logger.info(f"Test Balanced Accuracy: {test_acc:.4f}")
        logger.info(f"Best Parameters: {grid_search.best_params_}")

        # Store training metrics
        self.training_metrics = TrainingMetrics(
            train_accuracy=train_acc,
            test_accuracy=test_acc,
            train_accuracy_ci=(train_acc - train_ci_interval, train_acc + train_ci_interval),
            test_accuracy_ci=(test_acc - test_margin, test_acc + test_margin),
            confusion_matrix=confusion_matrix,
            best_params=grid_search.best_params_,
            label_encoder_classes=list(self.label_encoder.classes_)
            if self.label_encoder
            else [],
        )

        return {
            "pipeline": self.pipeline,
            "best_params": grid_search.best_params_,
            "train_accuracy": train_acc,
            "test_accuracy": test_acc,
            "confusion_matrix": confusion_matrix,
        }

    def _predict_admixed(
        self,
        projected: pd.DataFrame,
        train_pca: pd.DataFrame,
    ) -> pd.DataFrame:
        """Identify samples with complex admixture history.

        Samples closest to the global centroid rather than any specific
        ancestry centroid are classified as CAH (Complex Admixture History).

        Args:
            projected: DataFrame with projected samples (PCs + predicted label).
            train_pca: DataFrame with training PCA data and labels.

        Returns:
            Updated projected DataFrame with CAH labels.
        """
        # Copy training data for admixture calculation
        train_pca_admix = train_pca.copy()

        # Handle CAS (Central Asian) bimodal distribution
        cas_train = train_pca_admix[train_pca_admix["label"] == "CAS"]
        other_train = train_pca_admix[train_pca_admix["label"] != "CAS"]

        if len(cas_train) > 0:
            cas_ids = cas_train[["FID", "IID"]].reset_index(drop=True)
            cas_labels = cas_train["label"].reset_index(drop=True)
            cas_pcs = cas_train.drop(
                columns=["FID", "IID", "label"], errors="ignore"
            ).reset_index(drop=True)

            # Cluster CAS into two groups based on first 3 PCs
            birch = Birch(n_clusters=self.config.admixed.n_clusters_cas)
            birch.fit(cas_pcs[["PC1", "PC2", "PC3"]])
            cas_clusters = birch.predict(cas_pcs[["PC1", "PC2", "PC3"]])

            # Reconstruct CAS training data with subclusters
            cas_train_new = pd.concat([cas_ids, cas_pcs, cas_labels], axis=1)
            cas_train_new["label"] = np.where(cas_clusters == 0, "CAS", "CAS2")
            train_pca_admix = pd.concat(
                [other_train, cas_train_new], axis=0, ignore_index=True
            )

        # Calculate centroids for each ancestry
        pc_cols = [c for c in train_pca_admix.columns if c.startswith("PC")]
        pc_centroids = pd.DataFrame()

        train_pcs = train_pca_admix[pc_cols]
        pc_centroids["ALL"] = train_pcs.mean(axis=0)

        for ancestry in train_pca_admix["label"].unique():
            ancestry_data = train_pca_admix[train_pca_admix["label"] == ancestry]
            ancestry_pcs = ancestry_data[pc_cols]
            pc_centroids[ancestry] = ancestry_pcs.mean(axis=0)

        # Calculate distance to each centroid for projected samples
        projected_copy = projected.copy()
        projected_ids = projected_copy[["FID", "IID", "label"]]
        projected_pcs = projected_copy[pc_cols].copy()

        for ancestry in pc_centroids.columns:
            centroid = pc_centroids[ancestry]
            distances = projected_pcs.apply(
                lambda row: np.sqrt(np.sum((row - centroid) ** 2)), axis=1
            )
            projected_pcs[ancestry] = distances

        # Find minimum distance and corresponding ancestry
        distance_cols = list(pc_centroids.columns)
        projected_pcs["min_distance"] = projected_pcs[distance_cols].min(axis=1)
        projected_pcs["min_distance_ancestry"] = (
            projected_pcs[distance_cols]
            .eq(projected_pcs["min_distance"], axis=0)
            .idxmax(axis=1)
        )

        # Update labels - samples closest to ALL centroid become CAH
        result = pd.concat([projected_ids, projected_pcs[pc_cols]], axis=1)
        result["label"] = np.where(
            projected_pcs["min_distance_ancestry"] == "ALL",
            "CAH",
            projected_ids["label"],
        )

        return result

    def fit(
        self,
        raw_ref_data: pd.DataFrame,
        labels: pd.Series,  # type: ignore[type-arg]
        out_path: Optional[Path] = None,
    ) -> "AncestryModel":
        """Fit the ancestry model on reference data.

        Args:
            raw_ref_data: DataFrame with reference genotype data.
                Must have FID, IID columns followed by SNP columns.
            labels: Series mapping IID to ancestry label.
            out_path: Optional path for saving intermediate files.

        Returns:
            Self for method chaining.

        Raises:
            AncestryError: If training fails.
        """
        logger.info("Fitting ancestry model on reference data")

        # Prepare training data
        split_data = self._prepare_training_data(raw_ref_data, labels)

        # Run PCA
        pca_result = self._run_pca(
            X_train=split_data["X_train"],
            X_test=split_data["X_test"],
            train_ids=split_data["train_ids"],
            test_ids=split_data["test_ids"],
            train_labels=split_data["train_labels"],
            test_labels=split_data["test_labels"],
        )

        # Save PCA outputs if path provided
        if out_path:
            pca_result["train_pca"].to_csv(
                f"{out_path}_labeled_train_pca.txt", sep="\t", index=False
            )
            pca_result["ref_pca"].to_csv(
                f"{out_path}_labeled_ref_pca.txt", sep="\t", index=False
            )
            if self.pca_reducer:
                self.pca_reducer.save_eigenvalues(
                    Path(f"{out_path}_pca_eigenvalues.txt")
                )

        # Train classifier
        classifier_result = self._train_classifier(
            X_train=pca_result["X_train"],
            X_test=pca_result["X_test"],
            y_train=split_data["y_train"],
            y_test=split_data["y_test"],
        )

        # Captured at fit rather than at save: these are the versions that
        # produced *this* fit, and they travel with the pickle so the
        # single-file format carries provenance too.
        self.versions = package_versions()

        self._is_fitted = True
        logger.info("Ancestry model fitted successfully")

        return self

    def predict(
        self,
        raw_geno_data: pd.DataFrame,
        geno_ids: pd.DataFrame,
        out_path: Optional[Path] = None,
        detect_admixed: bool = True,
    ) -> AncestryPredictions:
        """Predict ancestry for new samples.

        Args:
            raw_geno_data: DataFrame with genotype data (samples x SNPs).
            geno_ids: DataFrame with FID, IID columns.
            out_path: Optional path for saving outputs.
            detect_admixed: Whether to detect admixed samples. Default True.

        Returns:
            AncestryPredictions with predicted ancestries.

        Raises:
            AncestryError: If model not fitted or prediction fails.
        """
        if not self._is_fitted:
            raise AncestryError("Model must be fitted before prediction")

        if self.pca_reducer is None or self.pipeline is None:
            raise AncestryError("Model is in invalid state")

        if self.label_encoder is None:
            raise AncestryError("Label encoder not set")

        logger.info(f"Predicting ancestry for {len(raw_geno_data)} samples")

        # Get SNP columns (exclude ID columns)
        snp_cols = [c for c in raw_geno_data.columns if c not in ("FID", "IID", "label")]
        X_new = raw_geno_data[snp_cols].values

        # Project through PCA
        X_new_pca = self.pca_reducer.transform(X_new, return_dataframe=True)
        assert isinstance(X_new_pca, pd.DataFrame)

        # Predict with pipeline
        y_pred_encoded = self.pipeline.predict(X_new_pca.values)
        y_pred = self.label_encoder.inverse_transform(y_pred_encoded)

        # Build projected DataFrame for admixture detection
        projected = pd.concat(
            [geno_ids.reset_index(drop=True), X_new_pca], axis=1
        )
        projected["label"] = y_pred

        # Detect admixed samples if requested
        if detect_admixed and self._train_pca is not None:
            projected = self._predict_admixed(projected, self._train_pca)
            y_pred = projected["label"].values

        # Save outputs if path provided
        output_path = None
        if out_path:
            output_path = Path(f"{out_path}_umap_linearsvc_predicted_labels.txt")
            projected[["FID", "IID", "label"]].to_csv(
                output_path, sep="\t", index=False
            )
            projected.to_csv(
                f"{out_path}_projected_new_pca.txt", sep="\t", index=False
            )

        # Create predictions result
        predictions_df = geno_ids.copy()
        predictions_df["predicted_ancestry"] = y_pred

        logger.info(f"Predicted ancestry counts: {pd.Series(y_pred).value_counts().to_dict()}")

        return AncestryPredictions(
            predictions=predictions_df,
            output_path=output_path,
        )

    def save(self, path: Path) -> Path:
        """Save fitted model to a directory.

        Creates a directory at ``path`` containing:
        - ``pipeline.pkl``: Full pickled AncestryModel.
        - ``common_snps.txt``: Plain text rsIDs, one per line.
        - ``metadata.json``: Config, best_params, training_metrics, and the
          library versions the fit ran under.
        - ``requirements.txt``: those same versions as pip requirements, so the
          fitting environment can be recreated rather than merely identified.

        Args:
            path: Output directory path.

        Returns:
            Path to the saved directory.

        Raises:
            AncestryError: If model not fitted.
        """
        if not self._is_fitted:
            raise AncestryError("Cannot save unfitted model")

        path = Path(path)
        path.mkdir(parents=True, exist_ok=True)

        # Full pickle
        pkl_path = path / "pipeline.pkl"
        with open(pkl_path, "wb") as f:
            pickle.dump(self, f)

        # Common SNPs (plain text, one per line)
        if self.common_snps is not None:
            snps_path = path / "common_snps.txt"
            with open(snps_path, "w") as f:
                for snp in self.common_snps:
                    f.write(f"{snp}\n")

        # The environment this fit needs, as something pip can consume.
        # Knowing a model drifted is only half of it; this is the other half.
        if self.versions is not None:
            reqs_path = path / "requirements.txt"
            with open(reqs_path, "w") as f:
                f.write(
                    "# Environment this ancestry model was fitted under.\n"
                    "# Recreate it with: pip install -r requirements.txt\n"
                )
                for line in requirements_lines(self.versions):
                    f.write(f"{line}\n")

        # Human-readable metadata
        metadata = self._build_metadata()
        meta_path = path / "metadata.json"
        with open(meta_path, "w") as f:
            json.dump(metadata, f, indent=2, default=str)

        logger.info(f"Model saved to directory: {path}")
        return path

    def _build_metadata(self) -> Dict[str, Any]:
        """Build metadata dictionary for JSON serialization."""
        metadata: Dict[str, Any] = {}

        metadata["config"] = self.config.to_dict()

        if self.best_params is not None:
            metadata["best_params"] = self.best_params

        if self.training_metrics is not None:
            metadata["training_metrics"] = self.training_metrics.to_dict()

        if self.common_snps is not None:
            metadata["n_common_snps"] = len(self.common_snps)

        if self.versions is not None:
            metadata["versions"] = dict(self.versions)

        return metadata

    @classmethod
    def load(cls, path: Path) -> "AncestryModel":
        """Load fitted model from a directory or a single pickle file.

        Supports two formats:
        - **Directory format**: reads ``pipeline.pkl`` from the directory.
          If the pickle lacks ``common_snps``, falls back to
          ``common_snps.txt`` in the same directory.
        - **Single-file format**: reads a ``.pkl`` file directly.

        Both must unpickle to an ``AncestryModel``. Models written by
        GenoTools 1.x hold an ``sklearn.pipeline.Pipeline`` instead and are
        rejected - "single-file" is 2.0's own layout, not a 1.x model.

        Args:
            path: Path to a model directory or a single ``.pkl`` file.

        Returns:
            Loaded AncestryModel.

        Raises:
            FileNotFoundError: If path doesn't exist.
            AncestryError: If file contains invalid model.
        """
        path = Path(path)
        if not path.exists():
            raise FileNotFoundError(f"Model not found: {path}")

        if path.is_dir():
            # Directory format
            pkl_path = path / "pipeline.pkl"
            if not pkl_path.exists():
                raise FileNotFoundError(
                    f"pipeline.pkl not found in model directory: {path}"
                )
            with open(pkl_path, "rb") as f:
                model = pickle.load(f)

            # Back-fill common_snps from text file if missing in pickle
            if getattr(model, "common_snps", None) is None:
                snps_path = path / "common_snps.txt"
                if snps_path.exists():
                    with open(snps_path, "r") as f:
                        model.common_snps = [
                            line.strip() for line in f if line.strip()
                        ]
        else:
            # Single-file format (2.0's own; a 1.x pickle fails the check below)
            with open(path, "rb") as f:
                model = pickle.load(f)

        if not isinstance(model, cls):
            # A 1.x model pickles an sklearn Pipeline. Naming that outright
            # saves the user guessing at a type mismatch, since --model is
            # exactly where a 1.x model gets pointed at 2.0.
            hint = ""
            if type(model).__module__.split(".")[0] == "sklearn":
                hint = (
                    " This looks like a GenoTools 1.x model, which 2.0 cannot "
                    "load. Pass a model directory written by 2.0, or retrain "
                    "by dropping --model and passing --ref-panel/--ref-labels."
                )
            raise AncestryError(
                f"Invalid model file: expected AncestryModel, got "
                f"{type(model)}.{hint}"
            )

        _warn_on_version_drift(model, path)

        logger.info(f"Model loaded from: {path}")
        return model

    @staticmethod
    def load_common_snps(model_dir: Path) -> List[str]:
        """Load common SNPs list from a model directory.

        This reads only the plain-text ``common_snps.txt`` file without
        deserializing the full pickle, which avoids dependency issues.

        Args:
            model_dir: Path to model directory.

        Returns:
            List of rsID strings.

        Raises:
            FileNotFoundError: If ``common_snps.txt`` not found.
        """
        snps_path = Path(model_dir) / "common_snps.txt"
        if not snps_path.exists():
            raise FileNotFoundError(
                f"common_snps.txt not found in: {model_dir}"
            )
        with open(snps_path, "r") as f:
            return [line.strip() for line in f if line.strip()]


def load_trained_pipeline(
    model_path: Path,
    X_test: np.ndarray,  # type: ignore[type-arg]
    y_test: np.ndarray,  # type: ignore[type-arg]
) -> Dict[str, Any]:
    """Load a trained sklearn Pipeline from pickle.

    This is a compatibility function for loading old-style models
    that only contain the UMAP+XGBoost pipeline.

    Args:
        model_path: Path to pickled pipeline.
        X_test: Test features for evaluation.
        y_test: Test labels for evaluation.

    Returns:
        Dictionary with pipeline, accuracy, and parameters.
    """
    from sklearn import metrics

    with open(model_path, "rb") as f:
        pipeline = pickle.load(f)

    test_acc = float(pipeline.score(X_test, y_test))
    y_pred = pipeline.predict(X_test)
    confusion_matrix = metrics.confusion_matrix(y_test, y_pred)
    params = pipeline.get_params()

    logger.info(f"Loaded pipeline with test accuracy: {test_acc:.4f}")

    return {
        "pipeline": pipeline,
        "test_accuracy": test_acc,
        "confusion_matrix": confusion_matrix,
        "params": params,
    }
