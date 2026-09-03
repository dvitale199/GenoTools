"""Tests for the training path of AncestryModel.

The defect these guard: `ClassifierConfig.learning_rate` was declared,
documented and validated but never passed to `XGBClassifier`, so XGBoost used
its own gblinear default of 0.5 — at which the optimizer sits at the edge of
numerical divergence. gblinear's default `updater="shotgun"` is Hogwild, so
with the thread count unpinned a race decided per run whether the fit
converged. A diverged fit reaches |intercept| ~1e15, saturates the softmax,
predicts one label for every sample, and used to be pickled anyway.

These assert against the constructed estimator rather than against a helper,
because the bug lived in the constructor call.
"""

import numpy as np
import pandas as pd
import pytest
from sklearn.pipeline import Pipeline

from genotools.ancestry import model as model_module
from genotools.ancestry.config import AncestryConfig, ClassifierConfig
from genotools.ancestry.model import AncestryModel


class _StubFittedPipeline:
    """What the search hands back, so the code after it has something to score.

    The real Pipeline is what these tests inspect; it is never fitted here
    because one UMAP fit at the real fold shape is ~16 seconds.
    """

    def __init__(self, y):
        self._y = np.asarray(y)

    def score(self, X, y):
        return 0.98

    def predict(self, X):
        return self._y[: len(X)]


class _FakeGridSearch:
    """Stand-in for GridSearchCV that records the estimator it was given.

    A real search over the default grid is 1080 pipeline fits at ~16s of UMAP
    each; these tests only need to know what was handed to it.
    """

    captured: list = []

    def __init__(self, estimator, param_grid, **kwargs):
        self.estimator = estimator
        self.param_grid = param_grid
        self.kwargs = kwargs
        type(self).captured.append(self)

    def fit(self, X, y):
        self.best_estimator_ = _StubFittedPipeline(y)
        self.best_params_ = {
            key: values[0] for key, values in self.param_grid.items()
        }
        self.best_score_ = 0.95
        self.cv_results_ = {
            "rank_test_score": np.array([1]),
            "std_test_score": np.array([0.01]),
            "mean_test_score": np.array([0.95]),
        }
        return self


@pytest.fixture
def fake_grid_search(monkeypatch):
    """Swap GridSearchCV in the model module's namespace."""
    _FakeGridSearch.captured = []
    monkeypatch.setattr(model_module, "GridSearchCV", _FakeGridSearch)
    return _FakeGridSearch


def _tiny_training_arrays(n_classes: int = 3, per_class: int = 8):
    """PC-shaped arrays big enough for a Pipeline to be constructed over."""
    rng = np.random.default_rng(0)
    y = np.repeat(np.arange(n_classes), per_class)
    X = rng.normal(size=(len(y), 4)) + y[:, None] * 5.0
    return X, y


def _model_with(classifier: ClassifierConfig) -> AncestryModel:
    model = AncestryModel(config=AncestryConfig(classifier=classifier))
    model.label_encoder = model_module.preprocessing.LabelEncoder()
    model.label_encoder.fit(["A", "B", "C"])
    return model


def _built_xgb(fake_grid_search) -> object:
    assert len(fake_grid_search.captured) == 1
    pipeline = fake_grid_search.captured[0].estimator
    assert isinstance(pipeline, Pipeline)
    return pipeline.named_steps["xgb"]


class TestClassifierConfigReachesTheBooster:
    """Every ClassifierConfig field must arrive at XGBClassifier."""

    def test_learning_rate_is_wired(self, fake_grid_search) -> None:
        """learning_rate reaches the booster instead of XGBoost's 0.5.

        This is the root-cause fix. Left unwired, the booster descends at 0.5
        and diverges on dense reference-panel embeddings.
        """
        X, y = _tiny_training_arrays()
        _model_with(ClassifierConfig(learning_rate=0.1))._train_classifier(
            X, X, y, y
        )
        assert _built_xgb(fake_grid_search).get_params()["learning_rate"] == 0.1

    def test_a_custom_learning_rate_is_not_overwritten(
        self, fake_grid_search
    ) -> None:
        """A non-default learning_rate is honoured, not silently replaced."""
        X, y = _tiny_training_arrays()
        _model_with(ClassifierConfig(learning_rate=0.05))._train_classifier(
            X, X, y, y
        )
        assert _built_xgb(fake_grid_search).get_params()["learning_rate"] == 0.05

    def test_n_jobs_is_wired_and_defaults_to_one(self, fake_grid_search) -> None:
        """n_jobs=1 reaches the booster, taking it off the Hogwild updater.

        Unset, `nthread: 0` uses every core and gblinear's shotgun updater is
        nondeterministic regardless of random_state.
        """
        X, y = _tiny_training_arrays()
        _model_with(ClassifierConfig())._train_classifier(X, X, y, y)
        assert _built_xgb(fake_grid_search).get_params()["n_jobs"] == 1

    def test_n_estimators_and_booster_are_wired(self, fake_grid_search) -> None:
        """The remaining fields arrive too, so none reads as settled-but-dead."""
        X, y = _tiny_training_arrays()
        _model_with(
            ClassifierConfig(n_estimators=42, booster="gblinear")
        )._train_classifier(X, X, y, y)
        params = _built_xgb(fake_grid_search).get_params()
        assert params["n_estimators"] == 42
        assert params["booster"] == "gblinear"

    def test_random_state_is_wired(self, fake_grid_search) -> None:
        """random_state is only meaningful once n_jobs is pinned."""
        X, y = _tiny_training_arrays()
        _model_with(ClassifierConfig(random_state=7))._train_classifier(
            X, X, y, y
        )
        assert _built_xgb(fake_grid_search).get_params()["random_state"] == 7
