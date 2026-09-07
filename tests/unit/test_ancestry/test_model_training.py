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
from genotools.ancestry.config import (
    AncestryConfig,
    ClassifierConfig,
    GridSearchConfig,
    PCAConfig,
    TrainingConfig,
)
from genotools.ancestry.model import FIT_FALLBACK_LEARNING_RATES, AncestryModel
from genotools.core.exceptions import AncestryError


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
        # No `best_estimator_`: the production code sets `refit=False` so that
        # the winning candidate is fitted and validated outside the search.
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
    """Swap GridSearchCV, and stub the fit that follows it.

    These tests inspect the estimator handed *to* the search. Letting
    `_fit_and_validate_candidate` run afterwards would fit a real UMAP on
    every one of them for nothing; `TestRealFit` covers that path.
    """
    _FakeGridSearch.captured = []
    monkeypatch.setattr(model_module, "GridSearchCV", _FakeGridSearch)

    def stub(self, pipeline, best_params, X_train, y_train, cv_score):
        validation = model_module.validate_fit(
            _Fitted(y_train), y_train, y_train, n_classes=len(np.unique(y_train))
        )
        return {
            "pipeline": _Fitted(y_train),
            "params": best_params,
            "validation": validation,
            "attempts": [{"learning_rate": 0.1, "accepted": True}],
        }

    monkeypatch.setattr(model_module.AncestryModel, "_fit_and_validate_candidate", stub)
    return _FakeGridSearch


class _Fitted:
    """A model-shaped stand-in for the fit that follows the search."""

    coef_ = np.array([[1.0]])
    intercept_ = np.array([1.0])

    def __init__(self, y):
        self._y = np.asarray(y)

    def predict(self, X):
        return self._y[: len(X)]

    def score(self, X, y):
        return 0.98


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


class _StagedPipeline:
    """A pipeline whose fit outcome is dictated per learning rate.

    `set_params` is what the fallback loop uses to change the learning rate,
    so the stub keys its behaviour off that. Every rate maps to either a
    healthy fit or the collapsed one: |intercept| 3e15 predicting a single
    label, which is what the real defect produced.
    """

    def __init__(self, outcomes, y_train):
        self.outcomes = outcomes
        self._y = np.asarray(y_train)
        self.params = {}
        self.fits = []

    # sklearn.base.clone needs these two.
    def get_params(self, deep=True):
        return {"outcomes": self.outcomes, "y_train": self._y}

    def set_params(self, **params):
        self.params.update(params)
        return self

    @property
    def learning_rate(self):
        return self.params.get("xgb__learning_rate")

    def fit(self, X, y):
        self.fits.append(self.learning_rate)
        healthy = self.outcomes[self.learning_rate]
        if healthy:
            self.coef_ = np.array([[1.03, -0.4]])
            self.intercept_ = np.array([4.24])
        else:
            self.coef_ = np.array([[1727.0]])
            self.intercept_ = np.array([3.0e15])
        return self

    def predict(self, X):
        n = len(X)
        if self.outcomes[self.learning_rate]:
            return self._y[:n]
        return np.full(n, self._y[0])

    def score(self, X, y):
        return 0.98 if self.outcomes[self.learning_rate] else 0.15


def _staged_model(outcomes, y, **training) -> AncestryModel:
    model = AncestryModel(
        config=AncestryConfig(training=TrainingConfig(**training))
    )
    model.label_encoder = model_module.preprocessing.LabelEncoder()
    model.label_encoder.fit([f"L{i}" for i in range(10)])
    return model


class TestFitFallback:
    """A collapsed fit must never become `self.pipeline`."""

    def _run(self, outcomes, **training):
        """Drive `_fit_and_validate_candidate` through the staged pipeline."""
        y = np.repeat(np.arange(10), 10)
        model = _staged_model(outcomes, y, **training)
        pipeline = _StagedPipeline(outcomes, y)
        return model, pipeline, y, model._fit_and_validate_candidate(
            pipeline=pipeline,
            best_params={"umap__a": 0.75, "xgb__lambda": 0.001},
            X_train=np.zeros((len(y), 3)),
            y_train=y,
            cv_score=0.94,
        )

    def test_a_healthy_first_fit_is_accepted_without_fallback(self) -> None:
        _, _, _, accepted = self._run({0.1: True})
        assert len(accepted["attempts"]) == 1
        assert accepted["attempts"][0]["accepted"] is True
        assert accepted["params"]["xgb__learning_rate"] == 0.1

    def test_the_configured_rate_is_tried_first(self) -> None:
        _, _, _, accepted = self._run({0.1: True, 0.05: True})
        assert [a["learning_rate"] for a in accepted["attempts"]] == [0.1]

    def test_a_collapsed_fit_steps_the_learning_rate_down(self) -> None:
        """Retrying the same rate is deterministic now, so it must change."""
        _, _, _, accepted = self._run({0.1: False, 0.05: True})
        rates = [a["learning_rate"] for a in accepted["attempts"]]
        assert rates == [0.1, 0.05]
        assert accepted["attempts"][0]["accepted"] is False
        assert accepted["attempts"][-1]["accepted"] is True
        assert accepted["params"]["xgb__learning_rate"] == 0.05

    def test_the_accepted_pipeline_is_the_one_that_passed(self) -> None:
        _, _, _, accepted = self._run({0.1: False, 0.05: True})
        assert accepted["pipeline"].learning_rate == 0.05
        assert not accepted["validation"].collapsed

    def test_a_successful_fallback_is_loud(self, caplog) -> None:
        import logging

        with caplog.at_level(logging.WARNING):
            self._run({0.1: False, 0.05: True})
        text = caplog.text
        assert "unusable fit" in text
        assert "learning_rate=0.05" in text
        assert "0.94" in text

    def test_exhausting_every_candidate_raises(self) -> None:
        outcomes = {0.1: False, 0.05: False, 0.01: False, 0.005: False}
        with pytest.raises(AncestryError) as excinfo:
            self._run(outcomes)
        message = str(excinfo.value)
        assert "no model was saved" in message
        assert "learning_rate=0.05" in message
        assert "--ancestry-min-fit-accuracy 0" in message

    def test_zero_fallbacks_refuses_immediately(self) -> None:
        with pytest.raises(AncestryError):
            self._run({0.1: False, 0.05: True}, fit_fallbacks=0)

    def test_the_fallback_budget_is_honoured(self) -> None:
        outcomes = dict.fromkeys(FIT_FALLBACK_LEARNING_RATES, False)
        with pytest.raises(AncestryError):
            self._run(outcomes, fit_fallbacks=1)
        model = _staged_model(outcomes, np.repeat(np.arange(10), 10),
                              fit_fallbacks=1)
        assert model._fit_attempt_learning_rates() == [0.1, 0.05]

    def test_every_fallback_is_below_the_configured_rate(self) -> None:
        """Raising the rate is the direction that diverges."""
        model = AncestryModel(
            config=AncestryConfig(
                classifier=ClassifierConfig(learning_rate=0.05),
                training=TrainingConfig(fit_fallbacks=3),
            )
        )
        rates = model._fit_attempt_learning_rates()
        assert rates[0] == 0.05
        assert all(rate < 0.05 for rate in rates[1:])

    def test_a_zero_floor_still_refuses_a_constant_fit(self) -> None:
        """--ancestry-min-fit-accuracy 0 drops the accuracy gate only.

        A diverged, single-label model is refused whatever the floor says --
        there is no threshold at which shipping it is the right answer.
        """
        with pytest.raises(AncestryError):
            self._run({0.1: False}, min_fit_balanced_accuracy=0.0,
                      fit_fallbacks=0)

    def test_attempts_carry_the_measurements(self) -> None:
        _, _, _, accepted = self._run({0.1: False, 0.05: True})
        first = accepted["attempts"][0]
        assert first["diverged"] is True
        assert first["max_abs_intercept"] == pytest.approx(3.0e15)
        assert first["n_classes_predicted"] == 1


class TestFailedCandidateCounting:
    """A candidate that raises is scored NaN and ranked last, silently."""

    def test_nan_candidates_are_counted_and_reported(
        self, fake_grid_search, caplog
    ) -> None:
        import logging

        class _WithFailures(_FakeGridSearch):
            def fit(self, X, y):
                super().fit(X, y)
                self.cv_results_ = {
                    "rank_test_score": np.array([1, 2, 3, 4]),
                    "std_test_score": np.array([0.01, 0.01, 0.01, 0.01]),
                    "mean_test_score": np.array([0.95, 0.90, np.nan, np.nan]),
                }
                return self

        fake_grid_search.captured = []
        X, y = _tiny_training_arrays(n_classes=3, per_class=8)
        model = _model_with(ClassifierConfig())
        with caplog.at_level(logging.WARNING):
            monkey = _WithFailures
            model_module.GridSearchCV = monkey
            try:
                model._train_classifier(X, X, y, y)
            finally:
                model_module.GridSearchCV = _FakeGridSearch
        assert model.training_metrics.n_failed_candidates == 2
        assert model.training_metrics.n_grid_candidates == 4
        assert "2 of 4 grid candidates failed to fit" in caplog.text


class TestSearchDoesNotRefit:
    """The winner is fitted here so it can be refused before it is pickled."""

    def test_grid_search_is_constructed_with_refit_false(
        self, fake_grid_search
    ) -> None:
        X, y = _tiny_training_arrays()
        _model_with(ClassifierConfig())._train_classifier(X, X, y, y)
        assert fake_grid_search.captured[0].kwargs["refit"] is False

    def test_cv_results_are_kept_for_the_report(self, fake_grid_search) -> None:
        X, y = _tiny_training_arrays()
        model = _model_with(ClassifierConfig())
        model._train_classifier(X, X, y, y)
        assert isinstance(model._cv_results, pd.DataFrame)
        assert len(model._cv_results) == 1


def _synthetic_reference(
    n_per: int = 16,
    n_snps: int = 120,
    n_labels: int = 4,
    structured: bool = True,
):
    """Genotype-shaped reference data, with or without ancestry structure.

    When `structured`, allele frequencies are shifted per label over the first
    half of the SNPs, so PCA finds the labels and the classifier has something
    to learn. When not, every label draws from the same frequencies -- there is
    no signal, so a fit cannot beat chance and the floor must refuse it. That
    is the only way to provoke an unusable model on demand now that the
    training race is gone. `PCAReducer` needs 50 samples, so keep n_per *
    n_labels above that.
    """
    rng = np.random.default_rng(7)
    names = [f"L{index}" for index in range(n_labels)]
    rows, labels = [], []
    for index, name in enumerate(names):
        freqs = rng.uniform(0.05, 0.95, n_snps)
        if structured:
            freqs[: n_snps // 2] *= 0.3 + 0.25 * index
        for _ in range(n_per):
            rows.append(rng.binomial(2, np.clip(freqs, 0.01, 0.99)))
            labels.append(name)
    frame = pd.DataFrame(
        np.asarray(rows, dtype=float), columns=[f"snp{k}" for k in range(n_snps)]
    )
    frame.insert(0, "IID", [f"S{k:04d}" for k in range(len(frame))])
    frame.insert(0, "FID", frame["IID"])
    return frame, pd.Series(labels, index=frame["IID"], name="label")


def _diverging_config(**training):
    """The same tiny grid, with a learning rate high enough to blow the fit up.

    10 is far outside anything the CLI offers; it is the only way to provoke
    the real failure on demand now that `n_jobs=1` has removed the thread race
    that used to do it at random.
    """
    import dataclasses

    from genotools.ancestry.config import ClassifierConfig

    return dataclasses.replace(
        _tiny_config(**training), classifier=ClassifierConfig(learning_rate=10.0)
    )


@pytest.fixture(scope="module")
def fitted_model():
    """One real fit, shared: each UMAP fit here costs seconds."""
    reference, labels = _synthetic_reference()
    return AncestryModel(config=_tiny_config()).fit(reference, labels)



def _tiny_config(**training):
    """The full grid collapsed to one candidate: 1080 fits become 2."""
    from genotools.ancestry.config import AncestryConfig

    return AncestryConfig(
        pca=PCAConfig(n_components=6),
        grid_search=GridSearchConfig(
            umap_n_neighbors=(5,),
            umap_n_components=(2,),
            umap_a=(1.0,),
            umap_b=(0.5,),
            xgb_lambda=(1.0,),
            cv_folds=2,
        ),
        training=TrainingConfig(n_jobs=1, **training),
    )


class TestRealFit:
    """One real `fit()` end to end, on a grid shrunk to a single candidate.

    Everything above stubs the search. This is the only test that runs the
    actual path a production training run takes, which is where the pieces
    could disagree about their interfaces without any of them being wrong.
    """

    def test_a_healthy_fit_records_every_measurement(self, fitted_model) -> None:
        metrics = fitted_model.training_metrics
        assert metrics is not None
        assert metrics.test_balanced_accuracy is not None
        assert metrics.train_balanced_accuracy is not None
        assert metrics.cv_balanced_accuracy == metrics.train_accuracy
        assert set(metrics.baseline_scores) == {"15-NN", "nearest-centroid"}
        assert metrics.n_grid_candidates == 1
        assert metrics.n_failed_candidates == 0
        assert metrics.fit_validation["collapsed"] is False
        assert metrics.fit_validation["diverged"] is False
        assert len(metrics.fit_attempts) == 1

    def test_the_learning_rate_used_is_recorded_in_best_params(
        self, fitted_model
    ) -> None:
        """The grid does not tune it, so nothing else would say what it was."""
        assert fitted_model.best_params["xgb__learning_rate"] == 0.1

    def test_the_fitted_pipeline_predicts(self, fitted_model) -> None:
        """`refit=False` means this pipeline was fitted by our code, not the
        search -- so it has to actually work."""
        assert fitted_model.is_fitted
        assert len(fitted_model.pipeline.predict(np.zeros((3, 6)))) == 3

    def test_a_real_divergence_fails_the_run(self, tmp_path) -> None:
        """The whole point: an unusable model is refused, not saved.

        The race that produced the original collapse is gone, so the failure
        is provoked the other way -- by raising the learning rate until the
        same optimizer diverges. The signature it produces is the PPMI one:
        |intercept| ~1e17, one label for every sample, balanced accuracy
        exactly 1/n_labels.
        """
        reference, labels = _synthetic_reference()
        model = AncestryModel(config=_diverging_config(fit_fallbacks=0))

        with pytest.raises(AncestryError) as excinfo:
            model.fit(reference, labels)

        message = str(excinfo.value)
        assert "no model was saved" in message
        assert "predicts 1/4 labels" in message
        assert not model.is_fitted
        with pytest.raises(AncestryError, match="Cannot save unfitted model"):
            model.save(tmp_path / "model")

    def test_the_fallback_rescues_a_diverged_fit(self, caplog) -> None:
        """Stepping the learning rate down is what recovery actually is."""
        import logging

        reference, labels = _synthetic_reference()
        model = AncestryModel(config=_diverging_config(fit_fallbacks=1))
        with caplog.at_level(logging.WARNING):
            model.fit(reference, labels)

        attempts = model.training_metrics.fit_attempts
        assert [a["learning_rate"] for a in attempts] == [10.0, 0.1]
        assert attempts[0]["diverged"] is True
        assert attempts[0]["n_classes_predicted"] == 1
        assert attempts[-1]["accepted"] is True
        assert model.best_params["xgb__learning_rate"] == 0.1
        assert "unusable fit" in caplog.text

    def test_a_zero_floor_still_refuses_a_diverged_fit(self) -> None:
        """No threshold makes shipping a constant classifier the right answer."""
        reference, labels = _synthetic_reference()
        config = _diverging_config(fit_fallbacks=0, min_fit_balanced_accuracy=0.0)
        with pytest.raises(AncestryError):
            AncestryModel(config=config).fit(reference, labels)
