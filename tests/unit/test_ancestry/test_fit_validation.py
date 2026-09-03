"""Tests for training-path fit validation.

The failure being guarded: a diverged gblinear fit with |intercept| ~1e15
predicts one label for every sample, scores exactly chance, and used to be
pickled and shipped. These drive the measurements that catch it, plus the two
things that must *not* be mistaken for it -- a rare class going unpredicted,
and a cheap baseline beating the model.
"""

import numpy as np
import pandas as pd
import pytest

from genotools.ancestry.fit_validation import (
    BASELINE_N_NEIGHBORS,
    MAX_HEALTHY_COEFFICIENT,
    MIN_FIT_CHANCE_MULTIPLE,
    WARN_UNPREDICTED_SUPPORT,
    baseline_scores,
    coefficient_health,
    fit_accuracy_floor,
    fit_validation_warnings,
    pc_columns,
    validate_fit,
)


class _Estimator:
    """Minimal stand-in exposing what `coefficient_health` reads."""

    def __init__(self, coef=None, intercept=None):
        if coef is not None:
            self.coef_ = np.asarray(coef, dtype=float)
        if intercept is not None:
            self.intercept_ = np.asarray(intercept, dtype=float)


class _Pipeline:
    """A stand-in whose final step is what carries the coefficients."""

    def __init__(self, final):
        self.named_steps = {"umap": object(), "xgb": final}


def _healthy() -> _Estimator:
    """Roughly what the working GP2 parity model measures."""
    return _Estimator(coef=[[1.03, -0.4], [0.2, 0.9]], intercept=[4.24, -3.1])


def _diverged() -> _Estimator:
    """Roughly what the collapsed PPMI model measures."""
    return _Estimator(coef=[[1727.0, -3.0]], intercept=[3.0e15])


class TestFitAccuracyFloor:
    """The bar is relative to chance, because the vocabulary is user-supplied."""

    def test_ten_classes(self) -> None:
        assert fit_accuracy_floor(10) == (0.30, "derived")

    def test_twenty_five_classes(self) -> None:
        """A fixed 0.70 would wrongly kill a usable model on a big panel."""
        floor, source = fit_accuracy_floor(25)
        assert floor == pytest.approx(0.12)
        assert source == "derived"

    def test_the_multiple_is_what_scales_it(self) -> None:
        for n in (4, 10, 25):
            assert fit_accuracy_floor(n)[0] == pytest.approx(
                MIN_FIT_CHANCE_MULTIPLE / n
            )

    def test_override_wins(self) -> None:
        assert fit_accuracy_floor(10, 0.8) == (0.8, "override")

    def test_zero_override_wins_over_the_derived_floor(self) -> None:
        """0 must disable the raise rather than fall through to the default."""
        assert fit_accuracy_floor(10, 0.0) == (0.0, "override")

    def test_few_classes_cannot_demand_more_than_one(self) -> None:
        assert fit_accuracy_floor(2)[0] == 1.0

    def test_degenerate_class_count(self) -> None:
        assert fit_accuracy_floor(1) == (0.0, "derived")


class TestCoefficientHealth:
    """Numerical health is checked at its source, not through its symptom."""

    def test_healthy_model(self) -> None:
        coef, intercept = coefficient_health(_healthy())
        assert coef == pytest.approx(1.03)
        assert intercept == pytest.approx(4.24)

    def test_diverged_model(self) -> None:
        coef, intercept = coefficient_health(_diverged())
        assert coef == pytest.approx(1727.0)
        assert intercept == pytest.approx(3.0e15)

    def test_reads_through_a_pipeline(self) -> None:
        """The production estimator is a Pipeline ending in the booster."""
        assert coefficient_health(_Pipeline(_healthy()))[0] == pytest.approx(1.03)

    def test_absent_coefficients_are_not_health(self) -> None:
        """A tree booster exposes none; None must not read as converged."""
        assert coefficient_health(_Estimator()) == (None, None)


class TestValidateFit:
    """The three checks, and the two things that are not failures."""

    def test_a_healthy_fit_passes(self) -> None:
        y = np.repeat(np.arange(10), 10)
        result = validate_fit(_healthy(), y, y, n_classes=10)
        assert not result.diverged
        assert not result.collapsed
        assert result.n_classes_predicted == 10
        assert result.train_balanced_accuracy == pytest.approx(1.0)

    def test_divergence_is_caught_even_when_predictions_look_fine(self) -> None:
        """The numerical bound catches the cause, not only the symptom."""
        y = np.repeat(np.arange(10), 10)
        result = validate_fit(_diverged(), y, y, n_classes=10)
        assert result.diverged
        assert result.collapsed

    def test_the_bound_is_where_it_says_it_is(self) -> None:
        y = np.repeat(np.arange(4), 5)
        just_under = _Estimator(
            coef=[[MAX_HEALTHY_COEFFICIENT]], intercept=[0.0]
        )
        just_over = _Estimator(
            coef=[[MAX_HEALTHY_COEFFICIENT * 1.01]], intercept=[0.0]
        )
        assert not validate_fit(just_under, y, y, n_classes=4).diverged
        assert validate_fit(just_over, y, y, n_classes=4).diverged

    def test_a_constant_classifier_is_collapsed(self) -> None:
        """The PPMI signature: one label for everything, chance accuracy."""
        y = np.repeat(np.arange(10), 10)
        pred = np.full_like(y, 7)
        result = validate_fit(_healthy(), y, pred, n_classes=10)
        assert result.n_classes_predicted == 1
        assert result.collapsed
        assert result.train_balanced_accuracy == pytest.approx(0.1)
        assert result.chance == pytest.approx(0.1)

    def test_a_dropped_rare_class_is_not_collapsed(self) -> None:
        """AAC and FIN are legitimate zero-prediction candidates.

        Both sit inside another label's cluster and hold under 3% of the
        panel; refusing a model over them would refuse every GP2 model.
        """
        y = np.concatenate([np.repeat(np.arange(9), 100), np.full(20, 9)])
        pred = y.copy()
        pred[pred == 9] = 0
        result = validate_fit(_healthy(), y, pred, n_classes=10)
        assert result.n_classes_predicted == 9
        assert not result.collapsed
        assert result.unpredicted_labels == ()

    def test_a_well_supported_unpredicted_class_warns(self) -> None:
        """Above the support threshold it is worth saying, still not a raise."""
        y = np.repeat(np.arange(10), 100)
        pred = y.copy()
        pred[pred == 9] = 0
        result = validate_fit(
            _healthy(), y, pred, n_classes=10, labels=list("ABCDEFGHIJ")
        )
        assert result.unpredicted_labels == ("J",)
        assert not result.collapsed
        assert any("never predicts J" in w for w in fit_validation_warnings(result))

    def test_support_threshold_boundary(self) -> None:
        """The threshold is the constant, not a number baked into a branch."""
        n_big = 1000
        rare = int(n_big * WARN_UNPREDICTED_SUPPORT * 0.5)
        y = np.concatenate([np.full(n_big, 0), np.full(rare, 1)])
        pred = np.zeros_like(y)
        assert validate_fit(_healthy(), y, pred, n_classes=2).unpredicted_labels == ()

    def test_a_zero_floor_records_but_refuses_nothing(self) -> None:
        """--ancestry-min-fit-accuracy 0 keeps the measurements, drops the gate."""
        y = np.repeat(np.arange(10), 10)
        pred = np.full_like(y, 3)
        result = validate_fit(
            _healthy(), y, pred, n_classes=10, min_balanced_accuracy=0.0
        )
        assert result.floor == 0.0
        assert result.floor_source == "override"
        # Still collapsed: a single-label model is refused whatever the floor.
        assert result.collapsed
        assert result.train_balanced_accuracy == pytest.approx(0.1)

    def test_below_the_floor_but_predicting_many_labels(self) -> None:
        """Weak-but-not-constant is caught by the accuracy check alone."""
        rng = np.random.default_rng(0)
        y = np.repeat(np.arange(10), 40)
        pred = rng.permutation(y)
        result = validate_fit(_healthy(), y, pred, n_classes=10)
        assert result.n_classes_predicted == 10
        assert not result.diverged
        assert result.collapsed

    def test_n_classes_is_taken_from_the_caller_not_the_labels(self) -> None:
        """A label with no training rows is still a class the model claims."""
        y = np.repeat(np.arange(3), 10)
        result = validate_fit(_healthy(), y, y, n_classes=10)
        assert result.n_classes_expected == 10
        assert result.floor == pytest.approx(0.3)

    def test_to_dict_is_json_ready(self) -> None:
        import json

        y = np.repeat(np.arange(4), 5)
        payload = validate_fit(_healthy(), y, y, n_classes=4).to_dict()
        assert json.loads(json.dumps(payload))["n_classes_predicted"] == 4

    def test_to_dict_carries_no_infinities(self) -> None:
        """JSON has no inf; a diverged fit must still serialize."""
        y = np.repeat(np.arange(4), 5)
        wild = _Estimator(coef=[[np.inf]], intercept=[np.nan])
        payload = validate_fit(wild, y, y, n_classes=4).to_dict()
        assert payload["max_abs_coefficient"] is None
        assert payload["max_abs_intercept"] is None

    def test_format_summary_names_the_numbers(self) -> None:
        y = np.repeat(np.arange(10), 10)
        text = validate_fit(_diverged(), y, np.full_like(y, 1), 10).format_summary()
        assert "diverged" in text
        assert "1/10" in text


class TestFitValidationWarnings:
    """Worst first, and each says what the reader should do about it."""

    def test_divergence_leads(self) -> None:
        y = np.repeat(np.arange(10), 10)
        result = validate_fit(_diverged(), y, np.full_like(y, 1), 10)
        warnings = fit_validation_warnings(result)
        assert "diverged" in warnings[0]
        assert "saturates" in warnings[0]

    def test_a_healthy_fit_warns_about_nothing(self) -> None:
        y = np.repeat(np.arange(10), 10)
        assert fit_validation_warnings(validate_fit(_healthy(), y, y, 10)) == []

    def test_a_cv_to_test_gap_is_called_out(self) -> None:
        y = np.repeat(np.arange(10), 10)
        result = validate_fit(_healthy(), y, y, 10)
        warnings = fit_validation_warnings(result, cv_score=0.95, test_score=0.15)
        assert any("did not carry over" in w for w in warnings)

    def test_a_winning_baseline_is_a_warning_and_says_it_is_not_a_failure(
        self,
    ) -> None:
        """On the real panel k-NN scores 0.956 against a healthy 0.953."""
        y = np.repeat(np.arange(10), 10)
        result = validate_fit(_healthy(), y, y, 10)
        warnings = fit_validation_warnings(
            result, test_score=0.953, baseline={"15-NN": 0.956}
        )
        assert any("not a failure" in w for w in warnings)

    def test_a_losing_baseline_is_silent(self) -> None:
        y = np.repeat(np.arange(10), 10)
        result = validate_fit(_healthy(), y, y, 10)
        assert (
            fit_validation_warnings(
                result, test_score=0.98, baseline={"15-NN": 0.90}
            )
            == []
        )


class TestPcColumns:
    """Ordered numerically, because a capped selection must take PC1..PCn."""

    def test_numeric_order_beats_lexical(self) -> None:
        frame = pd.DataFrame(columns=["PC10", "PC2", "PC1", "label"])
        assert pc_columns(frame) == ["PC1", "PC2", "PC10"]

    def test_cap_takes_the_first_n_in_numeric_order(self) -> None:
        frame = pd.DataFrame(columns=["PC10", "PC2", "PC1"])
        assert pc_columns(frame, 2) == ["PC1", "PC2"]

    def test_non_pc_columns_are_left_out(self) -> None:
        frame = pd.DataFrame(columns=["FID", "IID", "PC1", "label"])
        assert pc_columns(frame) == ["PC1"]


def _labeled_pcs(n_per: int = 30, separation: float = 8.0) -> pd.DataFrame:
    """Three well-separated clusters in three PCs, the shape of ref PCs."""
    rng = np.random.default_rng(1)
    names = ["AAA", "BBB", "CCC"]
    centres = np.eye(3) * separation
    rows, labels = [], []
    for index, name in enumerate(names):
        rows.append(rng.normal(size=(n_per, 3)) + centres[index])
        labels += [name] * n_per
    frame = pd.DataFrame(np.vstack(rows), columns=["PC1", "PC2", "PC3"])
    frame["label"] = labels
    return frame


class TestBaselineScores:
    """A second opinion sharing none of the pipeline's failure modes."""

    def test_both_arms_score_separable_clusters(self) -> None:
        frame = _labeled_pcs()
        scores = baseline_scores(frame.iloc[::2], frame.iloc[1::2])
        assert set(scores) == {f"{BASELINE_N_NEIGHBORS}-NN", "nearest-centroid"}
        assert all(score > 0.9 for score in scores.values())

    def test_overlapping_clusters_score_near_chance(self) -> None:
        frame = _labeled_pcs(separation=0.0)
        scores = baseline_scores(frame.iloc[::2], frame.iloc[1::2])
        assert all(score < 0.6 for score in scores.values())

    def test_k_is_capped_by_the_training_rows(self) -> None:
        """Fewer samples than neighbours must not raise."""
        frame = _labeled_pcs(n_per=2)
        scores = baseline_scores(frame.iloc[::2], frame.iloc[1::2])
        assert any(key.endswith("-NN") for key in scores)

    def test_a_frame_without_pcs_scores_nothing(self) -> None:
        frame = pd.DataFrame({"FID": ["a"], "IID": ["a"], "label": ["AAA"]})
        assert baseline_scores(frame, frame) == {}

    def test_a_frame_without_labels_scores_nothing(self) -> None:
        frame = _labeled_pcs().drop(columns=["label"])
        assert baseline_scores(frame, frame) == {}

    def test_the_pc_cap_is_honoured(self) -> None:
        """Separation lives only in PC3, so capping at 2 must lose it."""
        rng = np.random.default_rng(2)
        frame = pd.DataFrame(rng.normal(size=(90, 3)), columns=["PC1", "PC2", "PC3"])
        frame["label"] = np.repeat(["AAA", "BBB", "CCC"], 30)
        frame["PC3"] = np.repeat([0.0, 20.0, 40.0], 30)
        capped = baseline_scores(frame.iloc[::2], frame.iloc[1::2], n_pcs=2)
        full = baseline_scores(frame.iloc[::2], frame.iloc[1::2], n_pcs=3)
        assert max(full.values()) > max(capped.values())
