"""Prediction-path diagnostics at the model layer.

The projection and the classifier are stubbed so the PC coordinates are inputs
to the test rather than an outcome of a fit: what is under test is what the CAH
override does with a position, what it records about the decision, and what the
self-test makes of the reference panel. Fitting a real model would take minutes
and would measure UMAP, not this.
"""

import logging

import numpy as np
import pandas as pd
import pytest
from sklearn.preprocessing import LabelEncoder

from genotools.ancestry.model import AncestryModel

# Three well-separated reference groups. Their global centroid sits near
# (0, 10/3, 0) -- nowhere near any of them, which is what makes CAH reachable.
CENTROIDS = {"EUR": (10.0, 0.0, 0.0), "AFR": (-10.0, 0.0, 0.0), "EAS": (0.0, 10.0, 0.0)}
GLOBAL_CENTROID = (0.0, 10.0 / 3.0, 0.0)


class _StubPCAReducer:
    """Treats the first three feature columns as PC1-PC3."""

    def transform(self, X, return_dataframe=False):
        pcs = pd.DataFrame(np.asarray(X, dtype=float)[:, :3], columns=["PC1", "PC2", "PC3"])
        return pcs if return_dataframe else pcs.values


class _StubPipeline:
    """Returns encoded labels chosen by the test."""

    def __init__(self, encoded):
        self.encoded = np.asarray(encoded)

    def predict(self, X):
        return self.encoded[: len(X)]


def _train_pca(per_group=40, spread=0.2, seed=0):
    rng = np.random.default_rng(seed)
    frames = []
    for label, centre in CENTROIDS.items():
        offsets = rng.normal(0, spread, (per_group, 3))
        frame = pd.DataFrame(offsets + np.array(centre), columns=["PC1", "PC2", "PC3"])
        frame.insert(0, "IID", [f"{label}_{i}" for i in range(per_group)])
        frame.insert(0, "FID", [f"{label}_{i}" for i in range(per_group)])
        frame["label"] = label
        frames.append(frame)
    return pd.concat(frames, ignore_index=True)


def _model(predicted_labels, train_pca=None):
    encoder = LabelEncoder().fit(sorted(CENTROIDS))
    model = AncestryModel()
    model.pca_reducer = _StubPCAReducer()
    model.pipeline = _StubPipeline(encoder.transform(predicted_labels))
    model.label_encoder = encoder
    model._train_pca = train_pca if train_pca is not None else _train_pca()
    model._is_fitted = True
    return model


def _cohort(positions):
    frame = pd.DataFrame(positions, columns=["snp1", "snp2", "snp3"])
    frame.insert(0, "IID", [f"S{i}" for i in range(len(frame))])
    frame.insert(0, "FID", [f"S{i}" for i in range(len(frame))])
    return frame


def _predict(model, cohort, **kwargs):
    return model.predict(cohort, cohort[["FID", "IID"]], **kwargs)


# --- the override, and its working --------------------------------------------


def test_a_cohort_sitting_on_an_ancestry_centroid_keeps_its_label():
    cohort = _cohort([CENTROIDS["EUR"]] * 5)
    result = _predict(_model(["EUR"] * 5), cohort, collect_diagnostics=False)
    assert list(result.predictions["predicted_ancestry"]) == ["EUR"] * 5
    assert (result.decisions["margin_to_all"] > 0).all()
    assert (result.decisions["label_pre_admixture"] == "EUR").all()
    assert (result.decisions["nearest_centroid"] == "EUR").all()


def test_a_cohort_at_the_global_centroid_becomes_cah_but_keeps_the_classifier_label():
    """The all-CAH failure, and the evidence that separates its two causes."""
    cohort = _cohort([GLOBAL_CENTROID] * 6)
    result = _predict(
        _model(["EUR", "AFR", "EAS", "EUR", "AFR", "EAS"]), cohort,
        collect_diagnostics=False,
    )
    assert set(result.predictions["predicted_ancestry"]) == {"CAH"}
    assert set(result.decisions["label_pre_admixture"]) == {"EUR", "AFR", "EAS"}
    assert (result.decisions["margin_to_all"] < 0).all()
    assert (result.decisions["nearest_centroid"] == "ALL").all()


def test_the_decision_frame_carries_a_distance_per_centroid():
    cohort = _cohort([CENTROIDS["AFR"]] * 2)
    result = _predict(_model(["AFR"] * 2), cohort, collect_diagnostics=False)
    decisions = result.decisions
    for centroid in ["ALL", *CENTROIDS]:
        assert f"dist_{centroid}" in decisions.columns
    assert decisions.loc[0, "dist_AFR"] < decisions.loc[0, "dist_EUR"]
    assert decisions.loc[0, "nearest_ancestry"] == "AFR"
    assert decisions.loc[0, "dist_to_nearest_ancestry"] == pytest.approx(
        decisions.loc[0, "dist_AFR"]
    )


def test_margin_is_exactly_the_rule_the_override_applies():
    """`margin_to_all < 0` and `label == CAH` must never disagree."""
    positions = [CENTROIDS["EUR"], GLOBAL_CENTROID, CENTROIDS["EAS"], GLOBAL_CENTROID]
    cohort = _cohort(positions)
    result = _predict(_model(["EUR", "EUR", "EAS", "EAS"]), cohort, collect_diagnostics=False)
    decisions = result.decisions
    assert list(decisions["margin_to_all"] < 0) == list(decisions["label"] == "CAH")


def test_detection_can_be_turned_off_to_see_the_classifier():
    cohort = _cohort([GLOBAL_CENTROID] * 3)
    result = _predict(
        _model(["EUR", "AFR", "EAS"]), cohort,
        detect_admixed=False, collect_diagnostics=False,
    )
    assert list(result.predictions["predicted_ancestry"]) == ["EUR", "AFR", "EAS"]
    assert result.decisions is None


# --- what predict() reports ---------------------------------------------------


def test_predict_reports_drift_the_override_and_writes_the_decisions(tmp_path, caplog):
    cohort = _cohort([GLOBAL_CENTROID] * 4)
    with caplog.at_level(logging.INFO, logger="genotools"):
        result = _predict(
            _model(["EUR", "AFR", "EAS", "EUR"]), cohort,
            diagnostics_prefix=str(tmp_path / "cohort"),
        )

    diagnostics = result.diagnostics
    assert list(diagnostics.pc_drift["pc"]) == ["PC1", "PC2", "PC3"]
    assert diagnostics.admixture["fraction_cah"] == 1.0
    assert diagnostics.admixture["n_distinct_pre_admixture"] == 3
    assert any("labeled CAH" in w for w in diagnostics.warnings)

    written = pd.read_csv(tmp_path / "cohort_decisions.txt", sep="\t")
    assert len(written) == 4
    assert {"label_pre_admixture", "label", "margin_to_all"} <= set(written.columns)


def test_a_collapsed_cohort_is_reported_as_such(tmp_path):
    """Every sample at one point: the constant-fill signature, in PC space."""
    cohort = _cohort([GLOBAL_CENTROID] * 20)
    result = _predict(_model(["EUR"] * 20), cohort, diagnostics_prefix=str(tmp_path / "c"))
    assert any("collapsed to a point" in w for w in result.diagnostics.warnings)


def test_diagnostics_can_be_declined_entirely():
    cohort = _cohort([CENTROIDS["EUR"]] * 2)
    result = _predict(_model(["EUR"] * 2), cohort, collect_diagnostics=False)
    assert result.diagnostics is None


# --- the reference-panel self-test -------------------------------------------


def _labeled_ref(train_pca):
    """The panel as `get_raw_files` hands it over: IDs, features, label."""
    ref = train_pca[["FID", "IID"]].copy()
    ref[["snp1", "snp2", "snp3"]] = train_pca[["PC1", "PC2", "PC3"]].values
    ref["label"] = train_pca["label"].values
    return ref


def test_self_test_recovers_a_working_model():
    train_pca = _train_pca()
    ref = _labeled_ref(train_pca)
    model = _model(list(train_pca["label"]), train_pca=train_pca)

    result = model.self_test(ref)
    assert result["n_samples"] == len(ref)
    assert result["balanced_accuracy"] == pytest.approx(1.0)
    assert result["cah_fraction"] == 0.0
    assert result["top_confusions"] == {}


def test_self_test_indicts_the_override_when_the_panel_comes_back_cah(caplog):
    """The panel defines the centroids, so a CAH panel is not the cohort's fault."""
    train_pca = _train_pca()
    ref = _labeled_ref(train_pca)
    # Every reference sample projected onto the global centroid.
    ref[["snp1", "snp2", "snp3"]] = np.tile(GLOBAL_CENTROID, (len(ref), 1))
    model = _model(list(train_pca["label"]), train_pca=train_pca)

    with caplog.at_level(logging.WARNING, logger="genotools"):
        result = model.self_test(ref)

    assert result["cah_fraction"] == 1.0
    assert result["balanced_accuracy"] == 0.0
    # The classifier was right all along; only the override moved.
    assert result["balanced_accuracy_pre_admixture"] == pytest.approx(1.0)
    assert "not the cohort" in caplog.text


def test_self_test_separates_a_broken_classifier_from_a_broken_override(caplog):
    train_pca = _train_pca()
    ref = _labeled_ref(train_pca)
    # Classifier says EUR for everything; positions are still correct, so the
    # override leaves the labels alone and the accuracy collapse is the model's.
    model = _model(["EUR"] * len(ref), train_pca=train_pca)

    with caplog.at_level(logging.WARNING, logger="genotools"):
        result = model.self_test(ref)

    assert result["cah_fraction"] == 0.0
    assert result["balanced_accuracy"] == pytest.approx(1 / 3)
    assert "not in the cohort" in caplog.text
    assert result["top_confusions"]


def test_self_test_ignores_unlabeled_reference_samples():
    train_pca = _train_pca()
    ref = _labeled_ref(train_pca)
    ref.loc[: len(ref) // 2, "label"] = np.nan
    model = _model(list(train_pca["label"]), train_pca=train_pca)

    result = model.self_test(ref)
    assert result["n_samples"] == int(ref["label"].notna().sum())


def test_self_test_on_an_unlabeled_panel_says_so_rather_than_dividing_by_zero():
    train_pca = _train_pca()
    ref = _labeled_ref(train_pca)
    ref["label"] = np.nan
    result = _model(list(train_pca["label"]), train_pca=train_pca).self_test(ref)
    assert result == {"n_samples": 0, "error": "no labeled reference samples"}


def test_self_test_reports_the_cah_rate_per_true_label():
    train_pca = _train_pca()
    ref = _labeled_ref(train_pca)
    eas = ref["label"] == "EAS"
    ref.loc[eas, ["snp1", "snp2", "snp3"]] = np.tile(GLOBAL_CENTROID, (int(eas.sum()), 1))
    model = _model(list(train_pca["label"]), train_pca=train_pca)

    by_label = model.self_test(ref)["cah_fraction_by_label"]
    assert by_label["EAS"] == 1.0
    assert by_label["EUR"] == 0.0
