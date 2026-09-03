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

"""Is the model that came out of training usable at all?

A PPMI/NAPU run labeled 636 of 644 cohort samples `SAS` and predicted `SAS`
for all 802 held-out reference samples. Its reported `test_accuracy` of 0.1496
is exactly SAS's prevalence in that split -- the classifier had become a
constant function. Nothing in the pipeline noticed, so it was pickled and used
to label the cohort with full confidence.

The cause was numerical: gblinear descending at XGBoost's own default learning
rate, on its nondeterministic Hogwild updater, diverged to coefficients around
1.7e3 and intercepts around 3.0e15, which saturates the softmax. Round 19
removed the race (see `ClassifierConfig`), but a numerical failure that depends
on the data cannot be argued away by settings alone, so this module measures
the fitted model directly and the training path refuses to ship one that fails.

Three checks, cheapest and most direct first:

- **numerical health** -- `|coef|` and `|intercept|` bounds. This is the check
  that cleanly separates the two real models: a healthy GP2 model measured
  `|coef| 1.03 / |intercept| 4.24`, the collapsed one `1727 / 3.0e15`.
- **distinct predictions** -- a model that emits one label is not a classifier.
- **balanced accuracy against chance** -- relative to `1/n_classes`, never a
  fixed number, because the label vocabulary is user-supplied.

All of it is measured on the **training** set, never the held-out split, so the
test score stays an honest estimate and a retry leaks nothing.

That makes the accuracy check *stricter* than it looks, not weaker. The
Pipeline passes `y` to `UMAP.fit_transform`, so the training embedding is
supervised -- it already encodes the labels it is about to be scored against.
A working model therefore scores near 1.0 here whatever the genotypes say, and
`train_balanced_accuracy` is not a measure of whether ancestry is learnable
from the data. What it does measure is whether the fit converged at all: a
diverged model scores exactly `1 / n_classes` even on an embedding that knows
the answer. Read it as a health check, never as a quality estimate -- that is
what the test split and the baselines below are for.

Pure functions over arrays and frames, per the `select_het_outliers` pattern,
so they can be tested against hand-built input. Kept apart from
`diagnostics.py`, which instruments the *prediction* path: that module asks
whether the data reaching a model is fit to predict on, this one asks whether
the model that came out of training is fit to predict with.
"""

from dataclasses import dataclass
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from genotools.core.logging import get_logger

logger = get_logger(__name__)


#: A fitted linear model whose coefficients or intercepts exceed this has
#: diverged rather than converged. Healthy fits measure |coef| ~1 and
#: |intercept| ~4; the collapsed model measured |coef| 1727 and |intercept|
#: 3.0e15. Three orders of magnitude above healthy leaves the bound nowhere
#: near a working model while still catching divergence at its source.
MAX_HEALTHY_COEFFICIENT = 1e3

#: A fit must beat chance by this multiple to be accepted. Expressed relative
#: to `1 / n_classes` rather than as a fixed threshold because the label
#: vocabulary is user-supplied: 0.70 suits the 10-class GP2 panel but would
#: wrongly kill a usable model on a 25-class one.
MIN_FIT_CHANCE_MULTIPLE = 3.0

#: Share of the training set a label must hold before never being predicted is
#: worth a warning. Two GP2 labels are legitimate zero-prediction candidates --
#: AAC (1.8% of the panel, admixed, sitting between AFR and EUR) and FIN (2.5%,
#: inside EUR) -- and the smallest label that should be flagged is MDE at 3.8%.
WARN_UNPREDICTED_SUPPORT = 0.03

#: Neighbours for the k-NN second opinion. Large enough that the score is not
#: one sample's noise, small enough to resolve a label holding 1.8% of a
#: 4,008-sample panel.
BASELINE_N_NEIGHBORS = 15


def _f(value: Any) -> Optional[float]:
    """A float the JSON report can hold, or None where there is no number."""
    if value is None:
        return None
    number = float(value)
    if not np.isfinite(number):
        return None
    return number


# ---------------------------------------------------------------------------
# 1. What counts as beating chance?
# ---------------------------------------------------------------------------


def fit_accuracy_floor(
    n_classes: int,
    override: Optional[float] = None,
) -> Tuple[float, str]:
    """The balanced accuracy a fit must reach, and where that number came from.

    Derived as `MIN_FIT_CHANCE_MULTIPLE / n_classes` so the bar scales with the
    panel: 0.30 for GP2's 10 labels, 0.12 for a 25-label panel. An explicit
    `override` wins, including 0, which records every measurement but refuses
    nothing.

    Args:
        n_classes: Labels the model was trained to distinguish.
        override: A caller-supplied floor in [0, 1], or None to derive one.

    Returns:
        `(floor, source)` where source is "override" or "derived".
    """
    if override is not None:
        return float(override), "override"
    if n_classes < 2:
        return 0.0, "derived"
    return min(1.0, MIN_FIT_CHANCE_MULTIPLE / n_classes), "derived"


# ---------------------------------------------------------------------------
# 2. Did the fit converge, and is it a classifier?
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class FitValidation:
    """The measured health of one fitted classifier.

    Attributes:
        train_balanced_accuracy: Balanced accuracy on the training set. A
            collapsed model scores ~1/n_classes on the data it was fitted to,
            which is what makes this measurable without touching the test
            split.
        n_classes_expected: Labels the model was trained on.
        n_classes_predicted: Distinct labels it actually emits. 1 is the
            signature of a saturated softmax.
        max_abs_coefficient: Largest `|coef|`, or None for a booster that
            exposes none (any tree booster).
        max_abs_intercept: Largest `|intercept|`, same caveat.
        unpredicted_labels: Labels above `WARN_UNPREDICTED_SUPPORT` that the
            model never predicts. A warning, never a rejection.
        diverged: Whether the numerical bounds were exceeded.
        collapsed: Whether the fit is unusable -- diverged, or emitting fewer
            than two labels, or below the floor.
        floor: The balanced accuracy required.
        floor_source: "override" or "derived".
    """

    train_balanced_accuracy: float
    n_classes_expected: int
    n_classes_predicted: int
    max_abs_coefficient: Optional[float]
    max_abs_intercept: Optional[float]
    unpredicted_labels: Tuple[str, ...]
    diverged: bool
    collapsed: bool
    floor: float
    floor_source: str

    @property
    def chance(self) -> float:
        """Balanced accuracy a constant classifier gets on this panel."""
        if self.n_classes_expected < 1:
            return 0.0
        return 1.0 / self.n_classes_expected

    def to_dict(self) -> Dict[str, Any]:
        """The measurements, for the JSON report."""
        return {
            "train_balanced_accuracy": _f(self.train_balanced_accuracy),
            "n_classes_expected": int(self.n_classes_expected),
            "n_classes_predicted": int(self.n_classes_predicted),
            "max_abs_coefficient": _f(self.max_abs_coefficient),
            "max_abs_intercept": _f(self.max_abs_intercept),
            "unpredicted_labels": list(self.unpredicted_labels),
            "diverged": bool(self.diverged),
            "collapsed": bool(self.collapsed),
            "min_balanced_accuracy": _f(self.floor),
            "min_balanced_accuracy_source": self.floor_source,
        }

    def format_summary(self) -> str:
        """One line for the log."""
        health = "diverged" if self.diverged else "converged"
        coef = (
            "unmeasurable"
            if self.max_abs_coefficient is None
            else f"|coef| {self.max_abs_coefficient:.3g} "
            f"|intercept| {self.max_abs_intercept:.3g}"
        )
        return (
            f"fit {health}: {coef}; "
            f"predicts {self.n_classes_predicted}/{self.n_classes_expected} "
            f"labels; train balanced accuracy "
            f"{self.train_balanced_accuracy:.4f} against a floor of "
            f"{self.floor:.4f} ({self.floor_source}, chance "
            f"{self.chance:.4f})"
        )


def coefficient_health(
    estimator: Any,
) -> Tuple[Optional[float], Optional[float]]:
    """Largest `|coef|` and `|intercept|` on a fitted estimator, if it has them.

    Returns `(None, None)` for an estimator exposing neither -- a tree booster,
    or a Pipeline whose final step is one. Absence is not health, so callers
    must not read it as convergence.
    """
    if hasattr(estimator, "named_steps"):
        estimator = list(estimator.named_steps.values())[-1]

    values: List[Optional[float]] = []
    for name in ("coef_", "intercept_"):
        raw = getattr(estimator, name, None)
        if raw is None:
            values.append(None)
            continue
        array = np.abs(np.asarray(raw, dtype=float))
        values.append(float(array.max()) if array.size else None)
    return values[0], values[1]


def validate_fit(
    estimator: Any,
    y_true: Sequence[Any],
    y_pred: Sequence[Any],
    n_classes: int,
    min_balanced_accuracy: Optional[float] = None,
    labels: Optional[Sequence[str]] = None,
) -> FitValidation:
    """Measure a fitted classifier against the training data it was fitted to.

    Args:
        estimator: The fitted estimator, or a Pipeline ending in one.
        y_true: Training labels, encoded or not.
        y_pred: The estimator's predictions on the same rows.
        n_classes: Labels the model was trained to distinguish. Passed rather
            than inferred from `y_true`, because a label with no training rows
            is still a class the model claims to predict.
        min_balanced_accuracy: Floor override in [0, 1]; 0 refuses nothing.
        labels: Class names in encoder order, for naming unpredicted labels.
            Falls back to the encoded values.

    Returns:
        A `FitValidation`. Nothing is raised here -- deciding what a failed
        fit means belongs to the caller.
    """
    from sklearn.metrics import balanced_accuracy_score

    true = np.asarray(y_true)
    pred = np.asarray(y_pred)
    balanced = float(balanced_accuracy_score(true, pred))

    max_coef, max_intercept = coefficient_health(estimator)
    diverged = any(
        value is not None and value > MAX_HEALTHY_COEFFICIENT
        for value in (max_coef, max_intercept)
    )

    predicted = set(np.unique(pred).tolist())
    floor, floor_source = fit_accuracy_floor(n_classes, min_balanced_accuracy)

    collapsed = bool(
        diverged or len(predicted) < 2 or balanced < floor
    )

    return FitValidation(
        train_balanced_accuracy=balanced,
        n_classes_expected=int(n_classes),
        n_classes_predicted=len(predicted),
        max_abs_coefficient=max_coef,
        max_abs_intercept=max_intercept,
        unpredicted_labels=_unpredicted_labels(true, predicted, labels),
        diverged=diverged,
        collapsed=collapsed,
        floor=floor,
        floor_source=floor_source,
    )


def _unpredicted_labels(
    y_true: np.ndarray,  # type: ignore[type-arg]
    predicted: set,
    labels: Optional[Sequence[str]],
) -> Tuple[str, ...]:
    """Well-supported training labels the model never emits.

    Thin support is a legitimate reason for a label to go unpredicted, so only
    labels above `WARN_UNPREDICTED_SUPPORT` are reported.
    """
    if y_true.size == 0:
        return ()
    values, counts = np.unique(y_true, return_counts=True)
    missing = []
    for value, count in zip(values, counts):
        if value in predicted:
            continue
        if count / y_true.size < WARN_UNPREDICTED_SUPPORT:
            continue
        name = value
        if labels is not None and isinstance(value, (int, np.integer)):
            if 0 <= int(value) < len(labels):
                name = labels[int(value)]
        missing.append(str(name))
    return tuple(missing)


def fit_validation_warnings(
    validation: FitValidation,
    cv_score: Optional[float] = None,
    test_score: Optional[float] = None,
    baseline: Optional[Dict[str, float]] = None,
) -> List[str]:
    """Human-readable concerns about a fit, worst first."""
    warnings: List[str] = []

    if validation.diverged:
        warnings.append(
            f"the fit diverged rather than converged: "
            f"|coef| {validation.max_abs_coefficient:.3g}, "
            f"|intercept| {validation.max_abs_intercept:.3g}, against "
            f"{MAX_HEALTHY_COEFFICIENT:.0e}. A diverged linear model saturates "
            f"the softmax and predicts one label for every sample."
        )
    if validation.n_classes_predicted < 2:
        warnings.append(
            f"the model predicts a single label for every training sample, so "
            f"it is a constant function and not a classifier "
            f"({validation.n_classes_expected} labels were expected)."
        )
    elif validation.train_balanced_accuracy < validation.floor:
        warnings.append(
            f"training balanced accuracy {validation.train_balanced_accuracy:.4f} "
            f"is below the floor of {validation.floor:.4f} "
            f"({validation.floor_source}); chance on "
            f"{validation.n_classes_expected} labels is "
            f"{validation.chance:.4f}."
        )

    if validation.unpredicted_labels:
        warnings.append(
            f"the model never predicts "
            f"{', '.join(validation.unpredicted_labels)}, each holding at "
            f"least {100 * WARN_UNPREDICTED_SUPPORT:.0f}% of the training "
            f"set. Samples of those ancestries will be given some other "
            f"label."
        )

    if cv_score is not None and test_score is not None:
        gap = float(cv_score) - float(test_score)
        if gap > 0.1:
            warnings.append(
                f"cross-validated balanced accuracy {float(cv_score):.4f} is "
                f"{gap:.4f} above the held-out score {float(test_score):.4f}, "
                f"so the search's estimate did not carry over to the test "
                f"split."
            )

    if baseline and test_score is not None:
        best_name, best_score = max(baseline.items(), key=lambda kv: kv[1])
        if best_score > float(test_score):
            warnings.append(
                f"a {best_name} baseline on the raw reference PCs scores "
                f"{best_score:.4f} against the trained model's "
                f"{float(test_score):.4f}. The baseline can legitimately win, "
                f"so this is not a failure -- but a large gap means the "
                f"UMAP+booster pipeline is not earning its cost."
            )

    return warnings


# ---------------------------------------------------------------------------
# 3. A second opinion that shares none of the model's failure modes
# ---------------------------------------------------------------------------


def pc_columns(
    frame: pd.DataFrame,
    n_pcs: Optional[int] = None,
) -> List[str]:
    """PC column names in numeric order (PC2 before PC10), capped at `n_pcs`.

    Selected by name rather than by dropping known ID columns: a frame carrying
    an extra diagnostic column must not have it silently treated as a PC. And
    ordered numerically, because `startswith("PC")` alone puts PC10 before PC2
    and a capped selection then takes the wrong PCs.
    """
    cols = [c for c in frame.columns if str(c).startswith("PC")]

    def _index(name: str) -> int:
        try:
            return int(str(name)[2:])
        except ValueError:
            return 1 << 30

    cols.sort(key=_index)
    return cols[:n_pcs] if n_pcs else cols


def baseline_scores(
    train_pca: pd.DataFrame,
    test_pca: pd.DataFrame,
    n_pcs: Optional[int] = None,
    n_neighbors: int = BASELINE_N_NEIGHBORS,
) -> Dict[str, float]:
    """Balanced accuracy of two cheap classifiers on the raw reference PCs.

    Milliseconds on a 4,008 x 50 panel, and neither shares the trained
    pipeline's failure modes -- no UMAP, no gradient descent, nothing to
    diverge. On the collapsed PPMI run these would have read `0.956` beside
    the model's `0.150` in the report itself.

    A warning, never a gate: on that same panel the k-NN baseline scores 0.956
    against 0.953 for a *healthy* trained model, so the baseline winning is
    not evidence of a broken fit.

    Args:
        train_pca: Labeled reference PCs to fit on (`PC*` plus `label`).
        test_pca: Labeled reference PCs to score on, same columns.
        n_pcs: Cap on PCs used; None uses all present.
        n_neighbors: Neighbours for the k-NN arm.

    Returns:
        `{"15-NN": score, "nearest-centroid": score}`, keyed by the k actually
        used. Empty if either frame has no PCs or no labels.
    """
    from sklearn.metrics import balanced_accuracy_score
    from sklearn.neighbors import KNeighborsClassifier, NearestCentroid

    cols = pc_columns(train_pca, n_pcs)
    if not cols or "label" not in train_pca or "label" not in test_pca:
        return {}
    if train_pca.empty or test_pca.empty:
        return {}

    X_train = train_pca[cols].to_numpy(dtype=float)
    y_train = train_pca["label"].to_numpy()
    X_test = test_pca[cols].to_numpy(dtype=float)
    y_test = test_pca["label"].to_numpy()

    k = int(min(n_neighbors, len(y_train)))
    scores: Dict[str, float] = {}
    for name, estimator in (
        (f"{k}-NN", KNeighborsClassifier(n_neighbors=k)),
        ("nearest-centroid", NearestCentroid()),
    ):
        try:
            estimator.fit(X_train, y_train)
            scores[name] = float(
                balanced_accuracy_score(y_test, estimator.predict(X_test))
            )
        except Exception as error:  # pragma: no cover - defensive
            logger.debug(f"{name} baseline could not be scored: {error}")
    return scores
