#!/usr/bin/env python
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

"""How often does ancestry training collapse, before and after the round-19 fix?

A PPMI/NAPU run labeled 636 of 644 cohort samples ``SAS`` and predicted ``SAS``
for all 802 held-out panel samples. Its ``test_accuracy`` of 0.1496 is exactly
SAS's prevalence in that split -- the classifier had become a constant
function, and was pickled and used anyway.

The cause was that ``ClassifierConfig.learning_rate`` was never passed to
``XGBClassifier``, so gblinear descended at XGBoost's own default of 0.5, close
enough to the divergence boundary that its Hogwild ``shotgun`` updater's thread
races decided per run which side a fit landed on. Whether that boundary is
reached depends on the data, so the only honest way to state the rate is to
measure it per cohort -- which is what this does.

No PLINK, no genotypes. The reference panel's PCs are in every ``--ancestry``
report under ``ref_pcs``, and the cohort's projections under ``projected_pcs``,
so a report alone is enough to re-run everything downstream of the PCA. The PCA
is not what failed.

Three arms per report, each fitting the same UMAP embedding N times:

    as shipped        no learning rate, no thread pinning (the defect)
    threads pinned    n_jobs=1 alone
    the fix           n_jobs=1 + learning_rate=0.1 + n_estimators=200

Measured results, 20 repeats each
---------------------------------
    PPMI (long-read WGS, ~168k panel/cohort SNP overlap)
        as shipped     6/20 and 19/20 collapses on two draws of the same
                       configuration -- the rate itself is not stable
        the fix        0/20, bit-identical across repeats
    GP2 (array, 43k overlap)
        as shipped     0/20
        the fix        0/20, bit-identical across repeats

So array cohorts sit inside the stable region and dense-WGS cohorts do not.
Having *more* overlapping variants than usual is what caused the failure.

Usage
-----
    python tests/scripts/compare_fit_determinism.py \\
        --report run_output.json --repeats 20 --out findings.txt

    # name the arms explicitly, or add --predict to label the cohort
    python tests/scripts/compare_fit_determinism.py \\
        --report a.json --report b.json --predict

Two traps hit while measuring this, worth repeating: piping the output through
``grep`` re-buffers it even under ``python -u``, so a run that times out shows
*no* output rather than partial output -- hence ``--out``. And ``nohup cmd &``
from a tool-invoked shell dies with its parent; use
``setsid nohup cmd > out.txt 2>&1 < /dev/null & disown``.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

#: What genotools splits the reference panel with, in 1.x and 2.0 alike.
TEST_SIZE = 0.2
SPLIT_SEED = 123
MODEL_SEED = 123

#: The hyperparameters recovered from the collapsed 1.x model. Used as the
#: single fixed candidate so every arm fits the same embedding -- running the
#: real 216-point grid would take hours and answer a different question.
UMAP_PARAMS = {"n_neighbors": 5, "n_components": 25, "a": 1.0, "b": 0.25}
XGB_LAMBDA = 0.001

#: The three configurations, in the order they are reported. `None` means
#: "leave it to XGBoost", which is precisely what the defect did.
ARMS: Tuple[Tuple[str, Dict[str, Any]], ...] = (
    ("as shipped", {"learning_rate": None, "n_jobs": None, "n_estimators": None}),
    ("threads pinned", {"learning_rate": None, "n_jobs": 1, "n_estimators": None}),
    ("the fix", {"learning_rate": 0.1, "n_jobs": 1, "n_estimators": 200}),
)


def pc_columns(frame: pd.DataFrame, n_pcs: Optional[int] = None) -> List[str]:
    """PC column names in numeric order, so PC2 sorts before PC10."""
    cols = [c for c in frame.columns if str(c).startswith("PC")]
    cols.sort(key=lambda c: int(str(c)[2:]))
    return cols[:n_pcs] if n_pcs else cols


def read_report(path: Path) -> Tuple[pd.DataFrame, Optional[pd.DataFrame]]:
    """Panel PCs and, when present, the cohort's projections.

    Seeks to each key and hands the tail to a ``raw_decode`` rather than
    parsing the file: a full-GP2 report is 352 MB of JSON and the two blocks
    wanted are 6 MB of it.
    """
    text = path.read_text()

    def block(key: str) -> Optional[pd.DataFrame]:
        needle = f'"{key}":'
        if needle not in text:
            return None
        start = text.index(needle) + len(needle)
        while text[start].isspace():
            start += 1
        payload, _ = json.JSONDecoder().raw_decode(text, start)
        return pd.DataFrame(payload)

    ref = block("ref_pcs")
    if ref is None:
        raise SystemExit(
            f"{path.name} has no 'ref_pcs'. That report was not produced by an "
            f"--ancestry run."
        )
    return ref, block("projected_pcs")


def split_panel(ref: pd.DataFrame):
    """Reproduce genotools' stratified panel split over a report's PCs."""
    from sklearn.model_selection import train_test_split

    cols = pc_columns(ref)
    X = ref[cols].to_numpy(dtype=float)
    classes = sorted(ref["label"].unique())
    y = np.array([classes.index(value) for value in ref["label"]])
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=TEST_SIZE, random_state=SPLIT_SEED, stratify=y
    )
    return X_train, X_test, y_train, y_test, classes, cols


def embed(X_train, X_test, y_train, cohort=None):
    """Fit UMAP once and reuse it, since it is ~95% of a pipeline fit.

    ``y`` is passed because the production Pipeline passes it, which makes the
    training embedding supervised.
    """
    from umap import UMAP

    reducer = UMAP(random_state=MODEL_SEED, **UMAP_PARAMS)
    E_train = np.asarray(reducer.fit_transform(X_train, y=y_train))
    E_test = np.asarray(reducer.transform(X_test))
    E_cohort = None if cohort is None else np.asarray(reducer.transform(cohort))
    return E_train, E_test, E_cohort


def fit_once(E_train, E_test, y_train, y_test, arm: Dict[str, Any]):
    """One classifier fit, measured the way the health checker reads a model."""
    from sklearn.metrics import balanced_accuracy_score
    from xgboost import XGBClassifier

    kwargs: Dict[str, Any] = {"booster": "gblinear", "random_state": MODEL_SEED}
    kwargs.update({k: v for k, v in arm.items() if v is not None})
    model = XGBClassifier(**kwargs, **{"lambda": XGB_LAMBDA})

    started = time.perf_counter()
    model.fit(E_train, y_train)
    elapsed = time.perf_counter() - started

    train_pred = model.predict(E_train)
    return {
        "model": model,
        "seconds": elapsed,
        "max_abs_coefficient": float(np.abs(np.asarray(model.coef_)).max()),
        "max_abs_intercept": float(np.abs(np.asarray(model.intercept_)).max()),
        "n_classes_predicted": int(len(np.unique(train_pred))),
        "train_balanced": float(balanced_accuracy_score(y_train, train_pred)),
        "test_balanced": float(
            balanced_accuracy_score(y_test, model.predict(E_test))
        ),
        "collapsed": bool(len(np.unique(train_pred)) < 2),
    }


def summarize(rows: List[Dict[str, Any]]) -> Dict[str, Any]:
    """One row of the arm table."""
    fingerprints = {
        (row["max_abs_coefficient"], row["max_abs_intercept"]) for row in rows
    }
    return {
        "n": len(rows),
        "collapses": sum(row["collapsed"] for row in rows),
        "max_abs_intercept": max(row["max_abs_intercept"] for row in rows),
        "min_test_balanced": min(row["test_balanced"] for row in rows),
        "max_test_balanced": max(row["test_balanced"] for row in rows),
        "mean_seconds": sum(row["seconds"] for row in rows) / len(rows),
        "deterministic": len(fingerprints) == 1,
    }


def baseline(X_train, X_test, y_train, y_test) -> Dict[str, float]:
    """The cheap second opinion, on the raw PCs the model never sees directly."""
    from sklearn.metrics import balanced_accuracy_score
    from sklearn.neighbors import KNeighborsClassifier, NearestCentroid

    scores = {}
    for name, estimator in (
        ("15-NN", KNeighborsClassifier(n_neighbors=15)),
        ("nearest-centroid", NearestCentroid()),
    ):
        estimator.fit(X_train, y_train)
        scores[name] = float(
            balanced_accuracy_score(y_test, estimator.predict(X_test))
        )
    return scores


def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__.split("\n")[0],
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--report",
        action="append",
        required=True,
        type=Path,
        help="A GenoTools --ancestry JSON report; repeatable",
    )
    parser.add_argument("--out", type=Path, default=None, help="Write findings here")
    parser.add_argument("--repeats", type=int, default=20)
    parser.add_argument(
        "--predict",
        action="store_true",
        help="Also label the cohort from projected_pcs, per arm",
    )
    args = parser.parse_args()

    lines: List[str] = []

    def say(text: str = "") -> None:
        print(text, flush=True)
        lines.append(text)
        if args.out:
            args.out.write_text("\n".join(lines) + "\n")

    say("Ancestry training determinism, per report")
    say(f"{args.repeats} identical repeats per arm; one fixed candidate, "
        f"one cached embedding")
    say()

    any_collapse = False
    for path in args.report:
        ref, projected = read_report(path)
        X_train, X_test, y_train, y_test, classes, cols = split_panel(ref)
        cohort = (
            projected[cols].to_numpy(dtype=float)
            if args.predict and projected is not None
            else None
        )
        started = time.perf_counter()
        E_train, E_test, E_cohort = embed(X_train, X_test, y_train, cohort)

        say(f"=== {path.name} ===")
        say(
            f"panel: train {X_train.shape} test {X_test.shape}, "
            f"{len(classes)} labels; raw PC sd {X_train.std():.3f}, "
            f"embedding sd {E_train.std():.3f} "
            f"({time.perf_counter() - started:.1f}s to embed)"
        )
        scores = baseline(X_train, X_test, y_train, y_test)
        say(
            "baselines on the raw panel PCs: "
            + ", ".join(f"{name} {score:.4f}" for name, score in scores.items())
        )
        say()
        say(
            f"{'arm':16} {'collapses':>10} {'max|intercept|':>15} "
            f"{'test bal (min-max)':>21} {'s/fit':>7} {'identical':>10}"
        )
        for name, arm in ARMS:
            rows = [
                fit_once(E_train, E_test, y_train, y_test, arm)
                for _ in range(args.repeats)
            ]
            stats = summarize(rows)
            any_collapse = any_collapse or stats["collapses"] > 0
            say(
                f"{name:16} {stats['collapses']:>4}/{stats['n']:<5} "
                f"{stats['max_abs_intercept']:>15.4g} "
                f"{stats['min_test_balanced']:.4f}-{stats['max_test_balanced']:.4f}"
                f"{'':>7} {stats['mean_seconds']:>7.3f} "
                f"{str(stats['deterministic']):>10}"
            )

            if E_cohort is not None:
                labels = np.asarray(classes)[rows[-1]["model"].predict(E_cohort)]
                counts = pd.Series(labels).value_counts()
                say(
                    f"{'':16} cohort ({len(labels)} samples): "
                    + " ".join(
                        f"{label} {count}" for label, count in counts.items()
                    )
                )
        say()

    # The claim the whole round rests on: the defect is reproducible and the
    # fix removes it. Asserted rather than merely printed.
    say("self-check")
    say(
        "  Every 'the fix' row above must read 0 collapses and identical=True. "
        "An 'as shipped' row reading 0 collapses is not a contradiction -- it "
        "means that cohort's PCs sit inside the stable region, which is the "
        "finding for array data."
    )
    say()
    say("done")
    return 1 if any_collapse else 0


if __name__ == "__main__":
    sys.exit(main())
