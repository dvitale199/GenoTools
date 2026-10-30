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

"""Read a finished run's JSON report into plot-ready frames.

Everything here is pure pandas over the report's own tables, so each frame can
be tested against a hand-built report without PLINK or a plotting backend.
"""

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, Optional

import pandas as pd

# Not every report has an `outcome` column: reports from before 2.0 carry only
# the boolean `pass`, which cannot tell a skip from a failure.
_LEGACY_OUTCOME = {True: "pass", False: "fail"}


@dataclass(frozen=True)
class Report:
    """The parts of a GenoTools JSON report the plots read.

    Attributes:
        path: Report the frames came from.
        raw: The parsed JSON, for plots that need a key this class does not
            unpack.
        qc: The long `QC` table, one row per (step, ancestry, metric), in
            pipeline order. Empty when the run did no QC.
        input_samples: The input `.psam` table, or None.
        ancestry_labels: Per-sample predicted labels, or None for a run
            without `--ancestry`.
    """

    path: Path
    raw: Dict[str, Any] = field(repr=False)
    qc: pd.DataFrame
    input_samples: Optional[pd.DataFrame] = None
    ancestry_labels: Optional[pd.DataFrame] = None


def _frame(raw: Dict[str, Any], key: str) -> Optional[pd.DataFrame]:
    """A report table, or None when the run did not produce it.

    Rows come back in the order the run wrote them - JSON objects keep key
    order - which for `QC` is pipeline order.
    """
    if key not in raw:
        return None
    return pd.DataFrame(raw[key])


def load_report(path: Path) -> Report:
    """Load a GenoTools JSON report.

    Args:
        path: The `{out}.json` a run wrote.

    Returns:
        The report, with its tables as DataFrames.

    Raises:
        FileNotFoundError: The report does not exist.
        ValueError: The file is not a JSON object.
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Report not found: {path}")
    try:
        with open(path) as f:
            raw = json.load(f)
    except json.JSONDecodeError as exc:
        raise ValueError(f"{path} is not valid JSON: {exc}") from exc
    if not isinstance(raw, dict):
        raise ValueError(f"{path} is not a GenoTools report (expected a JSON object)")

    qc = _frame(raw, "QC")
    if qc is None:
        qc = pd.DataFrame(
            columns=["step", "pruned_count", "metric", "ancestry", "level", "pass"]
        )
    if "outcome" not in qc.columns:
        qc["outcome"] = qc["pass"].map(_LEGACY_OUTCOME)
    if "reason" not in qc.columns:
        qc["reason"] = None

    return Report(
        path=path,
        raw=raw,
        qc=qc,
        input_samples=_frame(raw, "input_samples"),
        ancestry_labels=_frame(raw, "ancestry_labels"),
    )


def _group_sizes(report: Report) -> Dict[str, int]:
    """Samples entering QC, per ancestry group.

    A flat run starts from the whole input; an ancestry run QCs each predicted
    group separately, so each group starts from its own predicted count. Groups
    with no QC rows (dropped by `--min-samples`) never appear in the frames
    that use this, so they need no entry.
    """
    sizes: Dict[str, int] = {}
    if report.input_samples is not None:
        sizes["all"] = len(report.input_samples)
    labels = report.ancestry_labels
    if labels is not None and "label" in labels.columns:
        sizes.update(labels["label"].astype(str).value_counts().to_dict())
    return sizes


def attrition_frame(report: Report) -> pd.DataFrame:
    """What each QC step removed, and what was left, per ancestry and level.

    Steps that report more than one metric (`related_prune` counts related and
    duplicated samples separately) are summed: the metrics are disjoint, and the
    waterfall asks how many left at that step, not why.

    Args:
        report: A loaded report.

    Returns:
        One row per (ancestry, level, step) in pipeline order, with columns
        `ancestry`, `level`, `step`, `outcome`, `removed`, `start` and
        `remaining`. `outcome` is carried so a skipped or failed step can be
        told apart from one that ran and removed nothing.
        `start` and `remaining` are NaN where the report does not record a
        starting count - always, for variants, whose input count the report
        does not carry.
    """
    columns = ["ancestry", "level", "step", "outcome", "removed", "start", "remaining"]
    qc = report.qc
    if qc.empty:
        return pd.DataFrame(columns=columns)

    sizes = _group_sizes(report)
    rows = []
    for (ancestry, level), group in qc.groupby(["ancestry", "level"], sort=False):
        by_step = group.groupby("step", sort=False)
        removed = by_step["pruned_count"].sum().astype(int)
        outcomes = by_step["outcome"].first()
        start = sizes.get(str(ancestry)) if level == "sample" else None
        remaining = start
        for step, n_removed in removed.items():
            if remaining is not None:
                remaining -= n_removed
            rows.append(
                {
                    "ancestry": str(ancestry),
                    "level": level,
                    "step": step,
                    "outcome": outcomes[step],
                    "removed": int(n_removed),
                    "start": start,
                    "remaining": remaining,
                }
            )
    return pd.DataFrame(rows, columns=columns).astype(
        {"start": "float", "remaining": "float"}
    )


def outcome_frame(report: Report) -> pd.DataFrame:
    """Each QC step's outcome, per ancestry group.

    Args:
        report: A loaded report.

    Returns:
        One row per (step, ancestry), steps in pipeline order, with columns
        `step`, `ancestry`, `outcome` (`pass` / `fail` / `skipped`) and
        `reason`. A step a group never ran has no row.
    """
    columns = ["step", "ancestry", "outcome", "reason"]
    if report.qc.empty:
        return pd.DataFrame(columns=columns)

    # Every metric row of one step carries the same outcome; keep the first.
    first = report.qc.drop_duplicates(subset=["step", "ancestry"], keep="first")
    out = first[columns].copy()
    out["ancestry"] = out["ancestry"].astype(str)
    # Object dtype first: pandas 3's string dtype turns None back into NaN.
    reason = out["reason"].astype(object)
    out["reason"] = reason.where(reason.notna(), None)
    return out.reset_index(drop=True)
