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

"""Plot-ready tables a QC step keeps from the statistics it computes.

Each step computes a per-sample or per-variant statistic, applies a threshold,
and until now deleted the statistic. These helpers turn it into something
`genotools-plot` can draw, with a size that does not grow with the array:

* per-sample tables are kept whole - a cohort has thousands of samples, not
  millions;
* per-variant statistics are kept as a fixed-bin histogram, plus the rows that
  failed (capped), which is enough to draw a distribution or a QQ plot.

A step writes its tables beside its output as `{out}_plotdata_{kind}.tsv` and
returns them in `FilterResult.plot_data`; the pipeline runner moves them into
the run's `{out}_plot_data/` directory.
"""

from pathlib import Path
from typing import Dict, Optional

import numpy as np
import pandas as pd

# 0.005-wide bins resolve the usual --geno / --mind thresholds (0.01-0.05).
MISSINGNESS_EDGES = np.linspace(0.0, 1.0, 201)

# -log10(p) in 0.05-wide bins up to 50, then one open bin. p == 0 lands there.
NEG_LOG10_P_EDGES = np.append(np.linspace(0.0, 50.0, 1001), np.inf)

# Failing variants kept per table, most extreme first. A badly genotyped array
# can fail --geno on a large fraction of its variants; the histogram still
# counts all of them.
MAX_FAILING_ROWS = 10_000


def plot_data_path(out_path: Path, kind: str) -> Path:
    """Where a step writes one plot table, beside its output prefix."""
    return out_path.parent / f"{out_path.name}_plotdata_{kind}.tsv"


def histogram(values: pd.Series, edges: np.ndarray) -> pd.DataFrame:
    """Counts of `values` in `edges`, closed on the left, last bin closed.

    Returns:
        One row per bin: `bin_lower`, `bin_upper`, `count`. Every bin is
        present, empty or not, so tables from different groups line up.
        NaNs are dropped.
    """
    values = pd.to_numeric(values, errors="coerce").dropna()
    counts, _ = np.histogram(values, bins=edges)
    return pd.DataFrame(
        {"bin_lower": edges[:-1], "bin_upper": edges[1:], "count": counts}
    )


def neg_log10_p(p: pd.Series) -> pd.Series:
    """-log10(p), with p == 0 mapped to +inf rather than raising."""
    p = pd.to_numeric(p, errors="coerce")
    with np.errstate(divide="ignore"):
        return -np.log10(p)


def pvalue_histogram(p: pd.Series, by: Optional[pd.Series] = None) -> pd.DataFrame:
    """Histogram of -log10(p), optionally one per group (`by`).

    Returns:
        `bin_lower`, `bin_upper`, `count` in -log10(p) units, with a leading
        group column named after `by` when it is given.
    """
    scores = neg_log10_p(p)
    if by is None:
        return histogram(scores, NEG_LOG10_P_EDGES)
    frames = []
    for group, idx in scores.groupby(by, sort=True).groups.items():
        frame = histogram(scores.loc[idx], NEG_LOG10_P_EDGES)
        frame.insert(0, by.name, group)
        frames.append(frame)
    return pd.concat(frames, ignore_index=True)


def failing_rows(
    frame: pd.DataFrame,
    mask: pd.Series,
    sort_by: str,
    ascending: bool,
    limit: int = MAX_FAILING_ROWS,
) -> pd.DataFrame:
    """The rows that failed a threshold, most extreme first, capped at `limit`."""
    return frame[mask].sort_values(sort_by, ascending=ascending).head(limit)


def write_plot_data(
    tables: Dict[str, pd.DataFrame], out_path: Path
) -> Dict[str, Path]:
    """Write each table as a TSV beside `out_path`.

    Returns:
        `kind -> path`, ready for `FilterResult.plot_data`.
    """
    written: Dict[str, Path] = {}
    for kind, frame in tables.items():
        path = plot_data_path(out_path, kind)
        frame.to_csv(path, sep="\t", index=False)
        written[kind] = path
    return written
