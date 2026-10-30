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

"""Colours and layout helpers every plot module shares.

Every figure is drawn twice, as a matplotlib PNG and a plotly figure, and both
take their colours from here so the two never drift apart.
"""

import math
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

import pandas as pd

# Colours from the dataviz reference palette (light mode). One series colour
# for the quantity a plot is about; inputs, survivors and context are neutral.
SERIES = "#2a78d6"
REMOVED = SERIES
NEUTRAL = "#b4b2a9"
INK = "#0b0b0b"
INK_SECONDARY = "#52514e"
SURFACE = "#fcfcfb"
GRID = "#e6e5e0"

# Threshold lines and the band a filter keeps: drawn in ink, never in a series
# or status colour, so a cutoff never reads as data.
THRESHOLD = INK_SECONDARY

# Categorical slots in fixed order, validated for colour-vision deficiency
# against SURFACE. Three slots sit under 3:1 contrast, so a categorical mark
# always carries its label too (legend, hover or direct label), never colour
# alone. There is no ninth colour: see `label_colors`.
CATEGORICAL = [
    "#2a78d6",  # blue
    "#eb6834",  # orange
    "#1baf7a",  # aqua
    "#eda100",  # yellow
    "#e87ba4",  # magenta
    "#008300",  # green
    "#6250d6",  # violet
    "#e34948",  # red
]

# Status colours, reserved for outcome; never reused for a series. Every cell
# also carries its outcome as text, so state is never colour alone.
OUTCOME_COLORS = {
    "pass": "#0ca30c",
    "fail": "#d03b3b",
    "skipped": "#fab219",
    "not run": "#f0efec",
}
OUTCOME_ORDER = ["pass", "fail", "skipped", "not run"]

# Reported step names are long and mostly the same suffix.
_STEP_LABELS = {"case_control_missingness_prune": "case_control"}


def step_label(step: str) -> str:
    """Short display name for a reported step (`het_prune` -> `het`)."""
    return _STEP_LABELS.get(step, step.removesuffix("_prune"))


def label_colors(labels: Iterable[str]) -> Dict[str, str]:
    """A colour per label that follows the label, never its rank.

    Labels are sorted before slots are handed out, so `EUR` gets the same
    colour in every plot of a run whichever labels a given plot shows - pass
    the run's whole vocabulary, not just one plot's subset. Ancestry labels
    are user-supplied (`--ref-labels`, a pickled encoder, or the invented
    `CAH`), so nothing here names one.

    Past eight labels there are no more distinguishable hues; the rest share
    NEUTRAL and must be told apart by their label (hover, legend or facet).
    """
    ordered = sorted({str(label) for label in labels})
    return {
        label: CATEGORICAL[i] if i < len(CATEGORICAL) else NEUTRAL
        for i, label in enumerate(ordered)
    }


def groups(frame: pd.DataFrame, column: str = "ancestry") -> List[str]:
    """Ancestry groups in a stable order: `all` first, then alphabetical."""
    found = sorted(str(g) for g in frame[column].unique())
    return sorted(found, key=lambda g: g != "all")


def grid(n: int, max_cols: int = 3) -> Tuple[int, int]:
    """Rows and columns for `n` small multiples."""
    cols = max(1, min(n, max_cols))
    return math.ceil(n / cols), cols


def style_axis(ax) -> None:
    """Recessive axes: no top/right spines, light grid behind the marks."""
    ax.set_facecolor(SURFACE)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(GRID)
    ax.tick_params(colors=INK_SECONDARY, labelsize="small")
    ax.yaxis.grid(True, color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)


def save_png(fig, plt, out_path: Path) -> Path:
    """Write a matplotlib figure on the shared surface and close it."""
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, facecolor=SURFACE)
    plt.close(fig)
    return out_path


def plotly_layout(fig, height: int, **overrides):
    """Apply the shared plotly surface, ink and margins; return the figure."""
    layout = {
        "height": height,
        "plot_bgcolor": SURFACE,
        "paper_bgcolor": SURFACE,
        "font": {"color": INK},
        "margin": {"t": 60, "l": 50, "r": 20, "b": 40},
    }
    layout.update(overrides)
    fig.update_layout(**layout)
    return fig
