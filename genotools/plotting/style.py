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
from dataclasses import dataclass
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

# Categorical slots, validated against SURFACE. Three sit under 3:1 contrast,
# so a categorical mark always carries its label too (legend, hover or direct
# label), never colour alone. The first eight are the dataviz reference
# palette; the last three were chosen to sit furthest from them.
CATEGORICAL = [
    "#2a78d6",  # blue
    "#eb6834",  # orange
    "#1baf7a",  # aqua
    "#eda100",  # yellow
    "#e87ba4",  # magenta
    "#008300",  # green
    "#6250d6",  # violet
    "#e34948",  # red
    "#14588f",  # navy
    "#883a2c",  # brick
    "#843753",  # wine
]

# Fixed colours for the default reference panel's labels, plus CAH (invented
# by admixture detection), so a group keeps its colour across runs and papers.
# This is a lookup with a fallback, not a vocabulary: `label_styles` handles
# any label. No eleven hues stay distinguishable in every pair, so the
# assignment gives look-alike colours to labels far apart in PC space; every
# pair of neighbouring clusters (the European/Middle-Eastern block, the
# African pair, the Asian groups, AMR and CAH against their sources) clears
# both the normal-vision and the CVD floor. Marker shape is the second cue.
PANEL_COLORS = {
    "AAC": "#008300",
    "AFR": "#e87ba4",
    "AJ": "#6250d6",
    "AMR": "#2a78d6",
    "CAH": "#843753",
    "CAS": "#14588f",
    "EAS": "#e34948",
    "EUR": "#1baf7a",
    "FIN": "#eb6834",
    "MDE": "#883a2c",
    "SAS": "#eda100",
}

# Marker shapes as (matplotlib, plotly) names: shape is the cue that survives
# where two colours do not. Boldest first; the star reads smallest.
MARKERS = [
    ("o", "circle"),
    ("s", "square"),
    ("D", "diamond"),
    ("^", "triangle-up"),
    ("h", "hexagon"),
    ("v", "triangle-down"),
    ("<", "triangle-left"),
    (">", "triangle-right"),
    ("P", "cross"),
    ("X", "x"),
    ("*", "star"),
]

# Fixed shapes for the panel labels: the boldest go to the crowded
# European/Middle-Eastern block, the star to CAH, which sits on its own.
PANEL_MARKERS = {
    "EUR": MARKERS[0],
    "FIN": MARKERS[1],
    "AJ": MARKERS[2],
    "MDE": MARKERS[3],
    "AFR": MARKERS[4],
    "AAC": MARKERS[5],
    "SAS": MARKERS[6],
    "AMR": MARKERS[7],
    "CAS": MARKERS[8],
    "EAS": MARKERS[9],
    "CAH": MARKERS[10],
}


@dataclass(frozen=True)
class LabelStyle:
    """How one ancestry label is drawn: colour plus a shape, never colour alone."""

    color: str
    marker: str  # matplotlib
    symbol: str  # plotly


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


def label_styles(labels: Iterable[str]) -> Dict[str, LabelStyle]:
    """A colour and marker per label that follow the label, never its rank.

    Panel labels take their `PANEL_COLORS` and `PANEL_MARKERS` entries. Any
    other label - the vocabulary is user-supplied via `--ref-labels` or a
    pickled encoder - takes the next colour and shape no label in this set
    already holds, in sorted order. Past the spare colours it falls back to
    NEUTRAL; past the spare shapes, shapes repeat. Either way the label itself
    (legend, hover) still names it. Pass the run's whole vocabulary, not one
    plot's subset, so every plot of a run agrees.
    """
    ordered = sorted({str(label) for label in labels})
    panel = [l for l in ordered if l in PANEL_COLORS]
    spare_colors = iter([c for c in CATEGORICAL if c not in {PANEL_COLORS[l] for l in panel}])
    spare_markers = [m for m in MARKERS if m not in {PANEL_MARKERS[l] for l in panel}] or MARKERS
    styles = {}
    n_custom = 0
    for label in ordered:
        if label in PANEL_COLORS:
            color, (marker, symbol) = PANEL_COLORS[label], PANEL_MARKERS[label]
        else:
            color = next(spare_colors, NEUTRAL)
            marker, symbol = spare_markers[n_custom % len(spare_markers)]
            n_custom += 1
        styles[label] = LabelStyle(color, marker, symbol)
    return styles


def label_colors(labels: Iterable[str]) -> Dict[str, str]:
    """Just the colours from `label_styles`."""
    return {label: style.color for label, style in label_styles(labels).items()}


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
