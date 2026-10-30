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

"""The figures, each drawn twice: a matplotlib PNG and a plotly figure.

Static PNGs go into papers and slides; the plotly figures go into the HTML
report, where hover carries the detail a PNG has no room for. Plotly could
write PNGs too, but only through kaleido and a headless Chrome, which is
exactly what is missing on a cluster node. Matplotlib needs neither.

Each figure reads one frame from `report.py` and nothing else.
"""

import math
from pathlib import Path
from typing import List, Tuple

import pandas as pd

# Colours from the dataviz reference palette. Removal is the one thing the
# attrition plots encode, so it gets the single series colour; inputs and
# survivors are neutral.
REMOVED = "#2a78d6"
NEUTRAL = "#b4b2a9"
INK = "#0b0b0b"
INK_SECONDARY = "#52514e"
SURFACE = "#fcfcfb"
GRID = "#e6e5e0"

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


def _attrition_labels(group: pd.DataFrame) -> List[str]:
    """Step labels, marking any step that did not run.

    A skipped or failed step removes nothing, which on its own reads the same
    as a step that ran and found nothing to remove.
    """
    return [
        step_label(step) if outcome == "pass" else f"{step_label(step)} ({outcome})"
        for step, outcome in zip(group["step"], group["outcome"])
    ]


def _grid(n: int, max_cols: int = 3) -> Tuple[int, int]:
    cols = min(n, max_cols)
    return math.ceil(n / cols), cols


def _style_axis(ax) -> None:
    ax.set_facecolor(SURFACE)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(GRID)
    ax.tick_params(colors=INK_SECONDARY, labelsize="small")
    ax.yaxis.grid(True, color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)


def _groups(frame: pd.DataFrame) -> List[str]:
    """Ancestry groups in a stable order: `all` first, then alphabetical."""
    groups = sorted(frame["ancestry"].unique())
    return sorted(groups, key=lambda g: g != "all")


# --- Attrition ---------------------------------------------------------------


def _waterfall_bars(group: pd.DataFrame) -> Tuple[List[str], List[float], List[float], List[str]]:
    """Labels, bottoms, heights and colours for one group's sample waterfall."""
    start = float(group["start"].iloc[0])
    labels = ["input"] + _attrition_labels(group) + ["kept"]
    bottoms = [0.0] + list(group["remaining"]) + [0.0]
    heights = [start] + list(group["removed"]) + [float(group["remaining"].iloc[-1])]
    colors = [NEUTRAL] + [REMOVED] * len(group) + [NEUTRAL]
    return labels, bottoms, heights, colors


def attrition_png(frame: pd.DataFrame, level: str, plt, out_path: Path) -> Path:
    """Small multiples, one panel per ancestry group.

    Samples are a waterfall from the input count down to what was kept.
    Variants are plain removal bars: the report does not record how many
    variants went in, so there is nothing to fall from.
    """
    frame = frame[frame["level"] == level]
    groups = _groups(frame)
    rows, cols = _grid(len(groups))
    fig, axes = plt.subplots(
        rows, cols, figsize=(4.6 * cols, 3.4 * rows), squeeze=False
    )
    fig.patch.set_facecolor(SURFACE)

    for ax, ancestry in zip(axes.flat, groups):
        group = frame[frame["ancestry"] == ancestry]
        _style_axis(ax)
        if level == "sample" and group["start"].notna().all():
            labels, bottoms, heights, colors = _waterfall_bars(group)
            start, kept = heights[0], heights[-1]
            title = f"{ancestry}: {start:,.0f} → {kept:,.0f}"
            if start:
                title += f" ({100 * kept / start:.1f}% kept)"
        else:
            labels = _attrition_labels(group)
            bottoms = [0.0] * len(group)
            heights = list(group["removed"])
            colors = [REMOVED] * len(group)
            title = f"{ancestry}: {int(group['removed'].sum()):,} {level}s removed"

        x = range(len(labels))
        ax.bar(x, heights, bottom=bottoms, color=colors, width=0.7)
        # Label the removals only; the axis already carries the totals.
        for i, (b, h, c) in enumerate(zip(bottoms, heights, colors)):
            if c == REMOVED and h > 0:
                ax.annotate(
                    f"−{h:,.0f}", (i, b + h), xytext=(0, 2),
                    textcoords="offset points", ha="center", va="bottom",
                    fontsize="x-small", color=INK_SECONDARY,
                )
        ax.set_xticks(list(x))
        ax.set_xticklabels(labels, rotation=35, ha="right")
        ax.set_title(title, loc="left", fontsize="medium", color=INK)
        ax.set_ylabel(f"{level}s", color=INK_SECONDARY, fontsize="small")

    for ax in list(axes.flat)[len(groups):]:
        ax.set_visible(False)

    fig.tight_layout()
    fig.savefig(out_path, dpi=150, facecolor=SURFACE)
    plt.close(fig)
    return out_path


def attrition_plotly(frame: pd.DataFrame, level: str):
    """The same small multiples as `attrition_png`, with counts on hover."""
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots

    frame = frame[frame["level"] == level]
    groups = _groups(frame)
    rows, cols = _grid(len(groups))
    fig = make_subplots(rows=rows, cols=cols, subplot_titles=groups)

    for i, ancestry in enumerate(groups):
        group = frame[frame["ancestry"] == ancestry]
        row, col = i // cols + 1, i % cols + 1
        steps = _attrition_labels(group)
        if level == "sample" and group["start"].notna().all():
            trace = go.Waterfall(
                x=["input"] + steps + ["kept"],
                y=[group["start"].iloc[0]] + [-r for r in group["removed"]] + [0],
                measure=["absolute"] + ["relative"] * len(group) + ["total"],
                decreasing={"marker": {"color": REMOVED}},
                totals={"marker": {"color": NEUTRAL}},
                increasing={"marker": {"color": NEUTRAL}},
                connector={"line": {"color": GRID}},
                name=ancestry,
                hovertemplate="%{x}: %{y:,}<extra>" + ancestry + "</extra>",
            )
        else:
            trace = go.Bar(
                x=steps,
                y=group["removed"],
                marker_color=REMOVED,
                name=ancestry,
                hovertemplate="%{x}: %{y:,} removed<extra>" + ancestry + "</extra>",
            )
        fig.add_trace(trace, row=row, col=col)

    fig.update_layout(
        showlegend=False,
        height=300 * rows + 60,
        plot_bgcolor=SURFACE,
        paper_bgcolor=SURFACE,
        font={"color": INK},
        margin={"t": 60, "l": 50, "r": 20, "b": 40},
    )
    fig.update_yaxes(gridcolor=GRID, title_text=f"{level}s")
    return fig


# --- Outcome grid ------------------------------------------------------------


def _outcome_matrix(frame: pd.DataFrame) -> pd.DataFrame:
    """Step x ancestry matrix of outcomes; `not run` where a group has no row."""
    steps = list(dict.fromkeys(frame["step"]))
    matrix = frame.pivot(index="step", columns="ancestry", values="outcome")
    matrix = matrix.reindex(index=steps, columns=_groups(frame))
    return matrix.fillna("not run")


def outcome_png(frame: pd.DataFrame, plt, out_path: Path) -> Path:
    """Step x ancestry grid of pass / fail / skipped, outcome written in each cell."""
    from matplotlib.colors import ListedColormap

    matrix = _outcome_matrix(frame)
    codes = matrix.apply(lambda col: col.map(OUTCOME_ORDER.index))
    cmap = ListedColormap([OUTCOME_COLORS[o] for o in OUTCOME_ORDER])

    n_steps, n_groups = matrix.shape
    fig, ax = plt.subplots(
        figsize=(1.6 + 1.1 * n_groups, 0.9 + 0.42 * n_steps)
    )
    fig.patch.set_facecolor(SURFACE)
    ax.imshow(codes.values, cmap=cmap, vmin=-0.5, vmax=len(OUTCOME_ORDER) - 0.5, aspect="auto")

    for (r, c), outcome in pd.DataFrame(matrix.values).stack().items():
        ax.text(
            c, r, outcome, ha="center", va="center", fontsize="x-small",
            color="white" if outcome in ("pass", "fail") else INK,
        )
    # A 2px surface gap between cells.
    ax.set_xticks([x - 0.5 for x in range(1, n_groups)], minor=True)
    ax.set_yticks([y - 0.5 for y in range(1, n_steps)], minor=True)
    ax.grid(which="minor", color=SURFACE, linewidth=2)
    ax.tick_params(which="minor", length=0)

    ax.set_xticks(range(n_groups))
    ax.set_xticklabels(matrix.columns, fontsize="small", color=INK_SECONDARY)
    ax.set_yticks(range(n_steps))
    ax.set_yticklabels([step_label(s) for s in matrix.index], fontsize="small", color=INK_SECONDARY)
    ax.xaxis.tick_top()
    for spine in ax.spines.values():
        spine.set_visible(False)

    fig.tight_layout()
    fig.savefig(out_path, dpi=150, facecolor=SURFACE)
    plt.close(fig)
    return out_path


def outcome_plotly(frame: pd.DataFrame):
    """The outcome grid, with each skip or failure's reason on hover."""
    import plotly.graph_objects as go

    matrix = _outcome_matrix(frame)
    reasons = (
        frame.pivot(index="step", columns="ancestry", values="reason")
        .reindex(index=matrix.index, columns=matrix.columns)
    )
    codes = matrix.apply(lambda col: col.map(OUTCOME_ORDER.index))
    hover = [
        [
            f"{step_label(step)} / {group}: {matrix.at[step, group]}"
            + (f"<br>{reasons.at[step, group]}" if pd.notna(reasons.at[step, group]) else "")
            for group in matrix.columns
        ]
        for step in matrix.index
    ]

    # A discrete colorscale: each code owns an equal band.
    n = len(OUTCOME_ORDER)
    colorscale = []
    for i, outcome in enumerate(OUTCOME_ORDER):
        colorscale += [[i / n, OUTCOME_COLORS[outcome]], [(i + 1) / n, OUTCOME_COLORS[outcome]]]

    fig = go.Figure(
        go.Heatmap(
            z=codes.values,
            x=list(matrix.columns),
            y=[step_label(s) for s in matrix.index],
            text=matrix.values,
            texttemplate="%{text}",
            hovertext=hover,
            hoverinfo="text",
            colorscale=colorscale,
            zmin=-0.5,
            zmax=n - 0.5,
            showscale=False,
            xgap=2,
            ygap=2,
        )
    )
    fig.update_layout(
        height=120 + 36 * len(matrix.index),
        plot_bgcolor=SURFACE,
        paper_bgcolor=SURFACE,
        font={"color": INK},
        margin={"t": 40, "l": 120, "r": 20, "b": 20},
    )
    fig.update_xaxes(side="top")
    fig.update_yaxes(autorange="reversed")
    return fig
