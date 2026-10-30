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

"""Sample QC plots: missingness, sex check, heterozygosity, relatedness.

Each plot is one small multiple per ancestry group (`all` first). The pure
logic - which samples a threshold flags, which bounds a group used - lives in
frame builders that return `HistPanel`s, so it can be tested on a hand-built
frame; the matplotlib and plotly renderers only draw what a panel holds.

Thresholds come from the report's `parameters` table for that group, never a
default: a threshold the report does not record is simply not drawn.
"""

from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

from genotools.plotting.plots.base import Plot
from genotools.plotting.report import Report, parameter, plot_tables
from genotools.plotting.style import (
    CATEGORICAL,
    GRID,
    INK,
    INK_SECONDARY,
    NEUTRAL,
    SERIES,
    SURFACE,
    THRESHOLD,
    grid,
    plotly_layout,
    save_png,
    style_axis,
)

BINS = 40

# KING's kinship boundaries between degrees of relationship (Manichaikul et
# al. 2010). Definitional constants, not settings: the step's own cutoffs are
# read from the report and drawn separately.
KING_DEGREES: List[Tuple[float, str]] = [
    (0.354, "duplicate / MZ"),
    (0.177, "1st degree"),
    (0.0884, "2nd degree"),
    (0.0442, "3rd degree"),
]

SEX_NAMES = {1: "male", 2: "female", 0: "unknown"}
# Reported sex is an identity, so it takes categorical slots; unknown is
# context rather than a category of interest.
SEX_COLORS = {1: CATEGORICAL[0], 2: CATEGORICAL[1], 0: NEUTRAL}


# --- Parameters and frames ---------------------------------------------------


def _number(value) -> Optional[float]:
    """A parameter as a float, or None if absent or not numeric."""
    try:
        out = float(value)
    except (TypeError, ValueError):
        return None
    return None if np.isnan(out) else out


def _threshold(report: Report, step: str, name: str, ancestry: str) -> Optional[float]:
    return _number(parameter(report, step, name, ancestry))


def _ids(frame: pd.DataFrame, fid: str = "#FID", iid: str = "IID") -> pd.Series:
    """`FID IID` labels for hover."""
    return frame[fid].astype(str) + " " + frame[iid].astype(str)


def _has_table(name: str):
    """`available` for a plot: the report kept this table for some group."""

    def available(report: Report) -> bool:
        try:
            return bool(plot_tables(report, name))
        except OSError:
            return False

    return available


def above_threshold(frame: pd.DataFrame, column: str, bound: Optional[float]) -> pd.DataFrame:
    """Rows whose `column` is strictly above `bound`; none if there is no bound."""
    if bound is None:
        return frame.iloc[0:0]
    return frame[frame[column] > bound]


def outside_bounds(
    frame: pd.DataFrame, column: str, lower: Optional[float], upper: Optional[float]
) -> pd.DataFrame:
    """Rows whose `column` is outside the keep band; an absent side flags nothing."""
    mask = pd.Series(False, index=frame.index)
    if lower is not None:
        mask |= frame[column] < lower
    if upper is not None:
        mask |= frame[column] > upper
    return frame[mask]


def sex_pass(frame: pd.DataFrame) -> pd.DataFrame:
    """The one `--check-sex` pass the sex plot shows.

    `sex_check.tsv` holds both passes. The step prunes a sample that PLINK
    flags in either, but the plot shows `filtered_x`: it is computed on the
    non-pseudoautosomal X with variants filtered for MAF, missingness and HWE,
    so F is not blurred by PAR or poorly genotyped sites, and the two
    cutoffs separate the sexes most cleanly. `all_x` is the fallback for a
    table without it. Showing both would double every sample and put two F
    values on one hover.
    """
    if "check" not in frame.columns:
        return frame
    for name in ("filtered_x", "all_x"):
        chosen = frame[frame["check"] == name]
        if not chosen.empty:
            return chosen
    return frame


# --- Histogram panels --------------------------------------------------------


@dataclass
class HistSeries:
    name: str
    color: str
    values: np.ndarray


@dataclass
class HistPanel:
    """Everything one histogram small multiple draws, already decided.

    `lines` are `(x, label)` threshold lines, `band` an optional keep band, and
    `flagged` the samples the threshold marks, with a hover string each.
    """

    ancestry: str
    title: str
    series: List[HistSeries]
    lines: List[Tuple[float, str]] = field(default_factory=list)
    band: Optional[Tuple[float, float]] = None
    flagged_x: np.ndarray = field(default_factory=lambda: np.array([]))
    flagged_hover: List[str] = field(default_factory=list)
    flag_label: str = ""


def _flag_title(ancestry: str, n: int, flagged: int, phrase: Optional[str]) -> str:
    if phrase is None:
        return f"{ancestry}: {n:,} samples"
    return f"{ancestry}: {flagged:,} of {n:,} {phrase}"


def missingness_panels(report: Report) -> List[HistPanel]:
    """Sample F_MISS per group, with the call-rate cutoff and who exceeds it."""
    panels = []
    for ancestry, frame in plot_tables(report, "callrate_sample_missingness").items():
        mind = _threshold(report, "callrate_prune", "mind", ancestry)
        flagged = above_threshold(frame, "F_MISS", mind)
        phrase = None if mind is None else f"above mind {mind:g}"
        panels.append(
            HistPanel(
                ancestry=ancestry,
                title=_flag_title(ancestry, len(frame), len(flagged), phrase),
                series=[HistSeries("samples", SERIES, frame["F_MISS"].to_numpy(float))],
                lines=[] if mind is None else [(mind, f"mind {mind:g}")],
                flagged_x=flagged["F_MISS"].to_numpy(float),
                flagged_hover=[
                    f"{who}<br>F_MISS {x:.4f}"
                    for who, x in zip(_ids(flagged), flagged["F_MISS"])
                ],
                flag_label=f"above mind ({len(flagged):,})",
            )
        )
    return panels


def sex_panels(report: Report) -> List[HistPanel]:
    """X-chromosome F by reported sex, with PLINK's two cutoffs."""
    panels = []
    for ancestry, table in plot_tables(report, "sex_check").items():
        frame = sex_pass(table)
        female_max = _threshold(report, "sex_prune", "female_max_f", ancestry)
        male_min = _threshold(report, "sex_prune", "male_min_f", ancestry)
        lines = []
        if female_max is not None:
            lines.append((female_max, f"female max {female_max:g}"))
        if male_min is not None:
            lines.append((male_min, f"male min {male_min:g}"))

        series = []
        for code in (2, 1, 0):
            values = frame.loc[frame["PEDSEX"] == code, "F"].to_numpy(float)
            if len(values):
                series.append(HistSeries(SEX_NAMES[code], SEX_COLORS[code], values))
        # PLINK's own verdict on this pass: it applied the cutoffs above.
        flagged = frame[frame["STATUS"] == "PROBLEM"]
        panels.append(
            HistPanel(
                ancestry=ancestry,
                title=_flag_title(
                    ancestry, len(frame), len(flagged), "flagged by --check-sex"
                ),
                series=series,
                lines=lines,
                flagged_x=flagged["F"].to_numpy(float),
                flagged_hover=[
                    f"{who}<br>F {f:.3f}, reported {SEX_NAMES.get(int(p), p)}, "
                    f"genetic {SEX_NAMES.get(int(s), s)}"
                    for who, f, p, s in zip(
                        _ids(flagged), flagged["F"], flagged["PEDSEX"], flagged["SNPSEX"]
                    )
                ],
                flag_label=f"flagged ({len(flagged):,})",
            )
        )
    return panels


def het_bounds(report: Report, ancestry: str) -> Tuple[Optional[float], Optional[float]]:
    """The keep band `het_prune` applied to this group.

    `parameter` prefers the `resolved` value, which is what differs by group
    under `--het sd`; the `requested` `f_lower`/`f_upper` do not.
    """
    return (
        _threshold(report, "het_prune", "het_lower", ancestry),
        _threshold(report, "het_prune", "het_upper", ancestry),
    )


def het_panels(report: Report) -> List[HistPanel]:
    """Heterozygosity F per group, with that group's resolved keep band."""
    panels = []
    for ancestry, frame in plot_tables(report, "het_f").items():
        lower, upper = het_bounds(report, ancestry)
        outliers = outside_bounds(frame, "F", lower, upper)
        lines = []
        if lower is not None:
            lines.append((lower, f"lower {lower:.3g}"))
        if upper is not None:
            lines.append((upper, f"upper {upper:.3g}"))
        band = (lower, upper) if lower is not None and upper is not None else None
        phrase = None if lower is None and upper is None else "outside the keep band"
        panels.append(
            HistPanel(
                ancestry=ancestry,
                title=_flag_title(ancestry, len(frame), len(outliers), phrase),
                series=[HistSeries("samples", SERIES, frame["F"].to_numpy(float))],
                lines=lines,
                band=band,
                flagged_x=outliers["F"].to_numpy(float),
                flagged_hover=[
                    f"{who}<br>F {f:.4f}" for who, f in zip(_ids(outliers), outliers["F"])
                ],
                flag_label=f"outliers ({len(outliers):,})",
            )
        )
    return panels


def _edges(panel: HistPanel) -> np.ndarray:
    """Shared bin edges for every series of one panel, so overlays line up."""
    values = np.concatenate([s.values for s in panel.series]) if panel.series else np.array([])
    values = values[np.isfinite(values)]
    return np.histogram_bin_edges(values, bins=BINS)


def hist_png(panels: List[HistPanel], xlabel: str, plt, out_path: Path) -> Path:
    """Histogram small multiples: one panel per group."""
    rows, cols = grid(len(panels))
    fig, axes = plt.subplots(rows, cols, figsize=(4.6 * cols, 3.4 * rows), squeeze=False)
    fig.patch.set_facecolor(SURFACE)
    for ax, panel in zip(axes.flat, panels):
        style_axis(ax)
        edges = _edges(panel)
        if panel.band is not None:
            ax.axvspan(*panel.band, color=THRESHOLD, alpha=0.07, linewidth=0, zorder=0)
        multi = len(panel.series) > 1
        for series in panel.series:
            ax.hist(
                series.values, bins=edges, color=series.color,
                alpha=0.6 if multi else 1.0, label=series.name if multi else None,
                histtype="stepfilled" if multi else "bar", linewidth=1.0,
                edgecolor=series.color if multi else SURFACE,
            )
        ymax = ax.get_ylim()[1]
        for x, label in panel.lines:
            ax.axvline(x, color=THRESHOLD, linestyle="--", linewidth=1.2, zorder=3)
            ax.annotate(
                label, (x, ymax), xytext=(3, -3), textcoords="offset points",
                rotation=90, ha="left", va="top", fontsize="x-small",
                color=INK_SECONDARY,
            )
        if len(panel.flagged_x):
            # A rug along the baseline: marks which samples, without a second
            # colour competing with the quantity.
            ax.plot(
                panel.flagged_x, np.zeros(len(panel.flagged_x)), "|", color=INK,
                markersize=9, markeredgewidth=1.2, clip_on=False, zorder=4,
                label=panel.flag_label,
            )
        if multi or len(panel.flagged_x):
            # Below the axis: inside, it lands on a threshold label or a tall bar.
            ax.legend(
                frameon=False, fontsize="x-small", loc="upper center",
                bbox_to_anchor=(0.5, -0.2), ncol=3,
            )
        ax.set_title(panel.title, loc="left", fontsize="medium", color=INK)
        ax.set_xlabel(xlabel, color=INK_SECONDARY, fontsize="small")
        ax.set_ylabel("samples", color=INK_SECONDARY, fontsize="small")
    for ax in list(axes.flat)[len(panels):]:
        ax.set_visible(False)
    return save_png(fig, plt, out_path)


def hist_plotly(panels: List[HistPanel], xlabel: str):
    """The same small multiples; hovering a rug mark names the sample."""
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots

    rows, cols = grid(len(panels))
    fig = make_subplots(rows=rows, cols=cols, subplot_titles=[p.title for p in panels])
    for i, panel in enumerate(panels):
        row, col = i // cols + 1, i % cols + 1
        edges = _edges(panel)
        xbins = {"start": float(edges[0]), "end": float(edges[-1]), "size": float(edges[1] - edges[0])}
        multi = len(panel.series) > 1
        if panel.band is not None:
            fig.add_vrect(
                x0=panel.band[0], x1=panel.band[1], fillcolor=THRESHOLD, opacity=0.07,
                line_width=0, layer="below", row=row, col=col,
            )
        for series in panel.series:
            fig.add_trace(
                go.Histogram(
                    x=series.values, xbins=xbins, marker_color=series.color,
                    opacity=0.6 if multi else 1.0, name=series.name,
                    legendgroup=series.name, showlegend=multi and i == 0,
                    hovertemplate=f"{series.name} (%{{x}}): %{{y}} samples<extra>{panel.ancestry}</extra>",
                ),
                row=row, col=col,
            )
        for x, label in panel.lines:
            fig.add_vline(
                x=x, line_dash="dash", line_color=THRESHOLD, line_width=1.2,
                annotation_text=label, annotation_font_size=10,
                annotation_font_color=INK_SECONDARY, annotation_position="top",
                row=row, col=col,
            )
        if len(panel.flagged_x):
            fig.add_trace(
                go.Scatter(
                    x=panel.flagged_x, y=np.zeros(len(panel.flagged_x)), mode="markers",
                    marker={"symbol": "line-ns-open", "color": INK, "size": 12,
                            "line": {"width": 1.5, "color": INK}},
                    name="flagged", legendgroup="flagged", showlegend=i == 0,
                    hovertext=panel.flagged_hover, hoverinfo="text",
                ),
                row=row, col=col,
            )
    plotly_layout(
        fig, 320 * rows + 60, barmode="overlay", showlegend=True,
        legend={"orientation": "h", "y": -0.18}, hovermode="closest",
    )
    fig.update_xaxes(title_text=xlabel, gridcolor=GRID)
    fig.update_yaxes(title_text="samples", gridcolor=GRID)
    return fig


# --- Relatedness -------------------------------------------------------------


@dataclass
class KinshipPanel:
    ancestry: str
    title: str
    ibs0: np.ndarray
    kinship: np.ndarray
    hover: List[str]
    cutoffs: List[Tuple[float, str]]


def kinship_panels(report: Report) -> List[KinshipPanel]:
    """KING pairs per group, with the step's own cutoffs."""
    panels = []
    for ancestry, frame in plot_tables(report, "related_kinship").items():
        cutoffs = []
        related = _threshold(report, "related_prune", "related_cutoff", ancestry)
        duplicated = _threshold(report, "related_prune", "duplicated_cutoff", ancestry)
        if related is not None:
            cutoffs.append((related, f"related cutoff {related:g}"))
        if duplicated is not None:
            cutoffs.append((duplicated, f"duplicate cutoff {duplicated:g}"))
        hover = [
            f"{a} / {b}<br>kinship {k:.4f}, IBS0 {i:.4f}"
            for a, b, k, i in zip(
                _ids(frame, "#FID1", "IID1"), _ids(frame, "FID2", "IID2"),
                frame["KINSHIP"], frame["IBS0"],
            )
        ]
        panels.append(
            KinshipPanel(
                ancestry=ancestry,
                title=f"{ancestry}: {len(frame):,} related pairs",
                ibs0=frame["IBS0"].to_numpy(float),
                kinship=frame["KINSHIP"].to_numpy(float),
                hover=hover,
                cutoffs=cutoffs,
            )
        )
    return panels


def _kinship_ylim(panel: KinshipPanel) -> Tuple[float, float]:
    """Room for every degree band and any pair below the lowest one."""
    low = min([KING_DEGREES[-1][0] * 0.5] + list(panel.kinship))
    high = max([0.5] + list(panel.kinship))
    pad = 0.02 * (high - low)
    return low - pad, high + pad


def _ibs0_xlim(panel: KinshipPanel) -> Tuple[float, float]:
    """IBS0 from zero, so how close to zero the pairs sit is honest."""
    top = float(panel.ibs0.max()) * 1.1 if len(panel.ibs0) else 0.0
    return 0.0, top or 0.01


def _band_mid(i: int, top: float) -> float:
    """Middle of degree band `i` (0 = duplicate), the top band ending at `top`."""
    upper = top if i == 0 else KING_DEGREES[i - 1][0]
    return (KING_DEGREES[i][0] + upper) / 2


def kinship_png(panels: List[KinshipPanel], plt, out_path: Path) -> Path:
    """Kinship (y) against IBS0 (x) with the KING degree bands."""
    from matplotlib.lines import Line2D

    rows, cols = grid(len(panels))
    fig, axes = plt.subplots(rows, cols, figsize=(4.8 * cols, 3.8 * rows), squeeze=False)
    fig.patch.set_facecolor(SURFACE)
    for ax, panel in zip(axes.flat, panels):
        style_axis(ax)
        ax.xaxis.grid(True, color=GRID, linewidth=0.8)
        low, high = _kinship_ylim(panel)
        ax.set_ylim(low, high)
        ax.set_xlim(*_ibs0_xlim(panel))
        for i, (bound, name) in enumerate(KING_DEGREES):
            ax.axhline(bound, color=NEUTRAL, linewidth=0.9, linestyle=":", zorder=1)
            ax.text(
                0.99, _band_mid(i, high), name, transform=ax.get_yaxis_transform(),
                ha="right", va="center", fontsize="x-small", color=INK_SECONDARY,
            )
        for bound, _ in panel.cutoffs:
            ax.axhline(bound, color=THRESHOLD, linewidth=1.2, linestyle="--", zorder=2)
        ax.scatter(panel.ibs0, panel.kinship, s=22, color=SERIES, alpha=0.8,
                   edgecolors=SURFACE, linewidths=0.5, zorder=3)
        handles = [
            Line2D([], [], color=THRESHOLD, linestyle="--", linewidth=1.2, label=label)
            for _, label in panel.cutoffs
        ]
        if handles:
            ax.legend(handles=handles, frameon=False, fontsize="x-small", loc="center right",
                      bbox_to_anchor=(1.0, 0.62))
        ax.set_title(panel.title, loc="left", fontsize="medium", color=INK)
        ax.set_xlabel("IBS0", color=INK_SECONDARY, fontsize="small")
        ax.set_ylabel("KING kinship", color=INK_SECONDARY, fontsize="small")
    for ax in list(axes.flat)[len(panels):]:
        ax.set_visible(False)
    return save_png(fig, plt, out_path)


def kinship_plotly(panels: List[KinshipPanel]):
    """The same scatter; hovering a point names both samples of the pair."""
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots

    rows, cols = grid(len(panels))
    fig = make_subplots(rows=rows, cols=cols, subplot_titles=[p.title for p in panels])
    for i, panel in enumerate(panels):
        row, col = i // cols + 1, i % cols + 1
        low, high = _kinship_ylim(panel)
        for j, (bound, name) in enumerate(KING_DEGREES):
            fig.add_hline(
                y=bound, line_dash="dot", line_color=NEUTRAL, line_width=1,
                annotation_text=name, annotation_font_size=10,
                annotation_font_color=INK_SECONDARY, annotation_position="top right",
                row=row, col=col,
            )
        for bound, label in panel.cutoffs:
            fig.add_hline(
                y=bound, line_dash="dash", line_color=THRESHOLD, line_width=1.2,
                annotation_text=label, annotation_font_size=10,
                annotation_font_color=INK_SECONDARY, annotation_position="bottom left",
                row=row, col=col,
            )
        fig.add_trace(
            go.Scatter(
                x=panel.ibs0, y=panel.kinship, mode="markers",
                marker={"color": SERIES, "size": 8, "opacity": 0.8,
                        "line": {"width": 0.5, "color": SURFACE}},
                hovertext=panel.hover, hoverinfo="text", showlegend=False,
            ),
            row=row, col=col,
        )
        fig.update_yaxes(range=[low, high], row=row, col=col)
        fig.update_xaxes(range=list(_ibs0_xlim(panel)), row=row, col=col)
    plotly_layout(fig, 340 * rows + 60, hovermode="closest")
    fig.update_xaxes(title_text="IBS0", gridcolor=GRID)
    fig.update_yaxes(title_text="KING kinship", gridcolor=GRID)
    return fig


# --- Registry ----------------------------------------------------------------

PLOTS: List[Plot] = [
    Plot(
        name="sample_missingness",
        title="Sample missingness",
        description=(
            "Distribution of per-sample missing-call rate, with the call-rate "
            "cutoff (--mind). Ticks along the baseline mark samples above it."
        ),
        available=_has_table("callrate_sample_missingness"),
        png=lambda r, plt, p: hist_png(missingness_panels(r), "F_MISS", plt, p),
        plotly=lambda r: hist_plotly(missingness_panels(r), "F_MISS"),
    ),
    Plot(
        name="sex_check",
        title="Sex check",
        description=(
            "X-chromosome F by reported sex (the QC-filtered X pass), with the "
            "female and male cutoffs. Ticks mark samples PLINK flagged; hover "
            "for reported and genetic sex."
        ),
        available=_has_table("sex_check"),
        png=lambda r, plt, p: hist_png(sex_panels(r), "X-chromosome F", plt, p),
        plotly=lambda r: hist_plotly(sex_panels(r), "X-chromosome F"),
    ),
    Plot(
        name="heterozygosity",
        title="Heterozygosity",
        description=(
            "Per-sample heterozygosity F with the keep band this group was "
            "pruned against (these differ by group under --het sd). Ticks mark "
            "outliers."
        ),
        available=_has_table("het_f"),
        png=lambda r, plt, p: hist_png(het_panels(r), "heterozygosity F", plt, p),
        plotly=lambda r: hist_plotly(het_panels(r), "heterozygosity F"),
    ),
    Plot(
        name="relatedness",
        title="Relatedness",
        description=(
            "KING kinship against IBS0 for pairs at or above 3rd degree, with "
            "the degree boundaries and the step's own cutoffs."
        ),
        available=_has_table("related_kinship"),
        png=lambda r, plt, p: kinship_png(kinship_panels(r), plt, p),
        plotly=lambda r: kinship_plotly(kinship_panels(r)),
    ),
]
