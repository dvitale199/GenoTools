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

"""Ancestry plots: reference and cohort PCs, the confusion matrix, label
counts and PC drift.

Principal components only: the report's UMAP embedding has 15-25 dimensions
(the grid search chooses), so a 2D or 3D view of it would show an arbitrary
slice of what the classifier sees (#304).

Every label is drawn with both a colour and a marker shape from
`style.label_styles`, given the run's whole vocabulary (reference labels,
predicted labels and count labels together) so each plot of a run agrees and
nothing here names a label: the vocabulary is user-supplied, and admixture
detection invents `CAH` at runtime.
"""

import re
from pathlib import Path
from typing import Callable, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

from genotools.plotting.plots.base import Plot
from genotools.plotting.report import Report
from genotools.plotting.style import (
    GRID,
    INK,
    INK_SECONDARY,
    SERIES,
    SURFACE,
    THRESHOLD,
    LabelStyle,
    label_styles,
    plotly_layout,
    save_png,
    style_axis,
)

# Cohort rows with no predicted label (a report that has PCs but no
# `ancestry_labels`) are drawn under this name; it is a placeholder, not an
# ancestry label.
UNLABELLED = "unlabelled"

SOURCE_REFERENCE = "reference"
SOURCE_COHORT = "cohort"

# The leading PC pairs, where ancestry structure lives.
PC_PAIRS = [(1, 2), (1, 3), (2, 3)]

_PC_NAME = re.compile(r"^PC(\d+)$")


# --- Frames --------------------------------------------------------------------


def _pc_names(columns) -> List[str]:
    """`PC<n>` column names in numeric order."""
    found = [(int(m.group(1)), c) for c in columns if (m := _PC_NAME.match(str(c)))]
    return [name for _, name in sorted(found)]


def _table_columns(report: Report, key: str) -> List[str]:
    """Column names of a report table without building the frame (cheap)."""
    table = report.raw.get(key)
    return [str(c) for c in table] if isinstance(table, dict) else []


def shared_pcs(report: Report) -> List[str]:
    """PCs present in both the reference and the cohort tables, in order.

    Empty unless the reference table also carries its `label` column.
    """
    ref = _table_columns(report, "ref_pcs")
    cohort = _table_columns(report, "projected_pcs")
    if "label" not in ref:
        return []
    return [pc for pc in _pc_names(ref) if pc in set(cohort)]


def pc_frame(report: Report) -> pd.DataFrame:
    """Reference and cohort samples on one frame, ready to scatter.

    Args:
        report: A loaded report with `ref_pcs` and `projected_pcs`.

    Returns:
        One row per sample with `source` (`reference` / `cohort`), `IID`,
        `label` and the shared `PC` columns. A cohort sample's label is its
        predicted ancestry; a reference sample's is its panel label.
    """
    pcs = shared_pcs(report)
    ref = pd.DataFrame(report.raw["ref_pcs"])
    cohort = pd.DataFrame(report.raw["projected_pcs"])

    labels = report.ancestry_labels
    if labels is not None and "label" in labels.columns:
        keys = [k for k in ("FID", "IID") if k in labels.columns and k in cohort.columns]
        if keys:
            cohort = _merge_labels(cohort, labels, keys)
    if "label" not in cohort.columns:
        cohort["label"] = UNLABELLED
    cohort["label"] = cohort["label"].astype(object).where(cohort["label"].notna(), UNLABELLED)

    frames = []
    for source, frame in ((SOURCE_REFERENCE, ref), (SOURCE_COHORT, cohort)):
        part = frame[pcs + ["label"]].copy()
        part["IID"] = frame["IID"].astype(str) if "IID" in frame.columns else ""
        part["label"] = part["label"].astype(str)
        part.insert(0, "source", source)
        frames.append(part)
    return pd.concat(frames, ignore_index=True)[["source", "IID", "label"] + pcs]


def _merge_labels(cohort: pd.DataFrame, labels: pd.DataFrame, keys: List[str]) -> pd.DataFrame:
    """Attach predicted labels to the cohort PCs by sample ID."""
    lookup = labels[keys + ["label"]].copy()
    cohort = cohort.drop(columns=["label"], errors="ignore").copy()
    for key in keys:
        lookup[key] = lookup[key].astype(str)
        cohort[key] = cohort[key].astype(str)
    return cohort.merge(lookup.drop_duplicates(keys), on=keys, how="left")


def label_vocabulary(report: Report) -> List[str]:
    """Every label the run mentions: panel, predicted and counted.

    Styles are assigned over this whole set, so a label keeps one colour and
    marker in every plot of the run, including a plot that only draws some
    of them.
    """
    labels = set()
    ref = report.raw.get("ref_pcs")
    if isinstance(ref, dict) and "label" in ref:
        labels |= {str(v) for v in pd.DataFrame(ref)["label"]}
    if report.ancestry_labels is not None and "label" in report.ancestry_labels.columns:
        labels |= {str(v) for v in report.ancestry_labels["label"]}
    counts = report.raw.get("ancestry_counts")
    if isinstance(counts, dict) and "label" in counts:
        labels |= {str(v) for v in pd.DataFrame(counts)["label"]}
    cm = report.raw.get("confusion_matrix")
    if isinstance(cm, dict):
        labels |= {str(c) for c in cm}
    return sorted(labels)


def styles_for(report: Report, extra=()) -> Dict[str, LabelStyle]:
    """The run's label styles (colour + marker), covering `extra` as well."""
    return label_styles(list(label_vocabulary(report)) + list(extra))


def counts_frame(report: Report) -> pd.DataFrame:
    """Predicted samples per label: columns `label`, `count`."""
    frame = pd.DataFrame(report.raw["ancestry_counts"])[["label", "count"]].copy()
    frame["label"] = frame["label"].astype(str)
    frame["count"] = frame["count"].astype(int)
    return frame.sort_values(["count", "label"], ascending=[False, True], ignore_index=True)


def confusion_frames(report: Report) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """The reference test-set confusion matrix as counts and row-normalised.

    Rows are the true label, columns the predicted one (sklearn's layout, which
    the report keeps). A row with no samples stays all zero rather than NaN.

    Returns:
        `(counts, normalised)`, both labelled by the report's own labels.
    """
    counts = pd.DataFrame(report.raw["confusion_matrix"]).astype(float)
    counts.index = counts.index.astype(str)
    counts.columns = counts.columns.astype(str)
    totals = counts.sum(axis=1)
    normalised = counts.div(totals.where(totals > 0, np.nan), axis=0).fillna(0.0)
    return counts, normalised


def accuracy_note(report: Report) -> str:
    """Held-out accuracy figures for the confusion matrix title.

    Balanced accuracy lives under `ancestry_fit`, which only a training run has;
    plain `test_accuracy` is on every run that kept its training metrics.
    """
    parts = []
    fit = report.raw.get("ancestry_fit")
    if isinstance(fit, dict) and fit.get("test_balanced_accuracy") is not None:
        parts.append(f"balanced accuracy {100 * float(fit['test_balanced_accuracy']):.1f}%")
    if report.raw.get("test_accuracy") is not None:
        parts.append(f"accuracy {100 * float(report.raw['test_accuracy']):.1f}%")
    return ", ".join(parts)


def drift_frame(report: Report) -> pd.DataFrame:
    """The cohort's per-PC position against the reference, from the report."""
    diagnostics = report.raw.get("ancestry_diagnostics")
    drift = diagnostics.get("pc_drift") if isinstance(diagnostics, dict) else None
    return pd.DataFrame(drift) if drift else pd.DataFrame()


def drift_limits(report: Report) -> Optional[Tuple[float, Tuple[float, float]]]:
    """The thresholds the run warned at: `(max_shift_sd, (ratio_low, ratio_high))`.

    Read from `ancestry_diagnostics.pc_drift_limits`, where the run recorded the
    values its warnings used. None for a report that did not record them: the
    plot then draws no band rather than guess at one.
    """
    diagnostics = report.raw.get("ancestry_diagnostics")
    limits = diagnostics.get("pc_drift_limits") if isinstance(diagnostics, dict) else None
    try:
        shift = float(limits["max_shift_sd"])
        low, high = (float(v) for v in limits["sd_ratio_range"])
    except (TypeError, KeyError, ValueError):
        return None
    return shift, (low, high)


# --- PC scatters ---------------------------------------------------------------


def _pairs(pcs: List[str]) -> List[Tuple[str, str]]:
    return [(f"PC{x}", f"PC{y}") for x, y in PC_PAIRS if f"PC{x}" in pcs and f"PC{y}" in pcs]


def _legend_order(frame: pd.DataFrame) -> List[str]:
    return sorted(frame["label"].unique())


def pairs_png(frame: pd.DataFrame, styles: Dict[str, LabelStyle], plt, out_path: Path) -> Path:
    """PC1/PC2, PC1/PC3 and PC2/PC3: faded reference, solid cohort."""
    pcs = _pc_names(frame.columns)
    pairs = _pairs(pcs)
    fig, axes = plt.subplots(1, len(pairs), figsize=(4.8 * len(pairs) + 1.6, 4.6), squeeze=False)
    fig.patch.set_facecolor(SURFACE)
    for ax, (xcol, ycol) in zip(axes[0], pairs):
        style_axis(ax)
        ax.xaxis.grid(True, color=GRID, linewidth=0.8)
        _scatter_2d(ax, frame, styles, xcol, ycol)
        ax.set_xlabel(xcol, color=INK_SECONDARY, fontsize="small")
        ax.set_ylabel(ycol, color=INK_SECONDARY, fontsize="small")
    handles, names = axes[0][0].get_legend_handles_labels()
    fig.legend(handles, names, loc="center left", bbox_to_anchor=(1.0, 0.5),
               frameon=False, fontsize="small", title="label")
    fig.suptitle("Reference panel (faded) and cohort (solid), by label", x=0.01,
                 ha="left", fontsize="medium", color=INK)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(out_path, dpi=150, facecolor=SURFACE, bbox_inches="tight")
    plt.close(fig)
    return out_path


def _scatter_2d(ax, frame, styles, xcol, ycol) -> None:
    for label in _legend_order(frame):
        style = styles[label]
        ref = frame[(frame["source"] == SOURCE_REFERENCE) & (frame["label"] == label)]
        cohort = frame[(frame["source"] == SOURCE_COHORT) & (frame["label"] == label)]
        if len(ref):
            ax.scatter(ref[xcol], ref[ycol], s=9, alpha=0.22, color=style.color,
                       marker=style.marker, linewidths=0)
        if len(cohort):
            ax.scatter(cohort[xcol], cohort[ycol], s=26, alpha=0.9, color=style.color,
                       marker=style.marker, edgecolors=SURFACE, linewidths=0.4, label=label)
        elif len(ref):
            # Legend entry for a label only the panel has.
            ax.scatter([], [], s=26, color=style.color, marker=style.marker, label=label)


# plotly's scatter3d draws only these marker symbols; 2D scatter has the rest
# (triangles, hexagon, star), which `style.MARKERS` uses.
SCATTER3D_SYMBOLS = [
    "circle", "square", "diamond", "cross", "x",
    "circle-open", "square-open", "diamond-open",
]


def symbols_3d(styles: Dict[str, LabelStyle]) -> Dict[str, str]:
    """Each label's marker for a 3D scatter.

    A label whose style symbol scatter3d can draw keeps it. One it cannot gets
    a supported symbol no other label holds, in sorted order, and repeats once
    those run out. Colour and the hover name still identify the label.
    """
    mapping = {l: st.symbol for l, st in styles.items() if st.symbol in SCATTER3D_SYMBOLS}
    spare = [s for s in SCATTER3D_SYMBOLS if s not in mapping.values()] or SCATTER3D_SYMBOLS
    for i, label in enumerate(l for l in sorted(styles) if l not in mapping):
        mapping[label] = spare[i % len(spare)]
    return mapping


def _plotly_traces(frame, styles, three_d: bool):
    """One reference and one cohort trace per label per panel."""
    traces = []
    symbols = symbols_3d(styles) if three_d else {l: st.symbol for l, st in styles.items()}
    for label in _legend_order(frame):
        style = styles[label]
        for source in (SOURCE_REFERENCE, SOURCE_COHORT):
            part = frame[(frame["source"] == source) & (frame["label"] == label)]
            if part.empty:
                continue
            is_ref = source == SOURCE_REFERENCE
            marker = {
                "symbol": symbols[label],
                "color": style.color,
                "size": (3 if three_d else 5) if is_ref else (4 if three_d else 8),
                "opacity": 0.2 if is_ref else 0.9,
                "line": {"width": 0 if is_ref else 0.5, "color": SURFACE},
            }
            traces.append((label, source, part, marker))
    return traces


def _hover(part: pd.DataFrame, label: str, source: str) -> List[str]:
    return [f"{iid}<br>{label} ({source})" for iid in part["IID"]]


def pairs_plotly(frame: pd.DataFrame, styles: Dict[str, LabelStyle]):
    """The same three panels with the sample ID and label on hover."""
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots

    pairs = _pairs(_pc_names(frame.columns))
    fig = make_subplots(rows=1, cols=len(pairs), subplot_titles=[f"{x} vs {y}" for x, y in pairs])
    for col, (xcol, ycol) in enumerate(pairs, start=1):
        for label, source, part, marker in _plotly_traces(frame, styles, False):
            fig.add_trace(
                go.Scatter(
                    x=part[xcol], y=part[ycol], mode="markers", marker=marker,
                    name=label if source == SOURCE_COHORT else f"{label} (reference)",
                    legendgroup=f"{label}-{source}",
                    showlegend=(col == 1),
                    text=_hover(part, label, source), hoverinfo="text+x+y",
                ),
                row=1, col=col,
            )
        fig.update_xaxes(title_text=xcol, gridcolor=GRID, row=1, col=col)
        fig.update_yaxes(title_text=ycol, gridcolor=GRID, row=1, col=col)
    return plotly_layout(fig, 460, legend={"title": {"text": "label"}})


def cube_png(frame: pd.DataFrame, styles: Dict[str, LabelStyle], plt, out_path: Path) -> Path:
    """A static 3D view for the PNG output; the HTML version rotates."""
    fig = plt.figure(figsize=(7.2, 6.0))
    fig.patch.set_facecolor(SURFACE)
    ax = fig.add_subplot(projection="3d")
    ax.set_facecolor(SURFACE)
    for label in _legend_order(frame):
        style = styles[label]
        has_cohort = bool(((frame["source"] == SOURCE_COHORT) & (frame["label"] == label)).any())
        for source in (SOURCE_REFERENCE, SOURCE_COHORT):
            part = frame[(frame["source"] == source) & (frame["label"] == label)]
            if part.empty:
                continue
            is_ref = source == SOURCE_REFERENCE
            ax.scatter(
                part["PC1"], part["PC2"], part["PC3"],
                s=7 if is_ref else 22, alpha=0.2 if is_ref else 0.9,
                color=style.color, marker=style.marker,
                linewidths=0 if is_ref else 0.3, edgecolors=SURFACE,
                # A label only the panel has still needs a legend entry.
                label=None if is_ref and has_cohort else label,
            )
    ax.set_xlabel("PC1", color=INK_SECONDARY)
    ax.set_ylabel("PC2", color=INK_SECONDARY)
    ax.set_zlabel("PC3", color=INK_SECONDARY)
    ax.tick_params(colors=INK_SECONDARY, labelsize="x-small")
    handles, names = ax.get_legend_handles_labels()
    fig.legend(handles, names, loc="center left", bbox_to_anchor=(0.98, 0.5),
               frameon=False, fontsize="small", title="label")
    ax.set_title("PC1-3: reference (faded) and cohort (solid)", loc="left",
                 fontsize="medium", color=INK)
    fig.savefig(out_path, dpi=150, facecolor=SURFACE, bbox_inches="tight")
    plt.close(fig)
    return out_path


def cube_plotly(frame: pd.DataFrame, styles: Dict[str, LabelStyle]):
    """Rotatable PC1-3 with the sample ID and label on hover."""
    import plotly.graph_objects as go

    fig = go.Figure()
    for label, source, part, marker in _plotly_traces(frame, styles, True):
        fig.add_trace(
            go.Scatter3d(
                x=part["PC1"], y=part["PC2"], z=part["PC3"], mode="markers",
                marker={k: v for k, v in marker.items() if k != "line"} | {"line": {"width": 0}},
                name=label if source == SOURCE_COHORT else f"{label} (reference)",
                legendgroup=f"{label}-{source}",
                text=_hover(part, label, source), hoverinfo="text+x+y+z",
            )
        )
    plotly_layout(fig, 620, legend={"title": {"text": "label"}})
    fig.update_layout(scene={"xaxis_title": "PC1", "yaxis_title": "PC2", "zaxis_title": "PC3"})
    return fig


# --- Confusion matrix ----------------------------------------------------------


def _title_with(base: str, note: str) -> str:
    return f"{base}\n{note}" if note else base


def confusion_png(counts: pd.DataFrame, normalised: pd.DataFrame, note: str, plt, out_path: Path) -> Path:
    """Row-normalised confusion matrix, each cell annotated with its fraction."""
    from matplotlib.colors import LinearSegmentedColormap

    cmap = LinearSegmentedColormap.from_list("share", [SURFACE, SERIES])
    n = len(normalised)
    side = max(4.5, 0.62 * n + 2.4)
    fig, ax = plt.subplots(figsize=(side, side * 0.9))
    fig.patch.set_facecolor(SURFACE)
    ax.imshow(normalised.values, cmap=cmap, vmin=0.0, vmax=1.0, aspect="auto")
    for r in range(n):
        for c in range(normalised.shape[1]):
            value = normalised.values[r, c]
            ax.text(c, r, f"{value:.2f}" if value else "0", ha="center", va="center",
                    fontsize="x-small", color="white" if value > 0.55 else INK)
    ax.set_xticks(range(normalised.shape[1]))
    ax.set_xticklabels(normalised.columns, rotation=45, ha="right", fontsize="small", color=INK_SECONDARY)
    ax.set_yticks(range(n))
    ax.set_yticklabels(normalised.index, fontsize="small", color=INK_SECONDARY)
    ax.set_xlabel("predicted", color=INK_SECONDARY, fontsize="small")
    ax.set_ylabel("true (reference test set)", color=INK_SECONDARY, fontsize="small")
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_title(_title_with("Confusion matrix: share of each true label", note),
                 loc="left", fontsize="medium", color=INK)
    return save_png(fig, plt, out_path)


def confusion_plotly(counts: pd.DataFrame, normalised: pd.DataFrame, note: str):
    """The same matrix with the raw count behind each cell on hover."""
    import plotly.graph_objects as go

    hover = [
        [f"true {t}, predicted {p}<br>{normalised.at[t, p]:.1%} of {t} ({counts.at[t, p]:.0f} samples)"
         for p in normalised.columns]
        for t in normalised.index
    ]
    fig = go.Figure(
        go.Heatmap(
            z=normalised.values, x=list(normalised.columns), y=list(normalised.index),
            zmin=0, zmax=1, colorscale=[[0, SURFACE], [1, SERIES]],
            text=[[f"{v:.2f}" if v else "0" for v in row] for row in normalised.values],
            texttemplate="%{text}", hovertext=hover, hoverinfo="text",
            xgap=1, ygap=1, colorbar={"title": "share"},
        )
    )
    plotly_layout(fig, 140 + 46 * len(normalised), title={"text": _title_with("Row-normalised", note).replace("\n", "<br>"), "font": {"size": 14}})
    fig.update_xaxes(title_text="predicted")
    fig.update_yaxes(title_text="true (reference test set)", autorange="reversed")
    return fig


# --- Label counts --------------------------------------------------------------


def counts_png(counts: pd.DataFrame, styles: Dict[str, LabelStyle], plt, out_path: Path) -> Path:
    """Predicted samples per label, in the label's own colour."""
    fig, ax = plt.subplots(figsize=(max(4.6, 0.8 * len(counts) + 2), 3.8))
    fig.patch.set_facecolor(SURFACE)
    style_axis(ax)
    x = range(len(counts))
    ax.bar(x, counts["count"], color=[styles[l].color for l in counts["label"]], width=0.7)
    for i, n in enumerate(counts["count"]):
        ax.annotate(f"{n:,}", (i, n), xytext=(0, 2), textcoords="offset points",
                    ha="center", va="bottom", fontsize="x-small", color=INK_SECONDARY)
    ax.set_xticks(list(x))
    ax.set_xticklabels(counts["label"], rotation=35, ha="right")
    ax.set_ylabel("samples", color=INK_SECONDARY, fontsize="small")
    ax.set_title(f"Predicted ancestry: {int(counts['count'].sum()):,} samples",
                 loc="left", fontsize="medium", color=INK)
    return save_png(fig, plt, out_path)


def counts_plotly(counts: pd.DataFrame, styles: Dict[str, LabelStyle]):
    import plotly.graph_objects as go

    fig = go.Figure(
        go.Bar(
            x=counts["label"], y=counts["count"],
            marker_color=[styles[l].color for l in counts["label"]],
            text=[f"{n:,}" for n in counts["count"]], textposition="outside",
            hovertemplate="%{x}: %{y:,} samples<extra></extra>",
        )
    )
    plotly_layout(fig, 380, showlegend=False)
    fig.update_yaxes(title_text="samples", gridcolor=GRID)
    return fig


# --- PC drift ------------------------------------------------------------------


def _outside(values: pd.Series, low: float, high: float) -> pd.Series:
    return (values < low) | (values > high)


def drift_png(drift: pd.DataFrame, limits, plt, out_path: Path) -> Path:
    """Where the cohort sits against the reference, per PC, with the warn band.

    Left: the cohort centre in reference standard deviations. Right: the
    cohort's spread as a ratio of the reference's (log scale). The shaded band
    is what the run accepted without a warning; a point outside it is one the
    run warned about. Without recorded limits there is no band.
    """
    fig, (ax_shift, ax_ratio) = plt.subplots(1, 2, figsize=(10.4, 4.0))
    fig.patch.set_facecolor(SURFACE)
    x = np.arange(len(drift))
    panels = (
        (ax_shift, "shift_sd", "cohort centre minus reference centre (reference SD)",
         None if limits is None else (-limits[0], limits[0]), False),
        (ax_ratio, "sd_ratio", "cohort SD / reference SD",
         None if limits is None else limits[1], True),
    )
    for ax, column, ylabel, band, log in panels:
        style_axis(ax)
        values = pd.to_numeric(drift[column], errors="coerce")
        if log:
            ax.set_yscale("log")
        if band is not None:
            ax.axhspan(band[0], band[1], color=THRESHOLD, alpha=0.10, linewidth=0)
            for edge in band:
                ax.axhline(edge, color=THRESHOLD, linewidth=1, linestyle="--")
            flagged = _outside(values, *band)
        else:
            flagged = pd.Series(False, index=values.index)
        ax.scatter(x[~flagged.values], values[~flagged], s=40, color=SERIES, zorder=3)
        ax.scatter(x[flagged.values], values[flagged], s=70, color=SERIES,
                   edgecolors=INK, linewidths=1.5, zorder=3)
        for xi, v in zip(x[flagged.values], values[flagged]):
            ax.annotate(f"{v:.3g}", (xi, v), xytext=(0, 7), textcoords="offset points",
                        ha="center", fontsize="x-small", color=INK)
        ax.set_xticks(list(x))
        ax.set_xticklabels(drift["pc"], fontsize="small")
        ax.set_ylabel(ylabel, color=INK_SECONDARY, fontsize="small")
        if band is None:
            ax.set_title("warning limits not recorded in this report", loc="left",
                         fontsize="small", color=INK_SECONDARY)
        else:
            ax.set_title(f"warns outside {band[0]:g} to {band[1]:g}", loc="left",
                         fontsize="small", color=INK_SECONDARY)
    ax_shift.axhline(0, color=GRID, linewidth=1, zorder=1)
    fig.suptitle("Cohort position in reference PC space", x=0.01, ha="left",
                 fontsize="medium", color=INK)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(out_path, dpi=150, facecolor=SURFACE)
    plt.close(fig)
    return out_path


def drift_plotly(drift: pd.DataFrame, limits):
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots

    fig = make_subplots(rows=1, cols=2, subplot_titles=(
        "Centre shift (reference SD)", "Spread ratio (cohort / reference)"))
    specs = (
        (1, "shift_sd", None if limits is None else (-limits[0], limits[0])),
        (2, "sd_ratio", None if limits is None else limits[1]),
    )
    for col, column, band in specs:
        values = pd.to_numeric(drift[column], errors="coerce")
        if band is not None:
            fig.add_hrect(y0=band[0], y1=band[1], fillcolor=THRESHOLD, opacity=0.10,
                          line_width=0, row=1, col=col,
                          exclude_empty_subplots=False)
            for edge in band:
                fig.add_hline(y=edge, line={"color": THRESHOLD, "dash": "dash", "width": 1},
                              row=1, col=col, exclude_empty_subplots=False)
            flagged = _outside(values, *band)
        else:
            flagged = pd.Series(False, index=values.index)
        fig.add_trace(
            go.Scatter(
                x=drift["pc"], y=values, mode="markers", showlegend=False,
                marker={"color": SERIES, "size": [12 if f else 8 for f in flagged],
                        "line": {"color": INK, "width": [2 if f else 0 for f in flagged]}},
                hovertemplate="%{x}: %{y:.3g}<extra></extra>",
            ),
            row=1, col=col,
        )
    fig.update_yaxes(gridcolor=GRID)
    fig.update_yaxes(type="log", row=1, col=2)
    note = "" if limits is not None else "warning limits not recorded in this report"
    plotly_layout(fig, 400, title={"text": note, "font": {"size": 12}})
    return fig


# --- Registry ------------------------------------------------------------------


def _has_pcs(n: int) -> Callable[[Report], bool]:
    return lambda report: len(shared_pcs(report)) >= n


def _has_table(key: str, *columns: str) -> Callable[[Report], bool]:
    return lambda report: all(c in _table_columns(report, key) for c in columns) and bool(report.raw.get(key))


def _has_drift(report: Report) -> bool:
    diagnostics = report.raw.get("ancestry_diagnostics")
    return isinstance(diagnostics, dict) and bool(diagnostics.get("pc_drift"))


def _pc_args(report: Report):
    return pc_frame(report), styles_for(report)


PLOTS: List[Plot] = [
    Plot(
        name="ancestry_pc_pairs",
        title="Ancestry: PC pairs",
        description=(
            "PC1/PC2, PC1/PC3 and PC2/PC3. The reference panel is faded and "
            "coloured by its label; the cohort is solid and coloured by its "
            "predicted label. Every label has its own colour and marker shape. "
            "Hover a point for its sample ID."
        ),
        available=_has_pcs(2),
        png=lambda r, plt, p: pairs_png(*_pc_args(r), plt, p),
        plotly=lambda r: pairs_plotly(*_pc_args(r)),
    ),
    Plot(
        name="ancestry_pc_3d",
        title="Ancestry: PC1-3, rotatable",
        description=(
            "The same encoding in three dimensions. Drag to rotate; hover a point "
            "for its sample ID and label. The PNG is a single fixed view."
        ),
        available=_has_pcs(3),
        png=lambda r, plt, p: cube_png(*_pc_args(r), plt, p),
        plotly=lambda r: cube_plotly(*_pc_args(r)),
    ),
    Plot(
        name="ancestry_confusion_matrix",
        title="Ancestry: confusion matrix",
        description=(
            "How the classifier labelled the held-out reference samples: each "
            "row is a true label and shows the share predicted as each label. "
            "A property of the model on the reference panel, not of the cohort."
        ),
        available=_has_table("confusion_matrix"),
        png=lambda r, plt, p: confusion_png(*confusion_frames(r), accuracy_note(r), plt, p),
        plotly=lambda r: confusion_plotly(*confusion_frames(r), accuracy_note(r)),
    ),
    Plot(
        name="ancestry_label_counts",
        title="Ancestry: predicted labels",
        description="How many cohort samples were predicted as each ancestry.",
        available=_has_table("ancestry_counts", "label", "count"),
        png=lambda r, plt, p: counts_png(counts_frame(r), styles_for(r), plt, p),
        plotly=lambda r: counts_plotly(counts_frame(r), styles_for(r)),
    ),
    Plot(
        name="ancestry_pc_drift",
        title="Ancestry: PC drift",
        description=(
            "Where the cohort landed in the reference's PC space. The shaded "
            "band is what the run accepted without a warning, drawn from the "
            "limits it recorded; a point outside it is one the run warned about."
        ),
        available=_has_drift,
        png=lambda r, plt, p: drift_png(drift_frame(r), drift_limits(r), plt, p),
        plotly=lambda r: drift_plotly(drift_frame(r), drift_limits(r)),
    ),
]
