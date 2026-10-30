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

"""Variant QC plots: missingness, and QQ plots for the HWE, case-control
and haplotype tests.

Everything here is drawn from the binned tables the QC steps keep (see
`plot_data` in docs/json_output_overview.md), never from per-variant data, so
a plot costs the same on a million-variant array as on a thousand-variant one.
A QQ plot is rebuilt from a -log10(p) histogram plus the exact rows of the
tests that failed: the bulk is quantised to the histogram's 0.05-wide bins,
and the tail, where a QQ plot is read, is exact.
"""

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional

import numpy as np
import pandas as pd

from genotools.core.logging import get_logger
from genotools.plotting.plots.base import Plot
from genotools.plotting.report import Report, parameter, plot_tables
from genotools.plotting.style import (
    GRID,
    INK,
    INK_SECONDARY,
    SERIES,
    SURFACE,
    THRESHOLD,
    grid,
    plotly_layout,
    save_png,
    style_axis,
)

logger = get_logger(__name__)

# A p of exactly 0 is -log10(p) = inf, and the histogram's last bin is open
# ("50 and above"). Neither can be placed on an axis, so they are drawn here and
# flagged as lower bounds.
OPEN_BIN_FLOOR = 50.0


def _float(value: Any) -> Optional[float]:
    """A parameter as a float, or None when it is absent or not numeric."""
    if value is None or isinstance(value, bool):
        return None
    try:
        out = float(value)
    except (TypeError, ValueError):
        return None
    return None if np.isnan(out) else out


def _threshold(report: Report, step: str, name: str, ancestry: str) -> Optional[float]:
    """A step's threshold for one group, or None (and so no line) if unrecorded."""
    return _float(parameter(report, step, name, ancestry))


def _usable(tables: Dict[str, pd.DataFrame]) -> bool:
    return any(not frame.empty for frame in tables.values())


def removed_count(report: Report, step: str, ancestry: str) -> Optional[int]:
    """Variants a step removed from one group, as the step reported them.

    Read from `report.qc`, not counted from a plot table: a step can remove
    more than the table lists (haplotype also removes each significant test's
    flanking variants) and the tables are capped. A step with several metrics
    is summed. None when the report has no row for it.
    """
    qc = report.qc
    rows = qc[(qc["step"] == step) & (qc["ancestry"].astype(str) == str(ancestry))]
    if rows.empty:
        return None
    return int(pd.to_numeric(rows["pruned_count"], errors="coerce").fillna(0).sum())


def _removed_label(n: Optional[int], noun: str = "variants") -> str:
    return "" if n is None else f"{n:,} {noun} removed"


# --- Variant missingness -----------------------------------------------------


@dataclass(frozen=True)
class MissingnessPanel:
    ancestry: str
    bins: pd.DataFrame
    threshold: Optional[float]
    removed: Optional[int]


def missingness_panels(report: Report) -> List[MissingnessPanel]:
    """One panel per group that kept a variant-missingness histogram."""
    bins = plot_tables(report, "geno_missingness_bins")
    panels = []
    for ancestry, frame in bins.items():
        if frame.empty:
            continue
        panels.append(
            MissingnessPanel(
                ancestry=ancestry,
                bins=frame,
                threshold=_threshold(report, "geno_prune", "geno", ancestry),
                removed=removed_count(report, "geno_prune", ancestry),
            )
        )
    return panels


def _missingness_xmax(panel: MissingnessPanel) -> float:
    """Right edge: the last occupied bin, but never short of the threshold."""
    occupied = panel.bins.loc[panel.bins["count"] > 0, "bin_upper"]
    xmax = float(occupied.max()) if len(occupied) else float(panel.bins["bin_upper"].max())
    if panel.threshold is not None:
        xmax = max(xmax, panel.threshold * 1.15)
    return min(xmax, float(panel.bins["bin_upper"].max()))


def _missingness_title(panel: MissingnessPanel) -> str:
    removed = _removed_label(panel.removed)
    return f"{panel.ancestry}: {removed}" if removed else panel.ancestry


def missingness_png(report: Report, plt, out_path: Path) -> Path:
    """Variant call-rate histogram per group, count on a log axis.

    Log, because the failing tail is a few variants per bin against a peak of
    thousands, and the tail is what the threshold is about.
    """
    panels = missingness_panels(report)
    rows, cols = grid(len(panels))
    fig, axes = plt.subplots(rows, cols, figsize=(4.6 * cols, 3.4 * rows), squeeze=False)
    fig.patch.set_facecolor(SURFACE)

    for ax, panel in zip(axes.flat, panels):
        style_axis(ax)
        b = panel.bins
        ax.bar(
            b["bin_lower"], b["count"], width=b["bin_upper"] - b["bin_lower"],
            align="edge", color=SERIES,
        )
        ax.set_yscale("log")
        ax.set_ylim(0.8, max(float(b["count"].max()), 1.0) * 2)
        ax.set_xlim(0, _missingness_xmax(panel))
        if panel.threshold is not None:
            ax.axvline(panel.threshold, color=THRESHOLD, linestyle="--", linewidth=1.2)
            ax.annotate(
                f"--geno {panel.threshold:g}", (panel.threshold, 1), xycoords=("data", "axes fraction"),
                xytext=(4, -4), textcoords="offset points", va="top", fontsize="x-small",
                color=THRESHOLD,
            )
        ax.set_title(_missingness_title(panel), loc="left", fontsize="medium", color=INK)
        ax.set_xlabel("variant missingness (F_MISS)", color=INK_SECONDARY, fontsize="small")
        ax.set_ylabel("variants (log scale)", color=INK_SECONDARY, fontsize="small")

    for ax in list(axes.flat)[len(panels):]:
        ax.set_visible(False)
    return save_png(fig, plt, out_path)


def missingness_plotly(report: Report):
    """The same small multiples, with the bin range and count on hover."""
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots

    panels = missingness_panels(report)
    rows, cols = grid(len(panels))
    fig = make_subplots(rows=rows, cols=cols, subplot_titles=[_missingness_title(p) for p in panels])
    for i, panel in enumerate(panels):
        row, col = i // cols + 1, i % cols + 1
        b = panel.bins
        fig.add_trace(
            go.Bar(
                x=b["bin_lower"], y=b["count"], width=b["bin_upper"] - b["bin_lower"],
                offset=0, marker_color=SERIES, name=panel.ancestry,
                customdata=b["bin_upper"],
                hovertemplate="F_MISS %{x:.3f} to %{customdata:.3f}: %{y:,} variants<extra>"
                + panel.ancestry + "</extra>",
            ),
            row=row, col=col,
        )
        fig.update_xaxes(range=[0, _missingness_xmax(panel)], title_text="variant missingness (F_MISS)", row=row, col=col)
        fig.update_yaxes(type="log", gridcolor=GRID, title_text="variants (log scale)", row=row, col=col)
        if panel.threshold is not None:
            fig.add_vline(
                x=panel.threshold, line={"color": THRESHOLD, "dash": "dash", "width": 1.2},
                annotation_text=f"--geno {panel.threshold:g}", annotation_position="top right",
                row=row, col=col,
            )
    plotly_layout(fig, 300 * rows + 60, showlegend=False)
    return fig


# --- QQ plots ------------------------------------------------------------------

QQ_COLUMNS = [
    "kind", "observed", "expected", "count", "rank_first", "rank_last",
    "label", "bin_lower", "bin_upper", "lower_bound", "n",
]


def neg_log10(p: pd.Series) -> pd.Series:
    """-log10(p), with p == 0 as +inf."""
    p = pd.to_numeric(p, errors="coerce")
    with np.errstate(divide="ignore"):
        return pd.Series(-np.log10(p.to_numpy(dtype=float)), index=p.index)


def qq_points(bins: pd.DataFrame, failing: Optional[pd.DataFrame] = None) -> pd.DataFrame:
    """Observed vs expected -log10(p) for a QQ plot, from binned p-values.

    Rank 1 is the most significant test. With `n` tests in all, the test at
    rank r has expected -log10((r - 0.5) / n).

    Args:
        bins: A `*_pvalue_bins` table: `bin_lower`, `bin_upper` (-log10(p), the
            last bin open with `bin_upper` = inf) and `count`. All the tests
            must be counted here: `n` is the sum of the counts.
        failing: The exact tests that failed, with `P` and optionally `label`
            (hover text). They must be the most extreme tests, which is how
            the steps write them (sorted by p, a prefix of the ranking, capped).
            They take ranks 1..m with their exact p, and the bins keep only the
            ranks after them, so the tail is exact and no test is drawn twice.

    Returns:
        One row per point: `kind` (`exact` or `bin`), `observed`, `expected`,
        `count` (tests the point stands for), `rank_first`/`rank_last`, `label`,
        the bin's edges (NaN for exact rows), `n`, and `lower_bound`, True where
        `observed` is a floor, not a value: an open bin (the true value is
        beyond it) or an exact p of 0 (infinite, drawn at the highest finite
        value). A bin is drawn at the middle of its ranks, at its midpoint.
    """
    counts = bins["count"].astype(np.int64)
    n = int(counts.sum())
    if n == 0:
        return pd.DataFrame(columns=QQ_COLUMNS)

    def expected(rank):
        return -np.log10((np.asarray(rank, dtype=float) - 0.5) / n)

    # Exact rows: ranks 1..m. A failing table longer than n is inconsistent
    # with the bins; trust the bins and keep the most extreme n rows.
    exact = None
    m = 0
    if failing is not None and len(failing):
        scores = neg_log10(failing["P"]).dropna()
        order = scores.sort_values(ascending=False, kind="stable").index[:n]
        m = len(order)
        ranks = np.arange(1, m + 1)
        exact = pd.DataFrame({
            "kind": "exact",
            "observed": scores.loc[order].to_numpy(),
            "expected": expected(ranks),
            "count": 1,
            "rank_first": ranks,
            "rank_last": ranks,
            "label": (
                failing.loc[order, "label"].astype(str).to_numpy()
                if "label" in failing else ""
            ),
            "bin_lower": np.nan,
            "bin_upper": np.nan,
            "lower_bound": False,
            "n": n,
        })

    # Bins, most significant first; keep only the ranks the exact rows left.
    b = bins.assign(count=counts).sort_values("bin_lower", ascending=False, kind="stable")
    last = b["count"].cumsum().to_numpy()
    first = last - b["count"].to_numpy() + 1
    first = np.maximum(first, m + 1)
    keep = (b["count"].to_numpy() > 0) & (first <= last)
    b, first, last = b[keep], first[keep], last[keep]
    mid_rank = (first + last) / 2
    open_bin = ~np.isfinite(b["bin_upper"].to_numpy(dtype=float))
    midpoint = np.where(
        open_bin, b["bin_lower"].to_numpy(dtype=float),
        (b["bin_lower"].to_numpy(dtype=float) + b["bin_upper"].to_numpy(dtype=float)) / 2,
    )
    binned = pd.DataFrame({
        "kind": "bin",
        "observed": midpoint,
        "expected": expected(mid_rank),
        "count": last - first + 1,
        "rank_first": first,
        "rank_last": last,
        "label": "",
        "bin_lower": b["bin_lower"].to_numpy(dtype=float),
        "bin_upper": b["bin_upper"].to_numpy(dtype=float),
        "lower_bound": open_bin,
        "n": n,
    })

    frames = [f for f in (exact, binned) if f is not None and len(f)]
    points = pd.concat(frames, ignore_index=True)[QQ_COLUMNS] if frames else binned
    # p == 0: infinite on the axis. Pin it to the highest finite value.
    infinite = ~np.isfinite(points["observed"].to_numpy(dtype=float))
    if infinite.any():
        finite = points.loc[~infinite, "observed"]
        points.loc[infinite, "observed"] = max(float(finite.max()) if len(finite) else 0.0, OPEN_BIN_FLOOR)
        points.loc[infinite, "lower_bound"] = True
    return points.astype({"lower_bound": bool})


def hwe_filter_test(bins: pd.DataFrame) -> Optional[str]:
    """Which of the HWE table's tests the filter actually applied.

    PLINK 1.9's `--hwe` tests controls only on case/control data, so that is
    `UNAFF`; with no phenotype it tests everyone (`ALL`). This is the choice
    `hwe_plot_tables` makes when it writes `hwe_failing.tsv` - UNAFF if present,
    else the first test - so the QQ's bins and its exact tail describe the same
    test. We follow that rather than the `filter_controls` parameter: with
    `--filter-controls` PLINK drops the cases before testing, so the only test
    present is the one the filter used, and "UNAFF if present" still lands on it.
    """
    if "test" not in bins.columns:
        return None
    tests = set(bins["test"].astype(str))
    if "UNAFF" in tests:
        return "UNAFF"
    return "ALL" if "ALL" in tests else (sorted(tests)[0] if tests else None)


@dataclass(frozen=True)
class QQSpec:
    """Where one step's QQ data and threshold live in the report."""

    name: str
    title: str
    step: str
    threshold_parameter: str
    bins_table: str
    failing_table: str
    noun: str  # what one failing row is
    id_columns: List[str]


QQ_SPECS: Dict[str, QQSpec] = {
    "hwe": QQSpec("hwe", "HWE", "hwe_prune", "hwe_threshold", "hwe_pvalue_bins",
                  "hwe_failing", "variants", ["SNP"]),
    "case_control": QQSpec("case_control", "Case-control missingness",
                           "case_control_missingness_prune", "p_threshold",
                           "case_control_pvalue_bins", "case_control_failing",
                           "variants", ["SNP"]),
    "haplotype": QQSpec("haplotype", "Haplotype missingness", "haplotype_prune",
                        "p_threshold", "haplotype_pvalue_bins", "haplotype_failing",
                        "tests", ["SNP", "HAPLOTYPE"]),
}


@dataclass(frozen=True)
class QQPanel:
    ancestry: str
    points: pd.DataFrame
    threshold: Optional[float]  # p, not -log10
    removed: Optional[int]  # variants the step removed, from report.qc
    exact: Optional[int]  # failing rows drawn as exact points
    test: Optional[str]  # HWE only


def _failing_for_qq(spec: QQSpec, failing: Optional[pd.DataFrame], test: Optional[str]) -> Optional[pd.DataFrame]:
    if failing is None or failing.empty:
        return failing
    if test is not None and "test" in failing.columns:
        failing = failing[failing["test"].astype(str) == test]
    ids = [c for c in spec.id_columns if c in failing.columns]
    label = failing[ids].astype(str).agg(" ".join, axis=1) if ids else ""
    return failing.assign(label=label)


def qq_panels(report: Report, spec: QQSpec) -> List[QQPanel]:
    """One panel per group that kept this step's p-value histogram."""
    all_bins = plot_tables(report, spec.bins_table)
    all_failing = plot_tables(report, spec.failing_table)
    panels = []
    for ancestry, bins in all_bins.items():
        if bins.empty:
            continue
        test = hwe_filter_test(bins) if spec.name == "hwe" else None
        if test is not None:
            bins = bins[bins["test"].astype(str) == test]
        failing = _failing_for_qq(spec, all_failing.get(ancestry), test)
        panels.append(
            QQPanel(
                ancestry=ancestry,
                points=qq_points(bins, failing),
                threshold=_threshold(report, spec.step, spec.threshold_parameter, ancestry),
                removed=removed_count(report, spec.step, ancestry),
                exact=None if failing is None else len(failing),
                test=test,
            )
        )
    return panels


def _threshold_y(panel: QQPanel) -> Optional[float]:
    if panel.threshold is None or panel.threshold <= 0:
        return None
    return float(-np.log10(panel.threshold))


def _qq_title(panel: QQPanel, spec: QQSpec) -> str:
    detail = []
    if panel.test:
        detail.append(f"{panel.test} test")
    if spec.noun == "tests" and panel.exact is not None:
        # The step also removes each significant test's flanking variants, so
        # the tests and the variants removed are different numbers.
        detail.append(f"{panel.exact:,} significant tests")
    removed = _removed_label(panel.removed)
    if removed:
        detail.append(removed)
    return f"{panel.ancestry}: " + ", ".join(detail) if detail else panel.ancestry


def _qq_limits(panel: QQPanel):
    pts = panel.points
    xmax = float(pts["expected"].max()) * 1.05 if len(pts) else 1.0
    ymax = float(pts["observed"].max()) if len(pts) else 1.0
    ty = _threshold_y(panel)
    if ty is not None:
        ymax = max(ymax, ty)
    return xmax, ymax * 1.05


def qq_png(report: Report, spec: QQSpec, plt, out_path: Path) -> Path:
    """QQ small multiples: binned bulk, exact tail, y = x and the p threshold."""
    panels = qq_panels(report, spec)
    rows, cols = grid(len(panels))
    fig, axes = plt.subplots(rows, cols, figsize=(4.6 * cols, 4.0 * rows), squeeze=False)
    fig.patch.set_facecolor(SURFACE)

    for k, (ax, panel) in enumerate(zip(axes.flat, panels)):
        style_axis(ax)
        ax.xaxis.grid(True, color=GRID, linewidth=0.8)
        pts = panel.points
        xmax, ymax = _qq_limits(panel)
        ax.plot([0, xmax], [0, xmax], color=THRESHOLD, linestyle=":", linewidth=1, label="y = x")
        ty = _threshold_y(panel)
        if ty is not None:
            ax.axhline(ty, color=THRESHOLD, linestyle="--", linewidth=1.2, label=f"p = {panel.threshold:g}")
        binned = pts[(pts["kind"] == "bin") & ~pts["lower_bound"]]
        floor = pts[pts["lower_bound"]]
        exact = pts[(pts["kind"] == "exact") & ~pts["lower_bound"]]
        ax.scatter(binned["expected"], binned["observed"], s=9, color=SERIES, alpha=0.7,
                   linewidths=0, label="binned")
        ax.scatter(exact["expected"], exact["observed"], s=14, color=SERIES, marker="D",
                   linewidths=0, label="exact")
        ax.scatter(floor["expected"], floor["observed"], s=22, facecolors="none",
                   edgecolors=SERIES, marker="^", linewidths=1, label="at least (open bin or p = 0)")
        ax.set_xlim(0, xmax)
        ax.set_ylim(0, ymax)
        ax.set_title(_qq_title(panel, spec), loc="left", fontsize="medium", color=INK)
        ax.set_xlabel("expected -log10(p)", color=INK_SECONDARY, fontsize="small")
        ax.set_ylabel("observed -log10(p)", color=INK_SECONDARY, fontsize="small")
        if k == 0:
            ax.legend(fontsize="x-small", frameon=False, loc="best")

    for ax in list(axes.flat)[len(panels):]:
        ax.set_visible(False)
    return save_png(fig, plt, out_path)


def qq_plotly(report: Report, spec: QQSpec):
    """The same QQ small multiples; exact points name their variant on hover."""
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots

    panels = qq_panels(report, spec)
    rows, cols = grid(len(panels))
    fig = make_subplots(rows=rows, cols=cols, subplot_titles=[_qq_title(p, spec) for p in panels])
    for i, panel in enumerate(panels):
        row, col = i // cols + 1, i % cols + 1
        first = i == 0
        pts = panel.points
        xmax, ymax = _qq_limits(panel)
        fig.add_trace(
            go.Scatter(x=[0, xmax], y=[0, xmax], mode="lines", name="y = x", legendgroup="ref",
                       line={"color": THRESHOLD, "dash": "dot", "width": 1}, hoverinfo="skip",
                       showlegend=first),
            row=row, col=col,
        )
        binned = pts[(pts["kind"] == "bin") & ~pts["lower_bound"]]
        exact = pts[(pts["kind"] == "exact") & ~pts["lower_bound"]]
        floor = pts[pts["lower_bound"]]
        fig.add_trace(
            go.Scatter(
                x=binned["expected"], y=binned["observed"], mode="markers", name="binned",
                legendgroup="bin", showlegend=first,
                marker={"color": SERIES, "size": 5, "opacity": 0.7},
                customdata=np.stack([binned["bin_lower"], binned["bin_upper"], binned["count"]], axis=-1),
                hovertemplate="%{customdata[2]:,} tests with -log10(p) in [%{customdata[0]:.2f}, "
                "%{customdata[1]:.2f})<br>expected %{x:.2f}<extra>" + panel.ancestry + "</extra>",
            ),
            row=row, col=col,
        )
        fig.add_trace(
            go.Scatter(
                x=exact["expected"], y=exact["observed"], mode="markers", name="exact",
                legendgroup="exact", showlegend=first,
                marker={"color": SERIES, "size": 6, "symbol": "diamond"},
                text=exact["label"],
                hovertemplate="%{text}<br>observed %{y:.2f}, expected %{x:.2f}<extra>"
                + panel.ancestry + "</extra>",
            ),
            row=row, col=col,
        )
        fig.add_trace(
            go.Scatter(
                x=floor["expected"], y=floor["observed"], mode="markers",
                name="at least (open bin or p = 0)", legendgroup="floor", showlegend=first,
                marker={"color": SERIES, "size": 8, "symbol": "triangle-up-open"},
                text=np.where(floor["label"] != "", floor["label"], "open bin"),
                customdata=floor["count"],
                hovertemplate="%{text}: at least %{y:.2f} (%{customdata:,} tests)<extra>"
                + panel.ancestry + "</extra>",
            ),
            row=row, col=col,
        )
        ty = _threshold_y(panel)
        if ty is not None:
            fig.add_hline(
                y=ty, line={"color": THRESHOLD, "dash": "dash", "width": 1.2},
                annotation_text=f"p = {panel.threshold:g}", annotation_position="bottom right",
                row=row, col=col,
            )
        fig.update_xaxes(range=[0, xmax], title_text="expected -log10(p)", row=row, col=col)
        fig.update_yaxes(range=[0, ymax], gridcolor=GRID, title_text="observed -log10(p)", row=row, col=col)
    plotly_layout(fig, 360 * rows + 60, legend={"orientation": "h", "y": -0.1})
    return fig


# --- Registry ----------------------------------------------------------------


def _has_tables(*names: str) -> Callable[[Report], bool]:
    """Available when every named table exists for at least one group.

    Never raises: a report that cannot be read for tables is a skip.
    """

    def available(report: Report) -> bool:
        try:
            return all(_usable(plot_tables(report, name)) for name in names)
        except Exception:  # noqa: BLE001 - `available` must not raise
            logger.debug("could not look for %s", names, exc_info=True)
            return False

    return available


def _qq_plot(spec: QQSpec, description: str) -> Plot:
    return Plot(
        name=f"{spec.name}_qq",
        title=f"{spec.title} QQ",
        description=description,
        available=_has_tables(spec.bins_table),
        png=lambda r, plt, p: qq_png(r, spec, plt, p),
        plotly=lambda r: qq_plotly(r, spec),
    )


_QQ_NOTE = (
    " Binned bulk, exact tail: points past the threshold are the individual "
    "failing tests; a hollow triangle is a floor (open last bin, or p = 0)."
)

PLOTS: List[Plot] = [
    Plot(
        name="variant_missingness",
        title="Variant missingness",
        description=(
            "Distribution of per-variant missingness, per ancestry group, with "
            "the --geno threshold. Counts are on a log axis."
        ),
        available=_has_tables("geno_missingness_bins"),
        png=missingness_png,
        plotly=missingness_plotly,
    ),
    _qq_plot(
        QQ_SPECS["hwe"],
        "QQ plot of Hardy-Weinberg p-values, per ancestry group, for the test "
        "the filter used (controls on case/control data, else all samples)." + _QQ_NOTE,
    ),
    _qq_plot(
        QQ_SPECS["case_control"],
        "QQ plot of case-control missingness p-values, per ancestry group." + _QQ_NOTE,
    ),
    _qq_plot(
        QQ_SPECS["haplotype"],
        "QQ plot of haplotype-missingness (--test-mishap) p-values, per "
        "ancestry group. Each point is a haplotype test, not a variant." + _QQ_NOTE,
    ),
]
