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

"""GWAS plots: Manhattan, QQ, PCs by phenotype and scree.

The Manhattan and QQ plots read the `.glm` file the `assoc` step wrote (its
absolute path is the resolved `gwas_output_path` parameter), keeping only the
`ADD` test rows, so logistic and linear output draw alike. The PC plot joins
the GWAS eigenvectors to the phenotype in the report's `input_samples`. Scree
plots read the GWAS `.eigenval` and, for ancestry training runs, the reference
PCA's `{out}_plot_data/ancestry/pca_eigenvalues.tsv`.

Every `available` is a cheap existence check that never raises.
"""

from pathlib import Path
from typing import List, Optional, Tuple

import numpy as np
import pandas as pd

from genotools.core.logging import get_logger
from genotools.plotting.plots.base import Plot
from genotools.plotting.report import Report, parameter, plot_data_dir
from genotools.plotting.style import (
    CATEGORICAL,
    GRID,
    INK,
    INK_SECONDARY,
    NEUTRAL,
    SERIES,
    SURFACE,
    THRESHOLD,
    plotly_layout,
    save_png,
    style_axis,
)

logger = get_logger(__name__)

GENOME_WIDE = 5e-8
# The smallest p a double holds; a p of exactly 0 is drawn here, not at infinity.
P_FLOOR = 1e-300
# Plotly figures are thinned past this many points (the PNG keeps every point).
PLOTLY_MAX_POINTS = 60_000

_CHROM_ALIASES = {"X": 23, "Y": 24, "XY": 25, "MT": 26, "M": 26}
_CHROM_LABELS = {23: "X", 24: "Y", 25: "XY", 26: "MT"}


# --- Locating files ----------------------------------------------------------


def gwas_group(report: Report) -> Optional[str]:
    """The ancestry group whose GWAS the plots draw: `all`, else the first with one."""
    params = report.parameters
    if params is None or params.empty or "parameter" not in params.columns:
        return None
    rows = params[(params["step"] == "assoc") & (params["parameter"] == "gwas_output_path")]
    groups = sorted(rows["ancestry"].astype(str).unique())
    if not groups:
        return None
    return "all" if "all" in groups else groups[0]


def _assoc_file(report: Report, name: str) -> Optional[Path]:
    """An `assoc` output path from the report, or None if unrecorded or gone."""
    try:
        group = gwas_group(report)
        if group is None:
            return None
        value = parameter(report, "assoc", name, group)
        if not value:
            return None
        path = Path(str(value))
        return path if path.is_file() else None
    except Exception:  # available() must never raise
        return None


def _ancestry_eigenvalues_file(report: Report) -> Optional[Path]:
    try:
        root = plot_data_dir(report)
        if root is None:
            return None
        path = root / "ancestry" / "pca_eigenvalues.tsv"
        return path if path.is_file() else None
    except Exception:
        return None


def _has_assoc_file(name: str):
    return lambda report: _assoc_file(report, name) is not None


# --- Frames ------------------------------------------------------------------


def chrom_number(chrom: pd.Series) -> pd.Series:
    """Chromosome as a number: `chr` stripped, X/Y/XY/MT as 23-26, else NaN."""
    text = chrom.astype(str).str.strip().str.replace("^chr", "", regex=True, case=False)
    text = text.str.upper().replace(_CHROM_ALIASES)
    return pd.to_numeric(text, errors="coerce")


def read_glm(path: Path) -> pd.DataFrame:
    """The test rows of a PLINK 2 `.glm` file as `chrom`, `pos`, `id`, `p`.

    Handles `.glm.logistic[.hybrid]` and `.glm.linear` alike: they share the
    columns used here. Covariate rows (any `TEST` other than `ADD`) are dropped
    when a `TEST` column exists; rows without a usable p are dropped, and a p of
    exactly 0 is clipped to `P_FLOOR`.
    """
    header = pd.read_csv(path, sep="\t", nrows=0).columns
    want = ["#CHROM", "CHROM", "POS", "ID", "TEST", "P"]
    cols = [c for c in header if c in want]
    frame = pd.read_csv(path, sep="\t", usecols=cols, dtype={"ID": str, "#CHROM": str, "CHROM": str}, na_values=["NA", "."])
    frame = frame.rename(columns={"#CHROM": "CHROM"})
    if "TEST" in frame.columns:
        frame = frame[frame["TEST"].astype(str) == "ADD"]
    p = pd.to_numeric(frame["P"], errors="coerce")
    frame = frame.assign(P=p)
    frame = frame[frame["P"].notna() & (frame["P"] >= 0) & (frame["P"] <= 1)]
    out = pd.DataFrame(
        {
            "chrom": frame["CHROM"].astype(str).to_numpy() if "CHROM" in frame else "",
            "pos": pd.to_numeric(frame["POS"], errors="coerce").to_numpy() if "POS" in frame else np.nan,
            "id": frame["ID"].astype(str).to_numpy() if "ID" in frame else "",
            "p": np.clip(frame["P"].to_numpy(dtype=float), P_FLOOR, 1.0),
        }
    )
    return out.reset_index(drop=True)


def manhattan_frame(glm: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Points with a cumulative genome position, plus one tick per chromosome.

    Returns:
        `(points, ticks)`. `points` has `chrom` (number), `pos`, `id`, `p`,
        `neglogp`, `x` and `band` (0/1, alternating by chromosome, so adjacent
        chromosomes read apart). `ticks` has `chrom`, `label`, `x`. Variants on
        a chromosome that is not a number or a recognised X/Y/MT are dropped.
    """
    cols = ["chrom", "pos", "id", "p", "neglogp", "x", "band"]
    frame = glm.assign(chrom=chrom_number(glm["chrom"]), pos=pd.to_numeric(glm["pos"], errors="coerce"))
    frame = frame.dropna(subset=["chrom", "pos"])
    if frame.empty:
        return pd.DataFrame(columns=cols), pd.DataFrame(columns=["chrom", "label", "x"])
    frame = frame.assign(chrom=frame["chrom"].astype(int))
    frame = frame.sort_values(["chrom", "pos"], kind="stable").reset_index(drop=True)
    extent = frame.groupby("chrom", sort=True)["pos"].agg(["min", "max"])
    # Pad chromosomes so neighbours do not touch.
    gap = float(extent["max"].max()) * 0.01
    length = (extent["max"] - extent["min"]) + gap
    offset = length.cumsum() - length
    frame["x"] = frame["pos"] - frame["chrom"].map(extent["min"]) + frame["chrom"].map(offset)
    frame["neglogp"] = -np.log10(frame["p"])
    order = {c: i for i, c in enumerate(extent.index)}
    frame["band"] = frame["chrom"].map(order) % 2
    centre = offset + (extent["max"] - extent["min"]) / 2
    ticks = pd.DataFrame(
        {
            "chrom": extent.index,
            "label": [_CHROM_LABELS.get(int(c), str(int(c))) for c in extent.index],
            "x": centre.to_numpy(),
        }
    )
    return frame[cols], ticks


def qq_frame(glm: pd.DataFrame, max_points: Optional[int] = None) -> pd.DataFrame:
    """Observed against expected -log10(p), most significant last.

    Args:
        glm: Output of `read_glm`.
        max_points: When set and exceeded, the bulk (expected below 3, where
            points pile up) is thinned to one point per 0.01 x 0.01 cell; the
            tail, where a QQ plot is read, is always kept whole.

    Returns:
        Columns `expected`, `observed`.
    """
    p = np.sort(glm["p"].to_numpy(dtype=float), kind="stable")
    n = len(p)
    if n == 0:
        return pd.DataFrame(columns=["expected", "observed"])
    expected = -np.log10((np.arange(1, n + 1) - 0.5) / n)
    observed = -np.log10(p)
    # Ascending p is descending significance; flip so the tail is last.
    frame = pd.DataFrame({"expected": expected[::-1], "observed": observed[::-1]})
    if max_points is not None and n > max_points:
        bulk = frame["expected"] < 3
        cells = frame[bulk].round(2).drop_duplicates()
        frame = pd.concat([cells, frame[~bulk]], ignore_index=True)
    return frame.reset_index(drop=True)


LAMBDA_ROWS = [
    ("lambda", "\u03bb GC"),
    ("lambda1000", "\u03bb\u2081\u2080\u2080\u2080"),
    ("lambda_maf", "\u03bb GC, MAF-filtered"),
    ("lambda1000_maf", "\u03bb\u2081\u2080\u2080\u2080, MAF-filtered"),
]


def lambda_values(report: Report) -> List[Tuple[str, float]]:
    """The report's GWAS inflation values as `(label, value)`, present ones only."""
    block = report.raw.get("GWAS")
    if not isinstance(block, dict) or "metric" not in block or "value" not in block:
        return []
    try:
        frame = pd.DataFrame(block)
        if "ancestry" in frame.columns:
            group = gwas_group(report) or "all"
            frame = frame[frame["ancestry"].astype(str) == group]
        values = dict(zip(frame["metric"].astype(str), pd.to_numeric(frame["value"], errors="coerce")))
    except Exception:
        return []
    return [(label, float(values[key])) for key, label in LAMBDA_ROWS if key in values and pd.notna(values[key])]


def _lambda_text(report: Report) -> str:
    return "\n".join(f"{label} = {value:.3f}" for label, value in lambda_values(report))


def pheno_name(report: Report) -> str:
    """The phenotype column the GWAS used (`PHENO1` unless `--pheno-name`)."""
    group = gwas_group(report) or "all"
    value = parameter(report, "assoc", "gwas.pheno_name", group, source="requested")
    return str(value) if value else "PHENO1"


_CASE_CONTROL = {1: "control", 2: "case"}


def pc_pheno_frame(report: Report) -> pd.DataFrame:
    """GWAS PC1/PC2 joined to each sample's phenotype.

    Returns:
        Columns `IID`, `PC1`, `PC2`, `pheno` (the raw value) and `group`:
        `control` / `case` for PLINK's 1/2 coding, `missing` for 0/-9/blank, and
        for a quantitative phenotype the numeric value is kept in `pheno` and
        `group` is `quantitative`. Empty when the files or columns are missing.
    """
    columns = ["IID", "PC1", "PC2", "pheno", "group"]
    path = _assoc_file(report, "eigenvec_path")
    samples = report.input_samples
    if path is None or samples is None:
        return pd.DataFrame(columns=columns)
    name = pheno_name(report)
    if name not in samples.columns:
        return pd.DataFrame(columns=columns)
    vec = pd.read_csv(path, sep="\t", dtype={"IID": str, "#FID": str, "FID": str})
    vec = vec.rename(columns={"#FID": "FID"})
    if not {"IID", "PC1", "PC2"} <= set(vec.columns):
        return pd.DataFrame(columns=columns)
    fid = "#FID" if "#FID" in samples.columns else "FID"
    keep = ["IID", name] + ([fid] if fid in samples.columns else [])
    pheno = samples[keep].copy()
    pheno["IID"] = pheno["IID"].astype(str)
    keys = ["IID"]
    if fid in pheno.columns and "FID" in vec.columns:
        pheno = pheno.rename(columns={fid: "FID"})
        pheno["FID"] = pheno["FID"].astype(str)
        keys = ["FID", "IID"]
    merged = vec.merge(pheno[keys + [name]], on=keys, how="inner").rename(columns={name: "pheno"})
    value = pd.to_numeric(merged["pheno"], errors="coerce")
    merged["pheno"] = value
    observed = set(value.dropna().unique())
    if observed <= {0.0, 1.0, 2.0, -9.0}:
        group = value.map(_CASE_CONTROL).astype(object)
        merged["group"] = group.where(group.notna(), "missing")
    else:
        merged["group"] = np.where(value.isna() | (value == -9), "missing", "quantitative")
        merged.loc[merged["group"] == "missing", "pheno"] = np.nan
    return merged[columns].reset_index(drop=True)


def scree_frame(path: Path) -> pd.DataFrame:
    """Eigenvalues by PC from a GWAS `.eigenval` or the ancestry TSV.

    Returns:
        Columns `pc` (1-based), `label`, `eigenvalue`, and `ratio` (explained
        variance ratio; NaN for a `.eigenval`, which holds no total).
    """
    text = Path(path).read_text().splitlines()
    first = text[0].split("\t") if text else []
    if "eigenvalue" in first:
        frame = pd.read_csv(path, sep="\t")
        values = pd.to_numeric(frame["eigenvalue"], errors="coerce")
        labels = frame["PC"].astype(str) if "PC" in frame.columns else None
        ratio = (
            pd.to_numeric(frame["explained_variance_ratio"], errors="coerce")
            if "explained_variance_ratio" in frame.columns
            else pd.Series(np.nan, index=frame.index)
        )
    else:
        values = pd.to_numeric(pd.Series([line.strip() for line in text if line.strip()]), errors="coerce")
        labels = None
        ratio = pd.Series(np.nan, index=values.index)
    out = pd.DataFrame({"eigenvalue": values.to_numpy(), "ratio": ratio.to_numpy()})
    out["pc"] = np.arange(1, len(out) + 1)
    out["label"] = labels.to_numpy() if labels is not None else [f"PC{i}" for i in out["pc"]]
    return out.dropna(subset=["eigenvalue"])[["pc", "label", "eigenvalue", "ratio"]].reset_index(drop=True)


def _thin(points: pd.DataFrame, limit: int = PLOTLY_MAX_POINTS) -> pd.DataFrame:
    """Keep every point with -log10(p) >= 2 and a regular stride of the rest."""
    if len(points) <= limit:
        return points
    top = points["neglogp"] >= 2
    rest = points[~top]
    budget = max(limit - int(top.sum()), 1000)
    stride = max(int(np.ceil(len(rest) / budget)), 1)
    return pd.concat([points[top], rest.iloc[::stride]]).sort_values("x", kind="stable")


# --- Manhattan ---------------------------------------------------------------

_BAND_COLORS = (SERIES, INK_SECONDARY)


def manhattan_png(points: pd.DataFrame, ticks: pd.DataFrame, plt, out_path: Path) -> Path:
    fig, ax = plt.subplots(figsize=(11, 4.5))
    style_axis(ax)
    for band, color in enumerate(_BAND_COLORS):
        part = points[points["band"] == band]
        ax.scatter(part["x"], part["neglogp"], s=6, color=color, alpha=0.8, linewidths=0, rasterized=True)
    ax.axhline(-np.log10(GENOME_WIDE), color=THRESHOLD, linestyle="--", linewidth=1)
    ax.text(
        1.0, -np.log10(GENOME_WIDE), " p = 5e-8", transform=ax.get_yaxis_transform(),
        ha="right", va="bottom", fontsize="small", color=THRESHOLD,
    )
    ax.set_xticks(ticks["x"])
    ax.set_xticklabels(ticks["label"], fontsize="x-small")
    ax.set_xlim(points["x"].min(), points["x"].max())
    ax.set_ylim(bottom=0)
    ax.set_xlabel("Chromosome", color=INK_SECONDARY)
    ax.set_ylabel("-log10(p)", color=INK_SECONDARY)
    ax.set_title(f"Manhattan: {len(points):,} variants", color=INK, loc="left")
    return save_png(fig, plt, out_path)


def manhattan_plotly(points: pd.DataFrame, ticks: pd.DataFrame):
    import plotly.graph_objects as go

    points = _thin(points)
    fig = go.Figure()
    for band, color in enumerate(_BAND_COLORS):
        part = points[points["band"] == band]
        fig.add_trace(
            go.Scattergl(
                x=part["x"], y=part["neglogp"], mode="markers",
                marker={"color": color, "size": 4, "opacity": 0.8},
                customdata=np.stack([part["id"], part["chrom"], part["pos"], part["p"]], axis=-1),
                hovertemplate="%{customdata[0]}<br>chr%{customdata[1]}:%{customdata[2]}"
                "<br>p = %{customdata[3]:.3g}<extra></extra>",
                showlegend=False,
            )
        )
    fig.add_hline(y=-np.log10(GENOME_WIDE), line={"color": THRESHOLD, "dash": "dash", "width": 1},
                  annotation_text="p = 5e-8", annotation_position="top right")
    plotly_layout(fig, 450, title="Manhattan")
    fig.update_xaxes(tickmode="array", tickvals=list(ticks["x"]), ticktext=list(ticks["label"]),
                     title="Chromosome", showgrid=False)
    fig.update_yaxes(title="-log10(p)", rangemode="tozero", gridcolor=GRID)
    return fig


# --- QQ ----------------------------------------------------------------------


def qq_png(frame: pd.DataFrame, note: str, plt, out_path: Path) -> Path:
    fig, ax = plt.subplots(figsize=(5.5, 5.5))
    style_axis(ax)
    top = float(max(frame["expected"].max(), frame["observed"].max())) * 1.05
    ax.plot([0, top], [0, top], color=THRESHOLD, linewidth=1, linestyle="--")
    ax.scatter(frame["expected"], frame["observed"], s=8, color=SERIES, linewidths=0, rasterized=True)
    ax.set_xlim(0, top)
    ax.set_ylim(0, top)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("Expected -log10(p)", color=INK_SECONDARY)
    ax.set_ylabel("Observed -log10(p)", color=INK_SECONDARY)
    ax.set_title("GWAS QQ", color=INK, loc="left")
    if note:
        ax.text(0.04, 0.96, note, transform=ax.transAxes, ha="left", va="top",
                fontsize="small", color=INK, linespacing=1.4)
    return save_png(fig, plt, out_path)


def qq_plotly(frame: pd.DataFrame, note: str):
    import plotly.graph_objects as go

    top = float(max(frame["expected"].max(), frame["observed"].max())) * 1.05
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=[0, top], y=[0, top], mode="lines", showlegend=False, hoverinfo="skip",
                             line={"color": THRESHOLD, "dash": "dash", "width": 1}))
    fig.add_trace(go.Scattergl(x=frame["expected"], y=frame["observed"], mode="markers", showlegend=False,
                               marker={"color": SERIES, "size": 5},
                               hovertemplate="expected %{x:.2f}<br>observed %{y:.2f}<extra></extra>"))
    plotly_layout(fig, 520, title="GWAS QQ")
    fig.update_xaxes(title="Expected -log10(p)", range=[0, top], gridcolor=GRID)
    fig.update_yaxes(title="Observed -log10(p)", range=[0, top], gridcolor=GRID, scaleanchor="x")
    if note:
        fig.add_annotation(x=0.03, y=0.97, xref="paper", yref="paper", showarrow=False, align="left",
                           xanchor="left", yanchor="top", text=note.replace("\n", "<br>"))
    return fig


# --- PCs by phenotype --------------------------------------------------------

# Case/control keep colour *and* shape; missing is faded context.
_GROUP_STYLE = {
    "control": (CATEGORICAL[0], "o", "circle"),
    "case": (CATEGORICAL[1], "s", "square"),
    "missing": (NEUTRAL, "o", "circle"),
}


def pcs_pheno_png(frame: pd.DataFrame, name: str, plt, out_path: Path) -> Path:
    fig, ax = plt.subplots(figsize=(6.5, 5.5))
    style_axis(ax)
    ax.xaxis.grid(True, color=GRID, linewidth=0.8)
    if (frame["group"] == "quantitative").any():
        part = frame[frame["group"] == "quantitative"]
        sc = ax.scatter(part["PC1"], part["PC2"], c=part["pheno"], cmap="viridis", s=18, linewidths=0, alpha=0.85)
        fig.colorbar(sc, ax=ax, label=name)
        miss = frame[frame["group"] == "missing"]
        ax.scatter(miss["PC1"], miss["PC2"], s=14, color=NEUTRAL, linewidths=0, alpha=0.6)
    else:
        for group in ("missing", "control", "case"):
            part = frame[frame["group"] == group]
            if part.empty:
                continue
            color, marker, _ = _GROUP_STYLE[group]
            ax.scatter(part["PC1"], part["PC2"], s=22, color=color, marker=marker, linewidths=0,
                       alpha=0.6 if group == "missing" else 0.8,
                       label=f"{group} (n={len(part):,})")
        ax.legend(frameon=False, fontsize="small", title=name, title_fontsize="small")
    ax.set_xlabel("PC1", color=INK_SECONDARY)
    ax.set_ylabel("PC2", color=INK_SECONDARY)
    ax.set_title(f"GWAS PC1 vs PC2 by {name}", color=INK, loc="left")
    return save_png(fig, plt, out_path)


def pcs_pheno_plotly(frame: pd.DataFrame, name: str):
    import plotly.graph_objects as go

    fig = go.Figure()
    if (frame["group"] == "quantitative").any():
        part = frame[frame["group"] == "quantitative"]
        fig.add_trace(go.Scatter(x=part["PC1"], y=part["PC2"], mode="markers", text=part["IID"],
                                 marker={"color": part["pheno"], "colorscale": "Viridis", "size": 7,
                                         "colorbar": {"title": name}},
                                 hovertemplate="%{text}<br>" + name + " = %{marker.color}<extra></extra>",
                                 showlegend=False))
    else:
        for group in ("missing", "control", "case"):
            part = frame[frame["group"] == group]
            if part.empty:
                continue
            color, _, symbol = _GROUP_STYLE[group]
            fig.add_trace(go.Scatter(x=part["PC1"], y=part["PC2"], mode="markers", text=part["IID"],
                                     name=f"{group} (n={len(part):,})",
                                     marker={"color": color, "symbol": symbol, "size": 7, "opacity": 0.75},
                                     hovertemplate="%{text}<extra>" + group + "</extra>"))
    plotly_layout(fig, 520, title=f"GWAS PC1 vs PC2 by {name}")
    fig.update_xaxes(title="PC1", gridcolor=GRID)
    fig.update_yaxes(title="PC2", gridcolor=GRID)
    return fig


# --- Scree -------------------------------------------------------------------


def scree_png(frame: pd.DataFrame, title: str, plt, out_path: Path) -> Path:
    fig, ax = plt.subplots(figsize=(6, 4))
    style_axis(ax)
    ax.plot(frame["pc"], frame["eigenvalue"], color=SERIES, marker="o", linewidth=1.5)
    ax.set_xticks(frame["pc"])
    ax.set_xticklabels(frame["label"], fontsize="x-small", rotation=45 if len(frame) > 10 else 0)
    ax.set_ylim(bottom=0)
    ax.set_xlabel("Principal component", color=INK_SECONDARY)
    ax.set_ylabel("Eigenvalue", color=INK_SECONDARY)
    ax.set_title(title, color=INK, loc="left")
    return save_png(fig, plt, out_path)


def scree_plotly(frame: pd.DataFrame, title: str):
    import plotly.graph_objects as go

    ratio = frame["ratio"].to_numpy()
    hover = [
        f"{lab}<br>eigenvalue {ev:.4g}" + ("" if np.isnan(r) else f"<br>{r:.1%} of variance")
        for lab, ev, r in zip(frame["label"], frame["eigenvalue"], ratio)
    ]
    fig = go.Figure(go.Scatter(x=frame["pc"], y=frame["eigenvalue"], mode="lines+markers", text=hover,
                               hoverinfo="text", line={"color": SERIES}, marker={"color": SERIES}))
    plotly_layout(fig, 420, title=title, showlegend=False)
    fig.update_xaxes(title="Principal component", tickmode="array", tickvals=list(frame["pc"]),
                     ticktext=list(frame["label"]), gridcolor=GRID)
    fig.update_yaxes(title="Eigenvalue", rangemode="tozero", gridcolor=GRID)
    return fig


# --- Registry ----------------------------------------------------------------


def _glm(report: Report) -> pd.DataFrame:
    return read_glm(_assoc_file(report, "gwas_output_path"))


def _manhattan_png(report, plt, path):
    return manhattan_png(*manhattan_frame(_glm(report)), plt, path)


def _manhattan_plotly(report):
    return manhattan_plotly(*manhattan_frame(_glm(report)))


def _qq_png(report, plt, path):
    return qq_png(qq_frame(_glm(report)), _lambda_text(report), plt, path)


def _qq_plotly(report):
    return qq_plotly(qq_frame(_glm(report), max_points=PLOTLY_MAX_POINTS), _lambda_text(report))


def _has_pcs(report: Report) -> bool:
    try:
        return (
            _assoc_file(report, "eigenvec_path") is not None
            and report.input_samples is not None
            and pheno_name(report) in report.input_samples.columns
        )
    except Exception:
        return False


GWAS_SCREE_TITLE = "GWAS PCA scree"
ANCESTRY_SCREE_TITLE = "Ancestry reference PCA scree"

PLOTS: List[Plot] = [
    Plot(
        name="gwas_manhattan",
        title="GWAS: Manhattan",
        description=(
            "-log10(p) of each tested variant (the additive test) by genome "
            "position, chromosomes alternating colour. The dashed line is "
            "genome-wide significance, p = 5e-8. Hover a point for its ID."
        ),
        available=_has_assoc_file("gwas_output_path"),
        png=_manhattan_png,
        plotly=_manhattan_plotly,
    ),
    Plot(
        name="gwas_qq",
        title="GWAS: QQ",
        description=(
            "Observed against expected -log10(p) under the null; the dashed "
            "diagonal is no inflation. Annotated with the run's genomic "
            "inflation (lambda) from the report."
        ),
        available=_has_assoc_file("gwas_output_path"),
        png=_qq_png,
        plotly=_qq_plotly,
    ),
    Plot(
        name="gwas_pcs_by_phenotype",
        title="GWAS: PC1 vs PC2 by phenotype",
        description=(
            "The GWAS covariate PCs coloured by case/control status. Cases and "
            "controls separating along a PC is confounding the model must "
            "absorb; a mixed cloud is what a clean design looks like."
        ),
        available=_has_pcs,
        png=lambda r, plt, p: pcs_pheno_png(pc_pheno_frame(r), pheno_name(r), plt, p),
        plotly=lambda r: pcs_pheno_plotly(pc_pheno_frame(r), pheno_name(r)),
    ),
    Plot(
        name="gwas_pca_scree",
        title="GWAS: PCA scree",
        description="Eigenvalue of each principal component used as a GWAS covariate.",
        available=_has_assoc_file("eigenval_path"),
        png=lambda r, plt, p: scree_png(scree_frame(_assoc_file(r, "eigenval_path")), GWAS_SCREE_TITLE, plt, p),
        plotly=lambda r: scree_plotly(scree_frame(_assoc_file(r, "eigenval_path")), GWAS_SCREE_TITLE),
    ),
    Plot(
        name="ancestry_pca_scree",
        title="Ancestry: reference PCA scree",
        description=(
            "Eigenvalue of each principal component of the ancestry reference "
            "panel, with the share of variance it explains on hover. Present "
            "only for runs that trained the ancestry model."
        ),
        available=lambda r: _ancestry_eigenvalues_file(r) is not None,
        png=lambda r, plt, p: scree_png(scree_frame(_ancestry_eigenvalues_file(r)), ANCESTRY_SCREE_TITLE, plt, p),
        plotly=lambda r: scree_plotly(scree_frame(_ancestry_eigenvalues_file(r)), ANCESTRY_SCREE_TITLE),
    ),
]
