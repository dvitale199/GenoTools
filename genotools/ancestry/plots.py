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

"""Diagnostic plots for an ancestry run.

The numbers in `diagnostics.py` fire on their own and land in the log and the
report; these pictures are for the case where they have already told you
*something* is wrong and you want to see *what*. A cohort projected onto its
reference panel shows in one glance what a table of standard deviations only
implies: whether it lands inside the reference cloud, beside it, or collapsed
into a single point somewhere off the edge.

Opt-in via `--ancestry-plots`, because these are binary artifacts and every
other output of a run is text.
"""

from pathlib import Path
from typing import List, Optional

import pandas as pd

from genotools.core.logging import get_logger
from genotools.ancestry.diagnostics import pc_columns

logger = get_logger(__name__)

def _cohort_style(n_samples: int) -> dict:
    """Marker styling for the cohort overlay, thinned by cohort size.

    The overlay's job is to show *where* the cohort sits relative to the
    reference, which it stops doing once it paints over it: 10,000 points at a
    fixed alpha turn every dense region into one black blob and hide the
    coloured panel underneath. Alpha and size fall with n so a large cohort
    reads as shading over a visible reference rather than as a silhouette.
    """
    alpha = min(0.35, max(0.04, 1500 / max(n_samples, 1)))
    return dict(
        color="black",
        marker="x",
        s=12 if n_samples <= 2000 else 6,
        alpha=alpha,
        linewidths=0.5,
    )


def _pyplot():
    """Import pyplot against a non-interactive backend, or return None.

    Matplotlib is a hard dependency, but a run must not die because a plotting
    backend is unhappy on some cluster node -- the plots are the least load-
    bearing output of the pipeline.
    """
    try:
        import matplotlib

        matplotlib.use("Agg", force=False)
        import matplotlib.pyplot as plt

        return plt
    except Exception as exc:  # pragma: no cover - environment-specific
        logger.warning(f"Diagnostic plots skipped: matplotlib unavailable ({exc})")
        return None


def plot_pca_projection(
    train_pca: pd.DataFrame,
    projected: pd.DataFrame,
    out_path: Path,
    pc_pairs: Optional[List[tuple]] = None,
) -> Optional[Path]:
    """Reference PCA coloured by ancestry, with the cohort overlaid.

    Args:
        train_pca: Reference PCA frame with PC columns and a `label` column.
        projected: Cohort PCA frame with the same PC columns.
        out_path: PNG path to write.
        pc_pairs: `(x, y)` PC index pairs. Default is the three leading pairs,
            which is where ancestry structure lives.

    Returns:
        The path written, or None if nothing could be plotted.
    """
    plt = _pyplot()
    if plt is None:
        return None

    available = [c for c in pc_columns(train_pca) if c in projected.columns]
    if len(available) < 2:
        logger.warning("Diagnostic plots skipped: fewer than 2 shared PCs")
        return None

    pc_pairs = pc_pairs or [(1, 2), (2, 3), (3, 4)]
    pairs = [
        (f"PC{x}", f"PC{y}")
        for x, y in pc_pairs
        if f"PC{x}" in available and f"PC{y}" in available
    ]
    if not pairs:
        pairs = [(available[0], available[1])]

    fig, axes = plt.subplots(1, len(pairs), figsize=(5.2 * len(pairs), 4.8))
    axes = [axes] if len(pairs) == 1 else list(axes)

    labels = (
        sorted(train_pca["label"].astype(str).unique())
        if "label" in train_pca.columns
        else []
    )
    colors = plt.get_cmap("tab20")

    for ax, (xcol, ycol) in zip(axes, pairs):
        if labels:
            for i, label in enumerate(labels):
                subset = train_pca[train_pca["label"].astype(str) == label]
                ax.scatter(
                    subset[xcol], subset[ycol],
                    s=6, alpha=0.5, color=colors(i % 20), label=label,
                )
        else:
            ax.scatter(train_pca[xcol], train_pca[ycol], s=6, alpha=0.5, color="lightgray")
        ax.scatter(
            projected[xcol], projected[ycol], label="cohort",
            **_cohort_style(len(projected)),
        )
        ax.set_xlabel(xcol)
        ax.set_ylabel(ycol)

    axes[0].set_title("reference panel (colour) vs projected cohort (black x)", loc="left")
    handles, legend_labels = axes[0].get_legend_handles_labels()
    fig.legend(
        handles, legend_labels, loc="center left", bbox_to_anchor=(1.0, 0.5),
        frameon=False, fontsize="small",
    )
    fig.tight_layout()
    fig.savefig(out_path, dpi=120, bbox_inches="tight")
    plt.close(fig)
    return out_path


def plot_centroid_distances(
    decisions: pd.DataFrame,
    out_path: Path,
) -> Optional[Path]:
    """How close each sample came to being called CAH.

    `margin_to_all` is the distance to the global centroid minus the distance
    to the nearest ancestry centroid, so the CAH rule is exactly "margin below
    zero". Plotting the distribution against that line shows whether an all-CAH
    result was a near-miss (the cohort sitting on the boundary) or a rout (the
    whole cohort far to the left of it, which means the projection, not the
    threshold).

    Args:
        decisions: Per-sample frame from `AncestryModel._predict_admixed`.
        out_path: PNG path to write.

    Returns:
        The path written, or None if there was nothing to plot.
    """
    plt = _pyplot()
    if plt is None or "margin_to_all" not in decisions.columns:
        return None

    margin = pd.to_numeric(decisions["margin_to_all"], errors="coerce").dropna()
    if margin.empty:
        return None

    fig, (ax_margin, ax_counts) = plt.subplots(1, 2, figsize=(11, 4.2))

    ax_margin.hist(margin, bins=60, color="steelblue", alpha=0.85)
    ax_margin.axvline(0.0, color="crimson", linestyle="--", linewidth=1.2)
    ax_margin.set_xlabel("distance to ALL centroid - distance to nearest ancestry")
    ax_margin.set_ylabel("samples")
    n_cah = int((margin < 0).sum())
    ax_margin.set_title(
        f"left of the red line is CAH: {n_cah}/{len(margin)} samples "
        f"({100 * n_cah / len(margin):.1f}%)",
        loc="left", fontsize="medium",
    )

    counts = decisions["label"].astype(str).value_counts()
    ax_counts.barh(list(counts.index)[::-1], list(counts.values)[::-1], color="slategray")
    ax_counts.set_xlabel("samples")
    ax_counts.set_title("final labels", loc="left", fontsize="medium")

    fig.tight_layout()
    fig.savefig(out_path, dpi=120)
    plt.close(fig)
    return out_path


def plot_ancestry_diagnostics(
    train_pca: Optional[pd.DataFrame],
    projected: Optional[pd.DataFrame],
    out_prefix: str,
    decisions: Optional[pd.DataFrame] = None,
) -> List[str]:
    """Write every diagnostic plot available for a run.

    Args:
        train_pca: Reference PCA frame, or None if the run has none.
        projected: Cohort PCA frame, or None.
        out_prefix: Path prefix; each plot appends its own suffix.
        decisions: Per-sample admixture decisions, if detection ran.

    Returns:
        Paths written, as strings, for the report.
    """
    written: List[str] = []

    if train_pca is not None and projected is not None:
        path = plot_pca_projection(
            train_pca, projected, Path(f"{out_prefix}_pca_diagnostic.png")
        )
        if path is not None:
            written.append(str(path))

    if decisions is not None and not decisions.empty:
        path = plot_centroid_distances(
            decisions, Path(f"{out_prefix}_centroid_distance.png")
        )
        if path is not None:
            written.append(str(path))

    if written:
        logger.info(f"Diagnostic plots written: {', '.join(written)}")
    return written
