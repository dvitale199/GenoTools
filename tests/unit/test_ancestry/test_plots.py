"""Tests for the ancestry diagnostic plots.

Plots are the least load-bearing output of a run, so what matters here is that
they are written when asked for, degrade rather than raise when there is
nothing to draw, and never take a run down with them.
"""

from pathlib import Path

import numpy as np
import pandas as pd

from genotools.ancestry.plots import (
    plot_ancestry_diagnostics,
    plot_centroid_distances,
    plot_pca_projection,
)


def _train_pca(n: int = 60) -> pd.DataFrame:
    rng = np.random.default_rng(0)
    return pd.DataFrame({
        "FID": [f"R{i}" for i in range(n)],
        "IID": [f"R{i}" for i in range(n)],
        "PC1": rng.normal(0, 1, n),
        "PC2": rng.normal(0, 1, n),
        "PC3": rng.normal(0, 1, n),
        "PC4": rng.normal(0, 1, n),
        "label": ["EUR" if i % 2 else "AFR" for i in range(n)],
    })


def _projected(n: int = 10) -> pd.DataFrame:
    rng = np.random.default_rng(1)
    return pd.DataFrame({
        "FID": [f"S{i}" for i in range(n)],
        "IID": [f"S{i}" for i in range(n)],
        "PC1": rng.normal(5, 0.1, n),
        "PC2": rng.normal(0, 0.1, n),
        "PC3": rng.normal(0, 0.1, n),
        "PC4": rng.normal(0, 0.1, n),
        "label": ["CAH"] * n,
    })


def _decisions(n: int = 10) -> pd.DataFrame:
    return pd.DataFrame({
        "FID": [f"S{i}" for i in range(n)],
        "IID": [f"S{i}" for i in range(n)],
        "label_pre_admixture": ["EUR"] * n,
        "label": ["CAH"] * n,
        "margin_to_all": np.linspace(-2, 1, n),
    })


def test_both_plots_are_written(tmp_path: Path) -> None:
    written = plot_ancestry_diagnostics(
        _train_pca(), _projected(), str(tmp_path / "run"), _decisions()
    )
    assert [Path(p).name for p in written] == [
        "run_pca_diagnostic.png",
        "run_centroid_distance.png",
    ]
    assert all(Path(p).stat().st_size > 0 for p in written)


def test_the_pca_panel_survives_a_cohort_with_only_two_pcs(tmp_path: Path) -> None:
    train = _train_pca()[["FID", "IID", "PC1", "PC2", "label"]]
    projected = _projected()[["FID", "IID", "PC1", "PC2"]]
    path = plot_pca_projection(train, projected, tmp_path / "two.png")
    assert path is not None and path.exists()


def test_a_single_shared_pc_is_not_plotted(tmp_path: Path) -> None:
    train = _train_pca()[["FID", "IID", "PC1", "label"]]
    projected = _projected()[["FID", "IID", "PC1"]]
    assert plot_pca_projection(train, projected, tmp_path / "one.png") is None


def test_an_unlabeled_reference_still_plots(tmp_path: Path) -> None:
    train = _train_pca().drop(columns=["label"])
    path = plot_pca_projection(train, _projected(), tmp_path / "nolabel.png")
    assert path is not None and path.exists()


def test_decisions_without_a_margin_produce_no_distance_plot(tmp_path: Path) -> None:
    decisions = _decisions().drop(columns=["margin_to_all"])
    assert plot_centroid_distances(decisions, tmp_path / "none.png") is None


def test_nothing_to_plot_is_an_empty_list_not_an_error(tmp_path: Path) -> None:
    assert plot_ancestry_diagnostics(None, None, str(tmp_path / "run"), None) == []


def test_the_cohort_overlay_thins_as_the_cohort_grows() -> None:
    """A 10k-sample overlay at a fixed alpha hides the reference underneath it,
    which is the one thing the panel exists to show."""
    from genotools.ancestry.plots import _cohort_style

    small, large = _cohort_style(200), _cohort_style(10_000)
    assert small["alpha"] > large["alpha"]
    assert small["s"] > large["s"]
    assert _cohort_style(10_000_000)["alpha"] >= 0.04
    assert _cohort_style(1)["alpha"] <= 0.35
