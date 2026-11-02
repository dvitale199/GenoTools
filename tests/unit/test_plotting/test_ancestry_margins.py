"""The CAH margin histogram (#309), against a hand-built decisions frame.

The table the plot reads is `centroid_margin_table(decisions)`, so the fixture
builds a decisions frame in the shape `AncestryModel._predict_admixed` returns
and runs it through that function: a change to either side breaks these tests.
Only this module's plots are rendered; nothing here asserts over `ALL_PLOTS`.
"""

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from genotools.ancestry.diagnostics import CENTROID_MARGIN_COLUMNS, centroid_margin_table
from genotools.plotting import render_plots
from genotools.plotting.plots import ancestry as plots
from genotools.plotting.report import load_report
from genotools.plotting.style import NEUTRAL

from .test_ancestry_plots import write_ancestry_report

# NORDIC, PAPUAN and CAS2 are not panel labels: CAS2 is a centroid the model
# invents by splitting CAS, and CAH is invented by the override itself.
GROUPS = {"EUR": (60, 2.0), "AFR": (40, 1.0), "PAPUAN": (25, 0.3), "CAS2": (15, -0.2)}


def decisions_frame(seed: int = 3) -> pd.DataFrame:
    """Shaped like `AncestryModel._predict_admixed`'s second return value."""
    rng = np.random.default_rng(seed)
    rows = []
    for group, (n, centre) in GROUPS.items():
        for margin in rng.normal(centre, 0.8, n):
            rows.append((group, float(margin)))
    frame = pd.DataFrame(rows, columns=["nearest_ancestry", "margin_to_all"])
    ids = [f"s{i}" for i in range(len(frame))]
    frame.insert(0, "IID", ids)
    frame.insert(0, "FID", ids)
    cah = frame["margin_to_all"] <= 0
    frame["label_pre_admixture"] = frame["nearest_ancestry"]
    frame["label"] = np.where(cah, "CAH", frame["nearest_ancestry"])
    frame["nearest_centroid"] = np.where(cah, "ALL", frame["nearest_ancestry"])
    frame["dist_to_all"] = 5.0
    frame["dist_to_nearest_ancestry"] = 5.0 - frame["margin_to_all"]
    return frame


@pytest.fixture
def decisions() -> pd.DataFrame:
    return decisions_frame()


@pytest.fixture
def report(tmp_path: Path, decisions: pd.DataFrame):
    # `write_ancestry_report` keeps this module's own fixture table beside it.
    return load_report(write_ancestry_report(tmp_path / "cohort.json"))


class TestMarginTable:
    def test_keeps_the_columns_the_plot_reads(self, decisions) -> None:
        table = centroid_margin_table(decisions)
        assert list(table.columns) == CENTROID_MARGIN_COLUMNS
        assert len(table) == len(decisions)
        assert table["margin_to_all"].tolist() == decisions["margin_to_all"].tolist()

    def test_none_when_detection_did_not_run(self, decisions) -> None:
        assert centroid_margin_table(None) is None
        assert centroid_margin_table(decisions.iloc[0:0]) is None

    def test_none_when_there_is_no_margin_column(self, decisions) -> None:
        # A one-label panel has no ancestry centroid to take a margin against.
        assert centroid_margin_table(decisions.drop(columns=["margin_to_all"])) is None

    def test_cah_is_exactly_the_non_positive_margins_in_the_fixture(self, decisions) -> None:
        table = centroid_margin_table(decisions)
        assert ((table["label"] == "CAH") == (table["margin_to_all"] <= 0)).all()
        assert (table.loc[table["label"] == "CAH", "nearest_centroid"] == "ALL").all()


class TestMarginFrame:
    def test_title_counts_cah_from_the_labels(self, report, decisions) -> None:
        n_cah = int((decisions["label"] == "CAH").sum())
        assert 0 < n_cah < len(decisions)
        title = plots.margin_title(plots.margin_frame(report))
        assert f"{n_cah:,} of {len(decisions):,} samples" in title
        assert f"({100 * n_cah / len(decisions):.1f}%)" in title

    def test_panels_are_largest_group_first_and_cover_everyone(self, report, decisions) -> None:
        panels = plots.margin_panels(plots.margin_frame(report))
        assert [name for name, _ in panels] == ["EUR", "AFR", "PAPUAN", "CAS2"]
        assert sum(len(part) for _, part in panels) == len(decisions)

    def test_bin_edges_are_one_set_with_zero_as_an_exact_edge(self, report) -> None:
        frame = plots.margin_frame(report)
        edges = plots.margin_bins(frame)
        assert 0.0 in edges
        assert np.allclose(np.diff(edges), np.diff(edges)[0])
        assert edges[0] <= frame["margin_to_all"].min()
        assert edges[-1] >= frame["margin_to_all"].max()
        # All-positive and all-negative cohorts still get a zero edge.
        for shift in (10.0, -10.0):
            moved = frame.assign(margin_to_all=frame["margin_to_all"].abs() * np.sign(shift) + shift)
            assert 0.0 in plots.margin_bins(moved)

    def test_equal_sized_groups_tie_break_by_name(self) -> None:
        frame = pd.DataFrame(
            {"nearest_ancestry": ["b", "a", "b", "a"], "margin_to_all": [1.0] * 4,
             "label": ["b", "a", "b", "a"]}
        )
        assert [n for n, _ in plots.margin_panels(frame)] == ["a", "b"]

    def test_unavailable_without_the_table(self, tmp_path: Path) -> None:
        loaded = load_report(write_ancestry_report(tmp_path / "none.json", drop=["centroid_margins"]))
        assert not plots.PLOTS[-1].available(loaded)
        assert plots.margin_frame(loaded).empty

    def test_available_never_raises_on_a_malformed_table(self, tmp_path: Path) -> None:
        path = write_ancestry_report(tmp_path / "bad.json", drop=["centroid_margins"])
        dest = tmp_path / "bad_plot_data" / "ancestry"
        dest.mkdir(parents=True)
        (dest / "centroid_margins.tsv").write_text("FID\tIID\nA\tB\n")
        assert not plot_named().available(load_report(path))

    def test_styles_cover_the_vocabulary_and_do_not_disturb_it(self, report) -> None:
        frame = plots.margin_frame(report)
        styles = plots.margin_styles(report, frame["nearest_ancestry"])
        base = plots.styles_for(report)
        assert "CAS2" not in base and "CAS2" in styles
        # Every label the run knows keeps its style; the invented centroid
        # taking a colour must not shift any of them.
        assert all(styles[label] == style for label, style in base.items())
        assert styles["CAS2"].color == NEUTRAL


def plot_named():
    (plot,) = [p for p in plots.PLOTS if p.name == "ancestry_cah_margins"]
    return plot


class TestRender:
    def test_available_with_the_table(self, report) -> None:
        assert plot_named().available(report)

    def test_png_is_written(self, report, tmp_path: Path) -> None:
        plt = pytest.importorskip("matplotlib.pyplot")
        out = tmp_path / "margins.png"
        assert plot_named().png(report, plt, out) == out
        assert out.stat().st_size > 1000

    def test_plotly_has_a_trace_per_panel(self, report) -> None:
        pytest.importorskip("plotly")
        fig = plot_named().plotly(report)
        assert len(fig.data) == 1 + 4
        assert "of" in fig.layout.title.text and "CAH" in fig.layout.title.text
        # The dashed CAH line sits at zero in every panel.
        assert len(fig.layout.shapes) == 5
        assert all(s.x0 == 0.0 for s in fig.layout.shapes)

    def test_render_plots_writes_png_and_html_for_this_plot_only(self, report, tmp_path: Path) -> None:
        pytest.importorskip("matplotlib")
        pytest.importorskip("plotly")
        result = render_plots(report.path, tmp_path / "out", plots=[plot_named()])
        assert [p.name for p in result.pngs] == ["ancestry_cah_margins.png"]
        assert result.skipped == []
        assert 'id="ancestry_cah_margins"' in result.html.read_text()
