"""Ancestry plots, against a synthetic report built to the real shapes.

No reference panel is needed: `write_ancestry_report` assembles what a
`--ancestry` run puts in its JSON (`PipelineOutput.to_dict()`), with the PCs
clustered by label so the figures look like something. The vocabulary includes
labels outside the default panel, to prove nothing here names one.
"""

import json
from pathlib import Path
from typing import Any, Dict, List, Optional

import numpy as np
import pandas as pd
import pytest

from genotools.ancestry.diagnostics import (
    WARN_PC_SD_RATIO,
    WARN_PC_SHIFT_SD,
    AncestryDiagnostics,
    pc_drift,
)
from genotools.plotting import render_plots
from genotools.plotting.plots import ancestry as plots
from genotools.plotting.plots.ancestry import PLOTS
from genotools.plotting.report import load_report
from genotools.plotting.style import PANEL_COLORS, label_styles

from .conftest import qc_row, write_report

N_PCS = 5

# Cluster centres in PC1-3. NORDIC and PAPUAN are not in the default reference
# panel: a user-supplied --ref-labels vocabulary.
CENTRES = {
    "EUR": (0.0, 0.0, 0.0),
    "AFR": (6.0, 1.0, -1.0),
    "EAS": (-5.0, 5.0, 1.0),
    "AMR": (-2.0, -4.0, 3.0),
    "NORDIC": (1.0, 1.5, -2.5),
    "PAPUAN": (-6.0, -6.0, -4.0),
}
COHORT_LABELS = {"EUR": 40, "AFR": 25, "PAPUAN": 20, "NORDIC": 15, "AMR": 10}
# One cohort PC sits far from the reference, so the drift plot has a flag.
DRIFTED_PC = "PC5"


def _cluster(rng: np.random.Generator, label: str, n: int) -> np.ndarray:
    centre = np.array(list(CENTRES[label]) + [0.0] * (N_PCS - 3))
    return centre + rng.normal(0.0, 0.8, size=(n, N_PCS))


def build_tables(seed: int = 7) -> Dict[str, Any]:
    """Reference and cohort PC frames, plus the cohort's predicted labels."""
    rng = np.random.default_rng(seed)
    pcs = [f"PC{i}" for i in range(1, N_PCS + 1)]

    ref_rows, ref_labels = [], []
    for label in CENTRES:
        ref_rows.append(_cluster(rng, label, 60))
        ref_labels += [label] * 60
    ref = pd.DataFrame(np.vstack(ref_rows), columns=pcs)
    ref.insert(0, "IID", [f"ref{i}" for i in range(len(ref))])
    ref.insert(0, "FID", ref["IID"])
    ref["label"] = ref_labels

    cohort_rows, cohort_labels = [], []
    for label, n in COHORT_LABELS.items():
        cohort_rows.append(_cluster(rng, label, n))
        cohort_labels += [label] * n
    cohort = pd.DataFrame(np.vstack(cohort_rows), columns=pcs)
    cohort[DRIFTED_PC] += 7.0 * ref[DRIFTED_PC].std()
    cohort.insert(0, "IID", [f"s{i}" for i in range(len(cohort))])
    cohort.insert(0, "FID", cohort["IID"])
    predicted = pd.DataFrame(
        {"FID": cohort["FID"], "IID": cohort["IID"], "label": cohort_labels}
    )
    return {"ref": ref, "cohort": cohort, "predicted": predicted}


def write_ancestry_report(
    path: Path,
    *,
    shuffle_labels: bool = False,
    drop: Optional[List[str]] = None,
    limits: Optional[Dict[str, Any]] = "default",  # type: ignore[assignment]
) -> Path:
    """A report shaped like a training `--ancestry` run's.

    Args:
        path: Where to write.
        shuffle_labels: Reverse `ancestry_labels`, so a join by position
            instead of by sample ID mislabels the cohort.
        drop: Report keys to leave out.
        limits: `pc_drift_limits` to record; "default" records the run's own,
            None leaves them out (a report from before they were recorded).
    """
    tables = build_tables()
    ref, cohort, predicted = tables["ref"], tables["cohort"], tables["predicted"]

    labels = sorted(CENTRES)
    confusion = pd.DataFrame(
        [
            [18 if r == c else (1 if abs(i - j) == 1 else 0) for j, c in enumerate(labels)]
            for i, r in enumerate(labels)
        ],
        index=labels,
        columns=labels,
    )
    counts = predicted["label"].value_counts().reset_index()
    counts.columns = ["label", "count"]

    diagnostics = AncestryDiagnostics(pc_drift=pc_drift(ref, cohort)).to_dict()
    if limits is None:
        diagnostics.pop("pc_drift_limits")
    elif limits != "default":
        diagnostics["pc_drift_limits"] = limits

    if shuffle_labels:
        predicted = predicted.iloc[::-1].reset_index(drop=True)

    report: Dict[str, Any] = {
        "ancestry_counts": counts.to_dict(),
        "ancestry_labels": predicted.to_dict(),
        "confusion_matrix": confusion.to_dict(),
        "test_accuracy": 0.9444,
        "ancestry_fit": {
            "test_balanced_accuracy": 0.9371,
            "cv_balanced_accuracy": 0.95,
            "train_balanced_accuracy": 0.97,
        },
        "ancestry_diagnostics": diagnostics,
        "ref_pcs": ref.to_dict(),
        "projected_pcs": cohort.to_dict(),
    }
    for key in drop or []:
        report.pop(key, None)
    path.write_text(json.dumps(report))
    return path


@pytest.fixture
def report_path(tmp_path: Path) -> Path:
    return write_ancestry_report(tmp_path / "cohort.json")


@pytest.fixture
def report(report_path: Path):
    return load_report(report_path)


def plot_named(name: str):
    (plot,) = [p for p in PLOTS if p.name == name]
    return plot


# --- availability ---------------------------------------------------------------


class TestAvailability:
    def test_every_plot_is_available_on_a_full_report(self, report) -> None:
        assert all(p.available(report) for p in PLOTS)

    def test_a_qc_only_report_has_none_of_them(self, tmp_path: Path) -> None:
        qc_only = load_report(write_report(tmp_path / "qc.json", [qc_row("callrate_prune", 1, "m")], n_input=5))
        assert not any(p.available(qc_only) for p in PLOTS)

    @pytest.mark.parametrize(
        "dropped, plot_names",
        [
            (["ref_pcs"], ["ancestry_pc_pairs", "ancestry_pc_3d"]),
            (["projected_pcs"], ["ancestry_pc_pairs", "ancestry_pc_3d"]),
            (["confusion_matrix"], ["ancestry_confusion_matrix"]),
            (["ancestry_counts"], ["ancestry_label_counts"]),
            (["ancestry_diagnostics"], ["ancestry_pc_drift"]),
        ],
    )
    def test_a_plot_needs_its_own_key(self, tmp_path: Path, dropped, plot_names) -> None:
        loaded = load_report(write_ancestry_report(tmp_path / "r.json", drop=dropped))
        unavailable = {p.name for p in PLOTS if not p.available(loaded)}
        assert unavailable == set(plot_names)

    def test_three_dimensions_need_three_shared_pcs(self, tmp_path: Path) -> None:
        raw = json.loads(write_ancestry_report(tmp_path / "r.json").read_text())
        for key in ("ref_pcs", "projected_pcs"):
            for pc in ("PC3", "PC4", "PC5"):
                raw[key].pop(pc)
        path = tmp_path / "two.json"
        path.write_text(json.dumps(raw))
        loaded = load_report(path)
        assert plot_named("ancestry_pc_pairs").available(loaded)
        assert not plot_named("ancestry_pc_3d").available(loaded)

    def test_drift_without_a_pc_drift_table_is_unavailable(self, tmp_path: Path) -> None:
        raw = json.loads(write_ancestry_report(tmp_path / "r.json").read_text())
        raw["ancestry_diagnostics"] = {"warnings": []}
        path = tmp_path / "nodrift.json"
        path.write_text(json.dumps(raw))
        assert not plot_named("ancestry_pc_drift").available(load_report(path))


# --- frames ---------------------------------------------------------------------


class TestPcFrame:
    def test_reference_and_cohort_rows(self, report) -> None:
        frame = plots.pc_frame(report)
        assert (frame["source"] == "reference").sum() == 60 * len(CENTRES)
        assert (frame["source"] == "cohort").sum() == sum(COHORT_LABELS.values())
        assert {"IID", "label", "PC1", "PC2", "PC3"} <= set(frame.columns)

    def test_cohort_takes_its_predicted_label_by_sample_id(self, tmp_path: Path) -> None:
        """Labels are joined on ID: the table order in the report is not trusted."""
        loaded = load_report(write_ancestry_report(tmp_path / "r.json", shuffle_labels=True))
        frame = plots.pc_frame(loaded)
        cohort = frame[frame["source"] == "cohort"].set_index("IID")["label"]
        truth = build_tables()["predicted"].set_index("IID")["label"]
        assert (cohort.loc[truth.index] == truth).all()

    def test_a_cohort_without_labels_is_still_drawable(self, tmp_path: Path) -> None:
        loaded = load_report(write_ancestry_report(tmp_path / "r.json", drop=["ancestry_labels"]))
        frame = plots.pc_frame(loaded)
        assert set(frame.loc[frame["source"] == "cohort", "label"]) == {plots.UNLABELLED}

    def test_vocabulary_is_reference_and_predicted_together(self, report) -> None:
        assert plots.label_vocabulary(report) == sorted(CENTRES)


class TestConfusion:
    def test_rows_sum_to_one_and_keep_true_label_orientation(self, report) -> None:
        counts, normalised = plots.confusion_frames(report)
        assert np.allclose(normalised.sum(axis=1), 1.0)
        # Row = true label, column = predicted: row/col labels survive the JSON trip.
        assert list(normalised.index) == list(normalised.columns) == sorted(CENTRES)
        first, second = sorted(CENTRES)[:2]
        assert normalised.at[first, second] == pytest.approx(
            counts.at[first, second] / counts.loc[first].sum()
        )

    def test_a_label_with_no_samples_stays_zero_not_nan(self, tmp_path: Path) -> None:
        raw = json.loads(write_ancestry_report(tmp_path / "r.json").read_text())
        first = sorted(CENTRES)[0]
        for col in raw["confusion_matrix"].values():
            col[first] = 0
        path = tmp_path / "empty_row.json"
        path.write_text(json.dumps(raw))
        _, normalised = plots.confusion_frames(load_report(path))
        assert (normalised.loc[first] == 0).all()
        assert not normalised.isna().any().any()

    def test_title_carries_both_accuracies(self, report) -> None:
        note = plots.accuracy_note(report)
        assert "balanced accuracy 93.7%" in note
        assert "accuracy 94.4%" in note

    def test_a_model_run_has_no_balanced_figure(self, tmp_path: Path) -> None:
        """`ancestry_fit` is training-only; the plot must not invent it."""
        loaded = load_report(write_ancestry_report(tmp_path / "r.json", drop=["ancestry_fit"]))
        note = plots.accuracy_note(loaded)
        assert "balanced" not in note and "accuracy 94.4%" in note


def test_counts_are_sorted_largest_first(report) -> None:
    frame = plots.counts_frame(report)
    assert list(frame["count"]) == sorted(COHORT_LABELS.values(), reverse=True)
    assert frame["count"].sum() == sum(COHORT_LABELS.values())


# --- the drift band comes from the report ----------------------------------------


class TestDriftBand:
    def test_the_fixture_records_the_runs_own_limits(self, report) -> None:
        assert plots.drift_limits(report) == (WARN_PC_SHIFT_SD, WARN_PC_SD_RATIO)

    def test_a_recorded_limit_moves_the_band(self, tmp_path: Path) -> None:
        loaded = load_report(
            write_ancestry_report(
                tmp_path / "r.json", limits={"max_shift_sd": 2.5, "sd_ratio_range": [0.5, 2.0]}
            )
        )
        limits = plots.drift_limits(loaded)
        assert limits == (2.5, (0.5, 2.0))
        fig = plots.drift_plotly(plots.drift_frame(loaded), limits)
        rects = [s for s in fig.layout.shapes if s.type == "rect"]
        assert sorted((s.y0, s.y1) for s in rects) == [(-2.5, 2.5), (0.5, 2.0)]
        lines = sorted(s.y0 for s in fig.layout.shapes if s.type == "line")
        assert lines == [-2.5, 0.5, 2.0, 2.5]

    def test_unrecorded_limits_draw_no_band_rather_than_guess(self, tmp_path: Path) -> None:
        loaded = load_report(write_ancestry_report(tmp_path / "r.json", limits=None))
        assert plots.drift_limits(loaded) is None
        fig = plots.drift_plotly(plots.drift_frame(loaded), None)
        assert not fig.layout.shapes

    def test_the_drifted_pc_lies_outside_the_band(self, report) -> None:
        drift = plots.drift_frame(report).set_index("pc")
        assert abs(drift.at[DRIFTED_PC, "shift_sd"]) > WARN_PC_SHIFT_SD
        assert abs(drift.at["PC1", "shift_sd"]) < WARN_PC_SHIFT_SD


# --- one vocabulary, one colour and marker per label, everywhere ------------------


def _trace_styles(fig) -> Dict[str, tuple]:
    """label -> (colour, symbol) for each cohort trace of a scatter figure."""
    found: Dict[str, tuple] = {}
    for trace in fig.data:
        label = trace.name.removesuffix(" (reference)")
        found.setdefault(label, set()).add((trace.marker.color, trace.marker.symbol))
    return found


class TestLabelEncoding:
    def test_every_scatter_trace_uses_the_labels_own_colour_and_marker(self, report) -> None:
        frame, styles = plots.pc_frame(report), plots.styles_for(report)
        expected = label_styles(plots.label_vocabulary(report))
        assert styles == expected
        symbols_3d = plots.symbols_3d(styles)
        for fig, symbols in (
            (plots.pairs_plotly(frame, styles), {l: s.symbol for l, s in expected.items()}),
            (plots.cube_plotly(frame, styles), symbols_3d),
        ):
            seen = _trace_styles(fig)
            assert set(seen) == set(CENTRES)
            for label, pairs in seen.items():
                assert pairs == {(expected[label].color, symbols[label])}

    def test_bars_use_the_same_colours_as_the_scatters(self, report) -> None:
        frame = plots.counts_frame(report)
        fig = plots.counts_plotly(frame, plots.styles_for(report))
        expected = plots.styles_for(report)
        assert list(fig.data[0].marker.color) == [expected[l].color for l in frame["label"]]

    def test_a_label_outside_the_default_panel_gets_a_distinct_colour_and_marker(self, report) -> None:
        styles = plots.styles_for(report)
        assert "PAPUAN" not in PANEL_COLORS and "NORDIC" not in PANEL_COLORS
        assert len({s.color for s in styles.values()}) == len(styles)
        assert len({s.symbol for s in styles.values()}) == len(styles)

    def test_3d_marker_shapes_are_ones_scatter3d_can_draw_and_stay_distinct(self, report) -> None:
        styles = plots.styles_for(report)
        symbols = plots.symbols_3d(styles)
        assert set(symbols.values()) <= set(plots.SCATTER3D_SYMBOLS)
        assert len(set(symbols.values())) == len(styles)
        # A label whose own shape scatter3d supports keeps it.
        assert all(symbols[l] == s.symbol for l, s in styles.items() if s.symbol in plots.SCATTER3D_SYMBOLS)

    def test_a_label_has_one_colour_in_every_plot(self, report) -> None:
        """A plot of fewer labels must not shift anyone's colour.

        The counts plot has no EAS (reference-only), so if it styled its own
        subset, a label outside the panel would take a different spare slot
        there than in the PC scatter.
        """
        by_name = {p.name: p for p in PLOTS}
        bars = by_name["ancestry_label_counts"].plotly(report).data[0]
        bar_colour = dict(zip(bars.x, bars.marker.color))
        pairs = by_name["ancestry_pc_pairs"].plotly(report)
        cohort = [t for t in pairs.data if t.name in bar_colour]
        assert {t.name for t in cohort} == set(bar_colour)
        for trace in cohort:
            assert trace.marker.color == bar_colour[trace.name], trace.name

    def test_reference_is_faded_and_cohort_solid(self, report) -> None:
        fig = plots.pairs_plotly(plots.pc_frame(report), plots.styles_for(report))
        ref = [t for t in fig.data if t.name.endswith("(reference)")]
        cohort = [t for t in fig.data if not t.name.endswith("(reference)")]
        assert ref and cohort
        assert max(t.marker.opacity for t in ref) < min(t.marker.opacity for t in cohort)

    def test_hover_names_the_sample_and_label(self, report) -> None:
        fig = plots.cube_plotly(plots.pc_frame(report), plots.styles_for(report))
        cohort = next(t for t in fig.data if t.name == "PAPUAN")
        assert cohort.text[0].startswith("s") and "PAPUAN" in cohort.text[0]
        assert cohort.type == "scatter3d"

    def test_pairs_have_three_panels(self, report) -> None:
        fig = plots.pairs_plotly(plots.pc_frame(report), plots.styles_for(report))
        assert len({t.xaxis for t in fig.data}) == 3


# --- rendering --------------------------------------------------------------------


class TestRender:
    def test_every_plot_renders_png_and_html(self, report_path: Path, tmp_path: Path) -> None:
        result = render_plots(report_path, tmp_path / "figs", plots=PLOTS)
        assert result.skipped == []
        assert sorted(p.name for p in result.pngs) == sorted(f"{p.name}.png" for p in PLOTS)
        assert all(p.stat().st_size > 5_000 for p in result.pngs)
        html = result.html.read_text()
        assert all(f'id="{p.name}"' in html for p in PLOTS)

    def test_unrecorded_limits_still_render(self, tmp_path: Path) -> None:
        path = write_ancestry_report(tmp_path / "old.json", limits=None)
        result = render_plots(path, tmp_path / "figs", plots=PLOTS)
        assert result.skipped == []

    def test_a_qc_only_report_lists_them_as_not_plotted(self, tmp_path: Path) -> None:
        path = write_report(tmp_path / "qc.json", [qc_row("callrate_prune", 1, "m")], n_input=5)
        result = render_plots(path, tmp_path / "figs", plots=PLOTS)
        assert result.pngs == []
        assert len(result.skipped) == len(PLOTS)


def test_the_retired_ancestry_plot_module_is_gone() -> None:
    """The old PNGs were replaced by this module; nothing may still import them."""
    import importlib.util

    import genotools.ancestry as ancestry_pkg

    assert importlib.util.find_spec("genotools.ancestry.plots") is None
    assert not hasattr(ancestry_pkg, "plot_ancestry_diagnostics")
