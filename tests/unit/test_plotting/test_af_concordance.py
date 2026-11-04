"""Allele-frequency concordance plot, against a hand-built retained table."""

import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from genotools.plotting import render_plots
from genotools.plotting.plots import ancestry as plots
from genotools.plotting.plots.ancestry import PLOTS
from genotools.plotting.report import load_report

from .conftest import qc_row, write_report

PLOT = next(p for p in PLOTS if p.name == "ancestry_af_concordance")
N_CONCORDANT, N_FLIP, N_DISCORDANT = 80, 6, 4


def af_table() -> pd.DataFrame:
    rng = np.random.default_rng(3)
    ref = np.linspace(0.05, 0.95, N_CONCORDANT)
    cohort = np.clip(ref + rng.normal(0, 0.02, N_CONCORDANT), 0, 1)
    flip_ref = np.array([0.1, 0.2, 0.3, 0.7, 0.8, 0.9])
    disc_ref = np.array([0.2, 0.4, 0.6, 0.8])
    disc_cohort = np.array([0.5, 0.9, 0.05, 0.4])
    r = np.concatenate([ref, flip_ref, disc_ref])
    c = np.concatenate([cohort, 1 - flip_ref, disc_cohort])
    n = len(r)
    return pd.DataFrame({
        "snp": [f"rs{i}" for i in range(n)],
        "counted_allele": "A", "other_allele": "G", "cohort_a1": "A", "cohort_a2": "G",
        "ref_freq": r, "cohort_freq": c, "abs_delta": np.abs(c - r),
        "discordant": np.abs(c - r) > 0.2,
        "flip_signature": [False] * N_CONCORDANT + [True] * N_FLIP + [False] * N_DISCORDANT,
    })


def make_report(tmp_path: Path, table=None, diagnostics=None, text_bools=False) -> Path:
    table = af_table() if table is None else table
    path = write_report(tmp_path / "run.json", [qc_row("callrate_prune", 1, "outlier_count")])
    out = tmp_path / "run_plot_data" / "ancestry"
    out.mkdir(parents=True)
    if text_bools:
        table = table.astype({"discordant": str, "flip_signature": str})
    table.to_csv(out / "allele_frequencies.tsv", sep="\t", index=False)
    raw = json.loads(path.read_text())
    raw["plot_data"] = str(tmp_path / "run_plot_data")
    if diagnostics:
        raw["ancestry_diagnostics"] = diagnostics
    path.write_text(json.dumps(raw))
    return path


@pytest.mark.parametrize("text_bools", [False, True])
def test_frame_kinds(tmp_path, text_bools):
    report = load_report(make_report(tmp_path, text_bools=text_bools))
    frame = plots.af_frame(report)
    counts = frame["kind"].value_counts().to_dict()
    assert counts == {"concordant": N_CONCORDANT, "flip": N_FLIP, "discordant": N_DISCORDANT}


def test_summary_from_sample_matches_hand_computation(tmp_path):
    report = load_report(make_report(tmp_path))
    frame = plots.af_frame(report)
    s = plots.af_summary(report, frame)
    t = af_table()
    assert s["source"] == "sample"
    assert s["correlation"] == pytest.approx(np.corrcoef(t.ref_freq, t.cohort_freq)[0, 1])
    assert s["flip_fraction"] == pytest.approx(N_FLIP / len(t))


def test_summary_prefers_report_numbers(tmp_path):
    diag = {"allele_frequency": {"n_snps": 90000, "correlation": 0.987,
                                 "fraction_flip_signature": 0.0123}}
    report = load_report(make_report(tmp_path, diagnostics=diag))
    s = plots.af_summary(report, plots.af_frame(report))
    assert s["source"] == "report" and s["correlation"] == 0.987
    note = plots.af_note(s)
    assert "0.987" in note and "1.23%" in note and "90,000" in note and "retained" in note


def test_summary_constant_column_has_no_correlation(tmp_path):
    t = af_table()
    t["cohort_freq"] = 0.5
    report = load_report(make_report(tmp_path, table=t))
    assert plots.af_summary(report, plots.af_frame(report))["correlation"] is None


def test_unavailable_without_table_or_columns(tmp_path):
    bare = load_report(write_report(tmp_path / "bare.json", [qc_row("callrate_prune", 1, "outlier_count")]))
    assert not PLOT.available(bare)
    sub = tmp_path / "x"
    sub.mkdir()
    partial = af_table().drop(columns=["flip_signature"])
    assert not PLOT.available(load_report(make_report(sub, table=partial)))


def test_available_with_table(tmp_path):
    assert PLOT.available(load_report(make_report(tmp_path)))


def test_plotly_traces_and_hover(tmp_path):
    report = load_report(make_report(tmp_path))
    fig = PLOT.plotly(report)
    names = [t.name for t in fig.data]
    assert any(n.startswith("flip signature") for n in names)
    flip = next(t for t in fig.data if t.name.startswith("flip signature"))
    assert len(flip.x) == N_FLIP
    assert "rs80" in flip.text[0] and "ref freq" in flip.text[0] and "counted_allele" in flip.text[0]


def test_png_and_render(tmp_path):
    report = load_report(make_report(tmp_path))
    out = tmp_path / "af.png"
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    PLOT.png(report, plt, out)
    assert out.stat().st_size > 1000
