"""Sample QC plots: the pure panel logic on hand-built reports, then rendering."""

import json
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np
import pandas as pd
import pytest

from genotools.plotting import render_plots
from genotools.plotting.plots import sample_qc
from genotools.plotting.plots.sample_qc import (
    PLOTS,
    above_threshold,
    het_bounds,
    het_panels,
    kinship_panels,
    missingness_panels,
    outside_bounds,
    sex_flagged,
    sex_pass,
    sex_panels,
)
from genotools.plotting.report import load_report

from .conftest import qc_row, write_report

PLOT_NAMES = ["sample_missingness", "sex_check", "heterozygosity", "relatedness"]


def _param(step: str, name: str, value, ancestry: str = "all", source: str = "requested") -> dict:
    return {"step": step, "parameter": name, "value": value, "ancestry": ancestry, "source": source}


def _report(
    tmp_path: Path,
    tables: Dict[str, Dict[str, pd.DataFrame]],
    params: Optional[List[dict]] = None,
):
    """A report on disk with `{group: {table: frame}}` plot tables and parameters."""
    path = write_report(tmp_path / "r.json", [qc_row("callrate_prune", 0, "outlier_count")], n_input=10)
    raw = json.loads(path.read_text())
    root = tmp_path / "r_plot_data"
    for group, named in tables.items():
        (root / group).mkdir(parents=True, exist_ok=True)
        for name, frame in named.items():
            frame.to_csv(root / group / f"{name}.tsv", sep="\t", index=False)
    if tables:
        raw["plot_data"] = str(root)
    if params is not None:
        raw["parameters"] = pd.DataFrame(params).to_dict()
    path.write_text(json.dumps(raw))
    return load_report(path)


def _miss(values) -> pd.DataFrame:
    return pd.DataFrame(
        {"#FID": [f"f{i}" for i in range(len(values))],
         "IID": [f"i{i}" for i in range(len(values))], "F_MISS": values}
    )


def _het(values) -> pd.DataFrame:
    return _miss(values).rename(columns={"F_MISS": "F"})


class TestThresholdLogic:
    def test_above_is_strict_and_needs_a_bound(self) -> None:
        frame = _miss([0.01, 0.05, 0.06])
        assert list(above_threshold(frame, "F_MISS", 0.05)["F_MISS"]) == [0.06]
        assert above_threshold(frame, "F_MISS", None).empty

    def test_outside_bounds_flags_each_side_and_skips_absent_ones(self) -> None:
        frame = _het([-0.5, -0.15, 0.0, 0.15, 0.4])
        assert list(outside_bounds(frame, "F", -0.15, 0.15)["F"]) == [-0.5, 0.4]
        assert list(outside_bounds(frame, "F", None, 0.15)["F"]) == [0.4]
        assert outside_bounds(frame, "F", None, None).empty

    def test_sex_pass_prefers_filtered_x_and_falls_back(self) -> None:
        both = pd.DataFrame({"F": [1.0, 2.0, 3.0], "check": ["all_x", "filtered_x", "filtered_x"]})
        assert list(sex_pass(both)["F"]) == [2.0, 3.0]
        only_all = both[both["check"] == "all_x"]
        assert list(sex_pass(only_all)["F"]) == [1.0]


class TestMissingness:
    def test_flags_samples_above_the_reported_mind(self, tmp_path: Path) -> None:
        report = _report(
            tmp_path, {"all": {"callrate_sample_missingness": _miss([0.01, 0.02, 0.2, 0.3])}},
            [_param("callrate_prune", "mind", 0.1)],
        )
        (panel,) = missingness_panels(report)
        assert panel.lines == [(0.1, "mind 0.1")]
        assert sorted(panel.flagged_x) == [0.2, 0.3]
        assert panel.title == "all: 2 of 4 above mind 0.1"
        assert any("f2 i2" in h for h in panel.flagged_hover)

    def test_missing_parameter_draws_no_line_and_flags_nothing(self, tmp_path: Path) -> None:
        report = _report(
            tmp_path, {"all": {"callrate_sample_missingness": _miss([0.01, 0.9])}}, params=[]
        )
        (panel,) = missingness_panels(report)
        assert panel.lines == [] and len(panel.flagged_x) == 0
        assert panel.title == "all: 2 samples"

    def test_each_group_uses_its_own_threshold(self, tmp_path: Path) -> None:
        tables = {g: {"callrate_sample_missingness": _miss([0.03, 0.2])} for g in ("all", "EUR")}
        report = _report(
            tmp_path, tables,
            [_param("callrate_prune", "mind", 0.02, "all"), _param("callrate_prune", "mind", 0.5, "EUR")],
        )
        by_group = {p.ancestry: p for p in missingness_panels(report)}
        assert list(by_group) == ["all", "EUR"]
        assert len(by_group["all"].flagged_x) == 2
        assert len(by_group["EUR"].flagged_x) == 0


class TestHet:
    def test_bounds_are_resolved_per_group(self, tmp_path: Path) -> None:
        tables = {g: {"het_f": _het([-0.3, 0.0, 0.3])} for g in ("all", "EUR", "AFR")}
        params = [
            _param("het_prune", "het_lower", -0.15, "all", "requested"),
            _param("het_prune", "het_upper", 0.15, "all", "requested"),
            _param("het_prune", "het_lower", -0.1, "EUR", "resolved"),
            _param("het_prune", "het_upper", 0.1, "EUR", "resolved"),
            _param("het_prune", "het_lower", -0.5, "AFR", "resolved"),
            _param("het_prune", "het_upper", 0.5, "AFR", "resolved"),
            # A requested bound the resolved one overrides.
            _param("het_prune", "het_lower", -0.15, "AFR", "requested"),
        ]
        report = _report(tmp_path, tables, params)
        assert het_bounds(report, "AFR") == (-0.5, 0.5)
        panels = {p.ancestry: p for p in het_panels(report)}
        assert panels["EUR"].band == (-0.1, 0.1)
        assert len(panels["EUR"].flagged_x) == 2
        assert panels["AFR"].band == (-0.5, 0.5)
        assert len(panels["AFR"].flagged_x) == 0
        assert len(panels["all"].flagged_x) == 2

    def test_one_missing_bound_draws_one_line_and_no_band(self, tmp_path: Path) -> None:
        report = _report(
            tmp_path, {"all": {"het_f": _het([-0.9, 0.0, 0.9])}},
            [_param("het_prune", "het_upper", 0.15, source="resolved")],
        )
        (panel,) = het_panels(report)
        assert panel.band is None
        assert [x for x, _ in panel.lines] == [0.15]
        assert list(panel.flagged_x) == [0.9]


class TestSex:
    def _table(self) -> pd.DataFrame:
        rows = []
        for check in ("all_x", "filtered_x"):
            rows += [
                ("a", "a", 1, 1, "OK", 1.0, check),
                ("b", "b", 2, 2, "OK", 0.0, check),
                ("c", "c", 2, 1, "PROBLEM", 0.9, check),
            ]
        rows.append(("d", "d", 0, 0, "OK", 0.5, "filtered_x"))
        return pd.DataFrame(rows, columns=["#FID", "IID", "PEDSEX", "SNPSEX", "STATUS", "F", "check"])

    def test_one_pass_split_by_reported_sex_with_both_cutoffs(self, tmp_path: Path) -> None:
        report = _report(
            tmp_path, {"all": {"sex_check": self._table()}},
            [_param("sex_prune", "female_max_f", 0.25), _param("sex_prune", "male_min_f", 0.75)],
        )
        (panel,) = sex_panels(report)
        assert {s.name: len(s.values) for s in panel.series} == {"female": 2, "male": 1, "unknown": 1}
        assert panel.lines == [(0.25, "female max 0.25"), (0.75, "male min 0.75")]
        assert list(panel.flagged_x) == [0.9]
        assert "reported female" in panel.flagged_hover[0]
        assert "flagged in all_x, filtered_x" in panel.flagged_hover[0]
        assert panel.title.startswith("all: 1 of 4 ")  # filtered_x only: 4 samples

    def test_missing_cutoffs_draw_no_lines(self, tmp_path: Path) -> None:
        report = _report(tmp_path, {"all": {"sex_check": self._table()}}, params=[])
        (panel,) = sex_panels(report)
        assert panel.lines == []


class TestSexFlagged:
    COLS = ["#FID", "IID", "PEDSEX", "SNPSEX", "STATUS", "F", "check"]

    def test_union_of_both_passes_matches_the_step(self) -> None:
        table = pd.DataFrame(
            [
                ("a", "a", 2, 1, "PROBLEM", 0.9, "all_x"),
                ("a", "a", 2, 1, "PROBLEM", 0.95, "filtered_x"),  # both passes
                ("b", "b", 2, 1, "PROBLEM", 0.8, "all_x"),        # all_x only
                ("b", "b", 2, 2, "OK", 0.1, "filtered_x"),
                ("c", "c", 1, 2, "OK", 0.0, "all_x"),
                ("c", "c", 1, 2, "PROBLEM", 0.0, "filtered_x"),   # filtered_x only
                ("d", "d", 1, 1, "OK", 1.0, "all_x"),
                ("d", "d", 1, 1, "OK", 1.0, "filtered_x"),
            ],
            columns=self.COLS,
        )
        flagged = sex_flagged(table).set_index("IID")
        assert sorted(flagged.index) == ["a", "b", "c"]
        # Each sits at its F in the drawn (filtered_x) pass.
        assert flagged.loc["b", "F"] == 0.1
        assert flagged.loc["a", "F"] == 0.95
        assert flagged.loc["a", "flagged_in"] == "all_x, filtered_x"
        assert flagged.loc["b", "flagged_in"] == "all_x"
        assert flagged.loc["c", "flagged_in"] == "filtered_x"

    def test_sample_flagged_only_by_the_other_pass_still_appears(self) -> None:
        table = pd.DataFrame(
            [
                ("a", "a", 2, 1, "PROBLEM", 0.8, "all_x"),
                ("z", "z", 1, 1, "OK", 1.0, "filtered_x"),
            ],
            columns=self.COLS,
        )
        flagged = sex_flagged(table)
        assert list(flagged["IID"]) == ["a"] and list(flagged["F"]) == [0.8]


class TestKinship:
    def _pairs(self) -> pd.DataFrame:
        return pd.DataFrame(
            {"#FID1": ["a", "c"], "IID1": ["a", "c"], "FID2": ["b", "d"], "IID2": ["b", "d"],
             "IBS0": [0.002, 0.01], "KINSHIP": [0.49, 0.2]}
        )

    def test_cutoffs_come_from_the_report(self, tmp_path: Path) -> None:
        report = _report(
            tmp_path, {"all": {"related_kinship": self._pairs()}},
            [_param("related_prune", "related_cutoff", 0.125),
             _param("related_prune", "duplicated_cutoff", 0.3)],
        )
        (panel,) = kinship_panels(report)
        assert [c for c, _ in panel.cutoffs] == [0.125, 0.3]
        assert "a a / b b" in panel.hover[0]

    def test_missing_cutoffs_draw_none(self, tmp_path: Path) -> None:
        report = _report(tmp_path, {"all": {"related_kinship": self._pairs()}}, params=[])
        assert kinship_panels(report)[0].cutoffs == []

    def test_degree_bands_are_king_boundaries(self) -> None:
        assert [b for b, _ in sample_qc.KING_DEGREES] == [0.354, 0.177, 0.0884, 0.0442]


class TestAvailable:
    @pytest.mark.parametrize("plot", PLOTS, ids=lambda p: p.name)
    def test_false_without_the_table_and_never_raises(self, plot, ancestry_report: Path) -> None:
        assert plot.available(load_report(ancestry_report)) is False

    def test_true_when_the_table_exists(self, tmp_path: Path) -> None:
        report = _report(tmp_path, {"all": {"het_f": _het([0.0])}})
        avail = {p.name: p.available(report) for p in PLOTS}
        assert avail == {"sample_missingness": False, "sex_check": False,
                         "heterozygosity": True, "relatedness": False}

    def test_registers_all_four(self) -> None:
        assert [p.name for p in PLOTS] == PLOT_NAMES


class TestRendering:
    def test_multi_group_report_renders_png_and_html(self, tmp_path: Path) -> None:
        tables = {
            g: {
                "callrate_sample_missingness": _miss(list(np.linspace(0.0, 0.2, 30))),
                "het_f": _het(list(np.linspace(-0.4, 0.4, 30))),
                "related_kinship": TestKinship()._pairs(),
                "sex_check": TestSex()._table(),
            }
            for g in ("all", "EUR", "AMR", "CAH")
        }
        params = [
            _param(s, n, v, g)
            for g in tables
            for s, n, v in [("callrate_prune", "mind", 0.1), ("het_prune", "het_lower", -0.2),
                            ("het_prune", "het_upper", 0.2)]
        ]
        report = _report(tmp_path, tables, params)
        result = render_plots(report.path, tmp_path / "out")
        names = {p.stem for p in result.pngs}
        assert set(PLOT_NAMES) <= names
        html = result.html.read_text()
        for name in PLOT_NAMES:
            assert f'id="{name}"' in html

    def test_golden_report_renders_every_plot(self, golden_sample_report: Path, tmp_path: Path) -> None:
        result = render_plots(golden_sample_report, tmp_path / "out")
        names = {p.stem for p in result.pngs}
        assert set(PLOT_NAMES) <= names
        assert not [s for s in result.skipped if any(n in s.lower() for n in ("missing", "sex", "het", "relat"))]

    def test_golden_counts_match_the_tables(self, golden_sample_report: Path) -> None:
        report = load_report(golden_sample_report)
        (miss,) = missingness_panels(report)
        assert len(miss.flagged_x) == 12
        (sex,) = sex_panels(report)
        assert len(sex.flagged_x) == 17  # what filter_sex pruned
        (het,) = het_panels(report)
        assert len(het.flagged_x) == 8
        (kin,) = kinship_panels(report)
        assert [c for c, _ in kin.cutoffs] == [0.0884, 0.354]
