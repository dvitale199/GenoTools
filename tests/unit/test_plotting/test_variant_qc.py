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

"""Variant QC plots: the QQ reconstruction on hand-built bins, thresholds,
availability and the HWE test choice, then a render of the golden report."""

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import pytest  # noqa: E402

from genotools.plotting.plots import variant_qc as vq  # noqa: E402
from genotools.plotting.plots.variant_qc import (  # noqa: E402
    PLOTS,
    QQ_SPECS,
    hwe_filter_test,
    missingness_panels,
    qq_panels,
    qq_points,
)
from genotools.plotting.report import load_report  # noqa: E402

from .conftest import qc_row, write_report  # noqa: E402

NAMES = ["variant_missingness", "hwe_qq", "case_control_qq", "haplotype_qq"]


def bins(counts, edges=(0, 1, 2, 3, np.inf)):
    """Bins [0,1) [1,2) [2,3) [3,inf) by default, with the given counts."""
    return pd.DataFrame({
        "bin_lower": list(edges[:-1]),
        "bin_upper": list(edges[1:]),
        "count": counts,
    })


MISS_BINS = bins([4, 3, 2, 1], edges=(0, 0.005, 0.01, 0.015, 0.02))


def expected(rank, n):
    return -np.log10((rank - 0.5) / n)


class TestQQPoints:
    def test_n_ranks_and_expected_from_counts(self) -> None:
        # Most significant first: [3,inf) holds 2, [2,3) 3, [1,2) 5, [0,1) 10.
        pts = qq_points(bins([10, 5, 3, 2]))
        assert (pts["n"] == 20).all()
        by_lower = pts.set_index("bin_lower")
        assert (by_lower.loc[3.0, ["rank_first", "rank_last"]] == [1, 2]).all()
        assert (by_lower.loc[2.0, ["rank_first", "rank_last"]] == [3, 5]).all()
        assert (by_lower.loc[1.0, ["rank_first", "rank_last"]] == [6, 10]).all()
        assert (by_lower.loc[0.0, ["rank_first", "rank_last"]] == [11, 20]).all()
        # A bin sits at the middle of its ranks.
        assert by_lower.loc[2.0, "expected"] == pytest.approx(expected(4, 20))
        assert by_lower.loc[0.0, "expected"] == pytest.approx(expected(15.5, 20))
        assert by_lower.loc[2.0, "observed"] == 2.5
        assert by_lower.loc[2.0, "count"] == 3

    def test_empty_bins_are_dropped_but_ranks_skip_them(self) -> None:
        pts = qq_points(bins([6, 0, 4, 0]))
        assert list(pts["bin_lower"]) == [2.0, 0.0]
        assert list(pts["rank_first"]) == [1, 5]

    def test_open_last_bin_is_a_floor(self) -> None:
        pts = qq_points(bins([10, 5, 3, 2])).set_index("bin_lower")
        assert pts.loc[3.0, "lower_bound"]
        assert pts.loc[3.0, "observed"] == 3.0  # its lower edge, not a midpoint
        assert not pts.loc[2.0, "lower_bound"]

    def test_exact_rows_replace_the_bins_they_cover(self) -> None:
        failing = pd.DataFrame({"P": [1e-4, 1e-6], "label": ["rsB", "rsA"]})  # unsorted
        pts = qq_points(bins([10, 5, 3, 2]), failing)
        exact = pts[pts["kind"] == "exact"]
        # Sorted most extreme first, exact observed values, ranks 1 and 2.
        assert list(exact["label"]) == ["rsA", "rsB"]
        assert exact["observed"].tolist() == pytest.approx([6.0, 4.0])
        assert exact["expected"].tolist() == pytest.approx([expected(1, 20), expected(2, 20)])
        # They stand for the two tests in the open bin, which is now gone.
        assert 3.0 not in set(pts["bin_lower"].dropna())
        assert pts["count"].sum() == 20  # no test drawn twice, none lost
        assert not exact["lower_bound"].any()

    def test_partly_covered_bin_keeps_only_the_remaining_ranks(self) -> None:
        failing = pd.DataFrame({"P": [1e-9, 1e-8, 1e-7, 1e-4], "label": list("abcd")})
        pts = qq_points(bins([10, 5, 3, 2]), failing)  # 4 exact: ranks 1-4
        rest = pts[pts["kind"] == "bin"].set_index("bin_lower")
        assert (rest.loc[2.0, ["rank_first", "rank_last", "count"]] == [5, 5, 1]).all()
        assert pts["count"].sum() == 20

    def test_infinite_p_is_pinned_and_flagged(self) -> None:
        failing = pd.DataFrame({"P": [0.0, 1e-60], "label": ["zero", "tiny"]})
        pts = qq_points(bins([10, 5, 3, 2]), failing)
        exact = pts[pts["kind"] == "exact"].set_index("label")
        assert np.isfinite(pts["observed"]).all()
        assert exact.loc["zero", "lower_bound"]
        assert exact.loc["zero", "observed"] >= exact.loc["tiny", "observed"] == pytest.approx(60)
        assert not exact.loc["tiny", "lower_bound"]

    def test_failing_table_longer_than_the_data_is_clipped(self) -> None:
        failing = pd.DataFrame({"P": [1e-3] * 30, "label": ["x"] * 30})
        pts = qq_points(bins([1, 1, 1, 1]), failing)
        assert pts["count"].sum() == 4

    def test_ties_keep_input_order(self) -> None:
        failing = pd.DataFrame({"P": [1e-5] * 4, "label": list("abcd")})
        pts = qq_points(bins([10, 5, 3, 2]), failing)
        assert list(pts[pts["kind"] == "exact"]["label"]) == list("abcd")

    def test_no_tests_gives_an_empty_frame(self) -> None:
        assert qq_points(bins([0, 0, 0, 0])).empty


class TestHweTest:
    def hwe_bins(self, tests):
        return pd.concat([bins([4, 3, 2, 1]).assign(test=t) for t in tests], ignore_index=True)

    def test_controls_when_present(self) -> None:
        assert hwe_filter_test(self.hwe_bins(["ALL", "AFF", "UNAFF"])) == "UNAFF"

    def test_all_when_no_controls_test(self) -> None:
        assert hwe_filter_test(self.hwe_bins(["ALL"])) == "ALL"

    def test_panel_uses_only_that_tests_bins_and_failing_rows(self, tmp_path: Path) -> None:
        report = _write_tables(tmp_path, {
            "hwe_pvalue_bins": self.hwe_bins(["ALL", "AFF", "UNAFF"]).assign(
                count=lambda f: np.where(f["test"] == "UNAFF", f["count"], f["count"] * 100)
            ),
            "hwe_failing": pd.DataFrame(
                {"SNP": ["rs1", "rs1"], "test": ["UNAFF", "ALL"], "P": [1e-5, 1e-5]}
            ),
        }, params={("hwe_prune", "hwe_threshold"): 1e-4})
        (panel,) = qq_panels(load_report(report), QQ_SPECS["hwe"])
        assert panel.test == "UNAFF"
        assert panel.points["n"].iloc[0] == 10  # UNAFF's total, not ALL's 1000
        assert panel.removed == 1


def _write_tables(tmp_path: Path, tables, groups=("all",), params=None) -> Path:
    """A report with `{out}_plot_data/{group}/{name}.tsv` and a parameters table."""
    root = tmp_path / "run_plot_data"
    for group in groups:
        (root / group).mkdir(parents=True, exist_ok=True)
        for name, frame in tables.items():
            frame.to_csv(root / group / f"{name}.tsv", sep="\t", index=False)
    rows = [
        {"step": s, "parameter": p, "value": v, "ancestry": g, "source": "requested"}
        for (s, p), v in (params or {}).items() for g in groups
    ]
    report = {
        "QC": pd.DataFrame([qc_row("geno_prune", 1, "geno_removed_count", level="variant")]).to_dict(),
        "plot_data": str(root),
        "parameters": pd.DataFrame(rows).to_dict() if rows else {},
    }
    path = tmp_path / "run.json"
    path.write_text(json.dumps(report))
    return path


class TestThresholds:
    def test_each_group_gets_its_own_threshold(self, tmp_path: Path) -> None:
        report = _write_tables(tmp_path, {"hwe_pvalue_bins": bins([4, 3, 2, 1]).assign(test="ALL")},
                               groups=("all", "EUR"))
        raw = json.loads(report.read_text())
        raw["parameters"] = pd.DataFrame([
            {"step": "hwe_prune", "parameter": "hwe_threshold", "value": 1e-4, "ancestry": "all", "source": "requested"},
            {"step": "hwe_prune", "parameter": "hwe_threshold", "value": 1e-8, "ancestry": "EUR", "source": "requested"},
        ]).to_dict()
        report.write_text(json.dumps(raw))
        panels = {p.ancestry: p.threshold for p in qq_panels(load_report(report), QQ_SPECS["hwe"])}
        assert panels == {"all": 1e-4, "EUR": 1e-8}

    def test_missing_parameter_means_no_line(self, tmp_path: Path) -> None:
        report = _write_tables(tmp_path, {"geno_missingness_bins": MISS_BINS})
        (panel,) = missingness_panels(load_report(report))
        assert panel.threshold is None
        out = vq.missingness_png(load_report(report), plt, tmp_path / "m.png")
        assert out.stat().st_size > 0

    def test_geno_threshold_comes_from_the_report(self, tmp_path: Path) -> None:
        report = _write_tables(tmp_path, {"geno_missingness_bins": MISS_BINS},
                               params={("geno_prune", "geno"): 0.02})
        (panel,) = missingness_panels(load_report(report))
        assert panel.threshold == 0.02

    @pytest.mark.parametrize("name", ["case_control", "haplotype"])
    def test_step_specific_threshold_names(self, tmp_path: Path, name: str) -> None:
        spec = QQ_SPECS[name]
        report = _write_tables(tmp_path, {spec.bins_table: bins([4, 3, 2, 1])},
                               params={(spec.step, "p_threshold"): 0.003})
        (panel,) = qq_panels(load_report(report), spec)
        assert panel.threshold == 0.003


class TestAvailable:
    @pytest.mark.parametrize("name", NAMES)
    def test_false_without_tables(self, ancestry_report: Path, name: str) -> None:
        plot = next(p for p in PLOTS if p.name == name)
        assert plot.available(load_report(ancestry_report)) is False

    @pytest.mark.parametrize("name", NAMES)
    def test_false_when_plot_data_is_missing_on_disk(self, golden_variant_report: Path, name: str, tmp_path: Path) -> None:
        raw = json.loads(golden_variant_report.read_text())
        raw["plot_data"] = str(tmp_path / "gone")
        golden_variant_report.write_text(json.dumps(raw))
        for p in golden_variant_report.parent.glob("output_plot_data"):
            p.rename(p.with_name("moved"))
        plot = next(p for p in PLOTS if p.name == name)
        assert plot.available(load_report(golden_variant_report)) is False

    def test_a_step_that_did_not_run_is_skipped_alone(self, golden_variant_report: Path) -> None:
        for f in (golden_variant_report.parent / "output_plot_data" / "all").glob("hwe_*.tsv"):
            f.unlink()
        report = load_report(golden_variant_report)
        avail = {p.name: p.available(report) for p in PLOTS if p.name in NAMES}
        assert avail == {"variant_missingness": True, "hwe_qq": False,
                         "case_control_qq": True, "haplotype_qq": True}

    def test_registered(self) -> None:
        assert [p.name for p in PLOTS] == NAMES


class TestGolden:
    def test_hwe_exact_tail_matches_the_failing_table(self, golden_variant_report: Path) -> None:
        report = load_report(golden_variant_report)
        (panel,) = qq_panels(report, QQ_SPECS["hwe"])
        failing = pd.read_csv(golden_variant_report.parent / "output_plot_data/all/hwe_failing.tsv", sep="\t")
        pts = panel.points
        assert (pts["kind"] == "exact").sum() == len(failing) == panel.removed
        assert pts["count"].sum() == pts["n"].iloc[0]
        # Every point past the threshold is exact, none are binned.
        past = pts[pts["observed"] > -np.log10(panel.threshold) + 0.05]
        assert (past["kind"] == "exact").all()

    def test_missingness_panel_counts_removed_variants(self, golden_variant_report: Path) -> None:
        (panel,) = missingness_panels(load_report(golden_variant_report))
        assert panel.threshold == 0.05
        assert panel.removed == 225

    @pytest.mark.parametrize("name", NAMES)
    def test_png_and_plotly_render(self, golden_variant_report: Path, tmp_path: Path, name: str) -> None:
        report = load_report(golden_variant_report)
        plot = next(p for p in PLOTS if p.name == name)
        assert plot.available(report)
        assert plot.png(report, plt, tmp_path / f"{name}.png").stat().st_size > 0
        fig = plot.plotly(report)
        assert len(fig.data) > 0
