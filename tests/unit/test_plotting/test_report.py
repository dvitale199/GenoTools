"""Frames the plots read, built from hand-built reports."""

import json
from pathlib import Path

import pandas as pd
import pytest

from genotools.plotting.report import (
    assoc_prefix,
    attrition_frame,
    load_report,
    outcome_frame,
    parameter,
    plot_data_dir,
    plot_tables,
)

from .conftest import qc_row, write_report

GOLDEN_SAMPLE = Path(__file__).parents[2] / "regression/golden/all_sample/output.json"


class TestLoadReport:
    def test_missing_report(self, tmp_path: Path) -> None:
        with pytest.raises(FileNotFoundError):
            load_report(tmp_path / "nope.json")

    def test_not_json(self, tmp_path: Path) -> None:
        path = tmp_path / "bad.json"
        path.write_text("{not json")
        with pytest.raises(ValueError, match="not valid JSON"):
            load_report(path)

    def test_not_an_object(self, tmp_path: Path) -> None:
        path = tmp_path / "list.json"
        path.write_text("[1, 2]")
        with pytest.raises(ValueError, match="not a GenoTools report"):
            load_report(path)

    def test_qc_keeps_pipeline_order_past_ten_rows(self, tmp_path: Path) -> None:
        """String row indices must not be re-sorted ("10" before "2")."""
        qc = [qc_row(f"step{i}_prune", i, "outlier_count") for i in range(12)]
        report = load_report(write_report(tmp_path / "r.json", qc, n_input=50))
        assert list(report.qc["step"]) == [f"step{i}_prune" for i in range(12)]

    def test_legacy_report_without_outcome(self, tmp_path: Path) -> None:
        """Pre-2.0 reports carry only `pass`."""
        path = tmp_path / "old.json"
        qc = pd.DataFrame([
            {"step": "sex_prune", "pruned_count": 0, "metric": "outlier_count",
             "ancestry": "all", "level": "sample", "pass": True},
            {"step": "het_prune", "pruned_count": 0, "metric": "outlier_count",
             "ancestry": "all", "level": "sample", "pass": False},
        ])
        path.write_text(json.dumps({"QC": qc.to_dict()}))
        report = load_report(path)
        assert list(report.qc["outcome"]) == ["pass", "fail"]

    def test_report_without_qc(self, tmp_path: Path) -> None:
        path = tmp_path / "anc_only.json"
        path.write_text(json.dumps({"run_info": {"version": "2.2.1"}}))
        report = load_report(path)
        assert report.qc.empty
        assert attrition_frame(report).empty
        assert outcome_frame(report).empty


class TestAttritionFrame:
    def test_golden_flat_run(self) -> None:
        """A real report: 500 in, then callrate 12, sex 17, het 8, related 0+4."""
        frame = attrition_frame(load_report(GOLDEN_SAMPLE))
        assert list(frame["step"]) == [
            "callrate_prune", "sex_prune", "het_prune", "related_prune"
        ]
        assert list(frame["removed"]) == [12, 17, 8, 4]
        assert list(frame["remaining"]) == [488, 471, 463, 459]
        assert set(frame["start"]) == {500}

    def test_groups_start_from_their_predicted_count(self, ancestry_report: Path) -> None:
        frame = attrition_frame(load_report(ancestry_report))
        samples = frame[frame["level"] == "sample"].set_index(["ancestry", "step"])
        assert samples.loc[("EUR", "callrate_prune"), "start"] == 200
        assert samples.loc[("AMR", "callrate_prune"), "start"] == 40
        # related_prune's two metrics are summed: 3 related + 1 duplicated.
        assert samples.loc[("EUR", "related_prune"), "removed"] == 4
        assert samples.loc[("EUR", "related_prune"), "remaining"] == 200 - 5 - 2 - 4

    def test_skipped_step_keeps_its_outcome(self, ancestry_report: Path) -> None:
        frame = attrition_frame(load_report(ancestry_report))
        het = frame[(frame["ancestry"] == "AMR") & (frame["step"] == "het_prune")]
        assert het["outcome"].item() == "skipped"
        assert het["removed"].item() == 0

    def test_variants_have_no_start(self, ancestry_report: Path) -> None:
        """The report does not record the input variant count."""
        frame = attrition_frame(load_report(ancestry_report))
        variants = frame[frame["level"] == "variant"]
        assert variants["start"].isna().all()
        assert variants["remaining"].isna().all()
        assert sorted(variants["removed"]) == [20, 100]


class TestOutcomeFrame:
    def test_one_row_per_step_and_group(self, ancestry_report: Path) -> None:
        frame = outcome_frame(load_report(ancestry_report))
        assert len(frame) == 8  # related_prune's two metric rows collapse
        assert not frame.duplicated(["step", "ancestry"]).any()

    def test_reason_survives(self, ancestry_report: Path) -> None:
        frame = outcome_frame(load_report(ancestry_report)).set_index(["step", "ancestry"])
        assert frame.loc[("het_prune", "AMR"), "outcome"] == "skipped"
        assert "50-sample floor" in frame.loc[("het_prune", "AMR"), "reason"]
        assert frame.loc[("het_prune", "EUR"), "reason"] is None


def _params(rows):
    return pd.DataFrame(
        rows, columns=["step", "parameter", "value", "ancestry", "source"]
    ).to_dict()


class TestParameter:
    def test_reads_the_golden_thresholds(self, golden_sample_report: Path) -> None:
        report = load_report(golden_sample_report)
        assert parameter(report, "callrate_prune", "mind") == 0.05
        assert parameter(report, "sex_prune", "male_min_f") == 0.75
        assert parameter(report, "het_prune", "het_lower") == -0.15

    def test_prefers_resolved_over_requested(self, tmp_path: Path) -> None:
        path = tmp_path / "r.json"
        path.write_text(json.dumps({"parameters": _params([
            ("het_prune", "het_lower", -0.25, "EUR", "requested"),
            ("het_prune", "het_lower", -0.12, "EUR", "resolved"),
        ])}))
        report = load_report(path)
        assert parameter(report, "het_prune", "het_lower", "EUR") == -0.12
        assert parameter(report, "het_prune", "het_lower", "EUR", source="requested") == -0.25

    def test_is_per_ancestry(self, tmp_path: Path) -> None:
        path = tmp_path / "r.json"
        path.write_text(json.dumps({"parameters": _params([
            ("callrate_prune", "mind", 0.02, "EUR", "requested"),
            ("callrate_prune", "mind", 0.05, "AFR", "requested"),
        ])}))
        report = load_report(path)
        assert parameter(report, "callrate_prune", "mind", "AFR") == 0.05
        assert parameter(report, "callrate_prune", "mind", "all") is None

    def test_missing_is_none(self, golden_sample_report: Path, tmp_path: Path) -> None:
        assert parameter(load_report(golden_sample_report), "hwe_prune", "hwe_threshold") is None
        # Recorded as None (het_sd in fixed mode) also reads as None.
        assert parameter(load_report(golden_sample_report), "het_prune", "het_sd") is None
        path = tmp_path / "r.json"
        path.write_text("{}")
        assert parameter(load_report(path), "callrate_prune", "mind") is None


class TestPlotTables:
    def test_reads_the_golden_tables(self, golden_sample_report: Path) -> None:
        tables = plot_tables(load_report(golden_sample_report), "het_f")
        assert list(tables) == ["all"]
        assert list(tables["all"].columns) == ["#FID", "IID", "F"]
        assert len(tables["all"]) > 0

    def test_missing_table_is_empty(self, golden_sample_report: Path) -> None:
        assert plot_tables(load_report(golden_sample_report), "hwe_failing") == {}

    def test_groups_all_first_then_alphabetical(self, tmp_path: Path) -> None:
        root = tmp_path / "cohort_plot_data"
        for group in ("SAS", "AFR", "EUR"):
            (root / group).mkdir(parents=True)
            (root / group / "het_f.tsv").write_text("#FID\tIID\tF\na\ta\t0.1\n")
        (root / "CAH").mkdir()  # skipped het: no table, no entry
        report = tmp_path / "cohort.json"
        report.write_text(json.dumps({"plot_data": str(root)}))
        assert list(plot_tables(load_report(report), "het_f")) == ["AFR", "EUR", "SAS"]

    def test_falls_back_to_beside_the_report(self, tmp_path: Path) -> None:
        """The recorded path is relative to the run's cwd, which may not be ours."""
        (tmp_path / "run.v2_plot_data" / "all").mkdir(parents=True)
        report = tmp_path / "run.v2.json"
        report.write_text(json.dumps({"plot_data": "elsewhere/run.v2_plot_data"}))
        assert plot_data_dir(load_report(report)) == tmp_path / "run.v2_plot_data"

    def test_no_plot_data(self, tmp_path: Path) -> None:
        report = tmp_path / "r.json"
        report.write_text("{}")
        assert plot_data_dir(load_report(report)) is None
        assert plot_tables(load_report(report), "het_f") == {}


class TestAssocPrefix:
    def test_flat_run(self, tmp_path: Path) -> None:
        path = tmp_path / "r.json"
        path.write_text(json.dumps({"pass_fail": {"assoc": {"output": "out/run_assoc"}}}))
        assert assoc_prefix(load_report(path)) == "out/run_assoc"

    def test_per_ancestry(self, tmp_path: Path) -> None:
        path = tmp_path / "r.json"
        path.write_text(json.dumps({
            "EUR_pass_fail": {"assoc": {"output": "out/EUR_assoc"}},
            "AFR_pass_fail": {"geno": {"output": "out/AFR_geno"}},
        }))
        report = load_report(path)
        assert assoc_prefix(report, "EUR") == "out/EUR_assoc"
        assert assoc_prefix(report, "AFR") is None
        assert assoc_prefix(report) is None

    def test_golden_qc_run_has_no_assoc(self, golden_qc_report: Path) -> None:
        assert assoc_prefix(load_report(golden_qc_report)) is None
