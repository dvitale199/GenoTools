"""Frames the plots read, built from hand-built reports."""

import json
from pathlib import Path

import pandas as pd
import pytest

from genotools.plotting.report import attrition_frame, load_report, outcome_frame

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
