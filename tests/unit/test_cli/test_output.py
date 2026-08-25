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

"""Tests for CLI output formatting."""

import json
from pathlib import Path
from typing import Any, Dict, List

import numpy as np
import pandas as pd
import pytest

from genotools.cli.output import (
    QCMetrics,
    GWASMetrics,
    PipelineOutput,
    write_results,
    build_metrics_dataframe,
    build_gwas_dataframe,
)


class TestQCMetrics:
    """Tests for QCMetrics dataclass."""

    def test_create_metrics(self) -> None:
        """QCMetrics can be created."""
        metrics = QCMetrics(
            step="callrate_prune",
            count=10,
            metric="outlier_count",
            ancestry="EUR",
            level="sample",
            passed=True,
        )
        assert metrics.step == "callrate_prune"
        assert metrics.count == 10
        assert metrics.ancestry == "EUR"

    def test_default_values(self) -> None:
        """QCMetrics has correct defaults."""
        metrics = QCMetrics(
            step="test",
            count=5,
            metric="count",
        )
        assert metrics.ancestry == "all"
        assert metrics.level == "sample"
        assert metrics.passed is True


class TestGWASMetrics:
    """Tests for GWASMetrics dataclass."""

    def test_create_metrics(self) -> None:
        """GWASMetrics can be created."""
        metrics = GWASMetrics(
            value=1.02,
            metric="lambda",
            ancestry="EUR",
        )
        assert metrics.value == 1.02
        assert metrics.metric == "lambda"
        assert metrics.ancestry == "EUR"

    def test_default_ancestry(self) -> None:
        """GWASMetrics defaults to 'all' ancestry."""
        metrics = GWASMetrics(value=1.0, metric="lambda")
        assert metrics.ancestry == "all"


class TestPipelineOutput:
    """Tests for PipelineOutput class."""

    @pytest.fixture
    def sample_output(self) -> PipelineOutput:
        """Create sample pipeline output."""
        return PipelineOutput(
            input_samples=pd.DataFrame({
                "#FID": ["FAM1", "FAM2"],
                "IID": ["SAMP1", "SAMP2"],
                "SEX": [1, 2],
            }),
            ancestry_counts=pd.DataFrame({
                "label": ["EUR", "AFR"],
                "count": [100, 50],
            }),
            qc_metrics=[
                QCMetrics(
                    step="callrate_prune",
                    count=5,
                    metric="outlier_count",
                    ancestry="EUR",
                    level="sample",
                    passed=True,
                ),
                QCMetrics(
                    step="geno_prune",
                    count=100,
                    metric="variant_count",
                    ancestry="EUR",
                    level="variant",
                    passed=True,
                ),
            ],
            gwas_metrics=[
                GWASMetrics(value=1.02, metric="lambda", ancestry="EUR"),
            ],
            pass_fail={
                "pass_fail": {
                    "callrate": {"status": True, "input": "/in", "output": "/out"},
                }
            },
        )

    def test_to_dict_structure(self, sample_output: PipelineOutput) -> None:
        """to_dict returns expected structure."""
        result = sample_output.to_dict()

        assert "input_samples" in result
        assert "ancestry_counts" in result
        assert "QC" in result
        assert "GWAS" in result
        assert "pass_fail" in result

    def test_to_dict_qc_format(self, sample_output: PipelineOutput) -> None:
        """to_dict QC section has correct format."""
        result = sample_output.to_dict()
        qc = result["QC"]

        # Should be dict format from DataFrame
        assert "step" in qc
        # Key is "pruned_count" in the JSON, per the pre-refactor contract.
        assert "pruned_count" in qc
        assert "count" not in qc
        assert "metric" in qc
        assert "ancestry" in qc
        assert "level" in qc

    def test_to_dict_gwas_format(self, sample_output: PipelineOutput) -> None:
        """to_dict GWAS section has correct format."""
        result = sample_output.to_dict()
        gwas = result["GWAS"]

        assert "value" in gwas
        assert "metric" in gwas
        assert "ancestry" in gwas

    def test_to_json(self, sample_output: PipelineOutput) -> None:
        """to_json produces valid JSON."""
        json_str = sample_output.to_json()
        parsed = json.loads(json_str)

        assert "input_samples" in parsed
        assert "ancestry_counts" in parsed

    def test_to_json_with_indent(self, sample_output: PipelineOutput) -> None:
        """to_json with indent produces formatted JSON."""
        json_str = sample_output.to_json(indent=2)
        # Indented JSON has newlines
        assert "\n" in json_str

    def test_save(self, sample_output: PipelineOutput, tmp_path: Path) -> None:
        """save writes to file correctly."""
        output_path = tmp_path / "output.json"
        sample_output.save(output_path)

        assert output_path.exists()

        with open(output_path) as f:
            loaded = json.load(f)

        assert "input_samples" in loaded
        assert "ancestry_counts" in loaded

    def test_empty_output(self) -> None:
        """Empty output converts correctly."""
        output = PipelineOutput()
        result = output.to_dict()

        # Should not include empty sections
        assert "input_samples" not in result
        assert "ancestry_counts" not in result
        assert "QC" not in result
        assert "GWAS" not in result

    def test_with_pruned_samples(self) -> None:
        """Output with pruned samples includes them."""
        output = PipelineOutput(
            pruned_samples=pd.DataFrame({
                "FID": ["FAM1", "FAM2"],
                "IID": ["SAMP1", "SAMP2"],
                "step": ["callrate", "het"],
            }),
        )
        result = output.to_dict()

        assert "pruned_samples" in result
        pruned = result["pruned_samples"]
        assert "FID" in pruned
        assert "IID" in pruned
        assert "step" in pruned

    def test_with_related_samples(self) -> None:
        """Output with related samples includes them."""
        output = PipelineOutput(
            related_samples=pd.DataFrame({
                "ID1": ["SAMP1"],
                "ID2": ["SAMP2"],
                "KINSHIP": [0.25],
            }),
        )
        result = output.to_dict()

        assert "related_samples" in result

    def test_empty_dataframe_not_included(self) -> None:
        """Empty DataFrames are not included in output."""
        output = PipelineOutput(
            pruned_samples=pd.DataFrame(),  # Empty
            related_samples=pd.DataFrame(),  # Empty
        )
        result = output.to_dict()

        assert "pruned_samples" not in result
        assert "related_samples" not in result


class TestWriteResults:
    """Tests for write_results function."""

    def test_write_results(self, tmp_path: Path) -> None:
        """write_results creates JSON file."""
        output = PipelineOutput(
            qc_metrics=[
                QCMetrics(step="test", count=5, metric="count"),
            ],
        )
        out_path = tmp_path / "results"
        write_results(output, out_path)

        json_path = tmp_path / "results.json"
        assert json_path.exists()


class TestBuildMetricsDataframe:
    """Tests for build_metrics_dataframe function."""

    def test_empty_list(self) -> None:
        """Empty list returns empty DataFrame."""
        df = build_metrics_dataframe([])
        assert df.empty

    def test_single_metric(self) -> None:
        """Single metric creates correct DataFrame."""
        metrics = [
            QCMetrics(
                step="callrate_prune",
                count=10,
                metric="outlier_count",
                ancestry="EUR",
                level="sample",
                passed=True,
            ),
        ]
        df = build_metrics_dataframe(metrics)

        assert len(df) == 1
        assert df.iloc[0]["step"] == "callrate_prune"
        assert df.iloc[0]["pruned_count"] == 10
        assert df.iloc[0]["ancestry"] == "EUR"

    def test_multiple_metrics(self) -> None:
        """Multiple metrics create correct DataFrame."""
        metrics = [
            QCMetrics(step="callrate", count=5, metric="count"),
            QCMetrics(step="sex", count=3, metric="count"),
            QCMetrics(step="het", count=2, metric="count"),
        ]
        df = build_metrics_dataframe(metrics)

        assert len(df) == 3
        assert list(df["step"]) == ["callrate", "sex", "het"]


class TestBuildGwasDataframe:
    """Tests for build_gwas_dataframe function."""

    def test_empty_list(self) -> None:
        """Empty list returns empty DataFrame."""
        df = build_gwas_dataframe([])
        assert df.empty

    def test_single_metric(self) -> None:
        """Single metric creates correct DataFrame."""
        metrics = [
            GWASMetrics(value=1.02, metric="lambda", ancestry="EUR"),
        ]
        df = build_gwas_dataframe(metrics)

        assert len(df) == 1
        assert df.iloc[0]["value"] == 1.02
        assert df.iloc[0]["metric"] == "lambda"
        assert df.iloc[0]["ancestry"] == "EUR"

    def test_multiple_metrics(self) -> None:
        """Multiple metrics create correct DataFrame."""
        metrics = [
            GWASMetrics(value=1.02, metric="lambda", ancestry="EUR"),
            GWASMetrics(value=1.05, metric="lambda", ancestry="AFR"),
        ]
        df = build_gwas_dataframe(metrics)

        assert len(df) == 2
        assert list(df["ancestry"]) == ["EUR", "AFR"]


class TestQCMetricExtraction:
    """Regression cover for the JSON report schema.

    Two drifts the refactor introduced, both caught by comparing a real
    old-vs-new ancestry run: variant steps emitted a duplicate generic
    outlier_count row alongside their own metric (so summing pruned_count per
    step double-counted), and pruned_samples lost the ancestry label (so a
    pruned ID could not be traced back to its group).
    """

    @staticmethod
    def _extract(result_dict: dict, ancestry: str, tmp_path: Path):
        from genotools.cli.output import PipelineOutput

        class _Args:
            out_path = str(tmp_path / "out")

        out = PipelineOutput()
        metrics: list = []
        pruned: list = []
        out._extract_qc_metrics(
            result_dict=result_dict,
            ancestry=ancestry,
            sample_steps=["callrate", "sex"],
            variant_steps=["geno", "hwe"],
            metrics_list=metrics,
            pruned_dfs=pruned,
            related_dfs=[],
            args=_Args(),
        )
        return metrics, pruned

    def test_variant_step_emits_one_row_not_two(self, tmp_path: Path) -> None:
        """The generic outlier_count is dropped when the step has its own metric."""
        result = {
            "geno": {
                "pass": True,
                "step": "geno_prune",
                # What FilterResult.to_dict() hands over: the step's own count
                # plus a generic outlier_count for legacy in-memory consumers.
                "metrics": {"outlier_count": 20236, "geno_removed_count": 20236},
                "output": {},
            }
        }
        metrics, _ = self._extract(result, "EUR", tmp_path)

        assert [m.metric for m in metrics] == ["geno_removed_count"]
        assert metrics[0].count == 20236

    def test_variant_step_keeps_lone_outlier_count(self, tmp_path: Path) -> None:
        """A step reporting only outlier_count is still represented."""
        result = {
            "geno": {
                "pass": True,
                "step": "geno_prune",
                "metrics": {"outlier_count": 7},
                "output": {},
            }
        }
        metrics, _ = self._extract(result, "EUR", tmp_path)

        assert [m.metric for m in metrics] == ["outlier_count"]

    def test_sample_step_keeps_outlier_count(self, tmp_path: Path) -> None:
        """Sample steps report outlier_count and must not be filtered."""
        result = {
            "callrate": {
                "pass": True,
                "step": "callrate_prune",
                "metrics": {"outlier_count": 125},
                "output": {},
            }
        }
        metrics, _ = self._extract(result, "EUR", tmp_path)

        assert [m.metric for m in metrics] == ["outlier_count"]

    def _pruned_file(self, tmp_path: Path) -> str:
        path = tmp_path / "sex.outliers"
        path.write_text("#FID\tIID\nFAM1\tSAMP1\n")
        return str(path)

    def test_pruned_samples_labeled_in_ancestry_run(self, tmp_path: Path) -> None:
        result = {
            "sex": {
                "pass": True,
                "step": "sex_prune",
                "metrics": {"outlier_count": 1},
                "output": {"pruned_samples": self._pruned_file(tmp_path)},
            }
        }
        _, pruned = self._extract(result, "EUR", tmp_path)

        assert len(pruned) == 1
        assert list(pruned[0].columns) == ["#FID", "IID", "step", "label"]
        assert pruned[0]["label"].tolist() == ["EUR"]

    def test_pruned_samples_unlabeled_in_flat_run(self, tmp_path: Path) -> None:
        """A non-ancestry run has no group to label, and keeps the old shape."""
        result = {
            "sex": {
                "pass": True,
                "step": "sex_prune",
                "metrics": {"outlier_count": 1},
                "output": {"pruned_samples": self._pruned_file(tmp_path)},
            }
        }
        _, pruned = self._extract(result, "all", tmp_path)

        assert list(pruned[0].columns) == ["#FID", "IID", "step"]
