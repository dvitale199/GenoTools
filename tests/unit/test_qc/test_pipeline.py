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

"""Tests for QC pipeline orchestration."""

from pathlib import Path
from typing import Any
from unittest.mock import MagicMock, patch

import pytest

from genotools.core.exceptions import QCError
from genotools.qc.config import CallrateConfig, GenoConfig, SexConfig
from genotools.qc.pipeline import QCPipeline
from genotools.qc.results import FilterResult


class TestQCPipeline:
    """Tests for QCPipeline class."""

    @pytest.fixture
    def mock_genotype_data(self) -> MagicMock:
        """Create a mock GenotypeData object."""
        mock = MagicMock()
        mock.path = Path("/input/test")
        mock.sample_count = 500
        mock.variant_count = 10000
        mock.pgen = MagicMock()
        mock.pgen.exists.return_value = True
        return mock

    @pytest.fixture
    def mock_step_fn(self, mock_genotype_data: MagicMock) -> MagicMock:
        """Create a mock step function."""

        def create_output(sample_count: int, variant_count: int) -> MagicMock:
            output = MagicMock()
            output.path = Path("/output/step")
            output.sample_count = sample_count
            output.variant_count = variant_count
            return output

        def step_fn(
            data: Any, config: Any, out_path: Path
        ) -> FilterResult:
            output = create_output(
                data.sample_count - 10,  # Remove 10 samples
                data.variant_count,
            )
            return FilterResult(
                output=output,
                samples_removed=10,
                variants_removed=0,
                metrics={},
                log="",
                pruned_samples_file=None,
                step_name="test_step",
            )

        return MagicMock(side_effect=step_fn)

    def test_empty_pipeline(self, mock_genotype_data: MagicMock, tmp_path: Path) -> None:
        """Empty pipeline returns input data unchanged."""
        pipeline = QCPipeline(steps=[])
        result = pipeline.run(mock_genotype_data, tmp_path)

        assert result.input == mock_genotype_data
        assert result.output == mock_genotype_data
        assert len(result.step_results) == 0

    def test_single_step(
        self, mock_genotype_data: MagicMock, mock_step_fn: MagicMock, tmp_path: Path
    ) -> None:
        """Single step pipeline executes correctly."""
        pipeline = QCPipeline(
            steps=[("test_step", mock_step_fn, CallrateConfig())]
        )
        result = pipeline.run(mock_genotype_data, tmp_path)

        assert mock_step_fn.called
        assert len(result.step_results) == 1
        assert result.step_results[0][0] == "test_step"

    def test_multiple_steps(
        self, mock_genotype_data: MagicMock, tmp_path: Path
    ) -> None:
        """Multiple steps are chained correctly."""
        step_calls = []

        def make_step(name: str) -> MagicMock:
            def step_fn(data: Any, config: Any, out_path: Path) -> FilterResult:
                step_calls.append(name)
                output = MagicMock()
                output.path = out_path
                output.sample_count = data.sample_count - 5
                output.variant_count = data.variant_count
                return FilterResult(
                    output=output,
                    samples_removed=5,
                    variants_removed=0,
                    metrics={},
                    log="",
                    pruned_samples_file=None,
                    step_name=name,
                )

            return MagicMock(side_effect=step_fn)

        step1 = make_step("step1")
        step2 = make_step("step2")
        step3 = make_step("step3")

        pipeline = QCPipeline(
            steps=[
                ("step1", step1, CallrateConfig()),
                ("step2", step2, SexConfig()),
                ("step3", step3, GenoConfig()),
            ]
        )
        result = pipeline.run(mock_genotype_data, tmp_path)

        assert step_calls == ["step1", "step2", "step3"]
        assert len(result.step_results) == 3
        assert result.total_samples_removed == 15  # 5 + 5 + 5

    def test_warn_only_continues_on_error(
        self, mock_genotype_data: MagicMock, tmp_path: Path
    ) -> None:
        """warn_only=True continues pipeline after step failure."""
        step_calls = []

        def failing_step(data: Any, config: Any, out_path: Path) -> FilterResult:
            step_calls.append("failing")
            raise QCError("Step failed")

        def passing_step(data: Any, config: Any, out_path: Path) -> FilterResult:
            step_calls.append("passing")
            output = MagicMock()
            output.path = out_path
            output.sample_count = data.sample_count
            output.variant_count = data.variant_count
            return FilterResult(
                output=output,
                samples_removed=0,
                variants_removed=0,
                metrics={},
                log="",
                pruned_samples_file=None,
                step_name="passing",
            )

        pipeline = QCPipeline(
            steps=[
                ("failing", failing_step, CallrateConfig()),
                ("passing", passing_step, SexConfig()),
            ],
            warn_only=True,
        )
        result = pipeline.run(mock_genotype_data, tmp_path)

        assert step_calls == ["failing", "passing"]
        assert len(result.step_results) == 1  # Only passing step has result
        assert result.step_results[0][0] == "passing"

    def test_warn_only_false_raises_on_error(
        self, mock_genotype_data: MagicMock, tmp_path: Path
    ) -> None:
        """warn_only=False raises exception on step failure."""

        def failing_step(data: Any, config: Any, out_path: Path) -> FilterResult:
            raise QCError("Step failed")

        pipeline = QCPipeline(
            steps=[("failing", failing_step, CallrateConfig())],
            warn_only=False,
        )

        with pytest.raises(QCError, match="Step failed"):
            pipeline.run(mock_genotype_data, tmp_path)

    def test_raises_when_all_samples_removed(
        self, mock_genotype_data: MagicMock, tmp_path: Path
    ) -> None:
        """Raises QCError when all samples are removed."""

        def remove_all_step(data: Any, config: Any, out_path: Path) -> FilterResult:
            output = MagicMock()
            output.path = out_path
            output.sample_count = 0  # All samples removed
            output.variant_count = data.variant_count
            return FilterResult(
                output=output,
                samples_removed=data.sample_count,
                variants_removed=0,
                metrics={},
                log="",
                pruned_samples_file=None,
                step_name="remove_all",
            )

        pipeline = QCPipeline(
            steps=[("remove_all", remove_all_step, CallrateConfig())]
        )

        with pytest.raises(QCError, match="All samples removed"):
            pipeline.run(mock_genotype_data, tmp_path)

    def test_raises_when_all_variants_removed(
        self, mock_genotype_data: MagicMock, tmp_path: Path
    ) -> None:
        """Raises QCError when all variants are removed."""

        def remove_all_step(data: Any, config: Any, out_path: Path) -> FilterResult:
            output = MagicMock()
            output.path = out_path
            output.sample_count = data.sample_count
            output.variant_count = 0  # All variants removed
            return FilterResult(
                output=output,
                samples_removed=0,
                variants_removed=data.variant_count,
                metrics={},
                log="",
                pruned_samples_file=None,
                step_name="remove_all",
            )

        pipeline = QCPipeline(
            steps=[("remove_all", remove_all_step, GenoConfig())]
        )

        with pytest.raises(QCError, match="All variants removed"):
            pipeline.run(mock_genotype_data, tmp_path)

    def test_add_step_returns_new_pipeline(self) -> None:
        """add_step returns a new pipeline, doesn't modify original."""
        original = QCPipeline(steps=[])

        def dummy_step(data: Any, config: Any, out_path: Path) -> FilterResult:
            raise NotImplementedError

        new_pipeline = original.add_step("new_step", dummy_step, CallrateConfig())

        assert len(original.steps) == 0
        assert len(new_pipeline.steps) == 1
        assert new_pipeline is not original


class TestQCPipelineFactories:
    """Tests for QCPipeline factory methods."""

    def test_sample_qc_all_steps(self) -> None:
        """sample_qc creates pipeline with all steps when all configs provided."""
        pipeline = QCPipeline.sample_qc(
            callrate_config=CallrateConfig(),
            sex_config=SexConfig(),
            het_config=MagicMock(),
            related_config=MagicMock(),
        )

        assert len(pipeline.steps) == 4
        step_names = [name for name, _, _ in pipeline.steps]
        assert step_names == ["callrate", "sex", "het", "related"]

    def test_sample_qc_some_steps(self) -> None:
        """sample_qc creates pipeline with only specified steps."""
        pipeline = QCPipeline.sample_qc(
            callrate_config=CallrateConfig(),
            sex_config=None,
            het_config=MagicMock(),
            related_config=None,
        )

        assert len(pipeline.steps) == 2
        step_names = [name for name, _, _ in pipeline.steps]
        assert step_names == ["callrate", "het"]

    def test_sample_qc_warn_only(self) -> None:
        """sample_qc respects warn_only parameter."""
        pipeline = QCPipeline.sample_qc(
            callrate_config=CallrateConfig(),
            warn_only=True,
        )

        assert pipeline.warn_only is True

    def test_variant_qc_all_steps(self) -> None:
        """variant_qc creates pipeline with all steps when all configs provided."""
        pipeline = QCPipeline.variant_qc(
            geno_config=GenoConfig(),
            case_control_config=MagicMock(),
            haplotype_config=MagicMock(),
            hwe_config=MagicMock(),
            ld_config=MagicMock(),
        )

        assert len(pipeline.steps) == 5
        step_names = [name for name, _, _ in pipeline.steps]
        assert step_names == ["geno", "case_control", "haplotype", "hwe", "ld"]

    def test_variant_qc_some_steps(self) -> None:
        """variant_qc creates pipeline with only specified steps."""
        pipeline = QCPipeline.variant_qc(
            geno_config=GenoConfig(),
            case_control_config=None,
            haplotype_config=None,
            hwe_config=MagicMock(),
            ld_config=None,
        )

        assert len(pipeline.steps) == 2
        step_names = [name for name, _, _ in pipeline.steps]
        assert step_names == ["geno", "hwe"]
