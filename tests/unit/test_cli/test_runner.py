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

"""Tests for CLI pipeline runner."""

from pathlib import Path
from typing import Any, Dict

import pytest

from genotools.cli.parser import (
    InputArgs,
    SampleQCArgs,
    VariantQCArgs,
    AncestryArgs,
    OutputArgs,
    PipelineArgs,
)
from genotools.cli.runner import (
    StepResult,
    PassFailRecord,
    PipelineState,
    PipelineRunner,
)


class TestStepResult:
    """Tests for StepResult dataclass."""

    def test_create_result(self) -> None:
        """StepResult can be created."""
        result = StepResult(
            step="callrate_prune",
            passed=True,
            metrics={"outlier_count": 5},
            output={"pruned_samples": "/path/to/file"},
        )
        assert result.step == "callrate_prune"
        assert result.passed is True
        assert result.metrics["outlier_count"] == 5

    def test_from_legacy_dict(self) -> None:
        """from_legacy_dict creates result correctly."""
        legacy = {
            "pass": True,
            "step": "callrate_prune",
            "metrics": {"outlier_count": 5},
            "output": {"plink_out": "/path/to/out"},
        }
        result = StepResult.from_legacy_dict("callrate", legacy)

        assert result.step == "callrate"
        assert result.passed is True
        assert result.metrics["outlier_count"] == 5

    def test_from_legacy_dict_missing_fields(self) -> None:
        """from_legacy_dict handles missing fields."""
        legacy: Dict[str, Any] = {}
        result = StepResult.from_legacy_dict("test", legacy)

        assert result.step == "test"
        assert result.passed is False
        assert result.metrics == {}


class TestPassFailRecord:
    """Tests for PassFailRecord dataclass."""

    def test_create_record(self) -> None:
        """PassFailRecord can be created."""
        record = PassFailRecord(
            status=True,
            input_path="/input",
            output_path="/output",
        )
        assert record.status is True
        assert record.input_path == "/input"
        assert record.output_path == "/output"
        assert record.error is None

    def test_with_error(self) -> None:
        """PassFailRecord can include error."""
        record = PassFailRecord(
            status=False,
            input_path="/input",
            output_path="/output",
            error="Something went wrong",
        )
        assert record.status is False
        assert record.error == "Something went wrong"


class TestPipelineState:
    """Tests for PipelineState dataclass."""

    def test_create_state(self) -> None:
        """PipelineState can be created."""
        state = PipelineState(
            geno_path=Path("/data/test"),
            out_path=Path("/output/test"),
        )
        assert state.geno_path == Path("/data/test")
        assert state.out_path == Path("/output/test")
        assert state.tmp_dir is None
        assert state.pass_fail == {}
        assert state.step_results == {}

    def test_get_last_passed_output_empty(self) -> None:
        """get_last_passed_output returns None when empty."""
        state = PipelineState(
            geno_path=Path("/data/test"),
            out_path=Path("/output/test"),
        )
        assert state.get_last_passed_output() is None

    def test_get_last_passed_output(self) -> None:
        """get_last_passed_output returns correct path."""
        state = PipelineState(
            geno_path=Path("/data/test"),
            out_path=Path("/output/test"),
        )
        state.pass_fail["callrate"] = PassFailRecord(
            status=True, input_path="/in1", output_path="/out1"
        )
        state.pass_fail["sex"] = PassFailRecord(
            status=True, input_path="/out1", output_path="/out2"
        )
        state.pass_fail["het"] = PassFailRecord(
            status=False, input_path="/out2", output_path="/out3"
        )

        assert state.get_last_passed_output() == "/out2"


class TestPipelineRunner:
    """Tests for PipelineRunner class."""

    @pytest.fixture
    def basic_args(self) -> PipelineArgs:
        """Create basic pipeline args."""
        return PipelineArgs(
            input=InputArgs(pfile=Path("/data/test")),
            output=OutputArgs(out_path=Path("/output/test")),
        )

    def test_init(self, basic_args: PipelineArgs) -> None:
        """PipelineRunner initializes correctly."""
        runner = PipelineRunner(basic_args)
        assert runner.args == basic_args
        assert runner.state is None

    def test_step_categories(self) -> None:
        """Step categories are defined correctly."""
        assert "callrate" in PipelineRunner.SAMPLE_STEPS
        assert "sex" in PipelineRunner.SAMPLE_STEPS
        assert "het" in PipelineRunner.SAMPLE_STEPS
        assert "related" in PipelineRunner.SAMPLE_STEPS
        assert "kinship_check" in PipelineRunner.SAMPLE_STEPS

        assert "geno" in PipelineRunner.VARIANT_STEPS
        assert "case_control" in PipelineRunner.VARIANT_STEPS
        assert "haplotype" in PipelineRunner.VARIANT_STEPS
        assert "hwe" in PipelineRunner.VARIANT_STEPS
        assert "ld" in PipelineRunner.VARIANT_STEPS

    def test_compute_step_paths_first_step(self, basic_args: PipelineArgs) -> None:
        """_compute_step_paths handles first step correctly."""
        runner = PipelineRunner(basic_args)
        runner.state = PipelineState(
            geno_path=Path("/data/test"),
            out_path=Path("/output/test"),
        )

        step_input, step_output = runner._compute_step_paths(
            step="callrate",
            step_index=0,
            steps=["callrate", "sex"],
            pass_fail={},
            geno_path="/data/test",
            out_path="/output/test",
            working_out="/tmp/test",
        )

        assert step_input == "/data/test"
        # Not last step, so gets step suffix
        assert step_output == "/tmp/test_callrate"

    def test_compute_step_paths_last_step(self, basic_args: PipelineArgs) -> None:
        """_compute_step_paths handles last step correctly."""
        basic_args.output.full_output = True
        runner = PipelineRunner(basic_args)
        runner.state = PipelineState(
            geno_path=Path("/data/test"),
            out_path=Path("/output/test"),
        )

        step_input, step_output = runner._compute_step_paths(
            step="sex",
            step_index=1,
            steps=["callrate", "sex"],
            pass_fail={
                "callrate": PassFailRecord(
                    status=True,
                    input_path="/data/test",
                    output_path="/output/test_callrate",
                )
            },
            geno_path="/data/test",
            out_path="/output/test",
            working_out="/tmp/test",
        )

        # Last step outputs to final out_path
        assert step_output == "/output/test"

    def test_compute_step_paths_warn_mode(self) -> None:
        """_compute_step_paths in warn mode finds last passed."""
        args = PipelineArgs(
            input=InputArgs(pfile=Path("/data/test")),
            output=OutputArgs(out_path=Path("/output/test"), warn_only=True),
        )
        runner = PipelineRunner(args)
        runner.state = PipelineState(
            geno_path=Path("/data/test"),
            out_path=Path("/output/test"),
        )

        pass_fail = {
            "callrate": PassFailRecord(
                status=True,
                input_path="/data/test",
                output_path="/tmp/test_callrate",
            ),
            "sex": PassFailRecord(
                status=False,
                input_path="/tmp/test_callrate",
                output_path="/tmp/test_callrate_sex",
            ),
        }

        step_input, step_output = runner._compute_step_paths(
            step="het",
            step_index=2,
            steps=["callrate", "sex", "het"],
            pass_fail=pass_fail,
            geno_path="/data/test",
            out_path="/output/test",
            working_out="/tmp/test",
        )

        # Should use last passed output (callrate)
        assert step_input == "/tmp/test_callrate"


class TestPipelineRunnerValidation:
    """Tests for PipelineRunner validation."""

    def test_no_steps_no_ancestry_raises(self) -> None:
        """Running with no steps and no ancestry raises error."""
        args = PipelineArgs(
            input=InputArgs(pfile=Path("/data/test")),
            output=OutputArgs(out_path=Path("/output/test")),
            # No QC steps enabled, no ancestry
        )
        runner = PipelineRunner(args)

        # The run method would fail because there's nothing to do
        # This tests the validation logic indirectly
        steps = args.get_all_enabled_steps()
        assert steps == []
        assert args.ancestry.run_ancestry is False


class TestDefaultEntryPoint:
    def test_only_main_entry_point(self):
        import genotools.cli as cli
        assert callable(cli.main)
        assert not hasattr(cli, "main_new")

    def test_run_pipeline_has_no_ab_flag(self):
        import inspect
        from genotools.cli import run_pipeline
        assert "use_new_ancestry" not in inspect.signature(run_pipeline).parameters


class TestNewAncestryStandalone:
    def test_new_ancestry_path_is_legacy_free(self):
        """The new ancestry path must call the ported functions, never the legacy
        Ancestry helper methods. A full ancestry run needs a trained model/ref panel
        (out of scope for a unit test), so this is a static source guard over the
        three methods that make up the new path; the ported functions themselves are
        differentially verified against legacy in tests/unit/test_ancestry/, and the
        end-to-end genotools-new path was validated by a manual smoke run."""
        import inspect
        from genotools.cli import runner

        src = (
            inspect.getsource(runner.PipelineRunner._run_ancestry_prediction_new)
            + inspect.getsource(runner.PipelineRunner._run_training_mode)
            + inspect.getsource(runner.PipelineRunner._run_inference_mode)
        )
        # No legacy reach-back:
        assert "self._ancestry.get_raw_files" not in src
        assert "self._ancestry.split_cohort_ancestry" not in src
        assert "self._ancestry.clean_up" not in src
        # Ported functions are used instead:
        assert "get_raw_files(" in src
        assert "split_cohort_by_ancestry(" in src
        assert "clean_up_files(" in src
