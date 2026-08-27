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

"""Tests for QC result dataclasses."""

from pathlib import Path
from unittest.mock import MagicMock

import pytest

from genotools.qc.results import FilterResult


class TestFilterResult:
    """Tests for FilterResult dataclass."""

    @pytest.fixture
    def mock_genotype_data(self) -> MagicMock:
        """Create a mock GenotypeData object."""
        mock = MagicMock()
        mock.path = Path("/output/test")
        mock.sample_count = 100
        mock.variant_count = 1000
        return mock

    @pytest.fixture
    def sample_filter_result(self, mock_genotype_data: MagicMock) -> FilterResult:
        """Create a sample FilterResult for testing."""
        return FilterResult(
            output=mock_genotype_data,
            samples_removed=10,
            variants_removed=0,
            metrics={"threshold": 0.05},
            log="PLINK output here",
            pruned_samples_file=Path("/output/test.outliers"),
            step_name="callrate_prune",
        )

    def test_to_dict_format(self, sample_filter_result: FilterResult) -> None:
        """to_dict returns the correct legacy format."""
        result_dict = sample_filter_result.to_dict()

        assert result_dict["pass"] is True
        assert result_dict["step"] == "callrate_prune"
        assert "outlier_count" in result_dict["metrics"]
        assert result_dict["metrics"]["outlier_count"] == 10  # samples_removed
        assert result_dict["metrics"]["threshold"] == 0.05
        assert result_dict["output"]["pruned_samples"] == "/output/test.outliers"
        assert result_dict["output"]["plink_out"] == "/output/test"

    def test_to_dict_no_pruned_samples(self, mock_genotype_data: MagicMock) -> None:
        """to_dict handles None pruned_samples_file."""
        result = FilterResult(
            output=mock_genotype_data,
            samples_removed=0,
            variants_removed=0,
            metrics={},
            log="",
            pruned_samples_file=None,
            step_name="test_step",
        )
        result_dict = result.to_dict()

        assert result_dict["output"]["pruned_samples"] is None

    def test_to_dict_variant_removal(self, mock_genotype_data: MagicMock) -> None:
        """to_dict correctly sums sample and variant removals."""
        result = FilterResult(
            output=mock_genotype_data,
            samples_removed=5,
            variants_removed=100,
            metrics={},
            log="",
            pruned_samples_file=None,
            step_name="combined_prune",
        )
        result_dict = result.to_dict()

        assert result_dict["metrics"]["outlier_count"] == 105

    def test_frozen(self, sample_filter_result: FilterResult) -> None:
        """FilterResult is immutable."""
        with pytest.raises(AttributeError):
            sample_filter_result.samples_removed = 0  # type: ignore[misc]


class TestLegacyCompatibility:
    """Tests ensuring backward compatibility with old return format."""

    def test_filter_result_matches_old_format(self) -> None:
        """FilterResult.to_dict() matches the old SampleQC/VariantQC return format."""
        mock_output = MagicMock()
        mock_output.path = Path("/path/to/output")

        result = FilterResult(
            output=mock_output,
            samples_removed=5,
            variants_removed=0,
            metrics={"extra_metric": "value"},
            log="",
            pruned_samples_file=Path("/path/to/output.outliers"),
            step_name="callrate_prune",
        )

        result_dict = result.to_dict()

        # Check all required keys exist
        assert "pass" in result_dict
        assert "step" in result_dict
        assert "metrics" in result_dict
        assert "output" in result_dict

        # Check nested structure
        assert "outlier_count" in result_dict["metrics"]
        assert "pruned_samples" in result_dict["output"]
        assert "plink_out" in result_dict["output"]

        # Check types
        assert isinstance(result_dict["pass"], bool)
        assert isinstance(result_dict["step"], str)
        assert isinstance(result_dict["metrics"], dict)
        assert isinstance(result_dict["output"], dict)

        # Check custom metrics are preserved
        assert result_dict["metrics"]["extra_metric"] == "value"


class TestStepReportRegistry:
    """STEP_REPORT is the report identity of each step - the names used when a
    step produces no data (failed or skipped). Nothing else links it to the step
    modules, so these guards catch drift.
    """

    # Every step key the registry claims, paired with the function that owns it.
    STEP_FUNCTIONS = {
        "callrate": "filter_callrate",
        "sex": "filter_sex",
        "het": "filter_heterozygosity",
        "related": "filter_relatedness",
        "geno": "filter_variant_missingness",
        "case_control": "filter_case_control",
        "haplotype": "filter_haplotype",
        "hwe": "filter_hwe",
        "ld": "prune_ld",
    }

    def test_registry_covers_every_step_function(self) -> None:
        """A new QC step must register, or its failures vanish from the report."""
        from genotools.qc.results import STEP_REPORT

        assert set(STEP_REPORT) == set(self.STEP_FUNCTIONS)

    @pytest.mark.parametrize("step", sorted(STEP_FUNCTIONS))
    def test_registered_step_name_matches_the_module(self, step: str) -> None:
        """The registry's step name must match the one the module emits."""
        import inspect
        import re

        from genotools.qc import steps as step_module
        from genotools.qc.results import STEP_REPORT

        func = getattr(step_module, self.STEP_FUNCTIONS[step])
        source = inspect.getsource(func)
        found = re.findall(r'step_name = "([^"]+)"', source)
        assert found, f"no step_name literal in {self.STEP_FUNCTIONS[step]}"

        expected_name, _ = STEP_REPORT[step]
        assert expected_name in found, (
            f"STEP_REPORT['{step}'] reports as {expected_name!r} but "
            f"{self.STEP_FUNCTIONS[step]} emits {found!r}"
        )

    @pytest.mark.parametrize("step", sorted(STEP_FUNCTIONS))
    def test_registered_metrics_appear_in_the_module(self, step: str) -> None:
        """Every registered metric key must be one the step actually writes."""
        import inspect

        from genotools.qc import steps as step_module
        from genotools.qc.results import STEP_REPORT

        source = inspect.getsource(getattr(step_module, self.STEP_FUNCTIONS[step]))
        _, metric_keys = STEP_REPORT[step]
        for key in metric_keys:
            if key == "outlier_count":
                # Added by FilterResult.to_dict, not by the step body.
                continue
            assert f'"{key}"' in source, (
                f"STEP_REPORT['{step}'] claims metric {key!r}, "
                f"which {self.STEP_FUNCTIONS[step]} never writes"
            )


class TestUnrunResult:
    """A step that produced no data still has to appear in the report."""

    def test_failed_step_reports_zeroed_metrics(self) -> None:
        from genotools.qc.results import unrun_result

        result = unrun_result("hwe", "fail", "plink2 exited 13", "/tmp/out_hwe")
        assert result is not None
        assert result["pass"] is False
        assert result["outcome"] == "fail"
        assert result["reason"] == "plink2 exited 13"
        assert result["step"] == "hwe_prune"
        assert result["metrics"] == {"hwe_removed_count": 0}

    def test_skipped_step_is_distinguishable_from_a_failure(self) -> None:
        from genotools.qc.results import unrun_result

        skipped = unrun_result("het", "skipped", "too few samples", "/tmp/out")
        failed = unrun_result("het", "fail", "plink2 exited 13", "/tmp/out")
        assert skipped is not None and failed is not None
        # Both are pass=False; only outcome separates them.
        assert skipped["pass"] is failed["pass"] is False
        assert skipped["outcome"] == "skipped"
        assert failed["outcome"] == "fail"

    def test_relatedness_reports_both_counts(self) -> None:
        from genotools.qc.results import unrun_result

        result = unrun_result("related", "skipped", "why", "/tmp/out")
        assert result is not None
        assert result["metrics"] == {"related_count": 0, "duplicated_count": 0}

    def test_unregistered_step_reports_nothing(self) -> None:
        """assoc/kinship_check carry no QC metrics and are absent either way."""
        from genotools.qc.results import unrun_result

        assert unrun_result("assoc", "fail", "why", "/tmp/out") is None

    def test_passing_result_declares_its_outcome(self) -> None:
        """to_dict must set outcome too, so every result dict has the field."""
        from pathlib import Path

        from genotools.core.genotypes import GenotypeData
        from genotools.qc.results import FilterResult

        result = FilterResult(
            output=GenotypeData(
                path=Path("/tmp/out"), format="pfile", sample_count=10, variant_count=5
            ),
            samples_removed=1,
            variants_removed=0,
            step_name="callrate_prune",
        )
        assert result.to_dict()["outcome"] == "pass"
        assert result.to_dict()["reason"] is None
