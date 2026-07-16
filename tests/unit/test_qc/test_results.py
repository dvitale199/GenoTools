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
