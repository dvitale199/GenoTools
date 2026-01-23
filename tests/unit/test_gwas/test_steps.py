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

"""Tests for GWAS step functions."""

import numpy as np
import pytest

from genotools.gwas.steps.association import calculate_inflation


class TestCalculateInflation:
    """Tests for calculate_inflation function."""

    def test_uniform_pvalues(self) -> None:
        """Uniform p-values should give lambda ~ 1.0."""
        # Generate uniform p-values (no inflation)
        np.random.seed(42)
        pvalues = np.random.uniform(0, 1, 10000)

        metrics = calculate_inflation(pvalues)

        # Lambda should be close to 1.0 for uniform p-values
        assert 0.95 <= metrics.lambda_gc <= 1.05

    def test_inflated_pvalues(self) -> None:
        """Inflated p-values should give lambda > 1.0."""
        # Generate inflated p-values (more small values than expected)
        np.random.seed(42)
        # Use chi-squared with higher nc to simulate inflation
        from scipy.stats import ncx2

        chi2_stats = ncx2.rvs(1, nc=2, size=10000)  # Inflated
        pvalues = 1 - ncx2.cdf(chi2_stats, 1, nc=0)

        metrics = calculate_inflation(pvalues)

        # Lambda should be > 1.0 for inflated p-values
        assert metrics.lambda_gc > 1.5

    def test_with_case_control_normalization(self) -> None:
        """Lambda1000 should be calculated when cases/controls provided."""
        np.random.seed(42)
        pvalues = np.random.uniform(0, 1, 10000)

        metrics = calculate_inflation(pvalues, n_cases=500, n_controls=500)

        assert metrics.lambda_gc is not None
        assert metrics.lambda_1000 is not None
        assert metrics.n_cases == 500
        assert metrics.n_controls == 500

    def test_without_case_control(self) -> None:
        """Lambda1000 should be None when cases/controls not provided."""
        np.random.seed(42)
        pvalues = np.random.uniform(0, 1, 10000)

        metrics = calculate_inflation(pvalues)

        assert metrics.lambda_gc is not None
        assert metrics.lambda_1000 is None
        assert metrics.n_cases is None
        assert metrics.n_controls is None

    def test_empty_pvalues(self) -> None:
        """Empty p-values should return NaN."""
        pvalues = np.array([])

        metrics = calculate_inflation(pvalues)

        assert np.isnan(metrics.lambda_gc)

    def test_all_nan_pvalues(self) -> None:
        """All NaN p-values should return NaN."""
        pvalues = np.array([np.nan, np.nan, np.nan])

        metrics = calculate_inflation(pvalues)

        assert np.isnan(metrics.lambda_gc)

    def test_some_nan_pvalues(self) -> None:
        """NaN p-values should be filtered out."""
        np.random.seed(42)
        pvalues = np.random.uniform(0, 1, 1000)
        # Add some NaN values
        pvalues[0:100] = np.nan

        metrics = calculate_inflation(pvalues)

        # Should still calculate lambda (not NaN)
        assert not np.isnan(metrics.lambda_gc)

    def test_lambda1000_normalization_direction(self) -> None:
        """Lambda1000 should be closer to 1 for larger sample sizes."""
        np.random.seed(42)
        # Generate slightly inflated p-values
        from scipy.stats import ncx2

        chi2_stats = ncx2.rvs(1, nc=0.5, size=10000)
        pvalues = 1 - ncx2.cdf(chi2_stats, 1, nc=0)

        # Small sample
        metrics_small = calculate_inflation(pvalues, n_cases=100, n_controls=100)

        # Large sample
        metrics_large = calculate_inflation(pvalues, n_cases=10000, n_controls=10000)

        # Lambda1000 should be further from 1 for small samples
        # (normalization inflates the deviation)
        if metrics_small.lambda_1000 is not None and metrics_large.lambda_1000 is not None:
            small_dev = abs(metrics_small.lambda_1000 - 1)
            large_dev = abs(metrics_large.lambda_1000 - 1)
            assert small_dev >= large_dev

    def test_pvalues_edge_cases(self) -> None:
        """Test edge case p-values (0 and 1)."""
        # P-value of exactly 0 would cause issues, but near 0 should work
        pvalues = np.array([1e-10, 0.5, 0.999999])

        metrics = calculate_inflation(pvalues)

        # Should not raise error and produce a valid result
        assert not np.isnan(metrics.lambda_gc)


class TestInflationWithRealDistributions:
    """Tests using realistic GWAS p-value distributions."""

    def test_well_calibrated_gwas(self) -> None:
        """Simulated well-calibrated GWAS (no stratification)."""
        np.random.seed(42)

        # Simulate null distribution with a few true signals
        n_variants = 100000
        pvalues = np.random.uniform(0, 1, n_variants)

        # Add some true signals (small p-values)
        n_signals = 100
        pvalues[0:n_signals] = np.random.uniform(1e-10, 1e-5, n_signals)

        metrics = calculate_inflation(pvalues)

        # Should still be close to 1.0 (signals don't affect median much)
        assert 0.95 <= metrics.lambda_gc <= 1.10

    def test_stratified_gwas(self) -> None:
        """Simulated GWAS with population stratification."""
        np.random.seed(42)

        # Simulate inflated distribution
        from scipy.stats import ncx2

        n_variants = 100000
        # Use chi-squared with small non-centrality to simulate mild inflation
        chi2_stats = ncx2.rvs(1, nc=0.3, size=n_variants)
        pvalues = 1 - ncx2.cdf(chi2_stats, 1, nc=0)

        metrics = calculate_inflation(pvalues, n_cases=1000, n_controls=1000)

        # Lambda should be notably > 1.0
        assert metrics.lambda_gc > 1.1
        assert metrics.lambda_1000 is not None
