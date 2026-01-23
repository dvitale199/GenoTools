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

"""Tests for GWAS result classes."""

from pathlib import Path

import pytest

from genotools.gwas.results import (
    PruningResult,
    PCAResult,
    InflationMetrics,
    GWASResult,
    AssocResult,
)


class TestPruningResult:
    """Tests for PruningResult."""

    def test_basic_creation(self) -> None:
        """Basic creation with required fields."""
        result = PruningResult(
            output_path=Path("/path/to/output"),
            variants_before=100000,
            variants_after=50000,
        )
        assert result.output_path == Path("/path/to/output")
        assert result.variants_before == 100000
        assert result.variants_after == 50000
        assert result.success is True

    def test_variants_removed_property(self) -> None:
        """variants_removed property calculates correctly."""
        result = PruningResult(
            output_path=Path("/path/to/output"),
            variants_before=100000,
            variants_after=50000,
        )
        assert result.variants_removed == 50000

    def test_to_dict(self) -> None:
        """to_dict returns legacy format."""
        result = PruningResult(
            output_path=Path("/path/to/output"),
            variants_before=100000,
            variants_after=50000,
            success=True,
        )
        d = result.to_dict()
        assert d["pass"] is True
        assert d["step"] == "plink_pruning"
        assert d["output"] == "/path/to/output"

    def test_failed_result(self) -> None:
        """Failed result has success=False."""
        result = PruningResult(
            output_path=Path("/path/to/output"),
            variants_before=100000,
            variants_after=0,
            success=False,
            log="Error message",
        )
        assert result.success is False
        assert result.to_dict()["pass"] is False


class TestPCAResult:
    """Tests for PCAResult."""

    def test_successful_result(self) -> None:
        """Successful PCA result."""
        result = PCAResult(
            eigenvec_path=Path("/path/to/output.eigenvec"),
            eigenval_path=Path("/path/to/output.eigenval"),
            n_pcs=10,
            n_samples=1000,
            success=True,
        )
        assert result.eigenvec_path == Path("/path/to/output.eigenvec")
        assert result.eigenval_path == Path("/path/to/output.eigenval")
        assert result.n_pcs == 10
        assert result.n_samples == 1000
        assert result.success is True

    def test_to_dict(self) -> None:
        """to_dict returns legacy format."""
        result = PCAResult(
            eigenvec_path=Path("/path/to/output.eigenvec"),
            eigenval_path=Path("/path/to/output.eigenval"),
            n_pcs=10,
            success=True,
        )
        d = result.to_dict()
        assert d["pass"] is True
        assert d["step"] == "plink_pca"
        assert d["output"]["eigenvec"] == "/path/to/output.eigenvec"
        assert d["output"]["eigenval"] == "/path/to/output.eigenval"

    def test_failed_result(self) -> None:
        """Failed PCA result has None paths."""
        result = PCAResult(
            eigenvec_path=None,
            eigenval_path=None,
            n_pcs=10,
            success=False,
        )
        d = result.to_dict()
        assert d["pass"] is False
        assert d["output"]["eigenvec"] is None
        assert d["output"]["eigenval"] is None


class TestInflationMetrics:
    """Tests for InflationMetrics."""

    def test_basic_metrics(self) -> None:
        """Basic metrics creation."""
        metrics = InflationMetrics(
            lambda_gc=1.05,
            lambda_1000=1.02,
            n_cases=500,
            n_controls=1000,
        )
        assert metrics.lambda_gc == 1.05
        assert metrics.lambda_1000 == 1.02
        assert metrics.n_cases == 500
        assert metrics.n_controls == 1000

    def test_with_maf_metrics(self) -> None:
        """Metrics with MAF-filtered values."""
        metrics = InflationMetrics(
            lambda_gc=1.05,
            lambda_1000=1.02,
            lambda_gc_maf=1.04,
            lambda_1000_maf=1.01,
            n_cases=500,
            n_controls=1000,
        )
        assert metrics.lambda_gc_maf == 1.04
        assert metrics.lambda_1000_maf == 1.01

    def test_to_dict_basic(self) -> None:
        """to_dict returns expected keys for basic metrics."""
        metrics = InflationMetrics(
            lambda_gc=1.05,
            lambda_1000=1.02,
            n_cases=500,
            n_controls=1000,
        )
        d = metrics.to_dict()
        assert d["lambda"] == 1.05
        assert d["lambda1000"] == 1.02
        assert d["cases"] == 500
        assert d["controls"] == 1000
        assert "lambda_maf" not in d
        assert "lambda1000_maf" not in d

    def test_to_dict_with_maf(self) -> None:
        """to_dict includes MAF metrics when present."""
        metrics = InflationMetrics(
            lambda_gc=1.05,
            lambda_1000=1.02,
            lambda_gc_maf=1.04,
            lambda_1000_maf=1.01,
            n_cases=500,
            n_controls=1000,
        )
        d = metrics.to_dict()
        assert d["lambda_maf"] == 1.04
        assert d["lambda1000_maf"] == 1.01

    def test_quantitative_trait_metrics(self) -> None:
        """Metrics for quantitative traits (no case/control)."""
        metrics = InflationMetrics(
            lambda_gc=1.03,
        )
        assert metrics.lambda_gc == 1.03
        assert metrics.lambda_1000 is None
        assert metrics.n_cases is None
        assert metrics.n_controls is None


class TestGWASResult:
    """Tests for GWASResult."""

    def test_successful_logistic_result(self) -> None:
        """Successful logistic regression result."""
        inflation = InflationMetrics(
            lambda_gc=1.05,
            lambda_1000=1.02,
            n_cases=500,
            n_controls=1000,
        )
        result = GWASResult(
            output_path=Path("/path/to/output.PHENO1.glm.logistic.hybrid"),
            phenotype_type="binary",
            inflation=inflation,
            n_variants_tested=500000,
            success=True,
        )
        assert result.output_path == Path("/path/to/output.PHENO1.glm.logistic.hybrid")
        assert result.phenotype_type == "binary"
        assert result.inflation.lambda_gc == 1.05
        assert result.n_variants_tested == 500000
        assert result.success is True

    def test_successful_linear_result(self) -> None:
        """Successful linear regression result."""
        inflation = InflationMetrics(lambda_gc=1.03)
        result = GWASResult(
            output_path=Path("/path/to/output.PHENO1.glm.linear"),
            phenotype_type="quantitative",
            inflation=inflation,
            n_variants_tested=500000,
            success=True,
        )
        assert result.phenotype_type == "quantitative"
        assert result.inflation.lambda_1000 is None

    def test_to_dict(self) -> None:
        """to_dict returns legacy format."""
        inflation = InflationMetrics(
            lambda_gc=1.05,
            lambda_1000=1.02,
            n_cases=500,
            n_controls=1000,
        )
        result = GWASResult(
            output_path=Path("/path/to/output.PHENO1.glm.logistic.hybrid"),
            phenotype_type="binary",
            inflation=inflation,
            success=True,
        )
        d = result.to_dict()
        assert d["pass"] is True
        assert d["step"] == "GWAS"
        assert d["metrics"]["lambda"] == 1.05
        assert d["metrics"]["lambda1000"] == 1.02
        assert d["metrics"]["cases"] == 500
        assert d["metrics"]["controls"] == 1000
        assert d["output"]["gwas_output"] == "/path/to/output.PHENO1.glm.logistic.hybrid"

    def test_failed_result(self) -> None:
        """Failed GWAS result."""
        result = GWASResult(
            output_path=None,
            phenotype_type="unknown",
            inflation=None,
            success=False,
        )
        d = result.to_dict()
        assert d["pass"] is False
        assert d["metrics"]["lambda"] is None
        assert d["output"]["gwas_output"] is None


class TestAssocResult:
    """Tests for AssocResult."""

    def test_pca_only(self) -> None:
        """Result with only PCA."""
        pca = PCAResult(
            eigenvec_path=Path("/path/to/output.eigenvec"),
            eigenval_path=Path("/path/to/output.eigenval"),
            n_pcs=10,
            success=True,
        )
        result = AssocResult(pca_result=pca, success=True)
        assert result.pca_result is not None
        assert result.gwas_result is None
        assert result.success is True

    def test_gwas_only(self) -> None:
        """Result with only GWAS."""
        gwas = GWASResult(
            output_path=Path("/path/to/output.PHENO1.glm.logistic.hybrid"),
            phenotype_type="binary",
            inflation=InflationMetrics(lambda_gc=1.05),
            success=True,
        )
        result = AssocResult(gwas_result=gwas, success=True)
        assert result.pca_result is None
        assert result.gwas_result is not None
        assert result.success is True

    def test_both_pca_and_gwas(self) -> None:
        """Result with both PCA and GWAS."""
        pca = PCAResult(
            eigenvec_path=Path("/path/to/output.eigenvec"),
            eigenval_path=Path("/path/to/output.eigenval"),
            n_pcs=10,
            success=True,
        )
        gwas = GWASResult(
            output_path=Path("/path/to/output.PHENO1.glm.logistic.hybrid"),
            phenotype_type="binary",
            inflation=InflationMetrics(lambda_gc=1.05),
            success=True,
        )
        result = AssocResult(pca_result=pca, gwas_result=gwas, success=True)
        assert result.pca_result is not None
        assert result.gwas_result is not None

    def test_neither_raises_error(self) -> None:
        """Must have at least one result."""
        with pytest.raises(ValueError, match="At least one"):
            AssocResult(pca_result=None, gwas_result=None)

    def test_eigenvec_path_property(self) -> None:
        """eigenvec_path property returns correct path."""
        pca = PCAResult(
            eigenvec_path=Path("/path/to/output.eigenvec"),
            eigenval_path=Path("/path/to/output.eigenval"),
            n_pcs=10,
            success=True,
        )
        result = AssocResult(pca_result=pca)
        assert result.eigenvec_path == Path("/path/to/output.eigenvec")

    def test_eigenvec_path_property_no_pca(self) -> None:
        """eigenvec_path property returns None when no PCA."""
        gwas = GWASResult(
            output_path=Path("/path/to/output.PHENO1.glm.logistic.hybrid"),
            phenotype_type="binary",
            inflation=InflationMetrics(lambda_gc=1.05),
            success=True,
        )
        result = AssocResult(gwas_result=gwas)
        assert result.eigenvec_path is None

    def test_gwas_output_path_property(self) -> None:
        """gwas_output_path property returns correct path."""
        gwas = GWASResult(
            output_path=Path("/path/to/output.PHENO1.glm.logistic.hybrid"),
            phenotype_type="binary",
            inflation=InflationMetrics(lambda_gc=1.05),
            success=True,
        )
        result = AssocResult(gwas_result=gwas)
        assert result.gwas_output_path == Path("/path/to/output.PHENO1.glm.logistic.hybrid")

    def test_inflation_metrics_property(self) -> None:
        """inflation_metrics property returns correct metrics."""
        inflation = InflationMetrics(lambda_gc=1.05, lambda_1000=1.02)
        gwas = GWASResult(
            output_path=Path("/path/to/output.PHENO1.glm.logistic.hybrid"),
            phenotype_type="binary",
            inflation=inflation,
            success=True,
        )
        result = AssocResult(gwas_result=gwas)
        assert result.inflation_metrics is not None
        assert result.inflation_metrics.lambda_gc == 1.05

    def test_to_dict_pca_only(self) -> None:
        """to_dict with only PCA."""
        pca = PCAResult(
            eigenvec_path=Path("/path/to/output.eigenvec"),
            eigenval_path=Path("/path/to/output.eigenval"),
            n_pcs=10,
            success=True,
        )
        result = AssocResult(pca_result=pca, success=True)
        d = result.to_dict()
        assert d["pass"] is True
        assert "pca" in d
        assert "gwas" not in d

    def test_to_dict_gwas_only(self) -> None:
        """to_dict with only GWAS."""
        gwas = GWASResult(
            output_path=Path("/path/to/output.PHENO1.glm.logistic.hybrid"),
            phenotype_type="binary",
            inflation=InflationMetrics(lambda_gc=1.05),
            success=True,
        )
        result = AssocResult(gwas_result=gwas, success=True)
        d = result.to_dict()
        assert d["pass"] is True
        assert "gwas" in d
        assert "pca" not in d

    def test_to_dict_both(self) -> None:
        """to_dict with both PCA and GWAS."""
        pca = PCAResult(
            eigenvec_path=Path("/path/to/output.eigenvec"),
            eigenval_path=Path("/path/to/output.eigenval"),
            n_pcs=10,
            success=True,
        )
        gwas = GWASResult(
            output_path=Path("/path/to/output.PHENO1.glm.logistic.hybrid"),
            phenotype_type="binary",
            inflation=InflationMetrics(lambda_gc=1.05),
            success=True,
        )
        result = AssocResult(pca_result=pca, gwas_result=gwas, success=True)
        d = result.to_dict()
        assert d["pass"] is True
        assert "pca" in d
        assert "gwas" in d
