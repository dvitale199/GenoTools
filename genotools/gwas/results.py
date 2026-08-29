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

"""GWAS/Association result dataclasses.

This module provides immutable result classes for PCA and GWAS analyses.
Each result class includes a to_dict() method for backward compatibility
with the legacy pipeline output format.

Usage:
    from genotools.gwas.results import PCAResult, GWASResult, AssocResult

    # Results are returned by the step functions
    pca_result = run_pca(data, config, out_path)
    print(pca_result.eigenvec_path)
    print(pca_result.to_dict())  # Legacy format
"""

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, Optional


@dataclass(frozen=True)
class PruningResult:
    """Result of variant pruning before PCA.

    Attributes:
        output_path: Path to pruned pfiles (without extension).
        variants_before: Number of variants before pruning.
        variants_after: Number of variants after pruning.
        exclusion_file: Path to the exclusion regions file used.
        success: Whether pruning completed successfully.
        log: PLINK log output.
    """

    output_path: Path
    variants_before: int
    variants_after: int
    exclusion_file: Optional[Path] = None
    success: bool = True
    log: str = ""

    @property
    def variants_removed(self) -> int:
        """Number of variants removed by pruning."""
        return self.variants_before - self.variants_after

    def to_dict(self) -> Dict[str, Any]:
        """Convert to legacy dictionary format."""
        return {
            "pass": self.success,
            "step": "plink_pruning",
            "output": str(self.output_path),
        }


@dataclass(frozen=True)
class PCAResult:
    """Result of PCA computation.

    Attributes:
        eigenvec_path: Path to eigenvector file (.eigenvec).
        eigenval_path: Path to eigenvalue file (.eigenval).
        n_pcs: Number of principal components computed.
        n_samples: Number of samples in PCA output.
        pruning_result: Result of variant pruning step.
        success: Whether PCA completed successfully.
        log: PLINK log output.
    """

    eigenvec_path: Optional[Path]
    eigenval_path: Optional[Path]
    n_pcs: int
    n_samples: int = 0
    pruning_result: Optional[PruningResult] = None
    success: bool = True
    log: str = ""

    def to_dict(self) -> Dict[str, Any]:
        """Convert to legacy dictionary format.

        Returns format expected by existing pipeline:
        {
            'pass': bool,
            'step': 'plink_pca',
            'output': {
                'eigenvec': str|None,
                'eigenval': str|None
            }
        }
        """
        return {
            "pass": self.success,
            "step": "plink_pca",
            "output": {
                "eigenvec": str(self.eigenvec_path) if self.eigenvec_path else None,
                "eigenval": str(self.eigenval_path) if self.eigenval_path else None,
            },
        }


@dataclass(frozen=True)
class InflationMetrics:
    """Genomic inflation (lambda) metrics.

    Attributes:
        lambda_gc: Raw genomic inflation factor (lambda GC).
        lambda_1000: Sample-size-normalized lambda (to 1000 cases + 1000 controls).
            Only computed for case-control studies.
        lambda_gc_maf: Lambda GC restricted to variants with MAF >= threshold.
            Only computed if maf_lambdas is True.
        lambda_1000_maf: Normalized lambda for MAF-filtered variants.
            Only computed if maf_lambdas is True.
        n_cases: Number of cases (for case-control studies).
        n_controls: Number of controls (for case-control studies).
    """

    lambda_gc: float
    lambda_1000: Optional[float] = None
    lambda_gc_maf: Optional[float] = None
    lambda_1000_maf: Optional[float] = None
    n_cases: Optional[int] = None
    n_controls: Optional[int] = None

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for JSON serialization."""
        result: Dict[str, Any] = {
            "lambda": self.lambda_gc,
            "lambda1000": self.lambda_1000,
            "cases": self.n_cases,
            "controls": self.n_controls,
        }
        if self.lambda_gc_maf is not None:
            result["lambda_maf"] = self.lambda_gc_maf
        if self.lambda_1000_maf is not None:
            result["lambda1000_maf"] = self.lambda_1000_maf
        return result


@dataclass(frozen=True)
class GWASResult:
    """Result of GWAS analysis.

    Attributes:
        output_path: Path to GWAS results file (.glm.logistic.hybrid or .glm.linear).
        phenotype_type: Type of phenotype ('binary' or 'quantitative').
        inflation: Genomic inflation metrics.
        n_variants_tested: Number of variants tested.
        success: Whether GWAS completed successfully.
        log: PLINK log output.
    """

    output_path: Optional[Path]
    phenotype_type: str  # 'binary' or 'quantitative'
    inflation: Optional[InflationMetrics] = None
    n_variants_tested: int = 0
    success: bool = True
    log: str = ""

    def to_dict(self) -> Dict[str, Any]:
        """Convert to legacy dictionary format.

        Returns format expected by existing pipeline:
        {
            'pass': bool,
            'step': 'GWAS',
            'metrics': {
                'lambda': float,
                'lambda1000': float|None,
                'cases': int|None,
                'controls': int|None,
                ...
            },
            'output': {
                'gwas_output': str|None
            }
        }
        """
        metrics: Dict[str, Any]
        if self.inflation:
            metrics = self.inflation.to_dict()
        else:
            metrics = {
                "lambda": None,
                "lambda1000": None,
                "cases": None,
                "controls": None,
            }

        return {
            "pass": self.success,
            "step": "GWAS",
            "metrics": metrics,
            "output": {
                "gwas_output": str(self.output_path) if self.output_path else None,
            },
        }


@dataclass(frozen=True)
class AssocResult:
    """Combined result of association analysis (PCA + GWAS).

    Attributes:
        pca_result: Result of PCA computation, or None if PCA was not run.
        gwas_result: Result of GWAS analysis, or None if GWAS was not run.
        success: Overall success flag (True if all requested analyses succeeded).
    """

    pca_result: Optional[PCAResult] = None
    gwas_result: Optional[GWASResult] = None
    success: bool = True

    def __post_init__(self) -> None:
        # Validate that at least one result is present
        if self.pca_result is None and self.gwas_result is None:
            raise ValueError("At least one of pca_result or gwas_result must be provided")

    def to_dict(self) -> Dict[str, Any]:
        """Convert to legacy dictionary format.

        Returns format expected by existing pipeline:
        {
            'pass': bool,
            'pca': {...},  # if PCA was run
            'gwas': {...}  # if GWAS was run
        }
        """
        result: Dict[str, Any] = {"pass": self.success}

        if self.pca_result is not None:
            result["pca"] = self.pca_result.to_dict()

        if self.gwas_result is not None:
            result["gwas"] = self.gwas_result.to_dict()

        return result

    @property
    def eigenvec_path(self) -> Optional[Path]:
        """Path to eigenvector file, if PCA was run."""
        if self.pca_result:
            return self.pca_result.eigenvec_path
        return None

    @property
    def gwas_output_path(self) -> Optional[Path]:
        """Path to GWAS results file, if GWAS was run."""
        if self.gwas_result:
            return self.gwas_result.output_path
        return None

    @property
    def inflation_metrics(self) -> Optional[InflationMetrics]:
        """Inflation metrics, if GWAS was run."""
        if self.gwas_result:
            return self.gwas_result.inflation
        return None
