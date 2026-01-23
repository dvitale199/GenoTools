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

"""GWAS/Association analysis module.

This module provides tools for genome-wide association analysis:
- PCA computation for population structure
- GWAS execution with logistic/linear regression
- Genomic inflation (lambda) calculation

Architecture:
- **Config classes**: Typed, validated configuration (config.py)
- **Result classes**: Immutable result containers (results.py)
- **Step functions**: Pure functions for each operation (steps/)
- **Pipeline**: Orchestration of PCA + GWAS (pipeline.py)

Usage:
    from genotools.gwas import (
        AssocPipeline,
        AssocConfig,
        PCAConfig,
        GWASConfig,
        run_association,
    )
    from genotools.core import GenotypeData

    # Simple usage with defaults
    data = GenotypeData.from_path("/path/to/genotypes")
    result = run_association(data, Path("/path/to/output"))

    # Custom configuration
    config = AssocConfig(
        pca=PCAConfig(n_pcs=20, build='hg38'),
        gwas=GWASConfig(maf_lambdas=True),
        run_pca=True,
        run_gwas=True,
    )
    pipeline = AssocPipeline(config)
    result = pipeline.run(data, Path("/path/to/output"))

    # Access results
    print(f"Lambda: {result.inflation_metrics.lambda_gc}")
    print(f"Eigenvecs: {result.eigenvec_path}")

    # Legacy format for backward compatibility
    legacy_dict = result.to_dict()
"""

# Configuration classes
from genotools.gwas.config import (
    PCAPruningConfig,
    PCAConfig,
    GWASConfig,
    CovariateConfig,
    AssocConfig,
    PLINK2_GLM_OPTIONS,
    HG19_EXCLUSION_REGIONS,
    HG38_EXCLUSION_REGIONS,
    get_exclusion_regions,
)

# Result classes
from genotools.gwas.results import (
    PruningResult,
    PCAResult,
    InflationMetrics,
    GWASResult,
    AssocResult,
)

# Step functions
from genotools.gwas.steps import (
    run_pca_pruning,
    run_pca,
    run_gwas,
    calculate_inflation,
)

# Pipeline
from genotools.gwas.pipeline import (
    AssocPipeline,
    run_association,
)

__all__ = [
    # Config
    "PCAPruningConfig",
    "PCAConfig",
    "GWASConfig",
    "CovariateConfig",
    "AssocConfig",
    "PLINK2_GLM_OPTIONS",
    "HG19_EXCLUSION_REGIONS",
    "HG38_EXCLUSION_REGIONS",
    "get_exclusion_regions",
    # Results
    "PruningResult",
    "PCAResult",
    "InflationMetrics",
    "GWASResult",
    "AssocResult",
    # Steps
    "run_pca_pruning",
    "run_pca",
    "run_gwas",
    "calculate_inflation",
    # Pipeline
    "AssocPipeline",
    "run_association",
]
