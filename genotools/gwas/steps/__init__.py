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

"""GWAS/Association step functions.

This module exports the pure functions for PCA and GWAS analysis.

Usage:
    from genotools.gwas.steps import run_pca_pruning, run_pca, run_gwas

    # Run PCA pruning
    pruning_result = run_pca_pruning(data, config, out_path)

    # Run PCA
    pca_result = run_pca(data, config, out_path)

    # Run GWAS
    gwas_result = run_gwas(data, config, out_path, covar_path=eigenvec_path)
"""

from genotools.gwas.steps.pca import run_pca_pruning, run_pca
from genotools.gwas.steps.association import run_gwas, calculate_inflation

__all__ = [
    "run_pca_pruning",
    "run_pca",
    "run_gwas",
    "calculate_inflation",
]
