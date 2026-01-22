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

"""QC module for genotype quality control.

This module provides pure functions for sample and variant quality control,
replacing the mutable SampleQC/VariantQC class-based architecture.

Example usage:
    from genotools.qc import (
        QCPipeline,
        filter_callrate,
        CallrateConfig,
        FilterResult,
    )

    # Run a single step
    data = GenotypeData.from_path(Path("input"))
    result = filter_callrate(data, CallrateConfig(mind=0.05), Path("output"))

    # Or compose into a pipeline
    pipeline = QCPipeline(
        steps=[
            ("callrate", filter_callrate, CallrateConfig(mind=0.05)),
            ("sex", filter_sex, SexConfig()),
        ]
    )
    qc_result = pipeline.run(data, Path("output_dir"))
"""

from genotools.qc.config import (
    CallrateConfig,
    CaseControlConfig,
    GenoConfig,
    HaplotypeConfig,
    HetConfig,
    HWEConfig,
    LDConfig,
    RelatedConfig,
    SexConfig,
)
from genotools.qc.pipeline import QCPipeline
from genotools.qc.results import FilterResult, QCResult
from genotools.qc.steps import (
    filter_callrate,
    filter_case_control,
    filter_haplotype,
    filter_heterozygosity,
    filter_hwe,
    filter_relatedness,
    filter_sex,
    filter_variant_missingness,
    prune_ld,
    verify_kinship,
)

__all__ = [
    # Pipeline
    "QCPipeline",
    # Results
    "FilterResult",
    "QCResult",
    # Configs
    "CallrateConfig",
    "CaseControlConfig",
    "GenoConfig",
    "HaplotypeConfig",
    "HetConfig",
    "HWEConfig",
    "LDConfig",
    "RelatedConfig",
    "SexConfig",
    # Step functions
    "filter_callrate",
    "filter_case_control",
    "filter_haplotype",
    "filter_heterozygosity",
    "filter_hwe",
    "filter_relatedness",
    "filter_sex",
    "filter_variant_missingness",
    "prune_ld",
    "verify_kinship",
]
