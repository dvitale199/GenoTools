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

"""QC step functions.

Each step is a pure function that takes:
    - data: GenotypeData (immutable input)
    - config: Step-specific configuration dataclass
    - out_path: Path for output files

And returns:
    - FilterResult with filtered output and metrics

Example:
    from genotools.qc.steps import filter_callrate
    from genotools.qc.config import CallrateConfig
    from genotools.core.genotypes import GenotypeData

    data = GenotypeData.from_path(Path("input"))
    result = filter_callrate(data, CallrateConfig(mind=0.05), Path("output"))
    print(f"Removed {result.samples_removed} samples")
"""

from genotools.qc.steps.callrate import filter_callrate
from genotools.qc.steps.case_control import filter_case_control
from genotools.qc.steps.haplotype import filter_haplotype
from genotools.qc.steps.heterozygosity import filter_heterozygosity
from genotools.qc.steps.hwe import filter_hwe
from genotools.qc.steps.ld_prune import prune_ld
from genotools.qc.steps.relatedness import filter_relatedness, verify_kinship
from genotools.qc.steps.sex import filter_sex
from genotools.qc.steps.variant_missingness import filter_variant_missingness

__all__ = [
    # Sample QC steps
    "filter_callrate",
    "filter_sex",
    "filter_heterozygosity",
    "filter_relatedness",
    "verify_kinship",
    # Variant QC steps
    "filter_variant_missingness",
    "filter_case_control",
    "filter_haplotype",
    "filter_hwe",
    "prune_ld",
]
