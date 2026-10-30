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

"""One module per plot area, each exporting `PLOTS`.

`ALL_PLOTS` is their concatenation, in the order the HTML report shows them.
A new plot goes into its area's `PLOTS`; neither this file nor `render.py`
changes.
"""

from typing import List

from genotools.plotting.plots import ancestry, gwas, qc_summary, sample_qc, variant_qc
from genotools.plotting.plots.base import Plot

ALL_PLOTS: List[Plot] = [
    *qc_summary.PLOTS,
    *sample_qc.PLOTS,
    *variant_qc.PLOTS,
    *ancestry.PLOTS,
    *gwas.PLOTS,
]

__all__ = ["ALL_PLOTS", "Plot"]
