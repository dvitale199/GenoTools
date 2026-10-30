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

"""Plots rendered from a finished run's JSON report.

Downstream of the pipeline, never part of it: this package reads the report a
run wrote and imports nothing from `cli`, `qc`, `gwas` or `ancestry`, so no
plotting failure can reach a run. `genotools --plots` calls `render_plots`
after the report is written; `genotools-plot` does the same for an existing
report.
"""

from genotools.plotting.render import PLOTS, RenderResult, default_out_dir, render_plots
from genotools.plotting.report import Report, attrition_frame, load_report, outcome_frame

__all__ = [
    "PLOTS",
    "RenderResult",
    "Report",
    "attrition_frame",
    "default_out_dir",
    "load_report",
    "outcome_frame",
    "render_plots",
]
