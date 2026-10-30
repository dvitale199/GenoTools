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

"""The `Plot` record every plot module exports a list of."""

from dataclasses import dataclass
from typing import Callable

from genotools.plotting.report import Report


@dataclass(frozen=True)
class Plot:
    """One plot: how to tell if the report supports it, and how to draw it.

    `png` and `plotly` both take the loaded report; `png` also takes pyplot and
    the path to write, and returns that path. `available` must be cheap and
    must not raise: a report without the plot's data is a skip, not an error.
    """

    name: str
    title: str
    description: str
    available: Callable[[Report], bool]
    png: Callable
    plotly: Callable
