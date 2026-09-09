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

"""Resolving which GenoTools is actually running.

`genotools.__version__` answers a different question than the one users ask.
From any directory holding a `genotools/` package -- the repo root, most
obviously -- `import genotools` picks up the *source tree* rather than the
installed distribution, because the cwd leads `sys.path` under
`python -m genotools`. The attribute then reports whatever the checkout says,
which may be neither what pip installed nor what is on `PATH`.

The distribution's metadata answers the question people mean -- "which install
is this?" -- because it is written by pip at install time and looked up by
distribution name rather than by import path. So `--version` reports the
distribution, and falls back to the imported package only when none is
installed (a source checkout run in place), where it says so.

Two consequences worth knowing, both correct rather than bugs:

- Under an **editable install** the metadata records the version that was
  current when `pip install -e .` ran, so bumping `__version__` in the source
  does not move `--version` until the install is redone. That is the honest
  answer: nothing was reinstalled.
- A directory on `sys.path` holding a built `*.egg-info` or `*.dist-info` --
  the repo root after a `python -m build` -- is itself a discoverable
  distribution, so from there the lookup can find the checkout's metadata
  ahead of site-packages. Reporting a version from a directory you built in is
  the one case where cwd still matters.
"""

from __future__ import annotations

from importlib import metadata

DISTRIBUTION = "the_real_genotools"


def resolve_version() -> str:
    """Return the installed distribution's version.

    Returns:
        The version recorded by pip for `the_real_genotools`. If no such
        distribution is installed, the imported package's `__version__` with a
        marker saying it came from a source tree rather than an install.
    """
    try:
        return metadata.version(DISTRIBUTION)
    except metadata.PackageNotFoundError:
        from .. import __version__

        return f"{__version__} (source tree; {DISTRIBUTION} is not installed)"


def version_string() -> str:
    """Return the full line printed by `genotools --version`."""
    return f"genotools {resolve_version()}"
