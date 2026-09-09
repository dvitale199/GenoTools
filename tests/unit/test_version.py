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

"""Which GenoTools is running, and the release metadata that has to agree.

`--version` exists because every other way of asking answers a subtly different
question: `genotools.__version__` reports whichever `genotools/` package the
import found, which from a repo checkout is the source tree rather than the
install. These tests pin the distinction rather than the value.
"""

import re
from importlib import metadata
from pathlib import Path

import pytest

from genotools import __version__
from genotools.cli.parser import parse_args
from genotools.core.version import DISTRIBUTION, resolve_version, version_string

REPO_ROOT = Path(__file__).resolve().parents[2]


class TestResolveVersion:
    """`resolve_version` reads the install, not the import."""

    def test_reports_the_installed_distribution(self, monkeypatch):
        """The distribution's metadata wins over the imported package."""
        monkeypatch.setattr(metadata, "version", lambda name: "9.9.9")
        assert resolve_version() == "9.9.9"

    def test_asks_for_the_distribution_by_name(self, monkeypatch):
        """The name queried is the PyPI one, not the import name."""
        asked = []

        def record(name):
            asked.append(name)
            return "9.9.9"

        monkeypatch.setattr(metadata, "version", record)
        resolve_version()
        assert asked == ["the_real_genotools"]

    def test_falls_back_to_the_package_and_says_so(self, monkeypatch):
        """With nothing installed, the source tree answers -- labeled."""

        def missing(name):
            raise metadata.PackageNotFoundError(name)

        monkeypatch.setattr(metadata, "version", missing)
        answer = resolve_version()
        assert __version__ in answer
        assert "source tree" in answer

    def test_version_string_names_the_program(self):
        assert version_string().startswith("genotools ")


class TestVersionFlag:
    """The flag short-circuits the parser."""

    def test_prints_and_exits_zero(self, capsys):
        with pytest.raises(SystemExit) as exc:
            parse_args(["--version"])
        assert exc.value.code == 0
        assert capsys.readouterr().out.strip() == version_string()

    def test_does_not_require_the_normally_required_args(self):
        """`--out` is required for a run; asking the version is not a run."""
        with pytest.raises(SystemExit) as exc:
            parse_args(["--version"])
        assert exc.value.code == 0


class TestReleaseMetadata:
    """Metadata that has to be bumped together, and silently isn't."""

    def test_changelog_leads_with_the_current_version(self):
        """The top `## x.y.z` section matches `__version__`."""
        changelog = (REPO_ROOT / "CHANGELOG.md").read_text()
        first = re.search(r"^## (\d+\.\d+\.\d+)", changelog, re.MULTILINE)
        assert first is not None, "CHANGELOG.md has no version section"
        assert first.group(1) == __version__

    def test_python_requires_matches_the_ci_interpreter(self):
        """The declared floor is the interpreter CI actually tests."""
        setup = (REPO_ROOT / "setup.py").read_text()
        floor = re.search(r"python_requires='>=(\d+\.\d+)'", setup)
        assert floor is not None, "setup.py declares no python_requires"

        ci = (REPO_ROOT / ".github/workflows/ci.yml").read_text()
        tested = set(re.findall(r'python-version: "(\d+\.\d+)"', ci))
        assert tested, "CI declares no python-version"
        assert floor.group(1) == min(tested, key=lambda v: tuple(map(int, v.split("."))))
