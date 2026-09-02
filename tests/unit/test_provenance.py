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

"""What software produced a result: library versions and external executables.

Both halves exist because the answer used to be unrecoverable after the fact —
a model fit under one umap and loaded under another produces different ancestry
calls with no error, and nothing recorded which plink2 ran.
"""

import os
import stat
from pathlib import Path

import pytest

from genotools import __version__
from genotools.core.provenance import (
    MODEL_PACKAGES,
    NOT_INSTALLED,
    describe_tool,
    executable_folder,
    package_versions,
    requirements_lines,
    tool_lines,
    tool_version,
    version_drift,
)


def _fake_tool(directory: Path, name: str, says: str) -> Path:
    """Write an executable that answers --version with ``says``."""
    directory.mkdir(parents=True, exist_ok=True)
    path = directory / name
    path.write_text(f'#!/bin/sh\necho "{says}"\n')
    path.chmod(path.stat().st_mode | stat.S_IXUSR)
    return path


class TestPackageVersions:
    def test_records_the_libraries_that_shape_a_fit(self) -> None:
        versions = package_versions()
        # umap-learn is the one that motivated this: it is pinned precisely
        # because a different version embeds differently.
        assert "umap-learn" in versions
        for package in MODEL_PACKAGES:
            assert package in versions, f"{package} is not recorded"

    def test_always_carries_genotools_and_python(self) -> None:
        versions = package_versions(())
        assert versions["genotools"] == __version__
        assert versions["python"].count(".") == 2, versions["python"]

    def test_a_missing_package_is_recorded_not_dropped(self) -> None:
        """Dropping it would make "absent then, present now" invisible."""
        versions = package_versions(["definitely-not-a-real-distribution"])
        assert versions["definitely-not-a-real-distribution"] == NOT_INSTALLED


class TestVersionDrift:
    def test_reports_only_what_changed(self) -> None:
        recorded = {"umap-learn": "0.5.3", "numpy": "1.23.5"}
        current = {"umap-learn": "0.5.3", "numpy": "2.3.5"}
        assert version_drift(recorded, current) == [("numpy", "1.23.5", "2.3.5")]

    def test_identical_environments_have_no_drift(self) -> None:
        recorded = package_versions()
        assert version_drift(recorded) == []

    def test_a_key_missing_from_current_is_not_drift(self) -> None:
        """An explicit ``current`` that omits a key cannot speak to it."""
        recorded = {"umap-learn": "0.5.3", "retired-metric": "1.0"}
        drift = version_drift(recorded, {"umap-learn": "0.5.3"})
        assert drift == []

    def test_an_uninstalled_package_is_drift(self) -> None:
        """Recorded then, gone now, is a real change to the environment."""
        drift = version_drift({"definitely-not-a-real-distribution": "1.0"})
        assert drift == [
            ("definitely-not-a-real-distribution", "1.0", NOT_INSTALLED)
        ]

    def test_preserves_the_recorded_order(self) -> None:
        recorded = {"a": "1", "b": "1", "c": "1"}
        current = {"a": "2", "b": "2", "c": "2"}
        assert [name for name, _, _ in version_drift(recorded, current)] == [
            "a", "b", "c"
        ]


class TestRequirementsLines:
    """Detecting drift is half of it; being able to undo it is the other half."""

    def test_versions_become_installable_pins(self) -> None:
        lines = requirements_lines({"umap-learn": "0.5.12", "numpy": "2.3.5"})
        assert "umap-learn==0.5.12" in lines
        assert "numpy==2.3.5" in lines

    def test_python_is_a_comment_not_a_requirement(self) -> None:
        """pip cannot install an interpreter, but dropping it silently would
        hide a difference that does affect unpickling."""
        lines = requirements_lines({"python": "3.11.2", "numpy": "2.3.5"})
        assert lines[0] == "# python 3.11.2"
        assert not any(line.startswith("python==") for line in lines)

    def test_an_absent_package_is_a_comment(self) -> None:
        """"not installed" is not a version pip can resolve."""
        lines = requirements_lines({"scipy": NOT_INSTALLED})
        assert not any(line.startswith("scipy==") for line in lines)
        assert any("scipy" in line and line.startswith("#") for line in lines)

    def test_old_umap_drags_in_the_setuptools_pin(self) -> None:
        """Otherwise the file recreates an environment that cannot import umap
        at all: umap<0.5.5 needs pkg_resources, which setuptools 82 deleted.
        A reproducibility file that does not reproduce is worse than none."""
        lines = requirements_lines({"umap-learn": "0.5.3"})
        assert any(line.startswith("setuptools<82") for line in lines)

    def test_new_umap_does_not(self) -> None:
        assert not any(
            line.startswith("setuptools")
            for line in requirements_lines({"umap-learn": "0.5.5"})
        )

    def test_version_compare_is_numeric_not_lexicographic(self) -> None:
        """"0.5.12" < "0.5.5" as strings. Getting this wrong would pin
        setuptools for every modern umap."""
        assert not any(
            line.startswith("setuptools")
            for line in requirements_lines({"umap-learn": "0.5.12"})
        )

    def test_an_unparseable_version_adds_no_bogus_pin(self) -> None:
        """Guessing wrong here corrupts a requirements file; omitting a helpful
        pin merely leaves it as it was."""
        for odd in ("weird-build", "", NOT_INSTALLED):
            assert not any(
                line.startswith("setuptools")
                for line in requirements_lines({"umap-learn": odd})
            )

    def test_a_model_with_no_umap_recorded_is_left_alone(self) -> None:
        assert requirements_lines({"numpy": "2.3.5"}) == ["numpy==2.3.5"]


class TestExecutableFolder:
    def test_honours_the_dep_dir_override(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """Must match dependencies.__get_executable_folder, or provenance
        describes a different plink2 than the one that actually runs."""
        monkeypatch.setenv("GENOTOOLS_DEP_DIR", str(tmp_path))
        assert executable_folder() == tmp_path

    def test_defaults_to_the_genotools_home(
        self, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        monkeypatch.delenv("GENOTOOLS_DEP_DIR", raising=False)
        assert executable_folder() == (
            Path.home() / ".genotools" / "misc" / "executables"
        )

    def test_matches_what_dependencies_actually_resolves(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """Pin the two implementations together rather than trusting the copy.

        ``executable_folder`` mirrors ``dependencies.__get_executable_folder``
        so it can probe without downloading. If the two ever diverge, every run
        log confidently names a plink2 that is not the one that ran.
        """
        import importlib

        import genotools.dependencies as dependencies

        deps_dir = tmp_path / "deps"
        monkeypatch.setenv("GENOTOOLS_DEP_DIR", str(deps_dir))
        try:
            importlib.reload(dependencies)
            # Module-level double-underscore names are not mangled at module
            # scope, but a lookup from inside this method would be - hence
            # __dict__. A KeyError here means dependencies was restructured
            # and this test needs rewriting, not deleting.
            resolved = Path(dependencies.__dict__["__executable_folder"])
            assert resolved == deps_dir, "the reload did not pick up the override"
            assert resolved == executable_folder()

            # And the default branch, which is the one a real user hits.
            monkeypatch.delenv("GENOTOOLS_DEP_DIR", raising=False)
            importlib.reload(dependencies)
            default = Path(dependencies.__dict__["__executable_folder"])
            assert default == executable_folder()
        finally:
            monkeypatch.delenv("GENOTOOLS_DEP_DIR", raising=False)
            importlib.reload(dependencies)


class TestDescribeTool:
    def test_an_unfetched_tool_is_named_not_downloaded(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """Describing a tool must never trigger the download that resolving it
        would - a callrate-only run should not pull PLINK 1.9 to name it."""
        monkeypatch.setenv("GENOTOOLS_DEP_DIR", str(tmp_path))
        monkeypatch.setattr("shutil.which", lambda name: None)

        info = describe_tool("plink")
        assert info.path is None
        assert "not fetched" in info.version
        assert not any(tmp_path.iterdir()), "describe_tool downloaded something"

    def test_reports_the_resolved_path_and_version(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        binary = _fake_tool(tmp_path, "plink2", "PLINK v2.00a5.10LM 64-bit")
        monkeypatch.setenv("GENOTOOLS_DEP_DIR", str(tmp_path))
        monkeypatch.setattr("shutil.which", lambda name: None)

        info = describe_tool("plink2")
        assert info.path == binary
        assert info.version == "PLINK v2.00a5.10LM 64-bit"
        assert info.shadowed_by is None

    def test_a_different_tool_on_path_is_flagged(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """The point of item 12: GenoTools ignores PATH, and "I upgraded
        plink2 and nothing changed" is otherwise very hard to work out."""
        owned = tmp_path / "owned"
        _fake_tool(owned, "plink2", "PLINK v2.00a5.10LM 64-bit")
        other = _fake_tool(tmp_path / "elsewhere", "plink2", "PLINK v2.00a6LM")
        monkeypatch.setenv("GENOTOOLS_DEP_DIR", str(owned))
        monkeypatch.setattr("shutil.which", lambda name: str(other))

        info = describe_tool("plink2")
        assert info.shadowed_by == other
        assert info.version == "PLINK v2.00a5.10LM 64-bit", (
            "must report the build GenoTools runs, not the one on PATH"
        )

    def test_the_same_tool_on_path_is_not_flagged(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """A false 'a different plink2 is on PATH' would train users to ignore
        the note that matters."""
        binary = _fake_tool(tmp_path, "plink2", "PLINK v2.00a5.10LM 64-bit")
        monkeypatch.setenv("GENOTOOLS_DEP_DIR", str(tmp_path))
        monkeypatch.setattr("shutil.which", lambda name: str(binary))

        assert describe_tool("plink2").shadowed_by is None

    def test_a_tool_that_cannot_run_degrades_to_unknown(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """Provenance must never be the reason a run fails."""
        broken = tmp_path / "plink2"
        broken.write_text("not an executable")
        monkeypatch.setenv("GENOTOOLS_DEP_DIR", str(tmp_path))
        monkeypatch.setattr("shutil.which", lambda name: None)

        assert describe_tool("plink2").version == "unknown"

    def test_version_read_from_stderr_when_stdout_is_empty(
        self, tmp_path: Path
    ) -> None:
        path = tmp_path / "noisy"
        path.write_text('#!/bin/sh\necho "KING 2.3.2" >&2\n')
        path.chmod(path.stat().st_mode | stat.S_IXUSR)
        assert tool_version(path) == "KING 2.3.2"

    def test_only_the_first_line_is_kept(self, tmp_path: Path) -> None:
        """PLINK 1.9 prints a licence blurb after its version line."""
        path = _fake_tool(tmp_path, "plink", "PLINK v1.9.0-b.8\nGPLv3 blurb")
        assert tool_version(path) == "PLINK v1.9.0-b.8"


class TestToolLines:
    def test_one_line_per_tool_plus_a_shadow_note(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        owned = tmp_path / "owned"
        _fake_tool(owned, "plink2", "PLINK v2.00a5.10LM 64-bit")
        other = _fake_tool(tmp_path / "elsewhere", "plink2", "PLINK v2.00a6LM")
        monkeypatch.setenv("GENOTOOLS_DEP_DIR", str(owned))
        monkeypatch.setattr("shutil.which", lambda name: str(other))

        lines = tool_lines(["plink2"])
        assert len(lines) == 2
        assert lines[0].startswith("plink2: ")
        assert str(owned / "plink2") in lines[0]
        assert "PATH" in lines[1] and str(other) in lines[1]

    def test_covers_every_tool_the_executors_resolve(self) -> None:
        """A tool missing from the default list leaves a silent gap in the
        provenance of every run that uses it.

        The executors are the only place GenoTools resolves an external binary,
        so the set of ``check_*`` calls there is the set worth describing.
        """
        import re

        import genotools.core.executors as executors

        source = Path(executors.__file__).read_text()
        resolved = {
            name.lower()
            for name in re.findall(r"check_(\w+)\(\)", source)
        }
        assert resolved, "no check_* calls found; executors was restructured"

        described = {line.split(":", 1)[0] for line in tool_lines()}
        assert resolved <= described, (
            f"executors resolve {sorted(resolved - described)}, "
            f"which no run log would record"
        )
