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

"""Golden files must not carry per-run noise.

Regenerating the goldens has to produce a diff only when behaviour actually
changed. It did not: run logs churned on every regeneration (timestamps, the
random temp-directory name, hostname, home directory, PLINK's RNG seed and
detected RAM/threads), and the temp-directory name leaked into the JSON
report's pass_fail paths as well. A regeneration that always diffs is a
regeneration nobody reads, which is how the generator stayed broken for
several rounds without anyone noticing.
"""

import importlib.util
import subprocess
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
GENERATOR = REPO_ROOT / "tests" / "scripts" / "generate_golden.py"


def _load_generator():
    """Import the generator script, which is not part of the package."""
    spec = importlib.util.spec_from_file_location("generate_golden", GENERATOR)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class TestNoTrackedGoldenLogs:
    """Run logs are byproducts, not references."""

    def test_no_log_files_are_tracked_under_golden(self) -> None:
        """`git add -A` after a regeneration must not re-add them.

        .gitignore is the guard; this asserts it actually holds. Nothing reads
        these files - tests/regression/test_logging.py owns the log contract
        and asserts it against a fresh run in tmp_path.
        """
        tracked = subprocess.run(
            ["git", "ls-files", "tests/regression/golden/"],
            cwd=REPO_ROOT,
            capture_output=True,
            text=True,
        )
        if tracked.returncode != 0:
            pytest.skip("not a git checkout")

        logs = [f for f in tracked.stdout.split() if f.endswith(".log")]
        assert logs == [], (
            f"{len(logs)} golden run log(s) are tracked again; they rewrite "
            f"themselves on every regeneration. First few: {logs[:3]}"
        )


class TestRunPathNormalization:
    """The temp-directory name must not reach a committed golden."""

    def test_replaces_the_random_temp_directory_name(self, tmp_path: Path) -> None:
        """Revert check: drop the _normalize_run_paths call from the generator
        and a regeneration diffs output.json every time."""
        module = _load_generator()
        report = tmp_path / "output.json"
        report.write_text(
            '{"pass_fail": {"callrate": {"output": '
            '"golden/all_sample/.3ic3nzix_tmp/output_callrate"}}}'
        )

        module._normalize_run_paths(report)

        assert ".3ic3nzix_tmp" not in report.read_text()
        assert ".GOLDEN_tmp" in report.read_text()

    def test_is_idempotent(self, tmp_path: Path) -> None:
        """Normalizing an already-normalized report changes nothing, or the
        placeholder itself would churn."""
        module = _load_generator()
        report = tmp_path / "output.json"
        report.write_text('{"output": "golden/.GOLDEN_tmp/out"}')
        before = report.read_text()

        module._normalize_run_paths(report)

        assert report.read_text() == before

    def test_leaves_real_paths_alone(self, tmp_path: Path) -> None:
        """Only the 8-character temp name matches; a real output path that
        happens to contain "tmp" must survive."""
        module = _load_generator()
        report = tmp_path / "output.json"
        report.write_text('{"output": "/data/my_tmp_run/output"}')

        module._normalize_run_paths(report)

        assert report.read_text() == '{"output": "/data/my_tmp_run/output"}'

    def test_missing_report_is_not_an_error(self, tmp_path: Path) -> None:
        """A step that writes no JSON must not break the generator."""
        module = _load_generator()
        module._normalize_run_paths(tmp_path / "absent.json")
