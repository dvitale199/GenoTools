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

"""Consolidated run-log integration tests.

core/logging.py::setup_logging() was never called at runtime, so the root
logger's WARNING default dropped every step's logger.info/error and the
consolidated {out}_all_logs.log came out header-only. These tests run the real
CLI and assert the consolidated log is actually populated with structured,
step-tagged records.

They require PLINK2 and the synthetic test data; they skip cleanly otherwise.
"""

import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]


def _run_cli(args: list[str]) -> subprocess.CompletedProcess:
    return subprocess.run(
        [sys.executable, "-m", "genotools", *args],
        cwd=str(REPO_ROOT),
        capture_output=True,
        text=True,
        timeout=300,
    )


def test_run_log_contains_step_markers(test_geno_path: Path, tmp_path: Path) -> None:
    """The consolidated {out}_all_logs.log must contain the steps' structured logs."""
    out = tmp_path / "logtest"

    res = _run_cli(
        [
            "--pfile", str(test_geno_path),
            "--out", str(out),
            "--callrate",
            "--geno",
        ]
    )
    assert res.returncode == 0, (
        f"Pipeline failed.\nSTDOUT:\n{res.stdout}\nSTDERR:\n{res.stderr}"
    )

    all_logs = Path(f"{out}_all_logs.log")
    assert all_logs.exists(), "consolidated {out}_all_logs.log was not created"

    content = all_logs.read_text()

    # The banner header alone is short; a populated log is materially longer.
    header_only_len = 400
    assert len(content) > header_only_len, (
        "consolidated log looks header-only (structured step logs were dropped):\n"
        f"{content!r}"
    )

    # Structured step-context markers written by core.logging.step_context.
    assert "[callrate_prune]" in content, "no callrate step marker in run log"
    assert "[geno_prune]" in content, "no geno step marker in run log"

    # Actual step messages emitted via logger.info must be present.
    assert "filtering complete" in content.lower(), (
        "step completion messages missing from run log"
    )
