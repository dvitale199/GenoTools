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

"""Consolidated + per-step run-log integration tests (round 7 redesign).

Exercise the real CLI and assert the round-7 logging contract:
  * ``{out}_all_logs.log`` is a per-step-sectioned consolidated log: banner,
    ``===== step =====`` headers, structured ``[step]`` records, the verbatim
    harvested PLINK output inlined after each step's summary, and a final run
    summary table;
  * per-step raw ``{out}_{step}.log`` files carry the exact PLINK output and
    persist even without ``--full_output``;
  * the dead ``{out}_cleaned_logs.log`` is gone;
  * the console (stderr) shows the concise curated stream, not the raw PLINK dump.

Requires PLINK2 and the synthetic test data; skips cleanly otherwise.
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


def _run_callrate_geno(test_geno_path: Path, out: Path) -> subprocess.CompletedProcess:
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
    return res


def test_consolidated_log_is_sectioned_with_raw(test_geno_path: Path, tmp_path: Path) -> None:
    """{out}_all_logs.log: banner, per-step sections, structured lines + raw PLINK."""
    out = tmp_path / "logtest"
    _run_callrate_geno(test_geno_path, out)

    all_logs = Path(f"{out}_all_logs.log")
    assert all_logs.exists(), "consolidated {out}_all_logs.log was not created"
    content = all_logs.read_text()

    # Materially longer than a bare banner.
    assert len(content) > 400, f"consolidated log looks header-only:\n{content!r}"

    # Per-step section headers.
    assert "===== callrate =====" in content, "missing callrate section header"
    assert "===== geno =====" in content, "missing geno section header"

    # Structured, step-tagged records.
    assert "[callrate_prune]" in content, "no callrate step marker in run log"
    assert "[geno_prune]" in content, "no geno step marker in run log"
    assert "filtering complete" in content.lower(), "step completion messages missing"

    # Verbatim harvested PLINK output inlined (raw source of truth).
    assert "PLINK v" in content, "raw PLINK output not inlined in consolidated log"

    # For each step, its structured summary precedes its raw PLINK block.
    i_marker = content.index("[callrate_prune]")
    i_raw = content.index("PLINK v")
    assert i_marker < i_raw, "raw PLINK block should follow the structured summary"

    # End-of-run summary table.
    assert "run summary" in content.lower(), "missing run summary block"
    assert "PASS" in content, "summary table missing pass/fail status"


def test_per_step_raw_logs_persist(test_geno_path: Path, tmp_path: Path) -> None:
    """Per-step {out}_{step}.log files carry raw PLINK output and persist (no --full_output)."""
    out = tmp_path / "logtest"
    _run_callrate_geno(test_geno_path, out)

    callrate_raw = Path(f"{out}_callrate.log")
    geno_raw = Path(f"{out}_geno.log")
    assert callrate_raw.exists(), "per-step {out}_callrate.log missing"
    assert geno_raw.exists(), "per-step {out}_geno.log missing"

    callrate_text = callrate_raw.read_text()
    assert "PLINK v" in callrate_text, "per-step raw log lacks PLINK output"
    assert "--mind" in callrate_text, "per-step raw log lacks the invoked command"

    # The per-step raw file is distinct from the consolidated log.
    assert callrate_text != Path(f"{out}_all_logs.log").read_text()


def test_cleaned_logs_is_gone(test_geno_path: Path, tmp_path: Path) -> None:
    """The dead 0-byte {out}_cleaned_logs.log is no longer produced."""
    out = tmp_path / "logtest"
    _run_callrate_geno(test_geno_path, out)
    assert not Path(f"{out}_cleaned_logs.log").exists(), (
        "cleaned_logs.log should be removed in the round-7 redesign"
    )


def test_console_is_curated_not_raw(test_geno_path: Path, tmp_path: Path) -> None:
    """The console (stderr) shows the concise step-tagged stream, not raw PLINK."""
    out = tmp_path / "logtest"
    res = _run_callrate_geno(test_geno_path, out)

    # Curated, step-tagged progress reaches the terminal.
    assert "[callrate_prune]" in res.stderr, "curated step line missing from console"
    assert "Run summary" in res.stderr, "run summary missing from console"

    # The verbose raw PLINK dump stays file-only.
    assert "PLINK v" not in res.stderr, "raw PLINK output leaked to the console"
    # No full file-format timestamps on the console (concise formatter).
    assert "genotools.qc.steps" not in res.stderr, "console used the full file format"
