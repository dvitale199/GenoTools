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

"""CLI failure-boundary tests.

Expected, user-actionable failures (``GenoToolsError``) must reach the terminal
as a one-line ``ERROR: ...`` message with a non-zero exit code -- not a Python
traceback. Unexpected exceptions are bugs and keep their traceback.

Requires the synthetic test data; skips cleanly otherwise.
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


def test_rerun_guard_reports_cleanly(test_geno_path: Path, tmp_path: Path) -> None:
    """Re-running on a used prefix: clean message, exit 1, no traceback."""
    out = tmp_path / "guarded"
    first = _run_cli(["--pfile", str(test_geno_path), "--out", str(out), "--callrate"])
    assert first.returncode == 0, f"first run failed:\n{first.stderr}"

    second = _run_cli(["--pfile", str(test_geno_path), "--out", str(out), "--callrate"])
    assert second.returncode == 1, "a blocked prefix must exit non-zero"
    assert "Traceback" not in second.stderr, (
        f"expected a message, got a traceback:\n{second.stderr}"
    )
    assert "ERROR:" in second.stderr
    assert "previously been run" in second.stderr
    # The suggested flag must be the one the parser actually accepts: the old
    # message said --skip_fails, which argparse rejects.
    assert "--skip-fails" in second.stderr
    assert "--skip_fails" not in second.stderr
    # Short enough to actually read.
    assert len(second.stderr.strip().splitlines()) <= 4, second.stderr


def test_no_steps_requested_reports_cleanly(
    test_geno_path: Path, tmp_path: Path
) -> None:
    """Nothing to do is a user error, not a crash."""
    out = tmp_path / "nosteps"
    res = _run_cli(["--pfile", str(test_geno_path), "--out", str(out)])

    assert res.returncode == 1
    assert "Traceback" not in res.stderr, res.stderr
    assert "ERROR:" in res.stderr
    assert "No QC steps" in res.stderr


def test_missing_input_reports_cleanly(tmp_path: Path) -> None:
    """A nonexistent --pfile is reported as a message, not a traceback.

    Missing files are signalled codebase-wide with FileNotFoundError, so the
    boundary handler has to cover it -- it is the most common user error.
    """
    res = _run_cli(
        [
            "--pfile", str(tmp_path / "does_not_exist"),
            "--out", str(tmp_path / "out"),
            "--callrate",
        ]
    )

    assert res.returncode == 1
    assert "Traceback" not in res.stderr, res.stderr
    assert "ERROR:" in res.stderr
    assert "does not exist" in res.stderr


def test_out_of_range_argument_reports_cleanly(
    test_geno_path: Path, tmp_path: Path
) -> None:
    """Config validation (ValueError) is a user error, not a crash.

    The parser accepts --callrate 1.5; the range check lives in the step config.
    """
    res = _run_cli(
        [
            "--pfile", str(test_geno_path),
            "--out", str(tmp_path / "badval"),
            "--callrate", "1.5",
        ]
    )

    assert res.returncode == 1
    assert "Traceback" not in res.stderr, res.stderr
    assert "mind must be between" in res.stderr


def test_unsupported_inference_flags_report_cleanly(
    test_geno_path: Path, tmp_path: Path
) -> None:
    """--container/--singularity/--cloud must fail loudly, not run locally.

    No execution path reads these, so accepting one would silently give the
    user local prediction while they believed it ran in a container or on
    cloud. They are also raised inside ``parse_args`` rather than during the
    run, so this doubles as the guard that config-validation errors raised at
    parse time still get the one-line ``ERROR:`` treatment.
    """
    for flag in ("--container", "--singularity", "--cloud"):
        res = _run_cli(
            [
                "--pfile", str(test_geno_path),
                "--out", str(tmp_path / f"inert{flag.strip('-')}"),
                "--ancestry",
                flag,
            ]
        )

        assert res.returncode == 1, f"{flag} must exit non-zero, got {res.returncode}"
        assert "Traceback" not in res.stderr, f"{flag}:\n{res.stderr}"
        assert f"{flag} is not supported" in res.stderr, f"{flag}:\n{res.stderr}"


def test_debug_flag_restores_the_traceback(
    test_geno_path: Path, tmp_path: Path
) -> None:
    """--debug re-raises, so the traceback is always one flag away."""
    res = _run_cli(
        [
            "--pfile", str(test_geno_path),
            "--out", str(tmp_path / "badval_debug"),
            "--callrate", "1.5",
            "--debug",
        ]
    )

    assert res.returncode != 0
    assert "Traceback" in res.stderr, "--debug should surface the full traceback"
    assert "mind must be between" in res.stderr


def test_skip_fails_rerun_preserves_prior_log(
    test_geno_path: Path, tmp_path: Path
) -> None:
    """--skip-fails bypasses the guard, so the prior log is rotated, not lost."""
    out = tmp_path / "rerun"
    first = _run_cli(["--pfile", str(test_geno_path), "--out", str(out), "--callrate"])
    assert first.returncode == 0, f"first run failed:\n{first.stderr}"

    original = Path(f"{out}_all_logs.log").read_text()
    assert "===== callrate =====" in original

    second = _run_cli(
        [
            "--pfile", str(test_geno_path),
            "--out", str(out),
            "--geno",
            "--skip-fails",
        ]
    )
    assert second.returncode == 0, f"re-run failed:\n{second.stderr}"

    rotated = Path(f"{out}_all_logs.log.1")
    assert rotated.exists(), "prior run's consolidated log was destroyed"
    assert rotated.read_text() == original

    # The new log is this run's, and it says where the old one went.
    current = Path(f"{out}_all_logs.log").read_text()
    assert "===== geno =====" in current
    assert "rerun_all_logs.log.1" in current
    assert "rerun_all_logs.log.1" in second.stderr, (
        "the rotation should be announced on the console too"
    )


def test_preflight_failure_restores_the_rotated_log(
    test_geno_path: Path, tmp_path: Path
) -> None:
    """A run that dies in pre-flight leaves the directory as it found it.

    Setting up logging rotates any existing log to ``.1``. If pre-flight then
    fails, the fresh log is removed -- so without the matching
    ``RunLog.restore_rotated()`` the prior run's log is left stranded at ``.1``
    and the expected path is simply gone.

    This drives the real ``_teardown_logging(remove_logs=True)`` call site.
    The unit test in ``test_core.py`` covers ``restore_rotated`` itself but
    performs the removal by hand, so it passes even if the runner never calls
    it (REFACTOR item 17).
    """
    out = tmp_path / "preflight"

    first = _run_cli(["--pfile", str(test_geno_path), "--out", str(out), "--callrate"])
    assert first.returncode == 0, f"first run failed:\n{first.stderr}"

    log = Path(f"{out}_all_logs.log")
    original = log.read_text()

    # --skip-fails gets past the re-run guard; the missing input then fails
    # pre-flight, after the log has already been rotated aside.
    second = _run_cli([
        "--pfile", "/nonexistent/nope",
        "--out", str(out),
        "--callrate",
        "--skip-fails",
    ])
    assert second.returncode == 1, second.stderr
    assert "does not exist" in second.stderr

    assert log.exists(), "the prior run's log was removed and never restored"
    assert log.read_text() == original, "the restored log is not the prior one"
    assert not Path(f"{out}_all_logs.log.1").exists(), (
        "the rotated copy was left stranded beside the outputs"
    )


def test_no_stray_tool_log_beside_outputs(
    test_geno_path: Path, tmp_path: Path
) -> None:
    """PLINK's own {out}.log is harvested and removed, not left as a duplicate."""
    out = tmp_path / "nostray"
    res = _run_cli(
        ["--pfile", str(test_geno_path), "--out", str(out), "--callrate", "--geno"]
    )
    assert res.returncode == 0, f"Pipeline failed:\n{res.stderr}"

    assert not Path(f"{out}.log").exists(), (
        "PLINK's stray {out}.log should be harvested into the run log and removed"
    )
    # Its content is still on disk, under the per-step and consolidated names.
    assert "PLINK v" in Path(f"{out}_geno.log").read_text()
    assert "PLINK v" in Path(f"{out}_all_logs.log").read_text()
