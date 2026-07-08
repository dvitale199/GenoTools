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

"""Regression tests for CLI runner bugs found in the refactor hardening audit."""

from pathlib import Path
from typing import Any, Dict, List

import pytest

from genotools.core.exceptions import QCError
from genotools.cli.parser import InputArgs, OutputArgs, PipelineArgs
from genotools.cli.runner import PipelineRunner, PipelineState


def _touch_pfiles(prefix: Path) -> None:
    """Create empty placeholder pfiles so input-existence checks pass."""
    for ext in (".pgen", ".pvar", ".psam"):
        Path(f"{prefix}{ext}").write_text("")


def _make_runner(tmp_path: Path, warn_only: bool) -> tuple[PipelineRunner, Path, Path]:
    """Build a PipelineRunner over touched input pfiles in full-output mode."""
    geno = tmp_path / "geno"
    out = tmp_path / "out"
    _touch_pfiles(geno)

    args = PipelineArgs(
        input=InputArgs(pfile=geno),
        output=OutputArgs(out_path=out, warn_only=warn_only, full_output=True),
    )
    runner = PipelineRunner(args)
    runner.state = PipelineState(geno_path=geno, out_path=out, tmp_dir=None)
    return runner, geno, out


def _fake_step_factory(failing_step: str):
    """Return a _run_single_step replacement that raises for one step.

    Passing steps create their output pfiles so the next step's input-existence
    check succeeds; the failing step raises QCError like a real failed step.
    """
    def fake_single_step(
        step: str, step_input: str, step_output: str, legacy_args: Dict[str, Any]
    ) -> Dict[str, Any]:
        if step == failing_step:
            raise QCError(f"{step} failed")
        _touch_pfiles(Path(step_output))
        return {"pass": True, "step": step, "metrics": {}, "output": {}}

    return fake_single_step


class TestWarnModeContinuesAfterRaisingStep:
    """Regression: refactored steps RAISE on failure, but _run_qc_pipeline had no
    try/except, so under --warn a raising step aborted the entire run instead of
    being recorded as failed and skipped (a --warn regression vs the legacy pipeline).
    """

    def test_warn_only_records_failure_and_continues(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        runner, geno, out = _make_runner(tmp_path, warn_only=True)
        monkeypatch.setattr(runner, "_run_single_step", _fake_step_factory("sex"))

        out_dict = runner._run_qc_pipeline(
            steps=["callrate", "sex", "het"],
            geno_path=str(geno),
            out_path=str(out),
        )

        pf = out_dict["pass_fail"]
        assert pf["callrate"]["status"] is True
        assert pf["sex"]["status"] is False          # failure recorded, not raised
        assert pf["het"]["status"] is True           # pipeline continued past it

    def test_no_warn_still_fails_fast(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """Without --warn, a raising step must propagate (fail-fast), not be swallowed."""
        runner, geno, out = _make_runner(tmp_path, warn_only=False)
        monkeypatch.setattr(runner, "_run_single_step", _fake_step_factory("sex"))

        with pytest.raises(QCError):
            runner._run_qc_pipeline(
                steps=["callrate", "sex", "het"],
                geno_path=str(geno),
                out_path=str(out),
            )
