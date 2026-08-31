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

import logging
import shutil
from pathlib import Path
from typing import Any, Dict, List

import pytest

from genotools.core.exceptions import QCError, ValidationError
from genotools.cli.parser import (
    AncestryArgs,
    HetSpec,
    InputArgs,
    OutputArgs,
    PipelineArgs,
    SampleQCArgs,
)
from genotools.cli.runner import PipelineRunner, PipelineState

SYNTHETIC = (
    Path(__file__).resolve().parents[2] / "data" / "synthetic" / "genotools_test"
)


def _plink2_available() -> bool:
    if shutil.which("plink2"):
        return True
    return (Path.home() / ".genotools" / "misc" / "executables" / "plink2").exists()


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
        step: str,
        step_input: str,
        step_output: str,
        legacy_args: Dict[str, Any],
        **_: Any,
    ) -> Dict[str, Any]:
        if step == failing_step:
            raise QCError(f"{step} failed")
        _touch_pfiles(Path(step_output))
        return {"pass": True, "step": step, "metrics": {}, "output": {}}

    return fake_single_step


def _write_traceable_pfiles(prefix: Path) -> None:
    """Create pfiles whose contents encode their own prefix.

    This lets a test tell *which* step's output ended up promoted to the final
    output path (the bytes are copied verbatim by shutil.copy2).
    """
    for ext in (".pgen", ".pvar", ".psam"):
        Path(f"{prefix}{ext}").write_text(f"{prefix.name}{ext}")


def _traceable_step_factory(failing_step: str):
    """Like _fake_step_factory but passing steps write traceable pfiles."""
    def fake_single_step(
        step: str,
        step_input: str,
        step_output: str,
        legacy_args: Dict[str, Any],
        **_: Any,
    ) -> Dict[str, Any]:
        if step == failing_step:
            raise QCError(f"{step} failed")
        _write_traceable_pfiles(Path(step_output))
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


class TestWarnModeFinalStepFailurePromotesLastPassedOutput:
    """Regression: when the LAST step raises under --warn, the pipeline must still
    produce final output pfiles at {out}.pgen/.pvar/.psam by promoting the output
    of the last step that passed (handled by _handle_final_step_failure). Only a
    *middle* raising step was previously covered; a terminal failure left no final
    output at all.
    """

    def test_terminal_failure_promotes_last_passed_output(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        runner, geno, out = _make_runner(tmp_path, warn_only=True)
        # het is the LAST step and raises; callrate + sex pass.
        monkeypatch.setattr(
            runner, "_run_single_step", _traceable_step_factory("het")
        )

        out_dict = runner._run_qc_pipeline(
            steps=["callrate", "sex", "het"],
            geno_path=str(geno),
            out_path=str(out),
        )

        pf = out_dict["pass_fail"]
        assert pf["callrate"]["status"] is True
        assert pf["sex"]["status"] is True
        assert pf["het"]["status"] is False  # terminal failure recorded, not raised

        # Final output must exist despite the terminal failure...
        for ext in (".pgen", ".pvar", ".psam"):
            final = Path(f"{out}{ext}")
            assert final.exists(), f"final output {final.name} was not produced"

        # ...and it must be the LAST PASSED step's output (sex), not the raw input
        # or the callrate output. sex ran on the callrate output, so its prefix is
        # "{out}_callrate_sex".
        expected_prefix = f"{out}_callrate_sex"
        assert Path(f"{out}.pgen").read_text() == f"{Path(expected_prefix).name}.pgen"

    def test_terminal_failure_without_warn_raises(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """Sanity check: without --warn a terminal failure still fails fast."""
        runner, geno, out = _make_runner(tmp_path, warn_only=False)
        monkeypatch.setattr(
            runner, "_run_single_step", _traceable_step_factory("het")
        )

        with pytest.raises(QCError):
            runner._run_qc_pipeline(
                steps=["callrate", "sex", "het"],
                geno_path=str(geno),
                out_path=str(out),
            )


class _CapturingAssoc:
    """Stand-in for run_association that records the AssocConfig it received."""

    def __init__(self) -> None:
        self.config: Any = None

    def __call__(self, data: Any, out_path: Path, config: Any) -> Any:
        self.config = config

        class _Result:
            def to_dict(self_inner) -> Dict[str, Any]:
                return {"pass": True, "step": "assoc", "metrics": {}, "output": {}}

        return _Result()


def _gwas_legacy_args(**overrides: Any) -> Dict[str, Any]:
    """Legacy args dict with the GWAS-relevant keys (mirrors to_legacy_dict)."""
    base: Dict[str, Any] = {
        "pca": None,
        "build": "hg38",
        "gwas": False,
        "covars": None,
        "covar_names": None,
        "maf_lambdas": False,
    }
    base.update(overrides)
    return base


class TestAssocStepThreadsPcaAndCovariateArgs:
    """Regression: the assoc branch of _run_single_step threaded only `build` and
    `maf_lambdas` into the config. It silently dropped the requested PC count
    (`pca`), so PCA always ran the default 10 PCs, and dropped external covariates
    (`covars`/`covar_names`), so `--covars` was ignored and GWAS fell back to PCA
    covariates. These assert the args reach AssocConfig.
    """

    def test_n_pcs_reaches_pca_config(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        runner, geno, out = _make_runner(tmp_path, warn_only=True)
        runner._initialize_modules()
        capture = _CapturingAssoc()
        monkeypatch.setitem(runner._new_modules, "run_gwas_association", capture)

        runner._run_single_step(
            "assoc", str(geno), f"{out}_assoc", _gwas_legacy_args(pca=3)
        )

        assert capture.config.run_pca is True
        assert capture.config.pca is not None
        assert capture.config.pca.n_pcs == 3  # was silently 10

    def test_external_covariates_reach_config(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        runner, geno, out = _make_runner(tmp_path, warn_only=True)
        runner._initialize_modules()
        capture = _CapturingAssoc()
        monkeypatch.setitem(runner._new_modules, "run_gwas_association", capture)

        covar_file = tmp_path / "external.cov"
        covar_file.write_text("#FID\tIID\tAGE\tSEX\n")

        runner._run_single_step(
            "assoc",
            str(geno),
            f"{out}_assoc",
            _gwas_legacy_args(gwas=True, covars=str(covar_file), covar_names="AGE,SEX"),
        )

        assert capture.config.covariates.covar_path == str(covar_file)
        assert capture.config.covariates.covar_names == "AGE,SEX"
        assert capture.config.covariates.has_external_covariates() is True

    def test_build_still_threaded(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """Guard the pre-existing `build` threading isn't lost in the fix."""
        runner, geno, out = _make_runner(tmp_path, warn_only=True)
        runner._initialize_modules()
        capture = _CapturingAssoc()
        monkeypatch.setitem(runner._new_modules, "run_gwas_association", capture)

        runner._run_single_step(
            "assoc", str(geno), f"{out}_assoc", _gwas_legacy_args(pca=10, build="hg19")
        )

        assert capture.config.pca.build == "hg19"


def _count_pc_columns(eigenvec_path: Path) -> int:
    """Count the PC columns (PC1, PC2, ...) in a PLINK2 .eigenvec header."""
    with open(eigenvec_path) as f:
        header = f.readline().strip().split("\t")
    return sum(1 for col in header if col.upper().startswith("PC"))


class TestAssocStepPcaProducesRequestedPcs:
    """End-to-end regression: --pca 3 must produce a 3-PC eigenvec. Before the fix
    the runner ignored the requested count and PLINK2 always wrote 10 PCs.
    """

    @pytest.mark.skipif(not _plink2_available(), reason="plink2 not available")
    def test_pca_3_produces_3_pc_eigenvec(self, tmp_path: Path) -> None:
        if not SYNTHETIC.with_suffix(".pgen").exists():
            pytest.skip("synthetic test data not found")

        local = tmp_path / "cohort"
        for ext in (".pgen", ".pvar", ".psam"):
            shutil.copy2(SYNTHETIC.with_suffix(ext), local.with_suffix(ext))

        args = PipelineArgs(
            input=InputArgs(pfile=local),
            output=OutputArgs(out_path=tmp_path / "out", full_output=True),
        )
        runner = PipelineRunner(args)
        runner._initialize_modules()

        out_assoc = tmp_path / "out_assoc"
        result = runner._run_single_step(
            "assoc", str(local), str(out_assoc), _gwas_legacy_args(pca=3)
        )

        assert result is not None and result["pass"] is True
        eigenvec = out_assoc.with_suffix(".eigenvec")
        assert eigenvec.exists(), "PCA did not write an eigenvec file"
        assert _count_pc_columns(eigenvec) == 3  # was 10


from genotools.core.validation import ValidationDecisions


class TestValidationDecisionsApplied:
    """Data-driven step-skip decisions from validate_input must reach the
    runner (reported as skipped, with a reason; filter_controls forced off)."""

    def test_skip_reasons_carry_decided_steps(self, tmp_path: Path) -> None:
        runner, _, _ = _make_runner(tmp_path, warn_only=False)
        runner._validation_decisions = ValidationDecisions(
            skip_reasons={
                "sex": "no X chromosome data is available",
                "case_control": "only cases or controls are available, not both",
                "het": "12 samples is fewer than the 50 PLINK requires",
            }
        )
        reasons = runner._skip_reasons(
            ["callrate", "sex", "het", "case_control", "hwe"]
        )
        assert set(reasons) == {"sex", "het", "case_control"}
        assert reasons["sex"] == "no X chromosome data is available"

    def test_skip_reasons_ignore_unrequested_steps(self, tmp_path: Path) -> None:
        """A decision about a step this run isn't doing must not leak in."""
        runner, _, _ = _make_runner(tmp_path, warn_only=False)
        runner._validation_decisions = ValidationDecisions(
            skip_reasons={"sex": "no X chromosome data is available"}
        )
        assert runner._skip_reasons(["callrate", "hwe"]) == {}

    def test_skip_reasons_empty_by_default(self, tmp_path: Path) -> None:
        runner, _, _ = _make_runner(tmp_path, warn_only=False)
        assert runner._skip_reasons(["callrate", "sex", "het"]) == {}

    def test_disable_filter_controls_reaches_legacy_args(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        runner, geno, out = _make_runner(tmp_path, warn_only=False)
        # Args request filter_controls=True; the decision must force it False,
        # so the assertion is non-trivial (default is already False).
        runner.args.variant_qc.filter_controls = True
        runner._validation_decisions = ValidationDecisions(disable_filter_controls=True)
        captured: Dict[str, Any] = {}

        def fake_single_step(step, step_input, step_output, legacy_args, **_):
            captured["filter_controls"] = legacy_args["filter_controls"]
            _touch_pfiles(Path(step_output))
            return {"pass": True, "step": step, "metrics": {}, "output": {}}

        monkeypatch.setattr(runner, "_run_single_step", fake_single_step)
        runner._run_qc_pipeline(steps=["hwe"], geno_path=str(geno), out_path=str(out))
        assert captured["filter_controls"] is False


class TestSetupLoggingArtifacts:
    """Round 7: _setup_logging installs a RunLog and no longer creates the dead
    cleaned_logs.log."""

    def test_installs_runlog_and_no_cleaned_logs(self, tmp_path: Path) -> None:
        from genotools.core.logging import raw_sink

        runner, geno, out = _make_runner(tmp_path, warn_only=False)
        try:
            runner._setup_logging()
            assert runner._runlog is not None
            assert Path(f"{out}_all_logs.log").exists()
            # The banner header is written at install time.
            assert Path(f"{out}_all_logs.log").read_text().strip() != ""
            # The dead 0-byte cleaned_logs.log is gone for good.
            assert not Path(f"{out}_cleaned_logs.log").exists()
        finally:
            if runner._runlog is not None:
                runner._runlog.close()
            raw_sink.set(None)


class TestRerunGuardRunsBeforeLogging:
    """Round 7: the re-run guard is a standalone pre-check that runs BEFORE
    logging setup truncates/creates the consolidated log."""

    def test_run_raises_and_preserves_existing_log(self, tmp_path: Path) -> None:
        from genotools.core.exceptions import ValidationError

        geno = tmp_path / "geno"
        out = tmp_path / "out"
        _touch_pfiles(geno)
        # A prior run's consolidated log exists.
        Path(f"{out}_all_logs.log").write_text("PRIOR RUN MARKER\n")

        args = PipelineArgs(
            input=InputArgs(pfile=geno),
            output=OutputArgs(out_path=out, warn_only=False, full_output=True),
        )
        # Enable a step so run() has work; the guard should fire first regardless.
        args.sample_qc.run_callrate = True
        runner = PipelineRunner(args)

        with pytest.raises(ValidationError):
            runner.run()

        # Guard ran before install_run_logging => the existing log was NOT
        # truncated/overwritten with a fresh banner.
        assert "PRIOR RUN MARKER" in Path(f"{out}_all_logs.log").read_text()


class TestNoPrintInRunnerSource:
    """Round 7: the four runtime print() calls are migrated to logging."""

    def test_no_print_calls_left(self) -> None:
        import genotools.cli.runner as runner_mod

        source = Path(runner_mod.__file__).read_text()
        assert "print(" not in source, "runner.py must route all output through logging"


class TestRunSummary:
    """Round 7: end-of-run summary table reaches the console and consolidated log."""

    def test_summary_written_to_log_and_records(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        import logging as _logging
        from genotools.core.logging import raw_sink

        runner, geno, out = _make_runner(tmp_path, warn_only=False)
        runner._setup_logging()

        # Capture the curated console stream directly off the genotools logger
        # (install_run_logging sets propagate=False, so caplog's root handler
        # would miss these records).
        captured: List[_logging.LogRecord] = []

        class _Capture(_logging.Handler):
            def emit(self, record):
                captured.append(record)

        cap = _Capture()
        _logging.getLogger("genotools").addHandler(cap)

        def fake_single_step(step, step_input, step_output, legacy_args, **_):
            _touch_pfiles(Path(step_output))
            return {
                "pass": True, "step": step,
                "metrics": {"outlier_count": 7}, "output": {},
            }

        monkeypatch.setattr(runner, "_run_single_step", fake_single_step)
        try:
            runner._run_qc_pipeline(
                steps=["callrate", "geno"], geno_path=str(geno), out_path=str(out)
            )
            runner._emit_run_summary()
        finally:
            _logging.getLogger("genotools").removeHandler(cap)
            if runner._runlog is not None:
                runner._runlog.close()
            raw_sink.set(None)

        content = Path(f"{out}_all_logs.log").read_text()
        assert "run summary" in content.lower()
        assert "callrate" in content and "geno" in content
        assert "7 removed" in content and "PASS" in content
        # Also surfaced as log records (→ console).
        assert any("Run summary" in r.getMessage() for r in captured)

        # The summary is rendered ONCE on disk: the console rows are marked
        # console_only so they don't duplicate the RunLog's aligned table.
        assert content.count("7 removed") == 2, (  # one row per step, no duplicates
            f"summary rows duplicated in the consolidated log:\n{content}"
        )
        assert "Run summary:" not in content, (
            "console summary header leaked into the consolidated log"
        )

    def test_step_paths_are_logged_to_file_not_console(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """The verbose "Running: <step> with input ..." line is file-only.

        It carries absolute temp-dir paths that are essential when debugging but
        pure noise on the curated console.
        """
        import logging as _logging
        from genotools.core.logging import raw_sink

        runner, geno, out = _make_runner(tmp_path, warn_only=False)
        runner._setup_logging()

        console_messages: List[str] = []

        class _ConsoleCapture(_logging.Handler):
            """Mimics the console handler, including its file_only filter."""

            def emit(self, record):
                if not getattr(record, "file_only", False):
                    console_messages.append(record.getMessage())

        cap = _ConsoleCapture()
        _logging.getLogger("genotools").addHandler(cap)

        def fake_single_step(step, step_input, step_output, legacy_args, **_):
            _touch_pfiles(Path(step_output))
            return {"pass": True, "step": step, "metrics": {}, "output": {}}

        monkeypatch.setattr(runner, "_run_single_step", fake_single_step)
        try:
            runner._run_qc_pipeline(
                steps=["callrate"], geno_path=str(geno), out_path=str(out)
            )
        finally:
            _logging.getLogger("genotools").removeHandler(cap)
            if runner._runlog is not None:
                runner._runlog.close()
            raw_sink.set(None)

        assert not any("Running: callrate" in m for m in console_messages), (
            f"verbose step-path line leaked to the console: {console_messages}"
        )
        # But it is preserved in the consolidated log, with the full paths.
        content = Path(f"{out}_all_logs.log").read_text()
        assert "Running: callrate" in content
        assert str(out) in content


class TestPreflightFailureDoesNotPoisonPrefix:
    """Round 7 fix: logging is set up before validate_input so the breakdown is
    captured, but a pre-flight failure (bad input / no steps) must remove the
    freshly-created log so the output prefix isn't blocked on the next run."""

    def test_validation_failure_removes_log(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        from genotools.core.exceptions import ValidationError
        from genotools.core.logging import raw_sink

        geno = tmp_path / "geno"
        out = tmp_path / "out"
        _touch_pfiles(geno)
        args = PipelineArgs(
            input=InputArgs(pfile=geno),
            output=OutputArgs(out_path=out, full_output=True),
        )
        args.sample_qc.run_callrate = True
        runner = PipelineRunner(args)

        def boom() -> None:
            raise ValidationError("bad input")

        monkeypatch.setattr(runner, "_convert_input_format", boom)

        with pytest.raises(ValidationError):
            runner.run()

        # The consolidated log created just before validation must be cleaned up,
        # so a corrected re-run over the same --out is not blocked by the guard.
        assert not Path(f"{out}_all_logs.log").exists()
        assert runner._runlog is None
        assert raw_sink.get() is None

    def test_no_steps_failure_removes_log(self, tmp_path: Path) -> None:
        from genotools.core.logging import raw_sink

        geno = tmp_path / "geno"
        out = tmp_path / "out"
        _touch_pfiles(geno)
        args = PipelineArgs(
            input=InputArgs(pfile=geno),
            output=OutputArgs(out_path=out, full_output=True),
        )
        runner = PipelineRunner(args)
        # No QC steps and no ancestry requested; _convert_input_format is a no-op
        # for already-pfile input with a touched (empty) psam/pvar, so patch it
        # out to isolate the "no steps" raise path.
        runner._convert_input_format = lambda: None  # type: ignore[method-assign]

        # A GenoToolsError subclass, so main() renders it as a message not a
        # traceback (ValidationError is not a ValueError).
        with pytest.raises(ValidationError, match="No QC steps"):
            runner.run()

        assert not Path(f"{out}_all_logs.log").exists()
        assert raw_sink.get() is None

    def test_preflight_failure_restores_a_rotated_log(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """A failed re-run leaves the directory exactly as it found it.

        With --skip-fails the guard is bypassed, so setup rotates the prior log
        aside; if the run then dies in pre-flight, the fresh log is removed AND
        the rotation is undone -- otherwise the prior log would be orphaned under
        a .1 suffix.
        """
        from genotools.core.logging import raw_sink

        geno = tmp_path / "geno"
        out = tmp_path / "out"
        _touch_pfiles(geno)
        Path(f"{out}_all_logs.log").write_text("PRIOR RUN CONTENT\n")

        args = PipelineArgs(
            input=InputArgs(pfile=geno),
            output=OutputArgs(out_path=out, full_output=True, skip_fails=True),
        )
        runner = PipelineRunner(args)
        runner._convert_input_format = lambda: None  # type: ignore[method-assign]

        with pytest.raises(ValidationError, match="No QC steps"):
            runner.run()

        assert Path(f"{out}_all_logs.log").read_text() == "PRIOR RUN CONTENT\n"
        assert not Path(f"{out}_all_logs.log.1").exists()
        assert raw_sink.get() is None


class TestStepOutcomesReachTheReport:
    """pass/fail/skipped all travel the same path into the report.

    Before 2.0 a step that produced no data was simply absent from the QC
    section, so "not requested" and "requested but impossible" were
    indistinguishable to any consumer of the JSON.
    """

    def test_skipped_step_is_recorded(self, tmp_path: Path, monkeypatch) -> None:
        runner, geno, out = _make_runner(tmp_path, warn_only=False)
        monkeypatch.setattr(runner, "_run_single_step", _fake_step_factory(""))

        out_dict = runner._run_qc_pipeline(
            steps=["callrate", "het", "hwe"],
            geno_path=str(geno),
            out_path=str(out),
            skip_reasons={"het": "12 samples is fewer than the 50 PLINK requires"},
        )

        assert out_dict["het"]["outcome"] == "skipped"
        assert out_dict["het"]["pass"] is False
        assert "12 samples" in out_dict["het"]["reason"]
        assert out_dict["het"]["metrics"] == {"outlier_count": 0}
        assert out_dict["pass_fail"]["het"]["outcome"] == "skipped"
        assert out_dict["pass_fail"]["het"]["status"] is False

    def test_skipped_step_does_not_run(self, tmp_path: Path, monkeypatch) -> None:
        runner, geno, out = _make_runner(tmp_path, warn_only=False)
        ran: list[str] = []

        def recording_step(step, step_input, step_output, legacy_args, **_):
            ran.append(step)
            _touch_pfiles(Path(step_output))
            return {"pass": True, "step": step, "metrics": {}, "output": {}}

        monkeypatch.setattr(runner, "_run_single_step", recording_step)
        runner._run_qc_pipeline(
            steps=["callrate", "het", "hwe"],
            geno_path=str(geno),
            out_path=str(out),
            skip_reasons={"het": "too few samples"},
        )
        assert ran == ["callrate", "hwe"]

    def test_skip_passes_the_pfile_chain_through(
        self, tmp_path: Path, monkeypatch
    ) -> None:
        """The step after a skip must read what the skipped step would have."""
        runner, geno, out = _make_runner(tmp_path, warn_only=False)
        inputs: dict[str, str] = {}

        def recording_step(step, step_input, step_output, legacy_args, **_):
            inputs[step] = step_input
            _touch_pfiles(Path(step_output))
            return {"pass": True, "step": step, "metrics": {}, "output": {}}

        monkeypatch.setattr(runner, "_run_single_step", recording_step)
        out_dict = runner._run_qc_pipeline(
            steps=["callrate", "het", "hwe"],
            geno_path=str(geno),
            out_path=str(out),
            skip_reasons={"het": "too few samples"},
        )

        # hwe reads callrate's output, not a path het never wrote.
        assert inputs["hwe"] == out_dict["pass_fail"]["callrate"]["output"]
        assert Path(f"{inputs['hwe']}.pgen").exists()

    def test_failed_step_is_recorded(self, tmp_path: Path, monkeypatch) -> None:
        """Regression: the refactor dropped failed steps from the QC section."""
        runner, geno, out = _make_runner(tmp_path, warn_only=True)
        monkeypatch.setattr(runner, "_run_single_step", _fake_step_factory("het"))

        out_dict = runner._run_qc_pipeline(
            steps=["callrate", "het", "hwe"],
            geno_path=str(geno),
            out_path=str(out),
        )

        assert out_dict["het"]["outcome"] == "fail"
        assert out_dict["het"]["pass"] is False
        assert "het failed" in out_dict["het"]["reason"]
        assert out_dict["het"]["metrics"] == {"outlier_count": 0}

    def test_failure_reason_reaches_the_report(
        self, tmp_path: Path, monkeypatch
    ) -> None:
        """PassFailRecord captured the error, then dropped it on serialization."""
        runner, geno, out = _make_runner(tmp_path, warn_only=True)
        monkeypatch.setattr(runner, "_run_single_step", _fake_step_factory("het"))

        out_dict = runner._run_qc_pipeline(
            steps=["callrate", "het", "hwe"],
            geno_path=str(geno),
            out_path=str(out),
        )
        assert "het failed" in out_dict["pass_fail"]["het"]["reason"]

    def test_passing_steps_keep_their_outcome(
        self, tmp_path: Path, monkeypatch
    ) -> None:
        runner, geno, out = _make_runner(tmp_path, warn_only=False)
        monkeypatch.setattr(runner, "_run_single_step", _fake_step_factory(""))

        out_dict = runner._run_qc_pipeline(
            steps=["callrate", "hwe"],
            geno_path=str(geno),
            out_path=str(out),
        )
        assert out_dict["pass_fail"]["callrate"]["outcome"] == "pass"
        assert out_dict["pass_fail"]["callrate"]["reason"] is None


def _psam_with(prefix: Path, n_samples: int) -> None:
    """Write a .psam with n_samples rows (plus header) at prefix."""
    lines = ["#FID\tIID\tSEX\tPHENO1"]
    lines += [f"F{i}\tI{i}\t1\t1" for i in range(n_samples)]
    Path(f"{prefix}.psam").write_text("\n".join(lines) + "\n")


class TestPerGroupHetFloor:
    """The cohort-level het guard runs before the ancestry split, so a small
    group can still fall under PLINK's LD floor. Decide against the group file.
    """

    def test_small_group_skips_het(self, tmp_path: Path) -> None:
        runner, _, _ = _make_runner(tmp_path, warn_only=False)
        group = tmp_path / "out_ancestry_FIN"
        _psam_with(group, 12)

        reasons = runner._skip_reasons(["callrate", "het", "geno"], str(group))
        assert set(reasons) == {"het"}
        assert "12 samples" in reasons["het"]

    def test_large_group_runs_het(self, tmp_path: Path) -> None:
        runner, _, _ = _make_runner(tmp_path, warn_only=False)
        group = tmp_path / "out_ancestry_EUR"
        _psam_with(group, 6927)

        assert runner._skip_reasons(["callrate", "het", "geno"], str(group)) == {}

    def test_floor_is_inclusive(self, tmp_path: Path) -> None:
        """MIN_HET_SAMPLES itself is enough - PLINK's error is 'less than 50'."""
        from genotools.core.validation import MIN_HET_SAMPLES

        runner, _, _ = _make_runner(tmp_path, warn_only=False)
        at_floor = tmp_path / "at_floor"
        _psam_with(at_floor, MIN_HET_SAMPLES)
        below = tmp_path / "below"
        _psam_with(below, MIN_HET_SAMPLES - 1)

        assert runner._skip_reasons(["het"], str(at_floor)) == {}
        assert "het" in runner._skip_reasons(["het"], str(below))

    def test_cohort_decision_is_not_overridden_by_group_size(
        self, tmp_path: Path
    ) -> None:
        """A large group must not resurrect a het skip decided cohort-wide."""
        runner, _, _ = _make_runner(tmp_path, warn_only=False)
        runner._validation_decisions = ValidationDecisions(
            skip_reasons={"het": "cohort-level reason"}
        )
        group = tmp_path / "big"
        _psam_with(group, 5000)

        reasons = runner._skip_reasons(["het"], str(group))
        assert reasons["het"] == "cohort-level reason"

    def test_het_not_requested_is_not_a_skip(self, tmp_path: Path) -> None:
        """A step you did not ask for stays absent, not reported as skipped."""
        runner, _, _ = _make_runner(tmp_path, warn_only=False)
        group = tmp_path / "tiny"
        _psam_with(group, 3)

        assert runner._skip_reasons(["callrate", "geno"], str(group)) == {}


class TestSkipAtTheEndOfTheChain:
    """A skipped step writes no files, so it must not claim the final output
    path - the last step that actually runs does.
    """

    @pytest.mark.parametrize("warn_only", [True, False])
    def test_trailing_skip_still_produces_output(
        self, tmp_path: Path, monkeypatch, warn_only: bool
    ) -> None:
        runner, geno, out = _make_runner(tmp_path, warn_only=warn_only)
        monkeypatch.setattr(runner, "_run_single_step", _fake_step_factory(""))

        runner._run_qc_pipeline(
            steps=["callrate", "het"],
            geno_path=str(geno),
            out_path=str(out),
            skip_reasons={"het": "too few samples"},
        )

        assert Path(f"{out}.pgen").exists(), "trailing skip swallowed the output"

    def test_last_runnable_step_writes_the_output_path(
        self, tmp_path: Path, monkeypatch
    ) -> None:
        runner, geno, out = _make_runner(tmp_path, warn_only=False)
        outputs: dict[str, str] = {}

        def recording_step(step, step_input, step_output, legacy_args, **_):
            outputs[step] = step_output
            _touch_pfiles(Path(step_output))
            return {"pass": True, "step": step, "metrics": {}, "output": {}}

        monkeypatch.setattr(runner, "_run_single_step", recording_step)
        runner._run_qc_pipeline(
            steps=["callrate", "geno", "het"],
            geno_path=str(geno),
            out_path=str(out),
            skip_reasons={"het": "too few samples"},
        )
        assert outputs["geno"] == str(out)

    def test_all_steps_skipped_passes_input_through(
        self, tmp_path: Path, monkeypatch
    ) -> None:
        """A group where nothing can run still needs output pfiles."""
        runner, geno, out = _make_runner(tmp_path, warn_only=True)
        monkeypatch.setattr(runner, "_run_single_step", _fake_step_factory(""))

        out_dict = runner._run_qc_pipeline(
            steps=["het"],
            geno_path=str(geno),
            out_path=str(out),
            skip_reasons={"het": "too few samples"},
        )

        assert Path(f"{out}.pgen").exists()
        assert out_dict["het"]["outcome"] == "skipped"


class TestNothingCanRunMessage:
    """"No QC steps requested" is wrong when a step was requested and ruled
    out. Say which, and why.
    """

    def test_message_names_the_ruled_out_step(self, tmp_path: Path) -> None:
        runner, _, _ = _make_runner(tmp_path, warn_only=False)
        runner._validation_decisions = ValidationDecisions(
            skip_reasons={"het": "12 samples is fewer than the 50 PLINK requires"}
        )
        skips = runner._skip_reasons(["het"])

        message = str(runner._nothing_to_run_error(skips))
        assert "het" in message
        assert "12 samples" in message
        assert "requested" not in message

    def test_message_unchanged_when_nothing_was_requested(
        self, tmp_path: Path
    ) -> None:
        runner, _, _ = _make_runner(tmp_path, warn_only=False)
        message = str(runner._nothing_to_run_error({}))
        assert message == "No QC steps or ancestry prediction requested"


class _StubTmpDir:
    """Stands in for tempfile.TemporaryDirectory - only ``.name`` is read."""

    def __init__(self, path: Path) -> None:
        self.name = str(path)


def _ancestry_runner(
    tmp_path: Path, full_output: bool, run_ancestry: bool = True
) -> tuple[PipelineRunner, Path, Path]:
    """Build a runner in ancestry mode, with a stand-in temp working dir."""
    geno = tmp_path / "geno"
    out = tmp_path / "out"
    _touch_pfiles(geno)
    work = tmp_path / "tmpwork"
    work.mkdir(exist_ok=True)

    args = PipelineArgs(
        input=InputArgs(pfile=geno),
        output=OutputArgs(out_path=out, full_output=full_output),
        ancestry=AncestryArgs(run_ancestry=run_ancestry),
    )
    runner = PipelineRunner(args)
    runner.state = PipelineState(
        geno_path=geno, out_path=out, tmp_dir=_StubTmpDir(work)
    )
    return runner, out, work


class TestAncestryGroupPathsWithoutFullOutput:
    """Without --full-output the ancestry split writes into the temp working
    directory, so ``{out}_ancestry_{label}`` never exists on disk. Anything
    that reads a group's files has to resolve that, or the run dies before its
    first QC step - which is exactly what happened for every ancestry run.
    """

    def test_group_prefix_resolves_into_tmp_dir(self, tmp_path: Path) -> None:
        runner, out, work = _ancestry_runner(tmp_path, full_output=False)

        resolved = runner._resolve_existing_geno(f"{out}_ancestry_FIN")

        assert resolved == f"{work}/out_ancestry_FIN"

    def test_full_output_keeps_the_nominal_prefix(self, tmp_path: Path) -> None:
        runner, out, _ = _ancestry_runner(tmp_path, full_output=True)

        nominal = f"{out}_ancestry_FIN"
        assert runner._resolve_existing_geno(nominal) == nominal

    def test_flat_run_keeps_the_nominal_prefix(self, tmp_path: Path) -> None:
        """A run without --ancestry reads the user's own input, not a copy."""
        runner, _, _ = _ancestry_runner(
            tmp_path, full_output=False, run_ancestry=False
        )

        assert runner._resolve_existing_geno("/data/cohort") == "/data/cohort"

    def test_skip_reasons_reads_the_group_that_exists(self, tmp_path: Path) -> None:
        """The regression: the group's .psam lives only in the temp dir, so
        deciding skips against the nominal prefix raises FileNotFoundError.
        """
        runner, out, work = _ancestry_runner(tmp_path, full_output=False)
        # Only the temp-dir copy exists, exactly as after a real split.
        _psam_with(work / "out_ancestry_FIN", 12)
        assert not Path(f"{out}_ancestry_FIN.psam").exists()

        resolved = runner._resolve_existing_geno(f"{out}_ancestry_FIN")
        reasons = runner._skip_reasons(["callrate", "het", "geno"], resolved)

        assert set(reasons) == {"het"}
        assert "12 samples" in reasons["het"]

    def test_unresolved_prefix_is_what_used_to_break(self, tmp_path: Path) -> None:
        """Pin the failure mode, so a future refactor that drops the
        resolution fails here instead of at runtime on a real cohort.
        """
        runner, out, work = _ancestry_runner(tmp_path, full_output=False)
        _psam_with(work / "out_ancestry_FIN", 12)

        with pytest.raises(FileNotFoundError):
            runner._skip_reasons(["het"], f"{out}_ancestry_FIN")

    def test_ancestry_loop_decides_skips_against_the_real_group_files(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """The regression, at the call site: drive the real per-group loop with
        the group files only where a split would put them. Reading the nominal
        prefix here raised FileNotFoundError and killed every ancestry run
        before its first QC step.
        """
        runner, out, work = _ancestry_runner(tmp_path, full_output=False)
        _psam_with(work / "out_ancestry_EUR", 6927)
        _psam_with(work / "out_ancestry_FIN", 12)
        assert not Path(f"{out}_ancestry_FIN.psam").exists()

        monkeypatch.setattr(runner, "_begin_section", lambda *a, **k: None)
        monkeypatch.setattr(runner, "_end_section", lambda *a, **k: None)
        monkeypatch.setattr(
            runner,
            "_run_ancestry_prediction_new",
            lambda: {"data": {"labels_list": ["EUR", "FIN"]}},
        )
        monkeypatch.setattr(runner, "_build_output", lambda: None)

        seen: Dict[str, Dict[str, str]] = {}

        def fake_qc(**kwargs: Any) -> None:
            seen[kwargs["ancestry_label"]] = kwargs["skip_reasons"]

        monkeypatch.setattr(runner, "_run_qc_pipeline", fake_qc)

        runner._run_with_ancestry(["callrate", "het", "geno"])

        assert set(seen) == {"EUR", "FIN"}
        assert seen["EUR"] == {}, "6927 samples clears the het floor"
        assert set(seen["FIN"]) == {"het"}
        assert "12 samples" in seen["FIN"]["het"]


def _psam_with_phenos(prefix: Path, phenos: List[int]) -> None:
    """Write a .psam at prefix with one row per entry in phenos."""
    lines = ["#FID\tIID\tSEX\tPHENO1"]
    lines += [f"F{i}\tI{i}\t1\t{p}" for i, p in enumerate(phenos)]
    Path(f"{prefix}.psam").write_text("\n".join(lines) + "\n")


class TestPerGroupCaseControl:
    """--test-missing needs both cases and controls. The cohort can hold both
    while an ancestry group holds only one, so the decision has to be re-made
    per group - otherwise the step raises and the same data-driven decision is
    reported as outcome="fail" here but "skipped" cohort-wide.
    """

    def test_group_with_only_cases_skips_case_control(self, tmp_path: Path) -> None:
        runner, _, _ = _make_runner(tmp_path, warn_only=False)
        group = tmp_path / "out_ancestry_FIN"
        _psam_with_phenos(group, [2] * 80)

        reasons = runner._skip_reasons(["geno", "case_control", "hwe"], str(group))
        assert set(reasons) == {"case_control"}
        assert "not both" in reasons["case_control"]

    def test_group_with_only_controls_skips_case_control(self, tmp_path: Path) -> None:
        runner, _, _ = _make_runner(tmp_path, warn_only=False)
        group = tmp_path / "out_ancestry_AAC"
        _psam_with_phenos(group, [1] * 80)

        assert "case_control" in runner._skip_reasons(["case_control"], str(group))

    def test_group_with_both_runs_case_control(self, tmp_path: Path) -> None:
        runner, _, _ = _make_runner(tmp_path, warn_only=False)
        group = tmp_path / "out_ancestry_EUR"
        _psam_with_phenos(group, [1, 2] * 40)

        assert runner._skip_reasons(["geno", "case_control"], str(group)) == {}

    def test_cohort_decision_is_not_overridden_by_the_group(
        self, tmp_path: Path
    ) -> None:
        """A group holding both must not resurrect a cohort-level skip."""
        runner, _, _ = _make_runner(tmp_path, warn_only=False)
        runner._validation_decisions = ValidationDecisions(
            skip_reasons={"case_control": "cohort-level reason"}
        )
        group = tmp_path / "both"
        _psam_with_phenos(group, [1, 2] * 40)

        reasons = runner._skip_reasons(["case_control"], str(group))
        assert reasons["case_control"] == "cohort-level reason"

    def test_case_control_not_requested_is_not_a_skip(self, tmp_path: Path) -> None:
        runner, _, _ = _make_runner(tmp_path, warn_only=False)
        group = tmp_path / "cases_only"
        _psam_with_phenos(group, [2] * 80)

        assert runner._skip_reasons(["geno", "hwe"], str(group)) == {}

    def test_no_dataset_means_no_per_group_decision(self, tmp_path: Path) -> None:
        """Called without a dataset (the nothing-to-run pre-flight), only the
        cohort decisions apply - there is nothing to read.
        """
        runner, _, _ = _make_runner(tmp_path, warn_only=False)
        assert runner._skip_reasons(["case_control"]) == {}

    def test_skipped_case_control_reports_zeroed_variant_metrics(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """The skip has to land in the report as a variant-step row, not vanish."""
        runner, geno, out = _make_runner(tmp_path, warn_only=False)
        monkeypatch.setattr(runner, "_run_single_step", _fake_step_factory(""))
        _psam_with_phenos(geno, [2] * 80)

        steps = ["geno", "case_control"]
        out_dict = runner._run_qc_pipeline(
            steps=steps,
            geno_path=str(geno),
            out_path=str(out),
            skip_reasons=runner._skip_reasons(steps, str(geno)),
        )

        row = out_dict["case_control"]
        assert row["outcome"] == "skipped"
        assert row["step"] == "case_control_missingness_prune"
        assert row["metrics"] == {"mis_removed_count": 0}
        assert out_dict["pass_fail"]["case_control"]["outcome"] == "skipped"

    def test_ancestry_loop_decides_case_control_per_group(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """The fix at its call site: drive the real per-group loop over a cohort
        that holds both phenotypes but splits into a group that holds one.
        """
        runner, out, work = _ancestry_runner(tmp_path, full_output=False)
        _psam_with_phenos(work / "out_ancestry_EUR", [1, 2] * 40)
        _psam_with_phenos(work / "out_ancestry_FIN", [2] * 80)

        monkeypatch.setattr(runner, "_begin_section", lambda *a, **k: None)
        monkeypatch.setattr(runner, "_end_section", lambda *a, **k: None)
        monkeypatch.setattr(
            runner,
            "_run_ancestry_prediction_new",
            lambda: {"data": {"labels_list": ["EUR", "FIN"]}},
        )
        monkeypatch.setattr(runner, "_build_output", lambda: None)

        seen: Dict[str, Dict[str, str]] = {}

        def fake_qc(**kwargs: Any) -> None:
            seen[kwargs["ancestry_label"]] = kwargs["skip_reasons"]

        monkeypatch.setattr(runner, "_run_qc_pipeline", fake_qc)

        runner._run_with_ancestry(["geno", "case_control"])

        assert seen["EUR"] == {}, "both phenotypes present - case_control runs"
        assert set(seen["FIN"]) == {"case_control"}
        assert "not both" in seen["FIN"]["case_control"]


def _psam_with_sexes(prefix: Path, sexes: List[int]) -> None:
    """Write a .psam at prefix with one row per entry in sexes."""
    lines = ["#FID\tIID\tSEX\tPHENO1"]
    lines += [f"F{i}\tI{i}\t{s}\t{1 + (i % 2)}" for i, s in enumerate(sexes)]
    Path(f"{prefix}.psam").write_text("\n".join(lines) + "\n")


class TestPerGroupSex:
    """--check-sex needs recorded sample sex. The cohort can hold it while an
    ancestry group holds none, so the decision is re-made per group - the same
    shape as the het floor and case-control checks.
    """

    def test_group_without_recorded_sex_skips_sex(self, tmp_path: Path) -> None:
        runner, _, _ = _make_runner(tmp_path, warn_only=False)
        group = tmp_path / "out_ancestry_CAS"
        _psam_with_sexes(group, [0] * 80)

        reasons = runner._skip_reasons(["callrate", "sex", "geno"], str(group))
        assert set(reasons) == {"sex"}
        assert reasons["sex"] == "no sample sex data is available"

    def test_group_with_one_sex_still_runs(self, tmp_path: Path) -> None:
        """Males only is enough - the cohort rule skips only when neither
        males nor females are recorded, and the group rule must match it.
        """
        runner, _, _ = _make_runner(tmp_path, warn_only=False)
        group = tmp_path / "males_only"
        _psam_with_sexes(group, [1] * 80)

        assert runner._skip_reasons(["sex"], str(group)) == {}

    def test_group_with_both_sexes_runs(self, tmp_path: Path) -> None:
        runner, _, _ = _make_runner(tmp_path, warn_only=False)
        group = tmp_path / "mixed"
        _psam_with_sexes(group, [1, 2] * 40)

        assert runner._skip_reasons(["callrate", "sex"], str(group)) == {}

    def test_cohort_no_x_decision_survives_the_group_check(
        self, tmp_path: Path
    ) -> None:
        """A cohort with no X chromosome rules sex out for every group, and the
        group's own sex column must not overwrite that reason - the pvar is
        shared, so no group can have X data the cohort lacks.
        """
        runner, _, _ = _make_runner(tmp_path, warn_only=False)
        runner._validation_decisions = ValidationDecisions(
            skip_reasons={"sex": "no X chromosome data is available"}
        )
        group = tmp_path / "nox_group"
        _psam_with_sexes(group, [0] * 80)

        reasons = runner._skip_reasons(["sex"], str(group))
        assert reasons["sex"] == "no X chromosome data is available"

    def test_sex_not_requested_is_not_a_skip(self, tmp_path: Path) -> None:
        runner, _, _ = _make_runner(tmp_path, warn_only=False)
        group = tmp_path / "unsexed"
        _psam_with_sexes(group, [0] * 80)

        assert runner._skip_reasons(["callrate", "geno"], str(group)) == {}

    def test_ancestry_loop_decides_sex_per_group(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """The fix at its call site: one group with sex recorded, one without."""
        runner, out, work = _ancestry_runner(tmp_path, full_output=False)
        _psam_with_sexes(work / "out_ancestry_EUR", [1, 2] * 40)
        _psam_with_sexes(work / "out_ancestry_CAS", [0] * 80)

        monkeypatch.setattr(runner, "_begin_section", lambda *a, **k: None)
        monkeypatch.setattr(runner, "_end_section", lambda *a, **k: None)
        monkeypatch.setattr(
            runner,
            "_run_ancestry_prediction_new",
            lambda: {"data": {"labels_list": ["EUR", "CAS"]}},
        )
        monkeypatch.setattr(runner, "_build_output", lambda: None)

        seen: Dict[str, Dict[str, str]] = {}
        monkeypatch.setattr(
            runner,
            "_run_qc_pipeline",
            lambda **kw: seen.update({kw["ancestry_label"]: kw["skip_reasons"]}),
        )

        runner._run_with_ancestry(["callrate", "sex"])

        assert seen["EUR"] == {}
        assert set(seen["CAS"]) == {"sex"}

    def test_all_three_checks_decide_together(self, tmp_path: Path) -> None:
        """A group can fail more than one precondition at once, and each has to
        be reported with its own reason rather than the first one found.
        """
        runner, _, _ = _make_runner(tmp_path, warn_only=False)
        group = tmp_path / "out_ancestry_TINY"
        lines = ["#FID\tIID\tSEX\tPHENO1"]
        lines += [f"F{i}\tI{i}\t0\t2" for i in range(12)]
        Path(f"{group}.psam").write_text("\n".join(lines) + "\n")

        reasons = runner._skip_reasons(
            ["callrate", "sex", "het", "geno", "case_control"], str(group)
        )
        assert set(reasons) == {"sex", "het", "case_control"}
        assert "12 samples" in reasons["het"]
        assert reasons["sex"] == "no sample sex data is available"
        assert "not both" in reasons["case_control"]


class _StubResult:
    """Minimal stand-in for FilterResult, enough for _run_qc_pipeline."""

    def __init__(self, out_path: str) -> None:
        self.out_path = out_path

    def to_dict(self) -> Dict[str, Any]:
        return {
            "pass": True,
            "step": "het_prune",
            "metrics": {},
            "output": {"plink_out": self.out_path},
        }


class TestHetConfigReachesEveryRunShape:
    """Regression: --amr-het was read in exactly one place, inside
    ``_run_with_ancestry``. ``_run_qc_only`` never passed het values at all, so
    the flag was silently inert in a flat run - and the per-ancestry production
    workflow runs ancestry once, then QCs each group as a separate *flat* job.
    That is precisely the path where it did nothing, with no warning.

    The fix routes both paths through ``SampleQCArgs.het_config_for``, so these
    tests assert the resolved config *arrives at the step*, not merely that the
    parser built it.
    """

    def _flat_runner(self, tmp_path: Path, sample_qc: SampleQCArgs):
        geno = tmp_path / "geno"
        out = tmp_path / "out"
        _touch_pfiles(geno)
        # het is auto-skipped below PLINK's 50-sample LD floor, and a skipped
        # step never reaches the config at all.
        _psam_with(geno, 200)
        args = PipelineArgs(
            input=InputArgs(pfile=geno),
            output=OutputArgs(out_path=out, full_output=True),
            sample_qc=sample_qc,
        )
        runner = PipelineRunner(args)
        runner.state = PipelineState(geno_path=geno, out_path=out, tmp_dir=None)
        return runner, geno, out

    @staticmethod
    def _capture_het(runner, monkeypatch) -> Dict[str, Any]:
        """Record the HetConfig filter_heterozygosity is actually called with.

        Stubs at the module boundary rather than at ``_run_single_step`` so the
        assertion covers the real config-building branch, not the plumbing
        around it.
        """
        runner._initialize_modules()
        seen: Dict[str, Any] = {}

        class _StubGenotypeData:
            @staticmethod
            def from_path(path):
                return path

        def fake_filter(data, config, out_path):
            seen["config"] = config
            _touch_pfiles(Path(out_path))
            return _StubResult(str(out_path))

        monkeypatch.setitem(
            runner._new_modules, "GenotypeData", _StubGenotypeData
        )
        monkeypatch.setitem(
            runner._new_modules, "filter_heterozygosity", fake_filter
        )
        return seen

    def test_flat_run_reaches_the_step_with_derived_bounds(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """THE regression. A flat ``--het sd`` run must reach
        filter_heterozygosity with auto_detect set.

        Drives ``_run_qc_only`` - the entry point the bug lived in - rather
        than handing ``_run_qc_pipeline`` a config directly, which would pass
        against the broken code and pin nothing.

        Revert check: drop ``het_config=`` from ``_run_qc_only``'s call and
        this fails with ``auto_detect=False`` and the base fixed bounds,
        exactly as --amr-het did.
        """
        runner, _, _ = self._flat_runner(
            tmp_path, SampleQCArgs(run_het=True, het_auto=True, het_auto_sd=2.0)
        )
        seen = self._capture_het(runner, monkeypatch)
        monkeypatch.setattr(runner, "_build_output", lambda: None)

        runner._run_qc_only(["het"])

        assert seen["config"].auto_detect is True
        assert seen["config"].auto_sd == 2.0

    def test_run_qc_only_passes_a_resolved_config(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """The bug was in _run_qc_only specifically: it called the pipeline
        without any het argument. Pin the call itself."""
        runner, _, _ = self._flat_runner(
            tmp_path, SampleQCArgs(run_het=True, het_auto=True, het_auto_sd=2.5)
        )
        monkeypatch.setattr(runner, "_build_output", lambda: None)

        captured: Dict[str, Any] = {}

        def fake_pipeline(**kwargs: Any) -> None:
            captured.update(kwargs)

        monkeypatch.setattr(runner, "_run_qc_pipeline", fake_pipeline)
        runner._run_qc_only(["het"])

        assert "het_config" in captured, (
            "_run_qc_only must pass a resolved HetConfig; omitting it is the "
            "original bug"
        )
        assert captured["het_config"].auto_detect is True
        assert captured["het_config"].auto_sd == 2.5

    def test_fixed_base_still_reaches_the_step_unchanged(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """The default path must not have moved."""
        runner, _, _ = self._flat_runner(
            tmp_path, SampleQCArgs(run_het=True, het_lower=-0.2, het_upper=0.2)
        )
        seen = self._capture_het(runner, monkeypatch)
        monkeypatch.setattr(runner, "_build_output", lambda: None)

        runner._run_qc_only(["het"])

        assert seen["config"].auto_detect is False
        assert (seen["config"].f_lower, seen["config"].f_upper) == (-0.2, 0.2)

    def test_pipeline_without_het_config_falls_back_to_the_base(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """Any caller that forgets the argument still gets the base spec, not
        the defaults - the old failure mode was reverting to fixed bounds."""
        runner, geno, out = self._flat_runner(
            tmp_path, SampleQCArgs(run_het=True, het_auto=True, het_auto_sd=1.5)
        )
        seen = self._capture_het(runner, monkeypatch)

        runner._run_qc_pipeline(
            steps=["het"], geno_path=str(geno), out_path=str(out)
        )

        assert seen["config"].auto_detect is True
        assert seen["config"].auto_sd == 1.5


class TestPerGroupHetConfig:
    """--het-ancestry resolves per group across the ancestry loop, on any
    label. The retired --amr-het hardcoded ``label == "AMR"``, the only
    user-facing feature that assumed one reference panel's naming - so it could
    never reach CAH, the synthetic admixed label admixture detection invents.
    """

    def _run_loop(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch,
        sample_qc: SampleQCArgs, labels: List[str],
    ) -> Dict[str, Any]:
        geno = tmp_path / "geno"
        out = tmp_path / "out"
        _touch_pfiles(geno)
        work = tmp_path / "tmpwork"
        work.mkdir(exist_ok=True)

        args = PipelineArgs(
            input=InputArgs(pfile=geno),
            output=OutputArgs(out_path=out, full_output=False),
            sample_qc=sample_qc,
            ancestry=AncestryArgs(run_ancestry=True),
        )
        runner = PipelineRunner(args)
        runner.state = PipelineState(
            geno_path=geno, out_path=out, tmp_dir=_StubTmpDir(work)
        )
        for label in labels:
            _psam_with(work / f"out_ancestry_{label}", 200)

        monkeypatch.setattr(runner, "_begin_section", lambda *a, **k: None)
        monkeypatch.setattr(runner, "_end_section", lambda *a, **k: None)
        monkeypatch.setattr(
            runner,
            "_run_ancestry_prediction_new",
            lambda: {"data": {"labels_list": labels}},
        )
        monkeypatch.setattr(runner, "_build_output", lambda: None)

        seen: Dict[str, Any] = {}

        def fake_qc(**kwargs: Any) -> None:
            seen[kwargs["ancestry_label"]] = kwargs["het_config"]

        monkeypatch.setattr(runner, "_run_qc_pipeline", fake_qc)
        runner._run_with_ancestry(["het"])
        return seen

    def test_override_applies_only_to_its_group(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        seen = self._run_loop(
            tmp_path,
            monkeypatch,
            SampleQCArgs(
                run_het=True,
                het_lower=-0.2,
                het_upper=0.2,
                het_by_ancestry={"AMR": HetSpec(auto=True, sd=3.0)},
            ),
            ["AMR", "EUR"],
        )
        assert seen["AMR"].auto_detect is True
        assert seen["EUR"].auto_detect is False
        assert (seen["EUR"].f_lower, seen["EUR"].f_upper) == (-0.2, 0.2)

    def test_mixed_overrides_on_an_adaptive_base(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        seen = self._run_loop(
            tmp_path,
            monkeypatch,
            SampleQCArgs(
                run_het=True,
                het_auto=True,
                het_auto_sd=2.0,
                het_by_ancestry={
                    "AMR": HetSpec(auto=True, sd=1.5),
                    "CAH": HetSpec(lower=-0.3, upper=0.3),
                },
            ),
            ["AMR", "CAH", "EUR"],
        )
        assert seen["AMR"].auto_sd == 1.5
        assert seen["CAH"].auto_detect is False
        assert (seen["CAH"].f_lower, seen["CAH"].f_upper) == (-0.3, 0.3)
        assert seen["EUR"].auto_sd == 2.0, "non-overridden group keeps the base"

    def test_non_amr_label_is_reachable(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """CAH is invented by admixture detection and appears in no reference
        panel, so --amr-het could not name it."""
        seen = self._run_loop(
            tmp_path,
            monkeypatch,
            SampleQCArgs(
                run_het=True, het_by_ancestry={"CAH": HetSpec(auto=True)}
            ),
            ["CAH", "EUR"],
        )
        assert seen["CAH"].auto_detect is True
        assert seen["EUR"].auto_detect is False

    def test_unknown_label_warns_and_names_the_predicted_groups(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch, caplog
    ) -> None:
        """An override matching no predicted group would otherwise do nothing
        silently - --amr-het's failure on any panel not spelling it "AMR"."""
        with caplog.at_level(logging.WARNING, logger="genotools"):
            seen = self._run_loop(
                tmp_path,
                monkeypatch,
                SampleQCArgs(
                    run_het=True, het_by_ancestry={"AMRR": HetSpec(auto=True)}
                ),
                ["AMR", "EUR"],
            )
        assert "AMRR" in caplog.text
        assert "AMR" in caplog.text and "EUR" in caplog.text
        assert seen["AMR"].auto_detect is False, "the typo must not match AMR"

    def test_matching_label_does_not_warn(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch, caplog
    ) -> None:
        with caplog.at_level(logging.WARNING, logger="genotools"):
            self._run_loop(
                tmp_path,
                monkeypatch,
                SampleQCArgs(
                    run_het=True, het_by_ancestry={"AMR": HetSpec(auto=True)}
                ),
                ["AMR", "EUR"],
            )
        assert "no such ancestry group" not in caplog.text

