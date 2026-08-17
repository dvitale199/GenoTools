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

import shutil
from pathlib import Path
from typing import Any, Dict, List

import pytest

from genotools.core.exceptions import QCError, ValidationError
from genotools.cli.parser import InputArgs, OutputArgs, PipelineArgs
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
        step: str, step_input: str, step_output: str, legacy_args: Dict[str, Any]
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
        step: str, step_input: str, step_output: str, legacy_args: Dict[str, Any]
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
    """Round-6: data-driven step-skip decisions from validate_input must be
    applied by the runner (skipped steps dropped; filter_controls forced off)."""

    def test_filter_steps_drops_decided_steps(self, tmp_path: Path) -> None:
        runner, _, _ = _make_runner(tmp_path, warn_only=False)
        runner._validation_decisions = ValidationDecisions(
            skip_sex=True, skip_case_control=True, skip_het=True
        )
        kept = runner._filter_steps_by_decisions(
            ["callrate", "sex", "het", "case_control", "hwe"]
        )
        assert kept == ["callrate", "hwe"]

    def test_filter_steps_noop_by_default(self, tmp_path: Path) -> None:
        runner, _, _ = _make_runner(tmp_path, warn_only=False)
        steps = ["callrate", "sex", "het"]
        assert runner._filter_steps_by_decisions(steps) == steps

    def test_disable_filter_controls_reaches_legacy_args(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        runner, geno, out = _make_runner(tmp_path, warn_only=False)
        # Args request filter_controls=True; the decision must force it False,
        # so the assertion is non-trivial (default is already False).
        runner.args.variant_qc.filter_controls = True
        runner._validation_decisions = ValidationDecisions(disable_filter_controls=True)
        captured: Dict[str, Any] = {}

        def fake_single_step(step, step_input, step_output, legacy_args):
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

        def fake_single_step(step, step_input, step_output, legacy_args):
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

        def fake_single_step(step, step_input, step_output, legacy_args):
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
