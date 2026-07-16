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

from genotools.core.exceptions import QCError
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
