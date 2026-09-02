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

"""Pipeline runner for GenoTools.

This module provides the main pipeline orchestration logic, coordinating
QC steps, ancestry prediction, and GWAS analysis.
"""

from __future__ import annotations

import os
import pathlib
import platform
import dataclasses
import shutil
import tempfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING, Any, Dict, List, Optional

import pandas as pd

from .parser import PipelineArgs
from .output import ParameterRecord, PipelineOutput, write_results
from ..core.exceptions import GenoToolsError, ValidationError
from ..core.logging import (
    RunLog,
    get_logger,
    install_run_logging,
    step_context,
    summary_tag,
    summary_tally,
)
from ..core.provenance import package_versions, tool_lines
from ..core.validation import (
    ValidationDecisions,
    case_control_skip_reason,
    guard_output_not_exists,
    het_skip_reason,
    sex_skip_reason,
)
from ..qc.results import STEP_REPORT, unrun_result

if TYPE_CHECKING:  # QC config is imported lazily in _initialize_modules
    from ..ancestry.diagnostics import AncestryDiagnostics
    from ..qc.config import HetConfig

logger = get_logger(__name__)


def _flatten_config(config: Any) -> Dict[str, Any]:
    """Flatten a frozen config dataclass into JSON-safe scalar leaves.

    Most step configs are already flat. ``AssocConfig`` nests PCA, GWAS and
    covariate configs inside itself, and the parameters table is long form -
    one scalar per row - so nested members become dotted keys ("pca.n_pcs")
    rather than being dropped or dumped as a dict.

    Args:
        config: A dataclass instance, or None for a step that takes no config.

    Returns:
        Mapping of parameter name to a JSON-serializable value.
    """
    if config is None:
        return {}
    # asdict() has already recursed into nested dataclasses; this flattens the
    # nested dicts it produced.
    return _flatten_mapping(dataclasses.asdict(config), "")


def _flatten_mapping(mapping: Dict[str, Any], prefix: str) -> Dict[str, Any]:
    """Flatten nested dicts to dotted keys, coercing leaves to JSON-safe values."""
    flat: Dict[str, Any] = {}
    for name, value in mapping.items():
        key = f"{prefix}{name}"
        if isinstance(value, dict):
            flat.update(_flatten_mapping(value, f"{key}."))
        else:
            flat[key] = _json_safe(value)
    return flat


def _json_safe(value: Any) -> Any:
    """Coerce a config value to something ``json.dump`` accepts.

    Paths are the only non-scalar that reaches here (``covar_path``), and a
    Path is not JSON-serializable - the report would fail to write at the very
    end of a long run.
    """
    if isinstance(value, Path):
        return str(value)
    return value


@dataclass
class StepResult:
    """Result from a single pipeline step."""

    step: str
    passed: bool
    metrics: Dict[str, Any] = field(default_factory=dict)
    output: Dict[str, Any] = field(default_factory=dict)
    error: Optional[str] = None
    outcome: str = "pass"
    reason: Optional[str] = None

    @classmethod
    def from_legacy_dict(cls, step: str, legacy: Dict[str, Any]) -> "StepResult":
        """Create from legacy output dictionary format.

        ``step`` is the pipeline key ("callrate"); the result carries the name
        the step reports under ("callrate_prune"). Keep the reported name - the
        pre-refactor JSON used it, and the per-ancestry path still does.
        """
        return cls(
            step=legacy.get("step", step),
            passed=legacy.get("pass", False),
            metrics=legacy.get("metrics", {}),
            output=legacy.get("output", {}),
            outcome=legacy.get(
                "outcome", "pass" if legacy.get("pass", False) else "fail"
            ),
            reason=legacy.get("reason"),
        )


@dataclass
class PassFailRecord:
    """Record of a step's outcome.

    ``status`` stays a bool for the pre-refactor JSON contract; ``outcome``
    distinguishes the two false cases ("fail" vs "skipped"), and ``reason``
    explains it. Both are serialized - before 2.0 the cause of a failure was
    captured here and then dropped on the way to the report.
    """

    status: bool
    input_path: str
    output_path: str
    outcome: str = "pass"
    reason: Optional[str] = None


@dataclass
class PipelineState:
    """Tracks the current state of the pipeline execution."""

    geno_path: Path
    out_path: Path
    tmp_dir: Optional[tempfile.TemporaryDirectory] = None  # type: ignore[type-arg]
    current_input: Optional[str] = None
    pass_fail: Dict[str, PassFailRecord] = field(default_factory=dict)
    step_results: Dict[str, StepResult] = field(default_factory=dict)
    ancestry_result: Optional[Dict[str, Any]] = None
    labels_list: List[str] = field(default_factory=list)
    parameters: List[ParameterRecord] = field(default_factory=list)

    def get_last_passed_output(self) -> Optional[str]:
        """Get the output path of the last passed step."""
        last_passed = None
        for step_name, record in self.pass_fail.items():
            if record.status:
                last_passed = record.output_path
        return last_passed


class PipelineRunner:
    """Runs the GenoTools QC and analysis pipeline.

    This class orchestrates the execution of QC steps, ancestry prediction,
    and GWAS analysis using pure function modules.

    Attributes:
        args: Validated pipeline arguments.
        state: Current pipeline execution state.
    """

    # Step categories
    SAMPLE_STEPS = ["callrate", "sex", "het", "related", "kinship_check"]
    VARIANT_STEPS = ["geno", "case_control", "haplotype", "hwe", "ld"]

    def __init__(self, args: PipelineArgs) -> None:
        """Initialize the pipeline runner.

        Args:
            args: Validated pipeline arguments.
        """
        self.args = args
        self.state: Optional[PipelineState] = None
        self._validation_decisions = ValidationDecisions()
        self._runlog: Optional[RunLog] = None

        # These will be set up during run()
        self._new_modules: Dict[str, Any] = {}
        self._config_classes: Dict[str, Any] = {}

    def run(self) -> PipelineOutput:
        """Run the complete pipeline.

        Returns:
            PipelineOutput containing all results.

        Raises:
            ValidationError: If no steps are specified, the output prefix is
                already used, or the input fails validation.
        """
        # Initialize state
        self._initialize_state()

        # Initialize modules
        self._initialize_modules()

        # Re-run guard runs *before* logging setup, since setup creates the
        # consolidated log file the guard checks for. (Decoupled from
        # validate_input so validate_input can log into the fresh log.)
        guard_output_not_exists(self.args.out_path, self.args.output.skip_fails)

        # Set up logging (builds the RunLog, installs handlers, writes banner).
        self._setup_logging()

        # Record what software is about to produce the results, before it does.
        self._log_software_provenance()

        # Pre-flight: convert input + validate (logged into the consolidated log),
        # then confirm there is work to do. A failure here happens *before* any
        # output is produced, so tear down and REMOVE the freshly-created log so
        # the output prefix isn't poisoned for the next run (restores the
        # pre-redesign behavior where a failed validation left no log behind).
        try:
            self._begin_section("input_preparation")
            try:
                with step_context("input"):
                    self._convert_input_format()
            finally:
                self._end_section(f"{self.args.out_path}_input_preparation.log")

            steps = self.args.get_all_enabled_steps()
            skips = self._skip_reasons(steps, str(self.args.geno_path))
            runnable = [s for s in steps if s not in skips]
            if not runnable and not self.args.ancestry.run_ancestry:
                # A ValidationError (not a bare ValueError) so main() reports it
                # as a one-line message instead of a traceback. Say which steps
                # the data ruled out - "none requested" is wrong and unhelpful
                # when the user did request one.
                raise self._nothing_to_run_error(skips)
        except Exception:
            self._teardown_logging(remove_logs=True)
            if self.state and self.state.tmp_dir:
                self.state.tmp_dir.cleanup()
            raise

        # Run the pipeline. From here the consolidated log persists even on
        # failure (partial output may exist), matching the re-run guard's intent.
        try:
            if self.args.ancestry.run_ancestry:
                result = self._run_with_ancestry(steps)
            else:
                result = self._run_qc_only(steps)

            self._emit_run_summary()
            return result
        finally:
            # Cleanup temporary directory
            if self.state and self.state.tmp_dir:
                self.state.tmp_dir.cleanup()
            self._teardown_logging(remove_logs=False)

    def _initialize_state(self) -> None:
        """Initialize pipeline state."""
        out_dir = self.args.out_path.parent
        if not out_dir.exists():
            out_dir.mkdir(parents=True, exist_ok=True)

        # Create temporary directory
        tmp_dir = tempfile.TemporaryDirectory(
            suffix="_tmp", prefix=".", dir=str(out_dir)
        )

        self.state = PipelineState(
            geno_path=self.args.geno_path,
            out_path=self.args.out_path,
            tmp_dir=tmp_dir,
            current_input=str(self.args.geno_path),
        )

    def _initialize_modules(self) -> None:
        """Initialize QC, GWAS, and ancestry modules."""
        # Import QC modules
        from ..qc import (
            filter_callrate,
            filter_sex,
            filter_heterozygosity,
            filter_relatedness,
            verify_kinship,
            filter_case_control,
            filter_haplotype,
            filter_hwe,
            filter_variant_missingness,
            prune_ld,
            CallrateConfig,
            SexConfig,
            HetConfig,
            RelatedConfig,
            CaseControlConfig,
            HaplotypeConfig,
            HWEConfig,
            GenoConfig,
            LDConfig,
        )
        from ..core.genotypes import GenotypeData

        # Import GWAS module
        from ..gwas import (
            run_association as run_gwas_association,
            AssocConfig,
            PCAConfig,
            GWASConfig,
            CovariateConfig,
        )

        # Import Ancestry module
        from ..ancestry import AncestryModel, AncestryConfig, AncestryDiagnostics

        # Store module references
        self._new_modules = {
            "filter_callrate": filter_callrate,
            "filter_sex": filter_sex,
            "filter_heterozygosity": filter_heterozygosity,
            "filter_relatedness": filter_relatedness,
            "verify_kinship": verify_kinship,
            "filter_case_control": filter_case_control,
            "filter_haplotype": filter_haplotype,
            "filter_hwe": filter_hwe,
            "filter_variant_missingness": filter_variant_missingness,
            "prune_ld": prune_ld,
            "run_gwas_association": run_gwas_association,
            "GenotypeData": GenotypeData,
            "AncestryModel": AncestryModel,
            "AncestryDiagnostics": AncestryDiagnostics,
        }

        # Store config classes
        self._config_classes = {
            "CallrateConfig": CallrateConfig,
            "SexConfig": SexConfig,
            "HetConfig": HetConfig,
            "RelatedConfig": RelatedConfig,
            "CaseControlConfig": CaseControlConfig,
            "HaplotypeConfig": HaplotypeConfig,
            "HWEConfig": HWEConfig,
            "GenoConfig": GenoConfig,
            "LDConfig": LDConfig,
            "AssocConfig": AssocConfig,
            "PCAConfig": PCAConfig,
            "GWASConfig": GWASConfig,
            "CovariateConfig": CovariateConfig,
            "AncestryConfig": AncestryConfig,
        }

    def _setup_logging(self) -> None:
        """Build the consolidated run log and install logging handlers.

        Installs a :class:`RunLog` that owns ``{out}_all_logs.log`` (banner +
        per-step sections + summary) and registers it as the run-scoped raw sink
        so the executor harvests each PLINK/KING ``.log`` into it. Also attaches
        a concise console handler for the curated progress stream. Runs after the
        re-run guard (which checks the log does not already exist) and before any
        step, so every ``logger.info``/``warning``/``error`` — including the
        input breakdown — is captured.

        ``--debug`` widens both streams to DEBUG; ``--quiet`` drops the console
        handler (the consolidated + per-step log files are still written).
        """
        assert self.state is not None
        out_path = str(self.state.out_path)
        self._runlog = install_run_logging(
            out_path,
            level="DEBUG" if self.args.output.debug else "INFO",
            console=not self.args.output.quiet,
        )

    def _teardown_logging(self, remove_logs: bool = False) -> None:
        """Close the RunLog and unregister the raw sink.

        With ``remove_logs=True`` (pre-flight failure), also delete the
        freshly-created consolidated + input-prep logs so the output prefix is
        not left poisoned for the next run -- and put back any log this run
        rotated aside, so a failed re-run leaves the directory as it found it.
        """
        from ..core.logging import raw_sink

        if self._runlog is None:
            return
        runlog = self._runlog
        log_path = runlog.path
        runlog.close()
        self._runlog = None
        raw_sink.set(None)
        if remove_logs:
            for p in (log_path, Path(f"{self.args.out_path}_input_preparation.log")):
                try:
                    p.unlink(missing_ok=True)
                except OSError:
                    pass
            runlog.restore_rotated()

    def _log_software_provenance(self) -> None:
        """Record the interpreter, GenoTools, and the external tools resolved.

        GenoTools resolves plink/plink2/KING from its own executable folder and
        never consults ``PATH``, which is good for reproducibility and bad for
        working out after the fact *which* build produced a result — especially
        on a box that has a different, newer plink2 on ``PATH`` that the
        pipeline quietly ignores. Recording it costs one ``--version`` call per
        tool that is already on disk; a tool that has not been fetched yet is
        named as such rather than downloaded here, so a callrate-only run does
        not pull PLINK 1.9 just to describe it.

        File-only: the console stream is the curated progress view, and this
        belongs in the durable log. Best-effort throughout — provenance must
        never be the reason a run fails.
        """
        self._begin_section("software")
        try:
            with step_context("software"):
                versions = package_versions(())
                logger.info(
                    f"GenoTools {versions['genotools']} "
                    f"(python {versions['python']}, {platform.platform()})",
                    extra={"file_only": True},
                )
                for line in tool_lines():
                    logger.info(line, extra={"file_only": True})
        except Exception as e:
            # Best-effort by design: a run that cannot describe its own tools
            # is still a run worth finishing.
            logger.debug(f"Could not record software provenance: {e}")
        finally:
            self._end_section()

    def _begin_section(self, title: str) -> None:
        """Open a consolidated-log section (no-op when logging isn't installed)."""
        if self._runlog is not None:
            self._runlog.begin_section(title)

    def _end_section(self, raw_log_path: Optional[str] = None) -> None:
        """Flush a section's buffered raw output (no-op without a RunLog)."""
        if self._runlog is not None:
            self._runlog.end_section(raw_log_path)

    def _emit_run_summary(self) -> None:
        """Emit the end-of-run summary table to the console and consolidated log."""
        assert self.state is not None
        rows = self._collect_summary_rows()
        if not rows:
            return
        # Console-only records: the RunLog writes its own aligned table below, so
        # emitting these to the file too would duplicate the summary.
        logger.info("Run summary:", extra={"console_only": True})
        for label, removed, outcome in rows:
            logger.info(
                f"  {label}: {removed} removed [{summary_tag(outcome)}]",
                extra={"console_only": True},
            )
        tally = summary_tally(rows)
        if tally:
            logger.warning(f"  {tally}", extra={"console_only": True})
        if self._runlog is not None:
            self._runlog.write_summary(rows)

    def _collect_summary_rows(self) -> List[tuple]:
        """Assemble (label, outliers_removed, outcome) rows from pipeline state.

        Covers the plain QC path (``state.pass_fail`` + ``state.step_results``)
        and the per-ancestry path (``state.ancestry_results``), labeling
        per-ancestry rows ``{label}/{step}``.
        """
        assert self.state is not None

        def _removed(metrics: Dict[str, Any]) -> int:
            metrics = metrics or {}
            if "outlier_count" in metrics:
                return int(metrics["outlier_count"])
            return int(
                metrics.get("related_count", 0)
                + metrics.get("duplicated_count", 0)
            )

        rows: List[tuple] = []
        # Plain QC path.
        for step, record in self.state.pass_fail.items():
            result = self.state.step_results.get(step)
            removed = _removed(result.metrics) if result else 0
            rows.append((step, removed, record.outcome))

        # Per-ancestry path.
        ancestry_results = getattr(self.state, "ancestry_results", {})
        for label, out_dict in ancestry_results.items():
            pf = out_dict.get("pass_fail", {})
            for step, rec in pf.items():
                step_result = out_dict.get(step, {})
                removed = _removed(step_result.get("metrics", {}))
                outcome = rec.get(
                    "outcome", "pass" if rec.get("status", False) else "fail"
                )
                rows.append((f"{label}/{step}", removed, outcome))

        return rows

    def _convert_input_format(self) -> None:
        """Convert input format to pfiles if needed."""
        input_format = self.args.input.input_format

        from ..core import GenotypeData

        if input_format == "bfile":
            bfile = str(self.args.input.bfile)
            GenotypeData.from_path(bfile).to_pfile(bfile)
        elif input_format == "vcf":
            GenotypeData.from_vcf(str(self.args.input.vcf))

        # Run upfront validation
        from ..core.validation import validate_input

        self._validation_decisions = validate_input(
            self.args.geno_path,
            self.args.out_path,
            skip_fails=self.args.output.skip_fails,
            sex_requested=self.args.sample_qc.run_sex,
            het_requested=self.args.sample_qc.run_het,
            hwe_requested=self.args.variant_qc.run_hwe,
            filter_controls=self.args.variant_qc.filter_controls,
            case_control_requested=self.args.variant_qc.run_case_control,
        )

    def _run_with_ancestry(self, steps: List[str]) -> PipelineOutput:
        """Run pipeline with ancestry prediction.

        Args:
            steps: List of QC steps to run after ancestry split.

        Returns:
            PipelineOutput with results.
        """
        assert self.state is not None

        # Run ancestry prediction (sectioned; the ported preprocessing runs many
        # PLINK invocations, all harvested under the "ancestry" step context).
        self._begin_section("ancestry_prediction")
        try:
            with step_context("ancestry"):
                ancestry_result = self._run_ancestry_prediction_new()
        finally:
            self._end_section(f"{self.args.out_path}_ancestry_prediction.log")
        self.state.ancestry_result = ancestry_result
        self.state.labels_list = ancestry_result["data"]["labels_list"]

        # If no QC steps, just return ancestry results
        if not steps:
            return self._build_output()

        # An override naming a group that was never predicted would otherwise
        # do nothing at all, silently - --amr-het's failure mode on any panel
        # that does not spell the group "AMR".
        self._warn_unmatched_het_labels()

        for label in self.state.labels_list:
            # Set up paths for this ancestry
            geno_path = f"{self.args.out_path}_ancestry_{label}"
            out_path = f"{self.args.out_path}_{label}"

            # Decided per group: an ancestry group's sample count can differ
            # sharply from the cohort validate_input saw. Counting samples reads
            # the group's .psam, so this needs the on-disk prefix, not the
            # nominal one.
            group_skips = self._skip_reasons(
                steps, self._resolve_existing_geno(geno_path)
            )

            # Announce the group so the console (where per-step "Running:" lines
            # are file-only) still shows which ancestry the step lines belong to.
            n_running = len(steps) - len(group_skips)
            logger.info(f"Ancestry group {label}: running {n_running} QC step(s)")

            # Create modified args for this ancestry
            self._run_qc_pipeline(
                steps=steps,
                geno_path=geno_path,
                out_path=out_path,
                het_config=self.args.sample_qc.het_config_for(label),
                ancestry_label=label,
                skip_reasons=group_skips,
            )

        return self._build_output()

    def _warn_unmatched_het_labels(self) -> None:
        """Warn once for each --het-ancestry label that was never predicted.

        The label vocabulary is entirely user-supplied - from the --ref-labels
        TSV, or from a pickled model's encoder - so a typo, or a panel that
        spells the group differently, produces an override that matches
        nothing. Naming the predicted labels turns that into a fixable
        message rather than a silently ignored flag.
        """
        assert self.state is not None
        predicted = set(self.state.labels_list or [])
        unmatched = sorted(
            set(self.args.sample_qc.het_by_ancestry) - predicted
        )
        if not unmatched:
            return
        logger.warning(
            f"--het-ancestry {', '.join(unmatched)}: no such ancestry group in "
            f"this run, so the override does nothing. Predicted groups: "
            f"{', '.join(sorted(predicted))}."
        )

    def _resolve_existing_geno(self, geno_path: str) -> str:
        """Resolve the prefix whose pfiles are actually on disk.

        Under ``--ancestry`` without ``--full-output`` the cohort split writes
        into the temp working directory, so the nominal
        ``{out}_ancestry_{label}`` prefix never exists. Mirrors the working-path
        setup in ``_run_qc_pipeline`` so the two cannot disagree about where a
        group's data lives.
        """
        if self.args.full_output or not self.args.ancestry.run_ancestry:
            return geno_path
        assert self.state is not None
        return f"{self.state.tmp_dir.name}/{pathlib.PurePath(geno_path).name}"  # type: ignore[union-attr]

    @staticmethod
    def _nothing_to_run_error(skips: Dict[str, str]) -> ValidationError:
        """Explain why the run has no work to do.

        "No QC steps requested" is wrong when the user did request steps and the
        data ruled every one of them out, so name them and say why.
        """
        if skips:
            detail = "; ".join(
                f"{step} ({reason})" for step, reason in sorted(skips.items())
            )
            return ValidationError(f"No QC steps can run on this data: {detail}")
        return ValidationError("No QC steps or ancestry prediction requested")

    def _skip_reasons(
        self, steps: List[str], geno_path: Optional[str] = None
    ) -> Dict[str, str]:
        """Resolve which of ``steps`` the data rules out, and why.

        Combines the cohort-level decisions from ``validate_input`` with checks
        that can only be made against a specific dataset. Skipped steps are
        reported, not dropped - see ``_run_qc_pipeline``.

        Args:
            steps: Steps requested for this dataset.
            geno_path: The dataset being decided about. Under --ancestry this
                is one ancestry group, whose samples can differ sharply from the
                cohort ``validate_input`` saw. Omitted only where there is no
                single dataset to read.
        """
        reasons = {
            step: reason
            for step, reason in self._validation_decisions.skip_reasons.items()
            if step in steps
        }

        # Re-decide the sample-derived checks against this dataset: the cohort
        # can clear PLINK's LD floor, hold sample sex and hold both phenotypes
        # while an ancestry group does none of those. Without this the step
        # raises instead, reporting the same data-driven decision as
        # outcome="fail" per group but "skipped" cohort-wide. A cohort-level
        # decision already made stands - a group cannot resurrect a step the
        # whole cohort ruled out.
        if geno_path is not None:
            for step, decide in (
                ("het", het_skip_reason),
                ("sex", sex_skip_reason),
                ("case_control", case_control_skip_reason),
            ):
                if step in steps and step not in reasons:
                    reason = decide(geno_path)
                    if reason is not None:
                        reasons[step] = reason

        return reasons

    def _run_qc_only(self, steps: List[str]) -> PipelineOutput:
        """Run QC pipeline without ancestry prediction.

        Args:
            steps: List of QC steps to run.

        Returns:
            PipelineOutput with results.
        """
        assert self.state is not None

        # A flat run has no labels, so only the base spec can apply - but it
        # must still be *passed*. Omitting it here is what made --amr-het inert
        # in exactly the run shape the per-ancestry production workflow uses.
        self._run_qc_pipeline(
            steps=steps,
            geno_path=str(self.state.geno_path),
            out_path=str(self.state.out_path),
            het_config=self.args.sample_qc.het_config_for(None),
            skip_reasons=self._skip_reasons(steps, str(self.state.geno_path)),
        )

        return self._build_output()

    def _run_ancestry_prediction_new(self) -> Dict[str, Any]:
        """Run ancestry prediction using the new AncestryModel.

        Uses the ported ``get_raw_files`` for PLINK preprocessing, then
        delegates to new AncestryModel for the ML pipeline (PCA, UMAP,
        XGBoost, admixture detection). Cohort splitting uses the ported
        ``split_cohort_by_ancestry``.

        Supports two modes:
        - **Training** (default): Train a new model from reference panel.
        - **Inference** (``--model``): Load a pre-trained model directory.

        Returns:
            Ancestry result dictionary in legacy format.
        """
        assert self.state is not None

        AncestryModel = self._new_modules["AncestryModel"]

        # Determine output path (same logic as legacy method)
        has_qc_steps = self.args.has_any_qc_steps()
        if has_qc_steps:
            out_path = f"{self.state.out_path}_ancestry"
        else:
            out_path = str(self.state.out_path)

        # Use temporary directory if not full output
        if not self.args.full_output:
            out_name = pathlib.PurePath(out_path).name
            actual_out = f"{self.state.tmp_dir.name}/{out_name}"  # type: ignore[union-attr]
        else:
            actual_out = out_path

        # Determine training vs inference mode
        model_path = self.args.ancestry.model_path
        is_inference = model_path is not None

        # One collector for the whole run: preprocessing fills in the SNP
        # overlap and the allele frequencies, prediction fills in the PC drift
        # and the admixture decision, and the report gets all of it. Its
        # artifacts are written to `out_path` rather than `actual_out`, so the
        # evidence survives a run that discards its temp directory.
        diagnostics = self._new_modules["AncestryDiagnostics"]()

        if is_inference:
            model, predictions, ref_pca, raw = self._run_inference_mode(
                AncestryModel, model_path, actual_out, diagnostics, out_path
            )
        else:
            model, predictions, ref_pca, raw = self._run_training_mode(
                AncestryModel, actual_out, out_path, diagnostics
            )

        # --- Common post-processing (both modes) ---

        # Positive control. High training accuracy says the model learned the
        # panel; this says the *prediction path* can still recover it. When it
        # passes and a cohort still comes back one label, the cohort's own
        # preprocessing is the only thing left.
        if self.args.ancestry.self_test:
            diagnostics.self_test = model.self_test(
                raw["raw_ref"],
                detect_admixed=self.args.ancestry.detect_admixed,
            )

        projected_pca_path = f"{actual_out}_projected_new_pca.txt"
        projected_pca = pd.read_csv(projected_pca_path, sep="\t")

        if self.args.ancestry.write_plots:
            from ..ancestry.plots import plot_ancestry_diagnostics

            diagnostics.plots = plot_ancestry_diagnostics(
                train_pca=model._train_pca,
                projected=projected_pca,
                out_prefix=out_path,
                decisions=predictions.decisions,
            )

        # UMAP visualization
        from ..ancestry.reducers.umap_reducer import run_umap

        umap_result: Optional[Dict[str, Any]] = None
        if ref_pca is not None:
            umap_result = run_umap(
                ref_pca=ref_pca,
                new_pca=projected_pca,
                new_labels=predictions.predictions["predicted_ancestry"],
                params=model.best_params,
            )

        # Cohort splitting (ported; legacy-free)
        from ..ancestry.cohort import split_cohort_by_ancestry

        labels_path = f"{actual_out}_umap_linearsvc_predicted_labels.txt"
        ancestry_split = split_cohort_by_ancestry(
            labels_path=labels_path,
            geno_path=str(self.state.geno_path),
            out_path=actual_out,
            min_samples=self.args.ancestry.min_samples,
            subset=self.args.ancestry.subset_ancestry,
        )

        # Build legacy-format result dict
        pred_ids = predictions.predictions.rename(
            columns={"predicted_ancestry": "label"}
        )

        data_dict: Dict[str, Any] = {
            "predict_data": {
                "ids": pred_ids[["FID", "IID", "label"]],
                "y_pred": predictions.predictions["predicted_ancestry"].values,
                "label_encoder": model.label_encoder,
            },
            "confusion_matrix": (
                model.training_metrics.confusion_matrix
                if model.training_metrics
                else None
            ),
            "train_pcs": model._train_pca,
            "ref_pcs": ref_pca,
            "projected_pcs": projected_pca,
            "total_umap": umap_result["total_umap"] if umap_result else None,
            "ref_umap": umap_result["ref_umap"] if umap_result else None,
            "new_samples_umap": umap_result["new_umap"] if umap_result else None,
            "label_encoder": model.label_encoder,
            "labels_list": ancestry_split["labels"],
            "pruned_samples": ancestry_split["pruned_samples"],
        }

        if model.common_snps is not None:
            data_dict["common_snps"] = model.common_snps

        metrics_dict: Dict[str, Any] = {
            "predicted_counts": predictions.predictions[
                "predicted_ancestry"
            ].value_counts(),
            "test_accuracy": (
                model.training_metrics.test_accuracy
                if model.training_metrics
                else None
            ),
        }

        outfiles_dict: Dict[str, Any] = {
            "predicted_labels": {"labels_outpath": labels_path},
            "split_paths": ancestry_split["paths"],
        }

        out_dict: Dict[str, Any] = {
            "step": "predict_ancestry",
            "data": data_dict,
            "metrics": metrics_dict,
            "output": outfiles_dict,
            "diagnostics": diagnostics.to_dict(),
        }

        # Cleanup intermediate files
        files_to_clean = [
            raw["out_paths"].get("bed", ""),
            f"{actual_out}_common_snps",
        ]
        from ..ancestry.preprocessing import clean_up_files

        clean_up_files(files_to_clean)

        # Move files if ancestry-only and not full output
        if not has_qc_steps and not self.args.full_output:
            self._move_ancestry_files(out_dict)

        return out_dict

    def _run_training_mode(
        self,
        AncestryModel: type,
        actual_out: str,
        out_path: str,
        diagnostics: "AncestryDiagnostics",
    ) -> tuple:
        """Train a new AncestryModel from reference panel.

        Args:
            AncestryModel: The AncestryModel class.
            actual_out: Output path prefix (may be temp dir).
            out_path: Final output path (for saving model).
            diagnostics: Run-wide diagnostics collector.

        Returns:
            Tuple of (model, predictions, ref_pca, raw).
        """
        from ..ancestry.preprocessing import get_raw_files

        # PLINK preprocessing (ported; legacy-free)
        raw = get_raw_files(
            geno_path=str(self.state.geno_path),
            ref_panel=str(self.args.ancestry.ref_panel),
            ref_labels=str(self.args.ancestry.ref_labels),
            out_path=actual_out,
            train=True,
            fill_strategy=self.args.ancestry.missing_fill,
            max_missing_fraction=self.args.ancestry.max_missing_snps,
            diagnostics=diagnostics,
            diagnostics_prefix=out_path,
        )
        raw_ref = raw["raw_ref"]
        raw_geno = raw["raw_geno"]

        # Transform data for AncestryModel
        labels = pd.Series(
            raw_ref["label"].values,
            index=raw_ref["IID"].values,
            name="label",
        )
        ref_data = raw_ref.drop(columns=["label"])
        geno_data = raw_geno.drop(columns=["label"])
        geno_ids = raw_geno[["FID", "IID"]]

        # Capture common SNPs from ref column names
        snp_columns = [c for c in ref_data.columns if c not in ("FID", "IID")]

        # Fit model
        model = AncestryModel()
        model.fit(ref_data, labels, out_path=Path(actual_out))

        # Store common SNPs and save model directory to final output path
        model.common_snps = list(snp_columns)
        model_save_dir = Path(f"{out_path}_ancestry_model")
        model.save(model_save_dir)
        logger.info(f"Ancestry model saved to: {model_save_dir}")

        # Predict
        predictions = model.predict(
            geno_data,
            geno_ids,
            out_path=Path(actual_out),
            detect_admixed=self.args.ancestry.detect_admixed,
            diagnostics=diagnostics,
            diagnostics_prefix=out_path,
        )

        # Load ref PCA written by fit()
        ref_pca_path = f"{actual_out}_labeled_ref_pca.txt"
        ref_pca = pd.read_csv(ref_pca_path, sep="\t")

        return model, predictions, ref_pca, raw

    def _run_inference_mode(
        self,
        AncestryModel: type,
        model_path: Path,
        actual_out: str,
        diagnostics: "AncestryDiagnostics",
        out_path: str,
    ) -> tuple:
        """Load a pre-trained AncestryModel and predict.

        Args:
            AncestryModel: The AncestryModel class.
            model_path: Path to a model directory, or a single .pkl written
                by GenoTools 2.0.
            actual_out: Output path prefix.
            diagnostics: Run-wide diagnostics collector.
            out_path: Final output prefix, for artifacts the user keeps.

        Returns:
            Tuple of (model, predictions, ref_pca, raw).
        """
        model = AncestryModel.load(model_path)
        logger.info(f"Loaded pre-trained ancestry model from: {model_path}")

        # Write the model's common SNPs to a file for the ported preprocessing
        if model.common_snps is not None:
            common_snps_file = f"{actual_out}_loaded_model.common_snps"
            with open(common_snps_file, "w") as f:
                for snp in model.common_snps:
                    f.write(f"{snp}\n")
        else:
            model_dir = Path(model_path)
            snps_file = model_dir / "common_snps.txt" if model_dir.is_dir() else None
            if snps_file and snps_file.exists():
                common_snps_file = f"{actual_out}_loaded_model.common_snps"
                shutil.copy2(snps_file, common_snps_file)
            else:
                raise FileNotFoundError(f"No common_snps found in model: {model_path}")

        from ..ancestry.preprocessing import get_raw_files

        # PLINK preprocessing (ported; legacy-free)
        raw = get_raw_files(
            geno_path=str(self.state.geno_path),
            ref_panel=str(self.args.ancestry.ref_panel),
            ref_labels=str(self.args.ancestry.ref_labels),
            out_path=actual_out,
            train=False,
            common_snps_file=common_snps_file,
            fill_strategy=self.args.ancestry.missing_fill,
            max_missing_fraction=self.args.ancestry.max_missing_snps,
            diagnostics=diagnostics,
            diagnostics_prefix=out_path,
        )
        raw_geno = raw["raw_geno"]
        geno_data = raw_geno.drop(columns=["label"])
        geno_ids = raw_geno[["FID", "IID"]]

        # Skip fit(), go straight to predict()
        predictions = model.predict(
            geno_data,
            geno_ids,
            out_path=Path(actual_out),
            detect_admixed=self.args.ancestry.detect_admixed,
            diagnostics=diagnostics,
            diagnostics_prefix=out_path,
        )

        # Use training PCA as ref PCA for UMAP visualization
        ref_pca = model._train_pca

        return model, predictions, ref_pca, raw

    def _move_ancestry_files(self, ancestry_result: Dict[str, Any]) -> None:
        """Move ancestry output files from tmp to final location."""
        assert self.state is not None
        assert self.state.tmp_dir is not None

        out_dir = str(self.state.out_path.parent)
        out_name = self.state.out_path.name

        for label in ancestry_result["data"]["labels_list"]:
            pgen = f"{self.state.tmp_dir.name}/{out_name}_{label}.pgen"
            if os.path.isfile(pgen):
                os.rename(pgen, f"{out_dir}/{out_name}_{label}.pgen")
                os.rename(
                    f"{self.state.tmp_dir.name}/{out_name}_{label}.psam",
                    f"{out_dir}/{out_name}_{label}.psam",
                )
                os.rename(
                    f"{self.state.tmp_dir.name}/{out_name}_{label}.pvar",
                    f"{out_dir}/{out_name}_{label}.pvar",
                )

    def _run_qc_pipeline(
        self,
        steps: List[str],
        geno_path: str,
        out_path: str,
        het_config: Optional["HetConfig"] = None,
        ancestry_label: str = "all",
        skip_reasons: Optional[Dict[str, str]] = None,
    ) -> Dict[str, Any]:
        """Run QC pipeline for a single dataset.

        Args:
            steps: List of QC steps to run.
            geno_path: Input genotype path.
            out_path: Output path.
            het_config: Resolved HetConfig for this dataset. Callers resolve it
                through ``SampleQCArgs.het_config_for`` rather than passing raw
                thresholds, so the base/override precedence is decided in one
                place for every run shape.
            ancestry_label: Ancestry label for this run.
            skip_reasons: Steps the input data rules out, mapped to why. They
                stay in ``steps`` and are reported as skipped rather than
                dropped, so the report distinguishes "not requested" from
                "requested but impossible".

        Returns:
            Dictionary with step results and pass/fail status.
        """
        assert self.state is not None

        # Set up working paths
        if self.args.full_output:
            working_geno = geno_path
            working_out = out_path
        else:
            geno_name = pathlib.PurePath(geno_path).name
            out_name = pathlib.PurePath(out_path).name
            working_geno = f"{self.state.tmp_dir.name}/{geno_name}"  # type: ignore[union-attr]
            working_out = f"{self.state.tmp_dir.name}/{out_name}"  # type: ignore[union-attr]

        # Initialize tracking
        step_paths = [working_out]
        pass_fail: Dict[str, PassFailRecord] = {}
        out_dict: Dict[str, Any] = {}

        # Get legacy args for parameter passing
        legacy_args = self.args.to_legacy_dict()
        if het_config is None:
            het_config = self.args.sample_qc.het_config_for(None)
        if self._validation_decisions.disable_filter_controls:
            legacy_args["filter_controls"] = False

        skip_reasons = skip_reasons or {}

        # A skipped step produces no files, so it cannot be the one that writes
        # out_path - the last step that actually runs is.
        runnable = [i for i, step in enumerate(steps) if step not in skip_reasons]
        last_runnable = runnable[-1] if runnable else None

        if last_runnable is None and steps:
            # Every requested step was ruled out (e.g. --ancestry --het on a
            # group under the LD floor). Nothing will write out_path, so carry
            # the input through unchanged rather than leaving the group with no
            # output at all.
            logger.warning(
                "No QC steps could run; passing input through unchanged."
            )
            source = (
                geno_path
                if (self.args.full_output or not self.args.ancestry.run_ancestry)
                else working_geno
            )
            for ext in (".pgen", ".psam", ".pvar"):
                if os.path.isfile(f"{source}{ext}"):
                    shutil.copy2(f"{source}{ext}", f"{out_path}{ext}")

        # Run each step
        for i, step in enumerate(steps):
            step_input, step_output = self._compute_step_paths(
                step=step,
                step_index=i,
                steps=steps,
                pass_fail=pass_fail,
                geno_path=geno_path if (self.args.full_output or not self.args.ancestry.run_ancestry) else working_geno,
                out_path=out_path,
                working_out=working_out,
                is_last=(i == last_runnable),
            )

            # A data-driven skip records its outcome and passes the chain
            # through untouched: output_path is the input, so the next step
            # reads exactly what this one would have.
            if step in skip_reasons:
                reason = skip_reasons[step]
                logger.warning(f"Skipping {step}: {reason}.")
                pass_fail[step] = PassFailRecord(
                    status=False,
                    input_path=step_input,
                    output_path=step_input,
                    outcome="skipped",
                    reason=reason,
                )
                skipped = unrun_result(step, "skipped", reason, step_input)
                if skipped is not None:
                    out_dict[step] = skipped
                continue

            step_paths.append(step_output)

            # Open a consolidated-log section for this step. Raw PLINK output
            # harvested during the step is buffered and flushed (to the section
            # and to a per-step {out}_{step}.log) by _end_section. Logs always
            # persist at the stable out_path location regardless of full_output.
            section_title = step if ancestry_label == "all" else f"{ancestry_label} / {step}"
            raw_log_path = f"{out_path}_{step}.log"
            self._begin_section(section_title)
            try:
                # File-only: the absolute temp-dir paths are essential for
                # debugging but pure noise on the console, where the step's own
                # first log line (and, under --ancestry, the group header)
                # already announces what is running.
                logger.info(
                    f"Running: {step} with input {step_input} and output: {step_output}",
                    extra={"file_only": True},
                )

                # Check if input exists
                if self.args.warn_only and not os.path.isfile(f"{step_input}.pgen"):
                    logger.warning(
                        f"Step {step} cannot be run! "
                        "All samples or variants were pruned in a previous step!"
                    )
                    reason = (
                        "all samples or variants were pruned in a previous step"
                    )
                    pass_fail[step] = PassFailRecord(
                        status=False,
                        input_path=step_input,
                        output_path=step_output,
                        outcome="fail",
                        reason=reason,
                    )
                    failed = unrun_result(step, "fail", reason, step_output)
                    if failed is not None:
                        out_dict[step] = failed
                    continue

                # Run the step. Refactored steps signal failure by raising a
                # GenoToolsError (QCError/ExternalToolError/...). Under --warn we
                # record the failure and continue from the last passed output;
                # otherwise we fail fast.
                try:
                    result = self._run_single_step(
                        step=step,
                        step_input=step_input,
                        step_output=step_output,
                        legacy_args=legacy_args,
                        het_config=het_config,
                        ancestry_label=ancestry_label,
                    )
                except GenoToolsError as e:
                    logger.error(f"Step {step} failed: {e}")
                    if not self.args.warn_only:
                        raise
                    logger.warning(f"Step {step} failed but continuing (--warn): {e}")
                    pass_fail[step] = PassFailRecord(
                        status=False,
                        input_path=step_input,
                        output_path=step_output,
                        outcome="fail",
                        reason=str(e),
                    )
                    # Record the failure in the report too, not just pass_fail:
                    # an absent QC row is indistinguishable from a step that was
                    # never requested.
                    failed = unrun_result(step, "fail", str(e), step_output)
                    if failed is not None:
                        out_dict[step] = failed
                    continue

                if result is not None:
                    out_dict[step] = result
                    pass_fail[step] = PassFailRecord(
                        status=result.get("pass", False),
                        input_path=step_input,
                        output_path=step_output,
                    )

                    # Clean up intermediate files
                    self._cleanup_intermediate_files(
                        step=step,
                        pass_fail=pass_fail,
                        out_dict=out_dict,
                        out_path=out_path,
                        geno_path=geno_path,
                    )
            finally:
                self._end_section(raw_log_path)

        # Handle final step failure with --warn
        # Only steps that ran can have produced (or failed to produce) out_path;
        # passing the full list would make a trailing skip look like the
        # terminal failure and copy out_path onto itself.
        self._handle_final_step_failure(
            [step for step in steps if step not in skip_reasons],
            pass_fail,
            out_path,
            geno_path,
            working_geno,
        )

        out_dict["paths"] = step_paths
        out_dict["pass_fail"] = {
            k: {
                "status": v.status,
                "outcome": v.outcome,
                "reason": v.reason,
                "input": v.input_path,
                "output": v.output_path,
            }
            for k, v in pass_fail.items()
        }

        # Store results by ancestry label
        if ancestry_label == "all":
            self.state.step_results = {
                k: StepResult.from_legacy_dict(k, v)
                for k, v in out_dict.items()
                if k not in ("paths", "pass_fail")
            }
            self.state.pass_fail = pass_fail
        else:
            # Store per-ancestry results
            if not hasattr(self.state, "ancestry_results"):
                self.state.ancestry_results = {}  # type: ignore[attr-defined]
            self.state.ancestry_results[ancestry_label] = out_dict  # type: ignore[attr-defined]

        return out_dict

    def _compute_step_paths(
        self,
        step: str,
        step_index: int,
        steps: List[str],
        pass_fail: Dict[str, PassFailRecord],
        geno_path: str,
        out_path: str,
        working_out: str,
        is_last: bool,
    ) -> tuple[str, str]:
        """Compute input and output paths for a step.

        Args:
            step: Current step name.
            step_index: Index in steps list.
            steps: Full list of steps.
            pass_fail: Pass/fail tracking dictionary.
            geno_path: Original genotype path.
            out_path: Final output path.
            working_out: Working directory output path.
            is_last: Whether this is the last step that will actually run.
                Skipped steps write nothing, so a skip at the end of the list
                must not claim ``out_path`` - the last *runnable* step does.

        Returns:
            Tuple of (step_input, step_output).
        """

        if self.args.warn_only:
            # Find last passed step
            last_passed_output = None
            for completed_step, record in pass_fail.items():
                if record.status:
                    last_passed_output = record.output_path

            if last_passed_output:
                step_input = last_passed_output
                step_output = f"{last_passed_output}_{step}"
            else:
                step_input = geno_path
                step_output = (
                    f"{out_path}_{step}"
                    if self.args.full_output
                    else f"{working_out}_{step}"
                )

            if is_last:
                step_output = out_path
        else:
            # Sequential processing
            if step_index == 0:
                step_input = geno_path
            else:
                prev_step = steps[step_index - 1]
                if prev_step in pass_fail:
                    step_input = pass_fail[prev_step].output_path
                else:
                    step_input = (
                        f"{out_path}_{steps[step_index - 1]}"
                        if self.args.full_output
                        else f"{working_out}_{steps[step_index - 1]}"
                    )

            if is_last:
                step_output = out_path
            else:
                step_output = (
                    f"{out_path}_{step}"
                    if self.args.full_output
                    else f"{working_out}_{step}"
                )

        return step_input, step_output

    def _run_single_step(
        self,
        step: str,
        step_input: str,
        step_output: str,
        legacy_args: Dict[str, Any],
        het_config: Optional["HetConfig"] = None,
        ancestry_label: str = "all",
    ) -> Optional[Dict[str, Any]]:
        """Run a single QC step using pure function modules.

        Args:
            step: Step name.
            step_input: Input path.
            step_output: Output path.
            legacy_args: Legacy argument dictionary.
            het_config: Resolved HetConfig for the het step. Passed in already
                resolved rather than rebuilt from legacy_args, which cannot
                express the derived-bounds mode.
            ancestry_label: Group this step is running for, recorded alongside
                the step's parameters. A per-group setting is only legible in
                the report if the report says which group it applied to.

        Returns:
            Step result dictionary in legacy format, or None if step skipped.

        Every branch assigns ``config`` and ``result`` and falls through to the
        shared tail, which records the parameters. Returning from inside a
        branch would mean a new step silently reports no settings.
        """
        GenotypeData = self._new_modules["GenotypeData"]

        # Load input data
        data = GenotypeData.from_path(Path(step_input))

        config: Any = None
        result: Any = None

        # Map step names to functions and config creation
        if step == "callrate":
            config = self._config_classes["CallrateConfig"](mind=legacy_args["callrate"])
            result = self._new_modules["filter_callrate"](data, config, Path(step_output))

        elif step == "sex":
            sex_vals = legacy_args["sex"]
            config = self._config_classes["SexConfig"](female_max_f=sex_vals[0], male_min_f=sex_vals[1])
            result = self._new_modules["filter_sex"](data, config, Path(step_output))

        elif step == "het":
            if het_config is None:
                het_vals = legacy_args["het"]
                het_config = self._config_classes["HetConfig"](
                    f_lower=het_vals[0], f_upper=het_vals[1]
                )
            config = het_config
            result = self._new_modules["filter_heterozygosity"](
                data, het_config, Path(step_output)
            )

        elif step == "related":
            config = self._config_classes["RelatedConfig"](
                related_cutoff=legacy_args["related_cutoff"],
                duplicated_cutoff=legacy_args["duplicated_cutoff"],
                prune_related=legacy_args["prune_related"],
                prune_duplicated=legacy_args["prune_duplicated"],
            )
            result = self._new_modules["filter_relatedness"](data, config, Path(step_output))

        elif step == "kinship_check":
            if platform.system() != "Linux":
                logger.warning("Relatedness Assessment can only run on a Linux OS!")
                return None
            result = self._new_modules["verify_kinship"](data, Path(step_output))

        elif step == "case_control":
            config = self._config_classes["CaseControlConfig"](p_threshold=legacy_args["case_control"])
            result = self._new_modules["filter_case_control"](data, config, Path(step_output))

        elif step == "haplotype":
            config = self._config_classes["HaplotypeConfig"](p_threshold=legacy_args["haplotype"])
            result = self._new_modules["filter_haplotype"](data, config, Path(step_output))

        elif step == "hwe":
            config = self._config_classes["HWEConfig"](
                hwe_threshold=legacy_args["hwe"],
                filter_controls=legacy_args["filter_controls"],
            )
            result = self._new_modules["filter_hwe"](data, config, Path(step_output))

        elif step == "geno":
            config = self._config_classes["GenoConfig"](geno=legacy_args["geno"])
            result = self._new_modules["filter_variant_missingness"](data, config, Path(step_output))

        elif step == "ld":
            ld_vals = legacy_args["ld"]
            config = self._config_classes["LDConfig"](
                window_size=ld_vals[0],
                step_size=ld_vals[1],
                r2_threshold=ld_vals[2],
            )
            result = self._new_modules["prune_ld"](data, config, Path(step_output))

        elif step == "assoc":
            # legacy_args["pca"] is the requested PC count (or None); it doubles as
            # the run-PCA flag. Thread it into PCAConfig so a non-default --pca N is
            # honored instead of silently defaulting to 10.
            n_pcs = legacy_args.get("pca")
            covar_path = legacy_args.get("covars")
            covariates = (
                self._config_classes["CovariateConfig"](
                    covar_path=covar_path,
                    covar_names=legacy_args.get("covar_names"),
                )
                if covar_path
                else self._config_classes["CovariateConfig"]()
            )
            config = self._config_classes["AssocConfig"](
                pca=self._config_classes["PCAConfig"](
                    n_pcs=n_pcs,
                    build=legacy_args.get("build", "hg38"),
                ) if n_pcs else None,
                gwas=self._config_classes["GWASConfig"](
                    maf_lambdas=legacy_args.get("maf_lambdas", False),
                ) if legacy_args.get("gwas") else None,
                covariates=covariates,
                run_pca=bool(n_pcs),
                run_gwas=bool(legacy_args.get("gwas", False)),
            )
            result = self._new_modules["run_gwas_association"](data, Path(step_output), config)

        if result is None:
            return None

        out = result.to_dict()
        self._record_parameters(
            step=step,
            config=config,
            resolved=out.get("parameters", {}),
            ancestry_label=ancestry_label,
        )
        return out

    def _record_parameters(
        self,
        step: str,
        config: Any,
        resolved: Dict[str, Any],
        ancestry_label: str,
    ) -> None:
        """Record what this step was configured with, and what it worked out.

        The QC report is a long (step, metric, pruned_count) table, so it can
        say how many samples a step pruned and nothing about the threshold that
        did the pruning. With a per-group setting - ``--het sd 2`` resolving to
        different bounds in every ancestry group - two rows of that table are
        not comparable, and nothing in the file says why. This is where the
        settings go instead.

        Both halves are recorded because neither is recoverable from the other:
        the config is the request, and a step that derives its bounds from the
        data is the only thing that knows what the request became.

        Args:
            step: Pipeline step key ("het"). Recorded under the name the step
                reports as ("het_prune"), matching the QC table so the two
                sections join.
            config: The frozen config dataclass the step was handed, or None
                for a step that takes no config.
            resolved: Values the step worked out at runtime, from
                ``FilterResult.parameters``. Empty for steps that derive
                nothing.
            ancestry_label: Group this step ran for, or "all" in a flat run.
        """
        assert self.state is not None

        reported = STEP_REPORT.get(step, (step, ()))[0]
        for name, value in _flatten_config(config).items():
            self.state.parameters.append(
                ParameterRecord(
                    step=reported,
                    parameter=name,
                    value=value,
                    ancestry=ancestry_label,
                    source="requested",
                )
            )
        for name, value in resolved.items():
            self.state.parameters.append(
                ParameterRecord(
                    step=reported,
                    parameter=name,
                    value=_json_safe(value),
                    ancestry=ancestry_label,
                    source="resolved",
                )
            )

    def _cleanup_intermediate_files(
        self,
        step: str,
        pass_fail: Dict[str, PassFailRecord],
        out_dict: Dict[str, Any],
        out_path: str,
        geno_path: str,
    ) -> None:
        """Clean up intermediate pfiles if not needed.

        Args:
            step: Current step name.
            pass_fail: Pass/fail tracking dictionary.
            out_dict: Step output dictionary.
            out_path: Final output path.
            geno_path: Original genotype path.
        """
        if self.args.full_output:
            return
        if self.args.warn_only:
            return
        if step in ("assoc", "ancestry", "kinship_check"):
            return

        # Don't remove if step failed with warn mode
        if (
            self.args.warn_only
            and "pass" in out_dict.get(step, {})
            and not out_dict[step]["pass"]
        ):
            return

        remove_path = pass_fail[step].input_path
        if (
            os.path.isfile(f"{remove_path}.pgen")
            and remove_path != out_path
            and remove_path != geno_path
        ):
            os.remove(f"{remove_path}.pgen")
            os.remove(f"{remove_path}.psam")
            os.remove(f"{remove_path}.pvar")

    def _handle_final_step_failure(
        self,
        steps: List[str],
        pass_fail: Dict[str, PassFailRecord],
        out_path: str,
        geno_path: str,
        working_geno: str,
    ) -> None:
        """Handle case where final step fails with --warn.

        Args:
            steps: List of steps.
            pass_fail: Pass/fail tracking dictionary.
            out_path: Final output path.
            geno_path: Original genotype path.
            working_geno: Working directory genotype path.
        """
        if not self.args.warn_only:
            return
        if len(steps) <= 1:
            return
        if steps[-1] not in pass_fail:
            return
        if pass_fail[steps[-1]].status:
            return

        # Find last passed output
        last_passed_output = None
        last_passed_step_name = None
        for step_name in steps[:-1]:
            if step_name in pass_fail and pass_fail[step_name].status:
                last_passed_output = pass_fail[step_name].output_path
                last_passed_step_name = step_name

        if last_passed_output and os.path.isfile(f"{last_passed_output}.pgen"):
            shutil.copy2(f"{last_passed_output}.pgen", f"{out_path}.pgen")
            shutil.copy2(f"{last_passed_output}.psam", f"{out_path}.psam")
            shutil.copy2(f"{last_passed_output}.pvar", f"{out_path}.pvar")
        elif last_passed_output:
            # All samples/variants pruned
            input_path = pass_fail[last_passed_step_name].input_path
            if os.path.isfile(f"{input_path}.pgen"):
                shutil.copy2(f"{input_path}.pgen", f"{out_path}.pgen")
                shutil.copy2(f"{input_path}.psam", f"{out_path}.psam")
                shutil.copy2(f"{input_path}.pvar", f"{out_path}.pvar")
        else:
            # No steps passed - copy (never rename) to preserve originals
            move_path = (
                geno_path
                if self.args.full_output or not self.args.ancestry.run_ancestry
                else working_geno
            )
            if os.path.isfile(f"{move_path}.pgen"):
                shutil.copy2(f"{move_path}.pgen", f"{out_path}.pgen")
                shutil.copy2(f"{move_path}.psam", f"{out_path}.psam")
                shutil.copy2(f"{move_path}.pvar", f"{out_path}.pvar")

    def _build_output(self) -> PipelineOutput:
        """Build the final pipeline output.

        Returns:
            PipelineOutput containing all results.
        """
        assert self.state is not None

        return PipelineOutput.from_runner_state(
            args=self.args,
            state=self.state,
        )


def run_pipeline(args: PipelineArgs) -> PipelineOutput:
    """Run the GenoTools pipeline.

    This is the main entry point for running the pipeline programmatically.

    Args:
        args: Validated pipeline arguments.

    Returns:
        PipelineOutput containing all results.
    """
    runner = PipelineRunner(args)
    return runner.run()
