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

Supports both legacy class-based modules and new pure function architecture.
Set environment variable GENOTOOLS_USE_NEW_MODULES=1 to use the new architecture.
"""

from __future__ import annotations

import logging
import os
import pathlib
import platform
import tempfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Callable

from .parser import PipelineArgs
from .output import PipelineOutput, write_results

logger = logging.getLogger(__name__)

# Feature flag to enable new pure function modules
USE_NEW_MODULES = os.environ.get("GENOTOOLS_USE_NEW_MODULES", "0") == "1"


@dataclass
class StepResult:
    """Result from a single pipeline step."""

    step: str
    passed: bool
    metrics: Dict[str, Any] = field(default_factory=dict)
    output: Dict[str, Any] = field(default_factory=dict)
    error: Optional[str] = None

    @classmethod
    def from_legacy_dict(cls, step: str, legacy: Dict[str, Any]) -> "StepResult":
        """Create from legacy output dictionary format."""
        return cls(
            step=step,
            passed=legacy.get("pass", False),
            metrics=legacy.get("metrics", {}),
            output=legacy.get("output", {}),
        )


@dataclass
class PassFailRecord:
    """Record of step pass/fail status."""

    status: bool
    input_path: str
    output_path: str
    error: Optional[str] = None


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
    and GWAS analysis. It supports both the new modular architecture and
    backward compatibility with legacy code.

    Attributes:
        args: Validated pipeline arguments.
        state: Current pipeline execution state.
    """

    # Step categories
    SAMPLE_STEPS = ["callrate", "sex", "het", "related", "kinship_check"]
    VARIANT_STEPS = ["case_control", "haplotype", "hwe", "geno", "ld"]

    def __init__(self, args: PipelineArgs) -> None:
        """Initialize the pipeline runner.

        Args:
            args: Validated pipeline arguments.
        """
        self.args = args
        self.state: Optional[PipelineState] = None

        # These will be set up during run()
        self._samp_qc: Any = None
        self._var_qc: Any = None
        self._ancestry: Any = None
        self._assoc: Any = None
        self._steps_dict: Dict[str, Callable[..., Dict[str, Any]]] = {}

    def run(self) -> PipelineOutput:
        """Run the complete pipeline.

        Returns:
            PipelineOutput containing all results.

        Raises:
            ValueError: If no input files or steps are specified.
        """
        # Initialize state
        self._initialize_state()

        # Initialize QC classes (legacy)
        self._initialize_qc_classes()

        # Convert input format if needed (must happen before logging setup
        # because upfront_check validates that log files don't exist)
        self._convert_input_format()

        # Set up logging (after upfront_check so log file creation doesn't trigger error)
        self._setup_logging()

        # Validate we have something to do
        steps = self.args.get_all_enabled_steps()
        if not steps and not self.args.ancestry.run_ancestry:
            raise ValueError("No QC steps or ancestry prediction requested")

        # Run the pipeline
        try:
            if self.args.ancestry.run_ancestry:
                result = self._run_with_ancestry(steps)
            else:
                result = self._run_qc_only(steps)

            return result
        finally:
            # Cleanup temporary directory
            if self.state and self.state.tmp_dir:
                self.state.tmp_dir.cleanup()

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

    def _initialize_qc_classes(self) -> None:
        """Initialize QC modules.

        Uses new pure function modules if USE_NEW_MODULES is enabled,
        otherwise falls back to legacy class-based modules.
        """
        if USE_NEW_MODULES:
            self._initialize_new_modules()
        else:
            self._initialize_legacy_modules()

    def _initialize_legacy_modules(self) -> None:
        """Initialize legacy QC classes for backward compatibility."""
        # Import here to avoid circular imports and allow gradual migration
        # Use importlib to explicitly import the legacy module FILES (qc.py, ancestry.py, gwas.py)
        # rather than the new package directories (qc/, ancestry/, gwas/)
        import importlib.util
        import sys
        from pathlib import Path

        genotools_dir = Path(__file__).parent.parent

        # Load legacy qc.py module
        qc_spec = importlib.util.spec_from_file_location(
            "genotools.legacy_qc", genotools_dir / "qc.py"
        )
        legacy_qc = importlib.util.module_from_spec(qc_spec)  # type: ignore
        sys.modules["genotools.legacy_qc"] = legacy_qc
        qc_spec.loader.exec_module(legacy_qc)  # type: ignore
        SampleQC = legacy_qc.SampleQC
        VariantQC = legacy_qc.VariantQC

        # Load legacy ancestry.py module
        ancestry_spec = importlib.util.spec_from_file_location(
            "genotools.legacy_ancestry", genotools_dir / "ancestry.py"
        )
        legacy_ancestry = importlib.util.module_from_spec(ancestry_spec)  # type: ignore
        sys.modules["genotools.legacy_ancestry"] = legacy_ancestry
        ancestry_spec.loader.exec_module(legacy_ancestry)  # type: ignore
        Ancestry = legacy_ancestry.Ancestry

        # Load legacy gwas.py module
        gwas_spec = importlib.util.spec_from_file_location(
            "genotools.legacy_gwas", genotools_dir / "gwas.py"
        )
        legacy_gwas = importlib.util.module_from_spec(gwas_spec)  # type: ignore
        sys.modules["genotools.legacy_gwas"] = legacy_gwas
        gwas_spec.loader.exec_module(legacy_gwas)  # type: ignore
        Assoc = legacy_gwas.Assoc

        self._samp_qc = SampleQC()
        self._var_qc = VariantQC()
        self._ancestry = Ancestry()
        self._assoc = Assoc()

        # Build step dispatch dictionary
        self._steps_dict = {
            "callrate": self._samp_qc.run_callrate_prune,
            "sex": self._samp_qc.run_sex_prune,
            "het": self._samp_qc.run_het_prune,
            "related": self._samp_qc.run_related_prune,
            "kinship_check": self._samp_qc.run_confirming_kinship,
            "case_control": self._var_qc.run_case_control_prune,
            "haplotype": self._var_qc.run_haplotype_prune,
            "hwe": self._var_qc.run_hwe_prune,
            "geno": self._var_qc.run_geno_prune,
            "ld": self._var_qc.run_ld_prune,
            "assoc": self._assoc.run_association,
        }

    def _initialize_new_modules(self) -> None:
        """Initialize new pure function modules.

        Sets up the step dispatch dictionary to use the new architecture
        with pure functions and typed configs.
        """
        logger.info("Using new pure function modules (GENOTOOLS_USE_NEW_MODULES=1)")

        # Import new QC modules
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

        # Import new GWAS module
        from ..gwas import run_association as run_gwas_association, AssocConfig, PCAConfig, GWASConfig

        # Import new Ancestry module
        from ..ancestry import AncestryModel, ReferencePanel, AncestryConfig

        # Store module references for later use
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
            "ReferencePanel": ReferencePanel,
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
            "AncestryConfig": AncestryConfig,
        }

        # We still need legacy ancestry module for now (the new ancestry module
        # doesn't have the same interface yet)
        # Use importlib to load the legacy ancestry.py file
        import importlib.util
        import sys
        from pathlib import Path as PathLib

        genotools_dir = PathLib(__file__).parent.parent
        ancestry_spec = importlib.util.spec_from_file_location(
            "genotools.legacy_ancestry", genotools_dir / "ancestry.py"
        )
        legacy_ancestry = importlib.util.module_from_spec(ancestry_spec)  # type: ignore
        sys.modules["genotools.legacy_ancestry"] = legacy_ancestry
        ancestry_spec.loader.exec_module(legacy_ancestry)  # type: ignore
        Ancestry = legacy_ancestry.Ancestry
        self._ancestry = Ancestry()

        # _steps_dict is not used when USE_NEW_MODULES is enabled
        # Instead, we call the new functions directly in _run_single_step_new()
        self._steps_dict = {}

    def _setup_logging(self) -> None:
        """Set up log files."""
        assert self.state is not None
        out_path = str(self.state.out_path)

        # Clear existing log files
        all_logs = f"{out_path}_all_logs.log"
        cleaned_logs = f"{out_path}_cleaned_logs.log"

        if os.path.exists(all_logs):
            os.remove(all_logs)
        if os.path.exists(cleaned_logs):
            os.remove(cleaned_logs)

        # Create new log files with header
        from ..utils import gt_header

        header = gt_header()
        with open(all_logs, "w") as fp:
            fp.write(header)
            fp.write("\n")
        with open(cleaned_logs, "w") as fp:
            pass

    def _convert_input_format(self) -> None:
        """Convert input format to pfiles if needed."""
        input_format = self.args.input.input_format

        if input_format == "bfile":
            from ..utils import bfiles_to_pfiles

            bfiles_to_pfiles(bfile_path=str(self.args.input.bfile))
        elif input_format == "vcf":
            from ..utils import vcf_to_pfiles

            vcf_to_pfiles(vcf_path=str(self.args.input.vcf))

        # Run upfront validation
        if not self.args.output.skip_fails:
            from ..utils import upfront_check

            legacy_dict = self.args.to_legacy_dict()
            upfront_check(str(self.args.geno_path), legacy_dict)

    def _run_with_ancestry(self, steps: List[str]) -> PipelineOutput:
        """Run pipeline with ancestry prediction.

        Args:
            steps: List of QC steps to run after ancestry split.

        Returns:
            PipelineOutput with results.
        """
        assert self.state is not None

        # Run ancestry prediction
        ancestry_result = self._run_ancestry_prediction()
        self.state.ancestry_result = ancestry_result
        self.state.labels_list = ancestry_result["data"]["labels_list"]

        # If no QC steps, just return ancestry results
        if not steps:
            return self._build_output()

        # Run QC for each ancestry group
        het_value = [self.args.sample_qc.het_lower, self.args.sample_qc.het_upper]

        for label in self.state.labels_list:
            # Set up paths for this ancestry
            geno_path = f"{self.args.out_path}_ancestry_{label}"
            out_path = f"{self.args.out_path}_{label}"

            # Handle AMR heterozygosity
            current_het = het_value
            if self.args.sample_qc.amr_het and label == "AMR":
                current_het = [-1.0, -1.0]

            # Create modified args for this ancestry
            self._run_qc_pipeline(
                steps=steps,
                geno_path=geno_path,
                out_path=out_path,
                het_values=current_het,
                ancestry_label=label,
            )

        return self._build_output()

    def _run_qc_only(self, steps: List[str]) -> PipelineOutput:
        """Run QC pipeline without ancestry prediction.

        Args:
            steps: List of QC steps to run.

        Returns:
            PipelineOutput with results.
        """
        assert self.state is not None

        self._run_qc_pipeline(
            steps=steps,
            geno_path=str(self.state.geno_path),
            out_path=str(self.state.out_path),
        )

        return self._build_output()

    def _run_ancestry_prediction(self) -> Dict[str, Any]:
        """Run ancestry prediction.

        Returns:
            Ancestry result dictionary.
        """
        assert self.state is not None

        # Determine output path
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

        # Configure ancestry object
        self._ancestry.geno_path = str(self.state.geno_path)
        self._ancestry.out_path = actual_out
        self._ancestry.final_out_path = out_path
        self._ancestry.ref_panel = (
            str(self.args.ancestry.ref_panel)
            if self.args.ancestry.ref_panel
            else None
        )
        self._ancestry.ref_labels = (
            str(self.args.ancestry.ref_labels)
            if self.args.ancestry.ref_labels
            else None
        )
        self._ancestry.model_path = (
            str(self.args.ancestry.model_path)
            if self.args.ancestry.model_path
            else None
        )
        self._ancestry.containerized = self.args.ancestry.use_container
        self._ancestry.singularity = self.args.ancestry.use_singularity
        self._ancestry.subset = self.args.ancestry.subset_ancestry
        self._ancestry.min_samples = self.args.ancestry.min_samples

        # Run prediction
        result = self._ancestry.run_ancestry()

        # Move files if ancestry-only and not full output
        if not has_qc_steps and not self.args.full_output:
            self._move_ancestry_files(result)

        return result

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
        het_values: Optional[List[float]] = None,
        ancestry_label: str = "all",
    ) -> Dict[str, Any]:
        """Run QC pipeline for a single dataset.

        Args:
            steps: List of QC steps to run.
            geno_path: Input genotype path.
            out_path: Output path.
            het_values: Optional heterozygosity thresholds.
            ancestry_label: Ancestry label for this run.

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
        if het_values is not None:
            legacy_args["het"] = het_values

        # Run each step
        for i, step in enumerate(steps):
            step_input, step_output = self._compute_step_paths(
                step=step,
                step_index=i,
                steps=steps,
                pass_fail=pass_fail,
                geno_path=geno_path if self.args.full_output else working_geno,
                out_path=out_path,
                working_out=working_out,
            )

            step_paths.append(step_output)
            logger.info(
                f"Running: {step} with input {step_input} and output: {step_output}"
            )
            print(f"Running: {step} with input {step_input} and output: {step_output}")

            # Check if input exists
            if self.args.warn_only and not os.path.isfile(f"{step_input}.pgen"):
                print(
                    f"Step {step} cannot be run! "
                    "All samples or variants were pruned in a previous step!"
                )
                pass_fail[step] = PassFailRecord(
                    status=False, input_path=step_input, output_path=step_output
                )
                continue

            # Run the step
            result = self._run_single_step(
                step=step,
                step_input=step_input,
                step_output=step_output,
                legacy_args=legacy_args,
            )

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

        # Handle final step failure with --warn
        self._handle_final_step_failure(steps, pass_fail, out_path, geno_path, working_geno)

        out_dict["paths"] = step_paths
        out_dict["pass_fail"] = {
            k: {"status": v.status, "input": v.input_path, "output": v.output_path}
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

        Returns:
            Tuple of (step_input, step_output).
        """
        is_last = step_index == len(steps) - 1

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
    ) -> Optional[Dict[str, Any]]:
        """Run a single QC step.

        Args:
            step: Step name.
            step_input: Input path.
            step_output: Output path.
            legacy_args: Legacy argument dictionary.

        Returns:
            Step result dictionary or None if step skipped.
        """
        if USE_NEW_MODULES:
            return self._run_single_step_new(step, step_input, step_output, legacy_args)

        return self._run_single_step_legacy(step, step_input, step_output, legacy_args)

    def _run_single_step_new(
        self,
        step: str,
        step_input: str,
        step_output: str,
        legacy_args: Dict[str, Any],
    ) -> Optional[Dict[str, Any]]:
        """Run a single QC step using new pure function modules.

        Args:
            step: Step name.
            step_input: Input path.
            step_output: Output path.
            legacy_args: Legacy argument dictionary.

        Returns:
            Step result dictionary in legacy format, or None if step skipped.
        """
        GenotypeData = self._new_modules["GenotypeData"]

        # Load input data
        data = GenotypeData.from_path(Path(step_input))

        # Map step names to functions and config creation
        if step == "callrate":
            config = self._config_classes["CallrateConfig"](mind=legacy_args["callrate"])
            result = self._new_modules["filter_callrate"](data, config, Path(step_output))
            return result.to_dict()

        elif step == "sex":
            sex_vals = legacy_args["sex"]
            config = self._config_classes["SexConfig"](female_max_f=sex_vals[0], male_min_f=sex_vals[1])
            result = self._new_modules["filter_sex"](data, config, Path(step_output))
            return result.to_dict()

        elif step == "het":
            het_vals = legacy_args["het"]
            config = self._config_classes["HetConfig"](f_lower=het_vals[0], f_upper=het_vals[1])
            result = self._new_modules["filter_heterozygosity"](data, config, Path(step_output))
            return result.to_dict()

        elif step == "related":
            config = self._config_classes["RelatedConfig"](
                related_cutoff=legacy_args["related_cutoff"],
                duplicated_cutoff=legacy_args["duplicated_cutoff"],
                prune_related=legacy_args["prune_related"],
                prune_duplicated=legacy_args["prune_duplicated"],
            )
            result = self._new_modules["filter_relatedness"](data, config, Path(step_output))
            return result.to_dict()

        elif step == "kinship_check":
            if platform.system() != "Linux":
                print("Relatedness Assessment can only run on a Linux OS!")
                return None
            result = self._new_modules["verify_kinship"](data, Path(step_output))
            return result.to_dict()

        elif step == "case_control":
            config = self._config_classes["CaseControlConfig"](p_threshold=legacy_args["case_control"])
            result = self._new_modules["filter_case_control"](data, config, Path(step_output))
            return result.to_dict()

        elif step == "haplotype":
            config = self._config_classes["HaplotypeConfig"](p_threshold=legacy_args["haplotype"])
            result = self._new_modules["filter_haplotype"](data, config, Path(step_output))
            return result.to_dict()

        elif step == "hwe":
            config = self._config_classes["HWEConfig"](
                hwe_threshold=legacy_args["hwe"],
                filter_controls=legacy_args["filter_controls"],
            )
            result = self._new_modules["filter_hwe"](data, config, Path(step_output))
            return result.to_dict()

        elif step == "geno":
            config = self._config_classes["GenoConfig"](geno=legacy_args["geno"])
            result = self._new_modules["filter_variant_missingness"](data, config, Path(step_output))
            return result.to_dict()

        elif step == "ld":
            ld_vals = legacy_args["ld"]
            config = self._config_classes["LDConfig"](
                window_size=ld_vals[0],
                step_size=ld_vals[1],
                r2_threshold=ld_vals[2],
            )
            result = self._new_modules["prune_ld"](data, config, Path(step_output))
            return result.to_dict()

        elif step == "assoc":
            # Use the new GWAS module
            config = self._config_classes["AssocConfig"](
                pca=self._config_classes["PCAConfig"](
                    build=legacy_args.get("build", "hg38"),
                ) if legacy_args.get("pca") else None,
                gwas=self._config_classes["GWASConfig"](
                    maf_lambdas=legacy_args.get("maf_lambdas", False),
                ) if legacy_args.get("gwas") else None,
                run_pca=legacy_args.get("pca", False),
                run_gwas=legacy_args.get("gwas", False),
            )
            result = self._new_modules["run_gwas_association"](data, Path(step_output), config)
            return result.to_dict()

        return None

    def _run_single_step_legacy(
        self,
        step: str,
        step_input: str,
        step_output: str,
        legacy_args: Dict[str, Any],
    ) -> Optional[Dict[str, Any]]:
        """Run a single QC step using legacy class-based modules.

        Args:
            step: Step name.
            step_input: Input path.
            step_output: Output path.
            legacy_args: Legacy argument dictionary.

        Returns:
            Step result dictionary or None if step skipped.
        """
        if step in self.SAMPLE_STEPS:
            self._samp_qc.geno_path = step_input
            self._samp_qc.out_path = step_output

            if step == "related":
                return self._steps_dict[step](
                    related_cutoff=legacy_args["related_cutoff"],
                    duplicated_cutoff=legacy_args["duplicated_cutoff"],
                    prune_related=legacy_args["prune_related"],
                    prune_duplicated=legacy_args["prune_duplicated"],
                )
            elif step == "kinship_check":
                if platform.system() != "Linux":
                    print("Relatedness Assessment can only run on a Linux OS!")
                    return None
                return self._steps_dict[step]()
            else:
                return self._steps_dict[step](legacy_args[step])

        elif step in self.VARIANT_STEPS:
            self._var_qc.geno_path = step_input
            self._var_qc.out_path = step_output

            if step == "hwe":
                return self._steps_dict[step](
                    hwe_threshold=legacy_args["hwe"],
                    filter_controls=legacy_args["filter_controls"],
                )
            elif step == "ld":
                return self._steps_dict[step](
                    window_size=legacy_args["ld"][0],
                    step_size=legacy_args["ld"][1],
                    r2_threshold=legacy_args["ld"][2],
                )
            else:
                return self._steps_dict[step](legacy_args[step])

        elif step == "assoc":
            self._assoc.geno_path = step_input
            self._assoc.out_path = step_output
            self._assoc.pca = legacy_args["pca"]
            self._assoc.build = legacy_args["build"]
            self._assoc.gwas = legacy_args["gwas"]
            self._assoc.covar_path = legacy_args["covars"]
            self._assoc.covar_names = legacy_args["covar_names"]
            self._assoc.maf_lambdas = legacy_args["maf_lambdas"]
            return self._steps_dict[step]()

        return None

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
        for step_name in steps[:-1]:
            if step_name in pass_fail and pass_fail[step_name].status:
                last_passed_output = pass_fail[step_name].output_path

        if last_passed_output and os.path.isfile(f"{last_passed_output}.pgen"):
            os.rename(f"{last_passed_output}.pgen", f"{out_path}.pgen")
            os.rename(f"{last_passed_output}.psam", f"{out_path}.psam")
            os.rename(f"{last_passed_output}.pvar", f"{out_path}.pvar")
        elif last_passed_output:
            # All samples/variants pruned
            input_path = pass_fail[last_passed_output].input_path
            if os.path.isfile(f"{input_path}.pgen"):
                os.rename(f"{input_path}.pgen", f"{out_path}.pgen")
                os.rename(f"{input_path}.psam", f"{out_path}.psam")
                os.rename(f"{input_path}.pvar", f"{out_path}.pvar")
        else:
            # No steps passed - use original input
            move_path = (
                geno_path
                if self.args.full_output or not self.args.ancestry.run_ancestry
                else working_geno
            )
            if os.path.isfile(f"{move_path}.pgen"):
                os.rename(f"{move_path}.pgen", f"{out_path}.pgen")
                os.rename(f"{move_path}.psam", f"{out_path}.psam")
                os.rename(f"{move_path}.pvar", f"{out_path}.pvar")

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
