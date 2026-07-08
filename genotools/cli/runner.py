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

import logging
import os
import pathlib
import platform
import shutil
import tempfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional

import pandas as pd

from .parser import PipelineArgs
from .output import PipelineOutput, write_results
from ..core.exceptions import GenoToolsError

logger = logging.getLogger(__name__)


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
    and GWAS analysis using pure function modules.

    Attributes:
        args: Validated pipeline arguments.
        state: Current pipeline execution state.
    """

    # Step categories
    SAMPLE_STEPS = ["callrate", "sex", "het", "related", "kinship_check"]
    VARIANT_STEPS = ["geno", "case_control", "haplotype", "hwe", "ld"]

    def __init__(self, args: PipelineArgs, use_new_ancestry: bool = False) -> None:
        """Initialize the pipeline runner.

        Args:
            args: Validated pipeline arguments.
            use_new_ancestry: If True, use new AncestryModel instead of legacy.
        """
        self.args = args
        self.state: Optional[PipelineState] = None
        self._use_new_ancestry = use_new_ancestry

        # These will be set up during run()
        self._ancestry: Any = None
        self._new_modules: Dict[str, Any] = {}
        self._config_classes: Dict[str, Any] = {}

    def run(self) -> PipelineOutput:
        """Run the complete pipeline.

        Returns:
            PipelineOutput containing all results.

        Raises:
            ValueError: If no input files or steps are specified.
        """
        # Initialize state
        self._initialize_state()

        # Initialize modules
        self._initialize_modules()

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
        from ..gwas import run_association as run_gwas_association, AssocConfig, PCAConfig, GWASConfig

        # Import Ancestry module
        from ..ancestry import AncestryModel, ReferencePanel, AncestryConfig

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

        # Load legacy ancestry module (ancestry.py) via importlib
        # since genotools.ancestry is the new package directory
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
        if self._use_new_ancestry:
            ancestry_result = self._run_ancestry_prediction_new()
        else:
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

    def _run_ancestry_prediction_new(self) -> Dict[str, Any]:
        """Run ancestry prediction using the new AncestryModel.

        Uses legacy get_raw_files() for PLINK preprocessing, then
        delegates to new AncestryModel for the ML pipeline (PCA, UMAP,
        XGBoost, admixture detection). Cohort splitting reuses legacy
        split_cohort_ancestry().

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

        # Configure shared legacy Ancestry fields
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
        self._ancestry.containerized = False
        self._ancestry.singularity = False
        self._ancestry.subset = self.args.ancestry.subset_ancestry
        self._ancestry.min_samples = self.args.ancestry.min_samples

        if is_inference:
            model, predictions, ref_pca, raw = self._run_inference_mode(
                AncestryModel, model_path, actual_out
            )
        else:
            model, predictions, ref_pca, raw = self._run_training_mode(
                AncestryModel, actual_out, out_path
            )

        # --- Common post-processing (both modes) ---

        # UMAP visualization
        from ..ancestry.reducers.umap_reducer import run_umap

        projected_pca_path = f"{actual_out}_projected_new_pca.txt"
        projected_pca = pd.read_csv(projected_pca_path, sep="\t")

        umap_result: Optional[Dict[str, Any]] = None
        if ref_pca is not None:
            umap_result = run_umap(
                ref_pca=ref_pca,
                new_pca=projected_pca,
                new_labels=predictions.predictions["predicted_ancestry"],
                params=model.best_params,
            )

        # Cohort splitting (reuse legacy)
        labels_path = f"{actual_out}_umap_linearsvc_predicted_labels.txt"
        ancestry_split = self._ancestry.split_cohort_ancestry(
            labels_path=labels_path
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
        }

        # Cleanup intermediate files
        files_to_clean = [
            raw["out_paths"].get("bed", ""),
            f"{actual_out}_common_snps",
        ]
        self._ancestry.clean_up(files_to_clean)

        # Move files if ancestry-only and not full output
        if not has_qc_steps and not self.args.full_output:
            self._move_ancestry_files(out_dict)

        return out_dict

    def _run_training_mode(
        self,
        AncestryModel: type,
        actual_out: str,
        out_path: str,
    ) -> tuple:
        """Train a new AncestryModel from reference panel.

        Args:
            AncestryModel: The AncestryModel class.
            actual_out: Output path prefix (may be temp dir).
            out_path: Final output path (for saving model).

        Returns:
            Tuple of (model, predictions, ref_pca, raw).
        """
        self._ancestry.model_path = None
        self._ancestry.train = True

        # PLINK preprocessing
        raw = self._ancestry.get_raw_files()
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
            detect_admixed=True,
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
    ) -> tuple:
        """Load a pre-trained AncestryModel and predict.

        Args:
            AncestryModel: The AncestryModel class.
            model_path: Path to model directory (or legacy .pkl).
            actual_out: Output path prefix.

        Returns:
            Tuple of (model, predictions, ref_pca, raw).
        """
        model = AncestryModel.load(model_path)
        logger.info(f"Loaded pre-trained ancestry model from: {model_path}")

        # Write common SNPs to temp file for legacy get_raw_files(train=False)
        if model.common_snps is not None:
            common_snps_file = f"{actual_out}_loaded_model.common_snps"
            with open(common_snps_file, "w") as f:
                for snp in model.common_snps:
                    f.write(f"{snp}\n")
            # Legacy get_raw_files() derives .common_snps path from model_path
            self._ancestry.model_path = f"{actual_out}_loaded_model.pkl"
        else:
            # Fall back: try common_snps.txt in model directory
            model_dir = Path(model_path)
            if model_dir.is_dir():
                snps_file = model_dir / "common_snps.txt"
                if snps_file.exists():
                    common_snps_file = f"{actual_out}_loaded_model.common_snps"
                    shutil.copy2(snps_file, common_snps_file)
                    self._ancestry.model_path = f"{actual_out}_loaded_model.pkl"
                else:
                    raise FileNotFoundError(
                        f"No common_snps found in model: {model_path}"
                    )
            else:
                raise FileNotFoundError(
                    f"No common_snps in model and path is not a directory: "
                    f"{model_path}"
                )

        self._ancestry.train = False

        # PLINK preprocessing (uses common SNPs for feature alignment)
        raw = self._ancestry.get_raw_files()
        raw_geno = raw["raw_geno"]
        geno_data = raw_geno.drop(columns=["label"])
        geno_ids = raw_geno[["FID", "IID"]]

        # Skip fit(), go straight to predict()
        predictions = model.predict(
            geno_data,
            geno_ids,
            out_path=Path(actual_out),
            detect_admixed=True,
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
                geno_path=geno_path if (self.args.full_output or not self.args.ancestry.run_ancestry) else working_geno,
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
                )
            except GenoToolsError as e:
                logger.error(f"Step {step} failed: {e}")
                if not self.args.warn_only:
                    raise
                print(f"Step {step} failed but continuing (--warn): {e}")
                pass_fail[step] = PassFailRecord(
                    status=False,
                    input_path=step_input,
                    output_path=step_output,
                    error=str(e),
                )
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
        """Run a single QC step using pure function modules.

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


def run_pipeline(args: PipelineArgs, use_new_ancestry: bool = False) -> PipelineOutput:
    """Run the GenoTools pipeline.

    This is the main entry point for running the pipeline programmatically.

    Args:
        args: Validated pipeline arguments.
        use_new_ancestry: If True, use new AncestryModel instead of legacy.

    Returns:
        PipelineOutput containing all results.
    """
    runner = PipelineRunner(args, use_new_ancestry=use_new_ancestry)
    return runner.run()
