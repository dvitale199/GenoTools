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

"""Output formatting and writing for GenoTools pipeline.

This module handles the construction and serialization of pipeline results
to JSON and other output formats.
"""

from __future__ import annotations

import json
import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, TYPE_CHECKING

import pandas as pd

if TYPE_CHECKING:
    from .parser import PipelineArgs
    from .runner import PipelineState


@dataclass
class QCMetrics:
    """QC step metrics."""

    step: str
    count: int
    metric: str
    ancestry: str = "all"
    level: str = "sample"  # or "variant"
    passed: bool = True  # Internal field, serialized as "pass"


@dataclass
class GWASMetrics:
    """GWAS result metrics."""

    value: float
    metric: str
    ancestry: str = "all"


@dataclass
class PipelineOutput:
    """Complete pipeline output.

    This class aggregates all results from the pipeline and provides
    methods for serialization to various formats.

    Attributes:
        input_samples: DataFrame of input samples.
        ancestry_counts: Ancestry prediction counts.
        ancestry_labels: Per-sample ancestry labels.
        confusion_matrix: Ancestry model confusion matrix.
        test_accuracy: Ancestry model test accuracy.
        ref_pcs: Reference PCA coordinates.
        projected_pcs: Projected sample PCA coordinates.
        total_umap: Combined UMAP coordinates.
        ref_umap: Reference UMAP coordinates.
        new_samples_umap: New sample UMAP coordinates.
        qc_metrics: List of QC metrics.
        gwas_metrics: List of GWAS metrics.
        pruned_samples: DataFrame of pruned samples.
        related_samples: DataFrame of related samples.
        pass_fail: Step pass/fail status by ancestry.
    """

    input_samples: Optional[pd.DataFrame] = None
    ancestry_counts: Optional[pd.DataFrame] = None
    ancestry_labels: Optional[pd.DataFrame] = None
    confusion_matrix: Optional[pd.DataFrame] = None
    test_accuracy: Optional[float] = None
    ref_pcs: Optional[pd.DataFrame] = None
    projected_pcs: Optional[pd.DataFrame] = None
    total_umap: Optional[pd.DataFrame] = None
    ref_umap: Optional[pd.DataFrame] = None
    new_samples_umap: Optional[pd.DataFrame] = None
    qc_metrics: List[QCMetrics] = field(default_factory=list)
    gwas_metrics: List[GWASMetrics] = field(default_factory=list)
    pruned_samples: Optional[pd.DataFrame] = None
    related_samples: Optional[pd.DataFrame] = None
    pass_fail: Dict[str, Dict[str, Any]] = field(default_factory=dict)

    @property
    def success(self) -> bool:
        """Check if pipeline completed successfully.

        Returns True if all steps passed (or if no steps were run).
        """
        for key, value in self.pass_fail.items():
            if isinstance(value, dict):
                # Check nested pass_fail structure (e.g., "pass_fail" or "EUR_pass_fail")
                for step_name, step_info in value.items():
                    if isinstance(step_info, dict) and "status" in step_info:
                        if not step_info["status"]:
                            return False
        return True

    @classmethod
    def from_runner_state(
        cls,
        args: "PipelineArgs",
        state: "PipelineState",
    ) -> "PipelineOutput":
        """Create PipelineOutput from runner state.

        Args:
            args: Pipeline arguments.
            state: Pipeline execution state.

        Returns:
            Populated PipelineOutput instance.
        """
        output = cls()

        # Load input samples
        psam_path = f"{args.geno_path}.psam"
        if os.path.isfile(psam_path):
            output.input_samples = pd.read_csv(psam_path, sep=r"\s+")

        # Process ancestry results
        if state.ancestry_result is not None:
            output._process_ancestry_result(state.ancestry_result)

        # Process QC results
        output._process_qc_results(state, args)

        return output

    def _process_ancestry_result(self, ancestry_result: Dict[str, Any]) -> None:
        """Process ancestry prediction results.

        Args:
            ancestry_result: Raw ancestry result dictionary.
        """
        # Ancestry counts
        if "metrics" in ancestry_result and "predicted_counts" in ancestry_result["metrics"]:
            counts_df = pd.DataFrame(
                ancestry_result["metrics"]["predicted_counts"]
            ).reset_index()
            counts_df.columns = ["label", "count"]
            self.ancestry_counts = counts_df

        # Ancestry labels
        if "data" in ancestry_result:
            data = ancestry_result["data"]

            if "predict_data" in data and "ids" in data["predict_data"]:
                self.ancestry_labels = pd.DataFrame(data["predict_data"]["ids"])

            # Confusion matrix
            if "confusion_matrix" in data and "label_encoder" in data:
                le = data["label_encoder"]
                cm = pd.DataFrame(data["confusion_matrix"])
                labels = le.inverse_transform([i for i in range(len(cm))])
                cm.columns = labels
                cm.index = labels
                self.confusion_matrix = cm

            # PCA and UMAP data
            if "ref_pcs" in data:
                self.ref_pcs = data["ref_pcs"]
            if "projected_pcs" in data:
                self.projected_pcs = data["projected_pcs"]
            if "total_umap" in data:
                self.total_umap = data["total_umap"]
            if "ref_umap" in data:
                self.ref_umap = data["ref_umap"]
            if "new_samples_umap" in data:
                self.new_samples_umap = data["new_samples_umap"]

        # Test accuracy
        if "metrics" in ancestry_result and "test_accuracy" in ancestry_result["metrics"]:
            self.test_accuracy = ancestry_result["metrics"]["test_accuracy"]

    def _process_qc_results(
        self,
        state: "PipelineState",
        args: "PipelineArgs",
    ) -> None:
        """Process QC step results.

        Args:
            state: Pipeline execution state.
            args: Pipeline arguments.
        """
        metrics_list: List[QCMetrics] = []
        pruned_dfs: List[pd.DataFrame] = []
        related_dfs: List[pd.DataFrame] = []

        sample_steps = ["callrate", "sex", "het", "related"]
        variant_steps = ["case_control", "haplotype", "hwe", "geno", "ld"]

        # Process ancestry-specific results
        if hasattr(state, "ancestry_results"):
            ancestry_results = getattr(state, "ancestry_results")
            for ancestry, result_dict in ancestry_results.items():
                self._extract_qc_metrics(
                    result_dict,
                    ancestry,
                    sample_steps,
                    variant_steps,
                    metrics_list,
                    pruned_dfs,
                    related_dfs,
                    args,
                )
                self.pass_fail[f"{ancestry}_pass_fail"] = result_dict.get("pass_fail", {})
        else:
            # Single run (no ancestry split)
            result_dict = {
                k: {"pass": v.passed, "step": v.step, "metrics": v.metrics, "output": v.output}
                for k, v in state.step_results.items()
            }
            self._extract_qc_metrics(
                result_dict,
                "all",
                sample_steps,
                variant_steps,
                metrics_list,
                pruned_dfs,
                related_dfs,
                args,
            )
            self.pass_fail["pass_fail"] = {
                k: {"status": v.status, "input": v.input_path, "output": v.output_path}
                for k, v in state.pass_fail.items()
            }

        self.qc_metrics = metrics_list

        # Combine pruned samples
        if pruned_dfs:
            combined = pd.concat(pruned_dfs, ignore_index=True)
            combined = combined.drop_duplicates(subset=["#FID", "IID"], ignore_index=True)
            combined = combined.rename({"#FID": "FID"}, axis=1)
            self.pruned_samples = combined

        # Combine related samples
        if related_dfs:
            self.related_samples = pd.concat(related_dfs, ignore_index=True)

    def _extract_qc_metrics(
        self,
        result_dict: Dict[str, Any],
        ancestry: str,
        sample_steps: List[str],
        variant_steps: List[str],
        metrics_list: List[QCMetrics],
        pruned_dfs: List[pd.DataFrame],
        related_dfs: List[pd.DataFrame],
        args: "PipelineArgs",
    ) -> None:
        """Extract metrics from a QC result dictionary.

        Args:
            result_dict: Step result dictionary.
            ancestry: Ancestry label.
            sample_steps: List of sample-level steps.
            variant_steps: List of variant-level steps.
            metrics_list: List to append metrics to.
            pruned_dfs: List to append pruned sample DataFrames.
            related_dfs: List to append related sample DataFrames.
            args: Pipeline arguments.
        """
        for step in sample_steps + variant_steps:
            if step not in result_dict:
                continue

            step_result = result_dict[step]
            level = "sample" if step in sample_steps else "variant"
            passed = step_result.get("pass", False)

            # Extract metrics
            if "metrics" in step_result:
                for metric, value in step_result["metrics"].items():
                    metrics_list.append(
                        QCMetrics(
                            step=step_result.get("step", step),
                            count=value,
                            metric=metric,
                            ancestry=ancestry,
                            level=level,
                            passed=passed,
                        )
                    )

            # Extract pruned samples
            if step in sample_steps and "output" in step_result:
                samplefile = step_result["output"].get("pruned_samples")
                if samplefile and os.path.isfile(samplefile):
                    pruned = pd.read_csv(samplefile, sep="\t")
                    if len(pruned) > 0:
                        pruned["step"] = step
                        pruned_dfs.append(pruned[["#FID", "IID", "step"]])

            # Extract related samples
            if step == "related" and "output" in step_result:
                relatedfile = step_result["output"].get("related_samples")
                if relatedfile and os.path.isfile(relatedfile):
                    related = pd.read_csv(relatedfile, sep=",")
                    if len(related) > 0:
                        # Write per-ancestry related file
                        if ancestry == "all":
                            related_out = f"{args.out_path}.related"
                        else:
                            related_out = f"{args.out_path}_{ancestry}.related"
                        related.to_csv(related_out, index=False)
                        related["ancestry"] = ancestry
                        related_dfs.append(related)

        # Extract GWAS metrics
        if "assoc" in result_dict and "gwas" in result_dict["assoc"]:
            gwas_result = result_dict["assoc"]["gwas"]
            if "metrics" in gwas_result:
                for metric, value in gwas_result["metrics"].items():
                    self.gwas_metrics.append(
                        GWASMetrics(
                            value=value,
                            metric=metric,
                            ancestry=ancestry,
                        )
                    )

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for JSON serialization.

        Returns:
            Dictionary suitable for JSON serialization.
        """
        result: Dict[str, Any] = {}

        # Input samples
        if self.input_samples is not None:
            result["input_samples"] = self.input_samples.to_dict()

        # Ancestry results
        if self.ancestry_counts is not None:
            result["ancestry_counts"] = self.ancestry_counts.to_dict()
        if self.ancestry_labels is not None:
            result["ancestry_labels"] = self.ancestry_labels.to_dict()
        if self.confusion_matrix is not None:
            result["confusion_matrix"] = self.confusion_matrix.to_dict()
        if self.test_accuracy is not None:
            result["test_accuracy"] = self.test_accuracy

        # PCA/UMAP data
        if self.ref_pcs is not None:
            result["ref_pcs"] = self.ref_pcs.to_dict()
        if self.projected_pcs is not None:
            result["projected_pcs"] = self.projected_pcs.to_dict()
        if self.total_umap is not None:
            result["total_umap"] = self.total_umap.to_dict()
        if self.ref_umap is not None:
            result["ref_umap"] = self.ref_umap.to_dict()
        if self.new_samples_umap is not None:
            result["new_samples_umap"] = self.new_samples_umap.to_dict()

        # QC metrics
        if self.qc_metrics:
            qc_df = pd.DataFrame([
                {
                    "step": m.step,
                    "count": m.count,
                    "metric": m.metric,
                    "ancestry": m.ancestry,
                    "level": m.level,
                    "pass": m.passed,
                }
                for m in self.qc_metrics
            ])
            result["QC"] = qc_df.to_dict()

        # GWAS metrics
        if self.gwas_metrics:
            gwas_df = pd.DataFrame([
                {
                    "value": m.value,
                    "metric": m.metric,
                    "ancestry": m.ancestry,
                }
                for m in self.gwas_metrics
            ])
            result["GWAS"] = gwas_df.to_dict()

        # Pruned samples
        if self.pruned_samples is not None and not self.pruned_samples.empty:
            result["pruned_samples"] = self.pruned_samples.to_dict()

        # Related samples
        if self.related_samples is not None and not self.related_samples.empty:
            result["related_samples"] = self.related_samples.to_dict()

        # Pass/fail status
        result.update(self.pass_fail)

        return result

    def to_json(self, indent: Optional[int] = None) -> str:
        """Convert to JSON string.

        Args:
            indent: JSON indentation level.

        Returns:
            JSON string representation.
        """
        return json.dumps(self.to_dict(), indent=indent)

    def save(self, path: Path) -> None:
        """Save results to JSON file.

        Args:
            path: Output file path.
        """
        with open(path, "w") as f:
            json.dump(self.to_dict(), f)


def write_results(output: PipelineOutput, out_path: Path) -> None:
    """Write pipeline results to output file.

    Args:
        output: Pipeline output to write.
        out_path: Output path prefix.
    """
    json_path = Path(f"{out_path}.json")
    output.save(json_path)


def build_metrics_dataframe(metrics: List[QCMetrics]) -> pd.DataFrame:
    """Build metrics DataFrame from list of QCMetrics.

    Args:
        metrics: List of QC metrics.

    Returns:
        DataFrame with metrics.
    """
    if not metrics:
        return pd.DataFrame()

    return pd.DataFrame([
        {
            "step": m.step,
            "count": m.count,
            "metric": m.metric,
            "ancestry": m.ancestry,
            "level": m.level,
            "pass": m.passed,
        }
        for m in metrics
    ])


def build_gwas_dataframe(metrics: List[GWASMetrics]) -> pd.DataFrame:
    """Build GWAS metrics DataFrame.

    Args:
        metrics: List of GWAS metrics.

    Returns:
        DataFrame with GWAS metrics.
    """
    if not metrics:
        return pd.DataFrame()

    return pd.DataFrame([
        {
            "value": m.value,
            "metric": m.metric,
            "ancestry": m.ancestry,
        }
        for m in metrics
    ])
