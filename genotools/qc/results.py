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

"""QC step result dataclasses."""

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Optional

from genotools.core.genotypes import GenotypeData


@dataclass(frozen=True)
class FilterResult:
    """Result of a single QC filter step.

    This is an immutable dataclass that captures the output of a QC step,
    including the filtered genotype data, metrics, and any outlier files.

    Attributes:
        output: New GenotypeData pointing to filtered output files.
        samples_removed: Number of samples filtered out by this step.
        variants_removed: Number of variants filtered out by this step.
        metrics: Step-specific metrics (e.g., thresholds used, statistics).
        log: PLINK/tool log output for debugging.
        pruned_samples_file: Path to .outliers file with filtered sample IDs,
            or None if no samples were filtered.
        related_samples_file: Path to .related file with kinship pair info,
            or None if not a relatedness step.
        step_name: Name of the QC step (e.g., 'callrate_prune').
    """

    output: GenotypeData
    samples_removed: int
    variants_removed: int
    metrics: dict[str, Any] = field(default_factory=dict)
    log: str = ""
    pruned_samples_file: Optional[Path] = None
    related_samples_file: Optional[Path] = None
    step_name: str = "unknown"

    def to_dict(self) -> dict[str, Any]:
        """Convert to legacy dictionary format for backward compatibility.

        Returns the standard format expected by existing pipeline code:
        {
            'pass': bool,
            'step': str,
            'metrics': {'outlier_count': int, ...},
            'output': {'pruned_samples': str|None, 'plink_out': str}
        }

        Returns:
            Dictionary in legacy format matching SampleQC/VariantQC output.
        """
        output_dict = {
            "pruned_samples": (
                str(self.pruned_samples_file)
                if self.pruned_samples_file
                else None
            ),
            "plink_out": str(self.output.path),
        }
        # Include related_samples for relatedness step (backward compatibility)
        if self.related_samples_file is not None:
            output_dict["related_samples"] = str(self.related_samples_file)

        # Build metrics dict - relatedness step uses related_count/duplicated_count
        # instead of outlier_count for backward compatibility
        is_relatedness_step = "related" in self.step_name
        if is_relatedness_step:
            # Relatedness step: only include related_count and duplicated_count
            metrics_dict = dict(self.metrics)
        else:
            # Other steps: include outlier_count
            metrics_dict = {
                "outlier_count": self.samples_removed + self.variants_removed,
                **self.metrics,
            }

        return {
            "pass": True,  # If FilterResult exists, step passed
            "step": self.step_name,
            "metrics": metrics_dict,
            "output": output_dict,
        }


@dataclass(frozen=True)
class QCResult:
    """Result of running the full QC pipeline.

    This captures the complete result of running multiple QC steps,
    including the final filtered data and results from each step.

    Attributes:
        input: Original input GenotypeData before any filtering.
        output: Final filtered GenotypeData after all steps.
        step_results: List of (step_name, FilterResult) tuples in execution order.
    """

    input: GenotypeData
    output: GenotypeData
    step_results: list[tuple[str, FilterResult]]

    @property
    def total_samples_removed(self) -> int:
        """Total number of samples removed across all steps."""
        return sum(r.samples_removed for _, r in self.step_results)

    @property
    def total_variants_removed(self) -> int:
        """Total number of variants removed across all steps."""
        return sum(r.variants_removed for _, r in self.step_results)

    def to_legacy_dict(self) -> dict[str, dict[str, Any]]:
        """Convert to legacy pass_fail dictionary format.

        Returns:
            Dictionary mapping step names to their result dictionaries,
            matching the format used by the old pipeline code.
        """
        return {name: result.to_dict() for name, result in self.step_results}

    def get_step_result(self, step_name: str) -> Optional[FilterResult]:
        """Get the result for a specific step by name.

        Args:
            step_name: Name of the step to retrieve.

        Returns:
            FilterResult for the step, or None if not found.
        """
        for name, result in self.step_results:
            if name == step_name:
                return result
        return None
