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


# The report identity of each QC step: the step name it reports under, and the
# metric keys it emits. A step that never produced data - it failed, or a
# data-driven decision skipped it - still has to appear in the report, and this
# is where those names come from. Keep in sync with each step module's
# ``step_name`` / metrics; ``tests/unit/test_qc/test_results.py`` asserts it.
STEP_REPORT: dict[str, tuple[str, tuple[str, ...]]] = {
    # Sample-level steps
    "callrate": ("callrate_prune", ("outlier_count",)),
    "sex": ("sex_prune", ("outlier_count",)),
    "het": ("het_prune", ("outlier_count",)),
    "related": ("related_prune", ("related_count", "duplicated_count")),
    # Variant-level steps
    "geno": ("geno_prune", ("geno_removed_count",)),
    "case_control": ("case_control_missingness_prune", ("mis_removed_count",)),
    "haplotype": ("haplotype_prune", ("haplotype_removed_count",)),
    "hwe": ("hwe_prune", ("hwe_removed_count",)),
    "ld": ("ld_prune", ("ld_removed_count",)),
}


def unrun_result(
    step: str,
    outcome: str,
    reason: str,
    plink_out: str,
) -> Optional[dict[str, Any]]:
    """Build the result dict for a step that ran no analysis.

    Covers both outcomes that produce no data: ``"fail"`` (the step raised) and
    ``"skipped"`` (a data-driven decision ruled it out). Counts are zeroed, the
    same shape the pre-refactor CLI reported for a failed step, so the row is
    present rather than silently absent.

    Args:
        step: Pipeline step key (e.g. ``"het"``).
        outcome: Either ``"fail"`` or ``"skipped"``.
        reason: Human-readable explanation, surfaced in the report.
        plink_out: Path the step would have written.

    Returns:
        Result dict, or None for steps that carry no QC metrics (``assoc``,
        ``kinship_check``), which the report does not list either way.
    """
    entry = STEP_REPORT.get(step)
    if entry is None:
        return None

    step_name, metric_keys = entry
    return {
        "pass": False,
        "outcome": outcome,
        "reason": reason,
        "step": step_name,
        "metrics": {key: 0 for key in metric_keys},
        "output": {"pruned_samples": None, "plink_out": plink_out},
    }


@dataclass(frozen=True)
class FilterResult:
    """Result of a single QC filter step.

    This is an immutable dataclass that captures the output of a QC step,
    including the filtered genotype data, metrics, and any outlier files.

    Attributes:
        output: New GenotypeData pointing to filtered output files.
        samples_removed: Number of samples filtered out by this step.
        variants_removed: Number of variants filtered out by this step.
        metrics: Step-specific counts, and only counts - the JSON report
            renders each one as a row in a long (step, metric, pruned_count)
            table. Descriptive or derived values go in ``parameters``.
        parameters: Values the step *resolved* at runtime that the caller could
            not know from the config alone - derived bounds, the mode actually
            taken. Reported under the JSON "parameters" section with
            source="resolved", never as a metric.
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
    parameters: dict[str, Any] = field(default_factory=dict)
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
            "outcome": "pass",
            "reason": None,
            "step": self.step_name,
            "metrics": metrics_dict,
            # Separate from metrics on purpose: these are descriptive, and the
            # metrics table's value column means "how many were pruned".
            "parameters": dict(self.parameters),
            "output": output_dict,
        }
