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

"""QC pipeline orchestration.

Composes QC steps into a sequential pipeline with configurable
error handling.
"""

import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable, Protocol

from genotools.core.exceptions import QCError
from genotools.core.genotypes import GenotypeData
from genotools.core.logging import step_context
from genotools.qc.results import FilterResult, QCResult

logger = logging.getLogger(__name__)


class QCStepProtocol(Protocol):
    """Protocol for QC step functions.

    All QC step functions must accept:
        - data: GenotypeData (input)
        - config: Step-specific configuration
        - out_path: Path for output files

    And return a FilterResult.
    """

    def __call__(
        self,
        data: GenotypeData,
        config: Any,
        out_path: Path,
    ) -> FilterResult: ...


# Type alias for step tuple
QCStep = tuple[str, Callable[..., FilterResult], Any]


@dataclass
class QCPipeline:
    """Composes QC steps into a sequential pipeline.

    The pipeline executes steps in order, passing the output of each
    step as input to the next. Handles error recovery based on
    warn_only setting.

    Attributes:
        steps: List of (name, step_function, config) tuples.
        warn_only: If True, continue on step failure (--warn flag behavior).
            Failed steps are skipped and the previous data is passed forward.

    Example:
        from genotools.qc import (
            QCPipeline,
            filter_callrate,
            filter_sex,
            CallrateConfig,
            SexConfig,
        )

        pipeline = QCPipeline(
            steps=[
                ("callrate", filter_callrate, CallrateConfig(mind=0.05)),
                ("sex", filter_sex, SexConfig()),
            ],
            warn_only=False,
        )
        result = pipeline.run(data, Path("output"))
    """

    steps: list[QCStep]
    warn_only: bool = False

    def run(self, data: GenotypeData, out_dir: Path) -> QCResult:
        """Execute all steps sequentially.

        Args:
            data: Input genotype data.
            out_dir: Directory for step outputs. Each step creates a
                subdirectory, with the final output at out_dir/output.

        Returns:
            QCResult with all step results and final filtered data.

        Raises:
            QCError: If a step fails and warn_only is False, or if
                all samples/variants are removed.
        """
        out_dir.mkdir(parents=True, exist_ok=True)

        current = data
        step_results: list[tuple[str, FilterResult]] = []

        logger.info(
            f"Starting QC pipeline with {len(self.steps)} steps "
            f"on {data.sample_count} samples, {data.variant_count} variants"
        )

        for i, (name, step_fn, config) in enumerate(self.steps):
            is_last_step = i == len(self.steps) - 1

            with step_context(name):
                logger.info(f"Starting step {i + 1}/{len(self.steps)}: {name}")

                # Determine output path
                if is_last_step:
                    step_out = out_dir / "output"
                else:
                    step_out = out_dir / name / "output"
                step_out.parent.mkdir(parents=True, exist_ok=True)

                try:
                    result = step_fn(current, config, step_out)
                    step_results.append((name, result))

                    logger.info(
                        f"Completed {name}: "
                        f"-{result.samples_removed} samples, "
                        f"-{result.variants_removed} variants "
                        f"({result.output.sample_count} samples, "
                        f"{result.output.variant_count} variants remaining)"
                    )

                    current = result.output

                    # Check for empty data
                    if current.sample_count == 0:
                        raise QCError(f"All samples removed at step: {name}")
                    if current.variant_count == 0:
                        raise QCError(f"All variants removed at step: {name}")

                except QCError as e:
                    if self.warn_only:
                        logger.warning(f"Step {name} failed (--warn mode): {e}")
                        # Continue with previous data
                        continue
                    else:
                        raise

        total_samples_removed = data.sample_count - current.sample_count
        total_variants_removed = data.variant_count - current.variant_count

        logger.info(
            f"QC pipeline complete: "
            f"removed {total_samples_removed} samples, {total_variants_removed} variants "
            f"({current.sample_count} samples, {current.variant_count} variants remaining)"
        )

        return QCResult(
            input=data,
            output=current,
            step_results=step_results,
        )

    def add_step(
        self,
        name: str,
        step_fn: Callable[..., FilterResult],
        config: Any,
    ) -> "QCPipeline":
        """Add a step to the pipeline.

        Returns a new QCPipeline with the additional step.
        Does not modify the original pipeline.

        Args:
            name: Step name for logging and result tracking.
            step_fn: Step function matching QCStepProtocol.
            config: Configuration for the step.

        Returns:
            New QCPipeline with the step added.
        """
        return QCPipeline(
            steps=self.steps + [(name, step_fn, config)],
            warn_only=self.warn_only,
        )

    @classmethod
    def sample_qc(
        cls,
        callrate_config: Any = None,
        sex_config: Any = None,
        het_config: Any = None,
        related_config: Any = None,
        warn_only: bool = False,
    ) -> "QCPipeline":
        """Create a standard sample QC pipeline.

        Args:
            callrate_config: CallrateConfig, or None to skip.
            sex_config: SexConfig, or None to skip.
            het_config: HetConfig, or None to skip.
            related_config: RelatedConfig, or None to skip.
            warn_only: If True, continue on step failure.

        Returns:
            QCPipeline configured for sample QC.
        """
        from genotools.qc.steps import (
            filter_callrate,
            filter_heterozygosity,
            filter_relatedness,
            filter_sex,
        )

        steps: list[QCStep] = []

        if callrate_config is not None:
            steps.append(("callrate", filter_callrate, callrate_config))
        if sex_config is not None:
            steps.append(("sex", filter_sex, sex_config))
        if het_config is not None:
            steps.append(("het", filter_heterozygosity, het_config))
        if related_config is not None:
            steps.append(("related", filter_relatedness, related_config))

        return cls(steps=steps, warn_only=warn_only)

    @classmethod
    def variant_qc(
        cls,
        geno_config: Any = None,
        case_control_config: Any = None,
        haplotype_config: Any = None,
        hwe_config: Any = None,
        ld_config: Any = None,
        warn_only: bool = False,
    ) -> "QCPipeline":
        """Create a standard variant QC pipeline.

        Args:
            geno_config: GenoConfig, or None to skip.
            case_control_config: CaseControlConfig, or None to skip.
            haplotype_config: HaplotypeConfig, or None to skip.
            hwe_config: HWEConfig, or None to skip.
            ld_config: LDConfig, or None to skip.
            warn_only: If True, continue on step failure.

        Returns:
            QCPipeline configured for variant QC.
        """
        from genotools.qc.steps import (
            filter_case_control,
            filter_haplotype,
            filter_hwe,
            filter_variant_missingness,
            prune_ld,
        )

        steps: list[QCStep] = []

        if geno_config is not None:
            steps.append(("geno", filter_variant_missingness, geno_config))
        if case_control_config is not None:
            steps.append(("case_control", filter_case_control, case_control_config))
        if haplotype_config is not None:
            steps.append(("haplotype", filter_haplotype, haplotype_config))
        if hwe_config is not None:
            steps.append(("hwe", filter_hwe, hwe_config))
        if ld_config is not None:
            steps.append(("ld", prune_ld, ld_config))

        return cls(steps=steps, warn_only=warn_only)
