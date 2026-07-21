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

"""Command-line interface for GenoTools.

This module provides the CLI layer for the GenoTools pipeline, including:
- Typed argument parsing with validation
- Pipeline orchestration
- Result formatting and output

Example usage:
    >>> from genotools.cli import parse_args, run_pipeline
    >>> args = parse_args(["--pfile", "data/test", "--out", "results/output", "--all-sample"])
    >>> result = run_pipeline(args)
    >>> result.save("results/output.json")
"""

from .parser import (
    InputArgs,
    SampleQCArgs,
    VariantQCArgs,
    AncestryArgs,
    GWASArgs,
    OutputArgs,
    PipelineArgs,
    create_parser,
    parse_args,
)
from .runner import (
    StepResult,
    PassFailRecord,
    PipelineState,
    PipelineRunner,
    run_pipeline,
)
from .output import (
    QCMetrics,
    GWASMetrics,
    PipelineOutput,
    write_results,
    build_metrics_dataframe,
    build_gwas_dataframe,
)

__all__ = [
    # Parser
    "InputArgs",
    "SampleQCArgs",
    "VariantQCArgs",
    "AncestryArgs",
    "GWASArgs",
    "OutputArgs",
    "PipelineArgs",
    "create_parser",
    "parse_args",
    # Runner
    "StepResult",
    "PassFailRecord",
    "PipelineState",
    "PipelineRunner",
    "run_pipeline",
    # Output
    "QCMetrics",
    "GWASMetrics",
    "PipelineOutput",
    "write_results",
    "build_metrics_dataframe",
    "build_gwas_dataframe",
    # Entry points
    "main",
]


def main():
    """Main entry point for the genotools CLI."""
    import sys

    args = parse_args()
    result = run_pipeline(args)

    # Save output JSON
    out_path = args.output.out_path
    result.save(out_path.with_suffix(".json"))

    # Return appropriate exit code
    sys.exit(0 if result.success else 1)
