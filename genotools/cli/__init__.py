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
    """Main entry point for the genotools CLI.

    Failures a user can act on are reported as a one-line ``ERROR:`` message on
    stderr with exit 1, not a Python traceback -- the detail belongs in the
    consolidated run log. That covers ``GenoToolsError`` (bad input, a failing
    external tool, a blocked output prefix), ``FileNotFoundError`` (the
    codebase-wide signal for a missing file), and the ``ValueError``/``TypeError``
    raised by config validation for out-of-range arguments.

    Any other exception is a bug and keeps its traceback. ``--debug`` re-raises
    everything, so the traceback is always one flag away.
    """
    import sys

    from genotools.core.exceptions import GenoToolsError

    # argparse reports its own errors and exits via SystemExit, which passes
    # through untouched. But the config dataclasses validate in __post_init__,
    # inside parse_args, so their ValueError/TypeError has to be caught here too
    # -- otherwise a bad argument prints a traceback rather than the one-line
    # ERROR promised above. --debug is read from argv because args may not exist
    # yet at that point.
    debug = "--debug" in sys.argv

    try:
        args = parse_args()
        debug = args.output.debug

        result = run_pipeline(args)

        # Save output JSON
        out_path = args.output.out_path
        result.save(out_path.with_suffix(".json"))
    except (GenoToolsError, FileNotFoundError, ValueError, TypeError) as e:
        if debug:
            raise
        print(f"ERROR: {e}", file=sys.stderr)
        print("Re-run with --debug for the full traceback.", file=sys.stderr)
        sys.exit(1)
    except KeyboardInterrupt:
        print("Interrupted.", file=sys.stderr)
        sys.exit(130)

    # Return appropriate exit code
    sys.exit(0 if result.success else 1)
