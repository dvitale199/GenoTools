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


def _render_plots(out_path) -> None:
    """Render `{out}_plots/` from the report a run just wrote.

    Plots are the least load-bearing output a run has, so nothing here may fail
    it: by now the report and the pfiles are already on disk, and a broken
    figure is a warning, not a lost run. `genotools-plot` reruns this alone.
    """
    from pathlib import Path

    from genotools.core.logging import get_logger

    logger = get_logger(__name__)
    try:
        from genotools.plotting import render_plots

        result = render_plots(Path(f"{out_path}.json"))
    except Exception as exc:
        logger.warning(
            f"Plots not rendered ({type(exc).__name__}: {exc}); the run's "
            f"results are unaffected. Retry with: genotools-plot --report {out_path}.json"
        )
        return
    logger.info(f"Plots: {result.html}")


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

        write_results(result, args.output.out_path)
        if args.output.plots:
            _render_plots(args.output.out_path)
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
