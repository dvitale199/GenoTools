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

"""`genotools-plot`: render a finished run's plots from its JSON report."""

import argparse
import sys
from pathlib import Path
from typing import List, Optional

from genotools.core.logging import setup_logging


def create_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="genotools-plot",
        description=(
            "Render plots from a GenoTools JSON report. Writes PNGs to "
            "{out-dir}/static/ and a self-contained {out-dir}/report.html."
        ),
    )
    parser.add_argument(
        "--report",
        type=Path,
        required=True,
        help="The {out}.json report a genotools run wrote",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=None,
        help="Directory to write plots to (default: {out}_plots beside the report)",
    )
    parser.add_argument(
        "--debug",
        action="store_true",
        help="Show the full traceback on error",
    )
    return parser


def main(argv: Optional[List[str]] = None) -> None:
    """Entry point for `genotools-plot`.

    Mirrors `genotools`: a bad report path or a file that is not a report is a
    one-line `ERROR:` with exit 1, and anything else keeps its traceback.
    """
    args = create_parser().parse_args(argv)
    setup_logging(level="DEBUG" if args.debug else "INFO")

    from genotools.plotting.render import render_plots

    try:
        result = render_plots(args.report, args.out_dir)
    except (FileNotFoundError, ValueError) as e:
        if args.debug:
            raise
        print(f"ERROR: {e}", file=sys.stderr)
        sys.exit(1)

    print(result.html)
