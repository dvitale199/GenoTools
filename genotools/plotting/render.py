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

"""Render every plot a report supports into one output directory.

    {out}_plots/
        static/<plot>.png
        report.html          self-contained: plotly.js is inlined

A plot whose data the run did not produce is skipped with a note, and a plot
that raises is skipped with a warning: plots are the least load-bearing output
GenoTools has, and one bad figure must not cost the rest.
"""

import html
from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Optional

from genotools.core.logging import get_logger
from genotools.plotting.plots import ALL_PLOTS, Plot
from genotools.plotting.report import Report, load_report

logger = get_logger(__name__)

# Every plot, in report order. Plots live in `plots/<area>.py`; adding one
# never touches this file.
PLOTS: List[Plot] = ALL_PLOTS

__all__ = ["PLOTS", "Plot", "RenderResult", "default_out_dir", "render_plots"]


@dataclass
class RenderResult:
    """What a render wrote, and what it left out and why."""

    out_dir: Path
    pngs: List[Path] = field(default_factory=list)
    html: Optional[Path] = None
    skipped: List[str] = field(default_factory=list)


def default_out_dir(report_path: Path) -> Path:
    """`{out}.json` -> `{out}_plots`.

    Appends rather than using `with_suffix`/`.stem`, which would rewrite a
    dotted prefix like `cohort.v2` (#296).
    """
    name = str(report_path)
    if name.endswith(".json"):
        name = name[: -len(".json")]
    return Path(f"{name}_plots")


def _pyplot():
    """Pyplot against a non-interactive backend, or None if it will not load."""
    try:
        import matplotlib

        matplotlib.use("Agg", force=False)
        import matplotlib.pyplot as plt

        return plt
    except Exception as exc:  # pragma: no cover - environment-specific
        logger.warning(f"PNG plots skipped: matplotlib unavailable ({exc})")
        return None


_PAGE = """<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>GenoTools Report</title>
<style>
  :root {{
    --page: #f9f9f7; --surface: #fcfcfb; --ink: #0b0b0b;
    --ink-2: #52514e; --rule: #e6e5e0;
  }}
  body {{
    margin: 0; background: var(--page); color: var(--ink);
    font: 15px/1.5 system-ui, -apple-system, "Segoe UI", sans-serif;
  }}
  main {{ max-width: 1200px; margin: 0 auto; padding: 24px 16px 48px; }}
  h1 {{ font-size: 22px; margin: 0 0 4px; }}
  h2 {{ font-size: 17px; margin: 0 0 4px; }}
  .meta, .desc {{ color: var(--ink-2); font-size: 13px; margin: 0 0 12px; }}
  .meta code {{ word-break: break-all; }}
  section {{
    background: var(--surface); border: 1px solid var(--rule);
    border-radius: 8px; padding: 16px; margin-top: 16px; overflow-x: auto;
  }}
  ul {{ color: var(--ink-2); font-size: 13px; }}
</style>
</head>
<body>
<main>
<h1>GenoTools report</h1>
<p class="meta">{meta}</p>
{sections}
{skipped}
</main>
</body>
</html>
"""


def _write_html(report: Report, sections: List[str], skipped: List[str], path: Path) -> Path:
    run_info = report.raw.get("run_info", {})
    meta = [f"Report: <code>{html.escape(str(report.path))}</code>"]
    if run_info.get("version"):
        meta.append(f"GenoTools {html.escape(str(run_info['version']))}")
    if run_info.get("invocation"):
        meta.append(f"<code>{html.escape(run_info['invocation'])}</code>")
    skipped_html = ""
    if skipped:
        items = "".join(f"<li>{html.escape(s)}</li>" for s in skipped)
        skipped_html = f"<section><h2>Not plotted</h2><ul>{items}</ul></section>"
    path.write_text(
        _PAGE.format(
            meta="<br>".join(meta),
            sections="\n".join(sections),
            skipped=skipped_html,
        )
    )
    return path


def render_plots(
    report_path: Path,
    out_dir: Optional[Path] = None,
    plots: Optional[List[Plot]] = None,
) -> RenderResult:
    """Render every plot the report supports.

    Args:
        report_path: The `{out}.json` a run wrote.
        out_dir: Where to write. Defaults to `{out}_plots` beside the report.
        plots: Plots to attempt; all of `PLOTS` by default.

    Returns:
        What was written and what was skipped.

    Raises:
        FileNotFoundError: The report does not exist.
        ValueError: The report is not a GenoTools JSON report.
    """
    report = load_report(Path(report_path))
    out_dir = Path(out_dir) if out_dir is not None else default_out_dir(report.path)
    static_dir = out_dir / "static"
    static_dir.mkdir(parents=True, exist_ok=True)

    result = RenderResult(out_dir=out_dir)
    plt = _pyplot()
    sections: List[str] = []
    include_js = True  # inline plotly.js once, in the first figure

    for plot in plots if plots is not None else PLOTS:
        if not plot.available(report):
            result.skipped.append(f"{plot.title}: no data for it in this report")
            continue
        try:
            if plt is not None:
                result.pngs.append(plot.png(report, plt, static_dir / f"{plot.name}.png"))
            figure_html = plot.plotly(report).to_html(
                full_html=False, include_plotlyjs=include_js, config={"displaylogo": False}
            )
            include_js = False
        except Exception as exc:
            logger.warning(f"Plot {plot.name} skipped: {type(exc).__name__}: {exc}")
            result.skipped.append(f"{plot.title}: failed to render ({exc})")
            continue
        sections.append(
            f'<section id="{plot.name}"><h2>{html.escape(plot.title)}</h2>'
            f'<p class="desc">{html.escape(plot.description)}</p>{figure_html}</section>'
        )

    result.html = _write_html(report, sections, result.skipped, out_dir / "report.html")
    logger.info(
        f"Plots written to {out_dir} ({len(result.pngs)} PNG, "
        f"{len(result.skipped)} skipped)"
    )
    return result
