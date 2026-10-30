"""Rendering a report into {out}_plots/, and the commands that do it."""

import ast
import sys
from pathlib import Path

import pytest

import genotools
from genotools.plotting import render_plots
from genotools.plotting.render import PLOTS, Plot, default_out_dir

from .conftest import qc_row, write_report

PLOTTING_DIR = Path(genotools.__file__).resolve().parent / "plotting"


class TestRenderPlots:
    def test_writes_every_supported_plot(self, ancestry_report: Path) -> None:
        result = render_plots(ancestry_report)
        assert result.out_dir == ancestry_report.parent / "cohort_plots"
        assert sorted(p.name for p in result.pngs) == [
            "sample_attrition.png", "step_outcomes.png", "variant_attrition.png"
        ]
        assert all(p.parent == result.out_dir / "static" for p in result.pngs)
        assert all(p.stat().st_size > 0 for p in result.pngs)
        assert result.skipped == []

    def test_html_is_self_contained(self, ancestry_report: Path) -> None:
        html = render_plots(ancestry_report).html.read_text()
        assert html.startswith("<!doctype html>")
        # plotly.js inlined (several MB), never fetched from a CDN.
        assert "<script src=" not in html
        assert len(html) > 1_000_000
        assert html.count("* plotly.js v") == 1
        for plot in PLOTS:
            assert f'id="{plot.name}"' in html

    def test_plot_without_data_is_listed_not_drawn(self, tmp_path: Path) -> None:
        """A sample-only run has nothing for the variant plot."""
        report = write_report(
            tmp_path / "s.json", [qc_row("callrate_prune", 1, "outlier_count")], n_input=10
        )
        result = render_plots(report)
        assert "variant_attrition.png" not in {p.name for p in result.pngs}
        assert any(s.startswith("Variant attrition") for s in result.skipped)
        assert "Not plotted" in result.html.read_text()

    def test_one_failing_plot_does_not_cost_the_rest(self, ancestry_report: Path) -> None:
        def boom(*_):
            raise RuntimeError("bad figure")

        broken = Plot("broken", "Broken", "", lambda r: True, boom, boom)
        result = render_plots(ancestry_report, plots=[broken] + PLOTS)
        assert len(result.pngs) == len(PLOTS)
        assert any("bad figure" in s for s in result.skipped)

    def test_out_dir_override(self, ancestry_report: Path, tmp_path: Path) -> None:
        result = render_plots(ancestry_report, tmp_path / "figs")
        assert result.html == tmp_path / "figs" / "report.html"

    def test_missing_report_raises(self, tmp_path: Path) -> None:
        with pytest.raises(FileNotFoundError):
            render_plots(tmp_path / "nope.json")


def test_default_out_dir_keeps_dotted_prefix() -> None:
    """#296: a prefix is not a filename."""
    assert default_out_dir(Path("runs/cohort.v2.json")) == Path("runs/cohort.v2_plots")


def test_plotting_imports_no_pipeline_package() -> None:
    """Plotting reads the report; it must not reach back into the pipeline."""
    banned = ("genotools.cli", "genotools.qc", "genotools.gwas", "genotools.ancestry")
    paths = sorted(PLOTTING_DIR.rglob("*.py"))
    assert PLOTTING_DIR / "plots" / "sample_qc.py" in paths
    for path in paths:
        tree = ast.parse(path.read_text())
        for node in ast.walk(tree):
            if isinstance(node, ast.ImportFrom):
                assert node.level == 0, f"{path.name}: use absolute imports"
                names = [node.module or ""]
            elif isinstance(node, ast.Import):
                names = [alias.name for alias in node.names]
            else:
                continue
            for name in names:
                assert not name.startswith(banned), f"{path.name} imports {name}"


def test_render_order_follows_the_areas() -> None:
    """Each area owns its plots; render.py only concatenates them."""
    from genotools.plotting.plots import ancestry, gwas, qc_summary, sample_qc, variant_qc

    expected = [
        *qc_summary.PLOTS, *sample_qc.PLOTS, *variant_qc.PLOTS,
        *ancestry.PLOTS, *gwas.PLOTS,
    ]
    assert PLOTS == expected
    assert len({p.name for p in PLOTS}) == len(PLOTS), "plot names must be unique"


def test_every_plot_renders_the_golden_reports(golden_qc_report: Path, tmp_path: Path) -> None:
    """Real reports with real plot tables: nothing available may fail to draw."""
    result = render_plots(golden_qc_report, tmp_path / "figs")
    assert not [s for s in result.skipped if "failed to render" in s], result.skipped
    assert result.pngs


class TestGenotoolsPlotCommand:
    def test_renders_and_prints_html_path(self, ancestry_report: Path, capsys) -> None:
        from genotools.plotting.cli import main

        main(["--report", str(ancestry_report)])
        out = capsys.readouterr().out.strip()
        assert out == str(ancestry_report.parent / "cohort_plots" / "report.html")

    def test_missing_report_is_one_error_line(self, tmp_path: Path, capsys) -> None:
        from genotools.plotting.cli import main

        with pytest.raises(SystemExit) as exc:
            main(["--report", str(tmp_path / "nope.json")])
        assert exc.value.code == 1
        assert capsys.readouterr().err.startswith("ERROR: Report not found")


class TestPlotsFlag:
    BASE = ["genotools", "--pfile", "/data/test", "--out", "/tmp/out"]

    def test_parse_args(self) -> None:
        from genotools.cli import parse_args

        assert parse_args(self.BASE[1:]).output.plots is False
        assert parse_args(self.BASE[1:] + ["--plots"]).output.plots is True

    def _run_main(self, monkeypatch, argv, render_calls) -> None:
        """Drive main() with the pipeline stubbed out."""
        import genotools.cli as cli

        class Result:
            success = True

        monkeypatch.setattr(sys, "argv", argv)
        monkeypatch.setattr(cli, "run_pipeline", lambda args: Result())
        monkeypatch.setattr(cli, "write_results", lambda result, out: None)
        monkeypatch.setattr(cli, "_render_plots", render_calls.append)
        with pytest.raises(SystemExit) as exc:
            cli.main()
        assert exc.value.code == 0

    def test_main_renders_after_writing_the_report(self, monkeypatch) -> None:
        calls = []
        self._run_main(monkeypatch, self.BASE + ["--plots"], calls)
        assert calls == [Path("/tmp/out")]

    def test_main_does_not_render_without_the_flag(self, monkeypatch) -> None:
        calls = []
        self._run_main(monkeypatch, self.BASE, calls)
        assert calls == []

    def test_plotting_failure_never_fails_the_run(self, tmp_path: Path, caplog) -> None:
        """No report at the prefix: a warning, not an exception."""
        from genotools.cli import _render_plots

        _render_plots(tmp_path / "missing")
        assert "Plots not rendered" in caplog.text
