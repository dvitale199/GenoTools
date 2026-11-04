"""Hand-built reports for the plotting tests."""

import json
import shutil
from pathlib import Path
from typing import Any, Dict, List, Optional

import pandas as pd
import pytest


def qc_row(
    step: str,
    count: int,
    metric: str,
    ancestry: str = "all",
    level: str = "sample",
    outcome: str = "pass",
    reason: Optional[str] = None,
) -> Dict[str, Any]:
    return {
        "step": step,
        "pruned_count": count,
        "metric": metric,
        "ancestry": ancestry,
        "level": level,
        "pass": outcome == "pass",
        "outcome": outcome,
        "reason": reason,
    }


def write_report(path: Path, qc: List[Dict[str, Any]], labels: Optional[Dict[str, int]] = None,
                 n_input: Optional[int] = None) -> Path:
    """Write a report shaped like `PipelineOutput.to_dict()`.

    Tables go through `DataFrame.to_dict()`, as the real writer does, so their
    row indices come back as strings.
    """
    report: Dict[str, Any] = {"QC": pd.DataFrame(qc).to_dict()}
    if labels:
        ids = [f"s{i}" for i in range(sum(labels.values()))]
        label_col = [label for label, n in labels.items() for _ in range(n)]
        report["ancestry_labels"] = pd.DataFrame(
            {"FID": ids, "IID": ids, "label": label_col}
        ).to_dict()
        n_input = n_input or len(ids)
    if n_input:
        ids = [f"s{i}" for i in range(n_input)]
        report["input_samples"] = pd.DataFrame({"#FID": ids, "IID": ids}).to_dict()
    path.write_text(json.dumps(report))
    return path


@pytest.fixture
def ancestry_report(tmp_path: Path) -> Path:
    """Two groups; AMR skips het, EUR's related step reports two metrics."""
    qc = [
        qc_row("callrate_prune", 5, "outlier_count", "EUR"),
        qc_row("het_prune", 2, "outlier_count", "EUR"),
        qc_row("related_prune", 3, "related_count", "EUR"),
        qc_row("related_prune", 1, "duplicated_count", "EUR"),
        qc_row("geno_prune", 100, "geno_removed_count", "EUR", level="variant"),
        qc_row("callrate_prune", 1, "outlier_count", "AMR"),
        qc_row("het_prune", 0, "outlier_count", "AMR", outcome="skipped",
               reason="40 samples, below the 50-sample floor"),
        qc_row("related_prune", 0, "related_count", "AMR"),
        qc_row("related_prune", 0, "duplicated_count", "AMR"),
        qc_row("geno_prune", 20, "geno_removed_count", "AMR", level="variant"),
    ]
    return write_report(tmp_path / "cohort.json", qc, labels={"EUR": 200, "AMR": 40})


GOLDEN_DIR = Path(__file__).resolve().parents[2] / "regression" / "golden"


def _copy_golden(name: str, tmp_path: Path) -> Path:
    """Copy one golden run's report and plot tables into `tmp_path`.

    Plots default to writing beside the report, so tests never point at
    `tests/regression/golden/` itself. The copy's `plot_data` is repointed at
    the copied tables.
    """
    src = GOLDEN_DIR / name
    dest = tmp_path / name
    dest.mkdir()
    shutil.copytree(src / "output_plot_data", dest / "output_plot_data")
    raw = json.loads((src / "output.json").read_text())
    raw["plot_data"] = str(dest / "output_plot_data")
    report = dest / "output.json"
    report.write_text(json.dumps(raw))
    return report


@pytest.fixture
def golden_sample_report(tmp_path: Path) -> Path:
    """A real `--all-sample` report on the synthetic data, with its plot tables."""
    return _copy_golden("all_sample", tmp_path)


@pytest.fixture
def golden_variant_report(tmp_path: Path) -> Path:
    """A real `--all-variant` report on the synthetic data, with its plot tables."""
    return _copy_golden("all_variant", tmp_path)


@pytest.fixture(params=["all_sample", "all_variant"])
def golden_qc_report(request, tmp_path: Path) -> Path:
    """Each real QC report in turn: for checks every plot must pass on both."""
    return _copy_golden(request.param, tmp_path)
