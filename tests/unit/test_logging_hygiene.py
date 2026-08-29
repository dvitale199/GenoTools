"""Static hygiene checks for the round-7 logging redesign.

Locks two invariants across the pipeline-path modules:
  * they bind the shared ``get_logger`` (not the bare ``logging.getLogger``), so
    every record lands under the ``genotools`` namespace and reaches the RunLog;
  * they route user-facing output through logging, not ``print()``.
"""

from pathlib import Path

import pytest

import genotools

GENOTOOLS_DIR = Path(genotools.__file__).resolve().parent

# Modules on the pipeline path that were standardized in round 7.
PIPELINE_MODULES = [
    "cli/runner.py",
    "core/validation.py",
    "core/executors.py",
    "ancestry/model.py",
    "ancestry/preprocessing.py",
    "gwas/steps/association.py",
    "gwas/steps/pca.py",
    "gwas/pipeline.py",
    "qc/steps/heterozygosity.py",
    "qc/steps/case_control.py",
    "qc/steps/sex.py",
    "qc/steps/ld_prune.py",
    "qc/steps/hwe.py",
    "qc/steps/relatedness.py",
    "qc/steps/haplotype.py",
    "qc/steps/variant_missingness.py",
    "qc/steps/callrate.py",
]


@pytest.mark.parametrize("rel", PIPELINE_MODULES)
def test_no_bare_getlogger(rel: str) -> None:
    src = (GENOTOOLS_DIR / rel).read_text()
    assert "logging.getLogger(" not in src, (
        f"{rel} should bind get_logger(__name__), not logging.getLogger"
    )


@pytest.mark.parametrize("rel", ["cli/runner.py", "core/validation.py"])
def test_no_print_calls(rel: str) -> None:
    """The runner and validation route everything through logging now."""
    src = (GENOTOOLS_DIR / rel).read_text()
    assert "print(" not in src, f"{rel} must not use print()"


def test_preprocessing_binds_get_logger() -> None:
    src = (GENOTOOLS_DIR / "ancestry/preprocessing.py").read_text()
    assert "from genotools.core.logging import get_logger" in src
    assert "logger = get_logger(__name__)" in src
