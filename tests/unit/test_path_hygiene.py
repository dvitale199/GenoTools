"""Static check: pipeline code never calls ``Path.with_suffix``.

GenoTools passes PLINK *prefixes* around -- extensionless paths like
``/data/cohort`` that stand for cohort.pgen/.pvar/.psam. ``with_suffix``
replaces whatever follows the last dot, so on a prefix whose name contains a
period it silently rewrites the name: ``genotools.test.v2`` -> ``genotools.test.pgen``.
Every QC step then failed to find its input and, under the default --warn, the
run exited 0 with the raw input copied out as the "cleaned" result.

Build file names from a prefix by appending, the way PLINK does:
``Path(f"{prefix}.pgen")``.
"""

from pathlib import Path

import pytest

import genotools

GENOTOOLS_DIR = Path(genotools.__file__).resolve().parent

# Packages on the 2.0 pipeline path. rglob, not a fixed list, so new modules
# are covered without anyone remembering to add them.
PIPELINE_PACKAGES = ["cli", "core", "qc", "gwas", "ancestry"]

PIPELINE_SOURCES = sorted(
    str(p.relative_to(GENOTOOLS_DIR))
    for pkg in PIPELINE_PACKAGES
    for p in (GENOTOOLS_DIR / pkg).rglob("*.py")
)


def test_pipeline_sources_found() -> None:
    """Guard the guard: an empty parametrize list would pass vacuously."""
    assert len(PIPELINE_SOURCES) > 20


@pytest.mark.parametrize("rel", PIPELINE_SOURCES)
def test_no_with_suffix(rel: str) -> None:
    src = (GENOTOOLS_DIR / rel).read_text()
    assert ".with_suffix(" not in src, (
        f"{rel} calls .with_suffix(), which mangles a PLINK prefix whose name "
        f'contains a period. Append instead: Path(f"{{prefix}}.ext").'
    )
