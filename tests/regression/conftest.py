"""Pytest fixtures for regression tests."""
import pytest
from pathlib import Path
import tempfile
import shutil

# Directory paths
TESTS_DIR = Path(__file__).parent.parent
REGRESSION_DIR = TESTS_DIR / "regression"
GOLDEN_DIR = REGRESSION_DIR / "golden"
DATA_DIR = TESTS_DIR / "data"
SYNTHETIC_DIR = DATA_DIR / "synthetic"


@pytest.fixture
def test_geno_path() -> Path:
    """
    Path to synthetic test genotype data.

    Returns the prefix path (without extension) for the synthetic test data.
    Skips if test data not found.
    """
    path = SYNTHETIC_DIR / "genotools_test"
    if not path.with_suffix(".pgen").exists():
        pytest.skip(
            "Synthetic test data not found. Run: "
            "python tests/data/synthetic/generate_test_data.py"
        )
    return path


@pytest.fixture
def outlier_manifest() -> Path:
    """
    Path to the outlier manifest CSV.

    This file documents all intentionally injected outliers in the test data.
    """
    path = SYNTHETIC_DIR / "outlier_manifest.csv"
    if not path.exists():
        pytest.skip(
            "Outlier manifest not found. Run: "
            "python tests/data/synthetic/generate_test_data.py"
        )
    return path


@pytest.fixture
def golden_dir() -> Path:
    """
    Path to golden reference outputs.

    Skips if golden files not generated yet.
    """
    if not GOLDEN_DIR.exists() or not any(GOLDEN_DIR.iterdir()):
        pytest.skip(
            "Golden files not found. Run: "
            "python tests/scripts/generate_golden.py --geno tests/data/synthetic/genotools_test --out tests/regression/golden"
        )
    return GOLDEN_DIR


@pytest.fixture
def golden_callrate(golden_dir: Path) -> Path:
    """Path to golden callrate output."""
    path = golden_dir / "callrate"
    if not path.exists():
        pytest.skip("Golden callrate files not found")
    return path


@pytest.fixture
def golden_sex(golden_dir: Path) -> Path:
    """Path to golden sex check output."""
    path = golden_dir / "sex"
    if not path.exists():
        pytest.skip("Golden sex files not found")
    return path


@pytest.fixture
def golden_het(golden_dir: Path) -> Path:
    """Path to golden heterozygosity output."""
    path = golden_dir / "het"
    if not path.exists():
        pytest.skip("Golden het files not found")
    return path


@pytest.fixture
def golden_related(golden_dir: Path) -> Path:
    """Path to golden relatedness output."""
    path = golden_dir / "related"
    if not path.exists():
        pytest.skip("Golden related files not found")
    return path


# Variant QC golden fixtures
@pytest.fixture
def golden_geno(golden_dir: Path) -> Path:
    """Path to golden variant missingness (geno) output."""
    path = golden_dir / "geno"
    if not path.exists():
        pytest.skip("Golden geno files not found")
    return path


@pytest.fixture
def golden_hwe(golden_dir: Path) -> Path:
    """Path to golden HWE output."""
    path = golden_dir / "hwe"
    if not path.exists():
        pytest.skip("Golden hwe files not found")
    return path


@pytest.fixture
def golden_case_control(golden_dir: Path) -> Path:
    """Path to golden case-control output."""
    path = golden_dir / "case_control"
    if not path.exists():
        pytest.skip("Golden case_control files not found")
    return path


@pytest.fixture
def golden_haplotype(golden_dir: Path) -> Path:
    """Path to golden haplotype output."""
    path = golden_dir / "haplotype"
    if not path.exists():
        pytest.skip("Golden haplotype files not found")
    return path


@pytest.fixture
def golden_ld(golden_dir: Path) -> Path:
    """Path to golden LD pruning output."""
    path = golden_dir / "ld"
    if not path.exists():
        pytest.skip("Golden ld files not found")
    return path


# Pipeline golden fixtures
@pytest.fixture
def golden_all_sample(golden_dir: Path) -> Path:
    """Path to golden --all_sample pipeline output."""
    path = golden_dir / "all_sample"
    if not path.exists():
        pytest.skip("Golden all_sample files not found")
    return path


@pytest.fixture
def golden_all_variant(golden_dir: Path) -> Path:
    """Path to golden --all_variant pipeline output."""
    path = golden_dir / "all_variant"
    if not path.exists():
        pytest.skip("Golden all_variant files not found")
    return path


@pytest.fixture
def tmp_output_dir(tmp_path: Path) -> Path:
    """
    Temporary directory for test outputs.

    Uses pytest's built-in tmp_path fixture.
    """
    return tmp_path


@pytest.fixture
def clean_output_dir(tmp_path: Path):
    """
    Factory fixture for creating clean output directories.

    Returns a function that creates named subdirectories.
    """
    def _make_dir(name: str) -> Path:
        out_dir = tmp_path / name
        out_dir.mkdir(parents=True, exist_ok=True)
        return out_dir
    return _make_dir


@pytest.fixture(scope="session")
def stable_venv_genotools():
    """
    Check if stable venv is available for comparison testing.

    Returns the path to stable venv's genotools or None.
    """
    stable_venv = Path(__file__).parent.parent.parent / ".venv-stable"
    genotools_bin = stable_venv / "bin" / "genotools"

    if not genotools_bin.exists():
        return None
    return genotools_bin
