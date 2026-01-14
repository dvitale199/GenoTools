# GenoTools Refactor Plan

This document outlines the phased architectural refactor for GenoTools.

---

## Executive Summary

GenoTools currently suffers from tight coupling, mutable state, and mixed concerns that make it difficult to test, debug, and extend. This refactor introduces a clean separation between QC (filter pipeline) and Ancestry (ML pipeline), with immutable data flow and pure functions.

**Approach:** Build new architecture alongside existing code, migrate module-by-module, then cutover.

---

## Current Problems

### All Issues

| ID | Issue | Priority | Phase |
|----|-------|----------|-------|
| A1 | Mutable State Architecture | P0 | 1, 2 |
| A2 | QC and Ancestry share execution model | P0 | 2, 3 |
| A3 | No separation of concerns | P1 | 4 |
| P1 | Repeated format conversions | P0 | 1 |
| R1 | Inconsistent error handling | P0 | 1 |
| R2 | Fragile file cleanup | P2 | 1 |
| R3 | Unsafe subprocess handling | P2 | 1 |
| M1 | Hardcoded magic values | P1 | 2, 3 |
| M2 | Duplicated logic patterns | P2 | 2 |
| M3 | Monolithic main function | P2 | 4 |
| D1 | No type hints or validation | P1 | 1, 2, 3 |
| D2 | Module-level dependency init | P1 | 1 |
| D3 | No structured logging | P1 | 1 |
| C1 | Boolean string parsing | P2 | 4 |

---

## Current Conventions to Preserve

These patterns from the existing codebase (documented in CLAUDE.md) must be maintained during the refactor.

### Standard Return Dictionary Format

Every QC method currently returns this structure:
```python
{
    'pass': bool,           # True if step completed successfully
    'step': str,            # Step identifier (e.g., 'callrate_prune')
    'metrics': {
        'outlier_count': int,  # Number of samples/variants pruned
        # ... other step-specific metrics
    },
    'output': {
        'pruned_samples': str,  # Path to pruned sample IDs (or None)
        'plink_out': str,       # Path to output pfiles (without extension)
        # ... other output files
    }
}
```

The new `FilterResult` dataclass must provide a `.to_dict()` method for backward compatibility.

### PLINK2 psam Column Preservation

**Critical:** All `--make-pgen` commands must preserve sample metadata:
```bash
--make-pgen psam-cols=fid,parents,sex,pheno1,phenos
```

### File Path Conventions

- Input/output paths never include extensions: `/path/to/data` → `data.pgen`, `data.pvar`, `data.psam`
- Intermediate files use step suffix: `{out_path}_{step}`
- Outlier files use `.outliers` extension

### Outlier File Format

Outlier files must be tab-separated with `#FID` header:
```python
df = df.rename({'FID': '#FID'}, axis=1)
df.to_csv(outliers_out, sep='\t', header=True, index=False)
```

### Platform Constraints

- **KING**: Linux only. Must check `platform.system() == 'Linux'` before use.

### Pinned Dependencies

- `umap-learn==0.5.3` - Required for model compatibility

### Utility Functions

Current code uses `shell_do()` and `concat_logs()` from `genotools.utils`. The new `run_command()` in `core/executors.py` should wrap these patterns.

---

## Phase 0: Regression Testing Framework

### Goal
Build a testing harness that ensures new implementations produce identical results to the current pip-installable version.

### Approach
1. **Generate golden files** using current implementation
2. **Compare new implementations** against golden outputs
3. **Run comparisons per-step** to isolate regressions

### Deliverables

#### 0.1 Directory Structure
```
tests/
├── regression/
│   ├── conftest.py          # Pytest fixtures for test data paths
│   ├── compare.py           # Comparison utilities
│   ├── test_qc_steps.py     # QC step regression tests
│   ├── test_ancestry.py     # Ancestry regression tests
│   └── golden/              # Expected outputs (generated once)
│       ├── callrate/
│       ├── sex/
│       ├── het/
│       ├── related/
│       └── ancestry/
├── data/                    # Test genotype data (user-supplied)
│   └── README.md
└── scripts/
    └── generate_golden.py   # One-time golden file generator
```

#### 0.2 `tests/regression/compare.py`
```python
"""Utilities for comparing genomic outputs."""
from dataclasses import dataclass
from pathlib import Path
import pandas as pd

@dataclass
class ComparisonResult:
    """Result of comparing two outputs."""
    equal: bool
    sample_diff: int
    variant_diff: int
    mismatched_samples: list[str]
    mismatched_variants: list[str]
    message: str

def compare_sample_ids(expected: Path, actual: Path) -> ComparisonResult:
    """Compare sample IDs between two psam/fam files."""
    ext = expected.suffix
    if ext == ".psam":
        exp_df = pd.read_csv(expected, sep="\t", dtype=str)
        act_df = pd.read_csv(actual, sep="\t", dtype=str)
        exp_ids = set(exp_df["IID"])
        act_ids = set(act_df["IID"])
    else:  # .fam
        exp_df = pd.read_csv(expected, sep=r"\s+", header=None, dtype=str)
        act_df = pd.read_csv(actual, sep=r"\s+", header=None, dtype=str)
        exp_ids = set(exp_df[1])
        act_ids = set(act_df[1])

    missing = exp_ids - act_ids
    extra = act_ids - exp_ids

    return ComparisonResult(
        equal=len(missing) == 0 and len(extra) == 0,
        sample_diff=len(missing) + len(extra),
        variant_diff=0,
        mismatched_samples=list(missing | extra),
        mismatched_variants=[],
        message=f"Missing: {len(missing)}, Extra: {len(extra)}"
    )

def compare_variant_ids(expected: Path, actual: Path) -> ComparisonResult:
    """Compare variant IDs between two pvar/bim files."""
    ext = expected.suffix
    if ext == ".pvar":
        exp_df = pd.read_csv(expected, sep="\t", comment="#", dtype=str)
        act_df = pd.read_csv(actual, sep="\t", comment="#", dtype=str)
        exp_ids = set(exp_df["ID"])
        act_ids = set(act_df["ID"])
    else:  # .bim
        exp_df = pd.read_csv(expected, sep="\t", header=None, dtype=str)
        act_df = pd.read_csv(actual, sep="\t", header=None, dtype=str)
        exp_ids = set(exp_df[1])
        act_ids = set(act_df[1])

    missing = exp_ids - act_ids
    extra = act_ids - exp_ids

    return ComparisonResult(
        equal=len(missing) == 0 and len(extra) == 0,
        sample_diff=0,
        variant_diff=len(missing) + len(extra),
        mismatched_samples=[],
        mismatched_variants=list(missing | extra),
        message=f"Missing: {len(missing)}, Extra: {len(extra)}"
    )

def compare_pfiles(expected_prefix: Path, actual_prefix: Path) -> ComparisonResult:
    """Compare two pfile sets (samples + variants)."""
    sample_result = compare_sample_ids(
        expected_prefix.with_suffix(".psam"),
        actual_prefix.with_suffix(".psam")
    )
    variant_result = compare_variant_ids(
        expected_prefix.with_suffix(".pvar"),
        actual_prefix.with_suffix(".pvar")
    )

    return ComparisonResult(
        equal=sample_result.equal and variant_result.equal,
        sample_diff=sample_result.sample_diff,
        variant_diff=variant_result.variant_diff,
        mismatched_samples=sample_result.mismatched_samples,
        mismatched_variants=variant_result.mismatched_variants,
        message=f"Samples: {sample_result.message}; Variants: {variant_result.message}"
    )

def compare_ancestry_predictions(
    expected: Path,
    actual: Path,
    tolerance: float = 0.0
) -> ComparisonResult:
    """Compare ancestry prediction outputs."""
    exp_df = pd.read_csv(expected, sep="\t")
    act_df = pd.read_csv(actual, sep="\t")

    # Merge on sample ID
    merged = exp_df.merge(act_df, on="IID", suffixes=("_exp", "_act"))

    # Check predicted ancestry matches
    mismatched = merged[merged["predicted_ancestry_exp"] != merged["predicted_ancestry_act"]]

    return ComparisonResult(
        equal=len(mismatched) == 0,
        sample_diff=len(mismatched),
        variant_diff=0,
        mismatched_samples=list(mismatched["IID"]),
        mismatched_variants=[],
        message=f"{len(mismatched)} samples have different predictions"
    )
```

#### 0.3 `scripts/generate_golden.py`
```python
#!/usr/bin/env python3
"""Generate golden reference files using current (pip-installed) implementation.

Run this ONCE before starting the refactor to capture expected outputs.

Usage:
    python scripts/generate_golden.py --geno tests/data/test_geno --out tests/regression/golden
"""
import argparse
from pathlib import Path
import shutil

# Import OLD implementation
from genotools.qc import SampleQC, VariantQC
from genotools.ancestry import Ancestry


def generate_qc_golden(geno_path: Path, out_dir: Path):
    """Run each QC step and save outputs."""

    # Sample QC steps
    sample_qc = SampleQC(geno_path=str(geno_path), out_path=str(out_dir / "callrate" / "output"))
    sample_qc.run_callrate_prune(0.02)
    # Copy output files to golden dir

    # Continue for each step...
    steps = [
        ("callrate", lambda sq: sq.run_callrate_prune(0.02)),
        ("sex", lambda sq: sq.run_sex_prune()),
        ("het", lambda sq: sq.run_het_prune()),
        ("related", lambda sq: sq.run_related_prune(0.0884)),
    ]

    for step_name, step_fn in steps:
        step_out = out_dir / step_name
        step_out.mkdir(parents=True, exist_ok=True)

        sq = SampleQC()
        sq.geno_path = str(geno_path)
        sq.out_path = str(step_out / "output")

        result = step_fn(sq)

        # Save metrics
        with open(step_out / "metrics.json", "w") as f:
            json.dump(result, f, indent=2)


def generate_ancestry_golden(geno_path: Path, ref_path: Path, out_dir: Path):
    """Run ancestry prediction and save outputs."""
    ancestry_out = out_dir / "ancestry"
    ancestry_out.mkdir(parents=True, exist_ok=True)

    # Run old implementation
    anc = Ancestry(geno_path=str(geno_path), ref_panel=str(ref_path), out_path=str(ancestry_out / "output"))
    result = anc.predict()

    # Predictions are saved by the old code


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--geno", type=Path, required=True, help="Test genotype file prefix")
    parser.add_argument("--ref", type=Path, help="Reference panel for ancestry (optional)")
    parser.add_argument("--out", type=Path, required=True, help="Golden output directory")
    args = parser.parse_args()

    args.out.mkdir(parents=True, exist_ok=True)

    print(f"Generating golden files from: {args.geno}")
    print(f"Output directory: {args.out}")

    generate_qc_golden(args.geno, args.out)

    if args.ref:
        generate_ancestry_golden(args.geno, args.ref, args.out)

    print("Golden files generated successfully!")


if __name__ == "__main__":
    main()
```

#### 0.4 `tests/regression/conftest.py`
```python
"""Pytest fixtures for regression tests."""
import pytest
from pathlib import Path

TESTS_DIR = Path(__file__).parent.parent
GOLDEN_DIR = TESTS_DIR / "regression" / "golden"
DATA_DIR = TESTS_DIR / "data"


@pytest.fixture
def test_geno_path() -> Path:
    """Path to test genotype data."""
    path = DATA_DIR / "test_geno"
    if not path.with_suffix(".pgen").exists() and not path.with_suffix(".bed").exists():
        pytest.skip("Test data not found. Place test data in tests/data/test_geno.*")
    return path


@pytest.fixture
def golden_dir() -> Path:
    """Path to golden reference outputs."""
    if not GOLDEN_DIR.exists():
        pytest.skip("Golden files not found. Run scripts/generate_golden.py first.")
    return GOLDEN_DIR
```

#### 0.5 `tests/regression/test_qc_steps.py`
```python
"""Regression tests comparing new QC implementation to golden reference."""
import pytest
from pathlib import Path

from genotools.core.genotypes import GenotypeData
from genotools.qc.steps.callrate import filter_callrate
from genotools.qc.config import CallrateConfig
from .compare import compare_pfiles


class TestCallrateRegression:
    """Regression tests for callrate filtering."""

    def test_output_matches_golden(self, test_geno_path: Path, golden_dir: Path, tmp_path: Path):
        """New implementation produces same output as old."""
        # Run new implementation
        data = GenotypeData.from_path(test_geno_path)
        result = filter_callrate(data, CallrateConfig(), tmp_path / "output")

        # Compare to golden
        comparison = compare_pfiles(
            golden_dir / "callrate" / "output",
            result.output.path
        )

        assert comparison.equal, (
            f"Output differs from golden:\n"
            f"  Sample diff: {comparison.sample_diff}\n"
            f"  Variant diff: {comparison.variant_diff}\n"
            f"  {comparison.message}"
        )

    def test_metrics_match(self, test_geno_path: Path, golden_dir: Path, tmp_path: Path):
        """Sample/variant counts match golden."""
        import json

        # Load golden metrics
        with open(golden_dir / "callrate" / "metrics.json") as f:
            golden_metrics = json.load(f)

        # Run new implementation
        data = GenotypeData.from_path(test_geno_path)
        result = filter_callrate(data, CallrateConfig(), tmp_path / "output")

        assert result.samples_removed == golden_metrics["samples_removed"]
        assert result.variants_removed == golden_metrics["variants_removed"]


# Similar test classes for sex, het, related, etc.
```

#### 0.6 `tests/data/README.md`
```markdown
# Test Data

Place test genotype files here for regression testing.

## Required Files

For pfile format:
- `test_geno.pgen`
- `test_geno.pvar`
- `test_geno.psam`

OR for bfile format:
- `test_geno.bed`
- `test_geno.bim`
- `test_geno.fam`

## Recommendations

- Use a small dataset (100-500 samples, 10k-50k variants)
- Include samples that will be filtered by each QC step
- Include known ancestry labels if testing ancestry prediction

## Generating Golden Files

After placing test data, run:

```bash
python scripts/generate_golden.py --geno tests/data/test_geno --out tests/regression/golden
```

This uses the CURRENT pip-installed version to generate expected outputs.
```

### What to Compare Per Step

| Step | Compare | Ignore |
|------|---------|--------|
| callrate | Sample IDs, variant IDs, counts | Log files, ordering |
| sex | Removed sample IDs, F-statistics (within tolerance) | Log files |
| het | Removed sample IDs, het rates (within tolerance) | Log files |
| related | Related pairs, removed sample IDs | Log files, pair ordering |
| ancestry | Predicted labels, probabilities (within tolerance) | Log files |

### Success Criteria
- [ ] `tests/regression/compare.py` implemented
- [ ] `scripts/generate_golden.py` implemented
- [ ] Golden files generated from current pip version
- [ ] `pytest tests/regression/` passes (baseline: old vs old)
- [ ] Test data documented in `tests/data/README.md`

### Usage During Refactor

```bash
# 1. Install current version
pip install the-real-genotools

# 2. Generate golden files (once)
python scripts/generate_golden.py --geno tests/data/test_geno --out tests/regression/golden

# 3. Implement new code...

# 4. Run regression tests
pytest tests/regression/ -v

# 5. If tests fail, investigate differences
pytest tests/regression/test_qc_steps.py::TestCallrateRegression -v --tb=long
```

---

## Target Architecture

### Core Insight

**QC and Ancestry are fundamentally different:**

| Aspect | QC | Ancestry |
|--------|----|---------|
| Model | Filter pipeline | ML pipeline |
| Data flow | N samples → N' samples (subset) | Samples → predictions |
| State | Stateless filters | Fitted model (UMAP, XGBoost) |
| Composition | Steps are reorderable | Phases are sequential |

They should not share the same execution model.

### Final Directory Structure

```
genotools/
├── __init__.py
├── __main__.py           # Thin entry point
│
├── core/                 # Shared infrastructure (Phase 1)
│   ├── __init__.py
│   ├── genotypes.py      # GenotypeData - immutable, format-aware
│   ├── executors.py      # External tool wrapper (PLINK, KING)
│   ├── logging.py        # Structured logging setup
│   ├── config.py         # Base config classes
│   └── exceptions.py     # Custom exception hierarchy
│
├── qc/                   # QC domain (Phase 2)
│   ├── __init__.py
│   ├── pipeline.py       # QCPipeline - composes steps
│   ├── config.py         # QCConfig, step configs
│   ├── results.py        # FilterResult, QCResult dataclasses
│   └── steps/            # Pure filter functions
│       ├── __init__.py
│       ├── callrate.py
│       ├── sex.py
│       ├── heterozygosity.py
│       ├── relatedness.py
│       ├── variant_missingness.py
│       ├── case_control.py
│       ├── haplotype.py
│       ├── hwe.py
│       └── ld_prune.py
│
├── ancestry/             # Ancestry domain (Phase 3)
│   ├── __init__.py
│   ├── model.py          # AncestryModel - fit/predict interface
│   ├── reference.py      # ReferencePanel management
│   ├── config.py         # AncestryConfig
│   ├── results.py        # AncestryPredictions dataclass
│   └── reducers/         # Dimensionality reduction
│       ├── __init__.py
│       ├── pca.py
│       └── umap.py
│
├── gwas/                 # GWAS domain (Phase 5)
│   └── ...
│
└── cli/                  # CLI layer (Phase 4)
    ├── __init__.py
    ├── parser.py         # Argparse with proper types
    ├── runner.py         # Orchestration logic
    └── output.py         # Result formatting, JSON serialization
```

---

## Phase 1: Core Foundation

### Goal
Build shared infrastructure that all other phases depend on.

### Issues Resolved
| ID | Issue | Resolution |
|----|-------|------------|
| A1 | Mutable state (partial) | GenotypeData is immutable |
| P1 | Format conversions | GenotypeData manages format, converts once |
| R1 | Error handling | Custom exception hierarchy |
| R2 | File cleanup | Context managers in executors.py |
| R3 | Subprocess handling | shlex.split() or list commands |
| D1 | Type hints (partial) | Dataclasses with full type hints |
| D2 | Module-level init | Lazy initialization in executors.py |
| D3 | Structured logging | core/logging.py with levels and context |

### Deliverables

#### 1.1 `core/exceptions.py`
```python
class GenoToolsError(Exception):
    """Base exception for all GenoTools errors."""
    pass

class QCError(GenoToolsError):
    """QC-specific errors."""
    pass

class AncestryError(GenoToolsError):
    """Ancestry-specific errors."""
    pass

class ExternalToolError(GenoToolsError):
    """External tool (PLINK, KING) failures."""
    def __init__(self, tool: str, command: str, returncode: int, stderr: str):
        self.tool = tool
        self.command = command
        self.returncode = returncode
        self.stderr = stderr
        super().__init__(f"{tool} failed with code {returncode}: {stderr[:200]}")
```

#### 1.2 `core/logging.py`
```python
import logging
from contextvars import ContextVar
from pathlib import Path
from typing import Optional

current_step: ContextVar[str] = ContextVar("current_step", default="")

class StepContextFilter(logging.Filter):
    def filter(self, record):
        record.step = current_step.get()
        return True

def setup_logging(level: str = "INFO", log_file: Optional[Path] = None):
    """Configure structured logging with optional file output."""
    formatter = logging.Formatter(
        "%(asctime)s [%(levelname)s] [%(step)s] %(message)s"
    )
    ...
```

#### 1.3 `core/genotypes.py`
```python
from dataclasses import dataclass
from pathlib import Path
from typing import Literal

@dataclass(frozen=True)
class GenotypeData:
    """Immutable reference to genotype files."""
    path: Path
    format: Literal["pfile", "bfile"]
    sample_count: int
    variant_count: int

    @classmethod
    def from_path(cls, path: Path) -> "GenotypeData":
        """Detect format and count samples/variants."""
        ...

    def to_bfile(self, out_path: Path) -> "GenotypeData":
        """Convert to bfile format, return NEW GenotypeData."""
        ...

    def to_pfile(self, out_path: Path) -> "GenotypeData":
        """Convert to pfile format, return NEW GenotypeData."""
        ...

    def ensure_bfile(self, work_dir: Path) -> "GenotypeData":
        """Return self if bfile, or convert and return new."""
        if self.format == "bfile":
            return self
        return self.to_bfile(work_dir / f"{self.path.stem}_bfile")
```

#### 1.4 `core/executors.py`
```python
import shlex
import subprocess
import platform
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path
from typing import Optional
import logging

from .exceptions import ExternalToolError

logger = logging.getLogger(__name__)

# Lazy-loaded executables
_plink: Optional[Path] = None
_plink2: Optional[Path] = None
_king: Optional[Path] = None

# CRITICAL: Always use this flag to preserve sample metadata
PLINK2_PSAM_COLS = "psam-cols=fid,parents,sex,pheno1,phenos"

def get_plink() -> Path:
    global _plink
    if _plink is None:
        _plink = _find_or_download("plink")
    return _plink

def get_king() -> Optional[Path]:
    """Get KING executable. Returns None on non-Linux platforms."""
    global _king
    if platform.system() != "Linux":
        logger.warning("KING is only available on Linux")
        return None
    if _king is None:
        _king = _find_or_download("king")
    return _king

@dataclass
class CommandResult:
    returncode: int
    stdout: str
    stderr: str
    log_file: Optional[Path] = None

def run_command(
    command: list[str],
    tool_name: str = "command",
    check: bool = True
) -> CommandResult:
    """Run external command with proper error handling.

    This wraps the existing shell_do() pattern with structured return.
    """
    logger.debug(f"Running: {' '.join(command)}")

    result = subprocess.run(
        command,
        capture_output=True,
        text=True
    )

    cmd_result = CommandResult(
        returncode=result.returncode,
        stdout=result.stdout,
        stderr=result.stderr
    )

    if check and result.returncode != 0:
        raise ExternalToolError(
            tool=tool_name,
            command=' '.join(command),
            returncode=result.returncode,
            stderr=result.stderr
        )

    return cmd_result

def run_plink2_make_pgen(
    input_path: Path,
    output_path: Path,
    input_format: str = "pfile",
    extra_args: list[str] = None
) -> CommandResult:
    """Run PLINK2 --make-pgen with required psam column preservation."""
    input_flag = "--pfile" if input_format == "pfile" else "--bfile"
    command = [
        str(get_plink2()),
        input_flag, str(input_path),
        "--make-pgen", PLINK2_PSAM_COLS,
        "--out", str(output_path)
    ]
    if extra_args:
        # Insert extra args before --out
        command = command[:-2] + extra_args + command[-2:]
    return run_command(command, tool_name="plink2")

@contextmanager
def temp_files(*paths: Path):
    """Context manager for temporary file cleanup."""
    try:
        yield
    finally:
        for path in paths:
            for ext in [".bed", ".bim", ".fam", ".pgen", ".pvar", ".psam", ".log"]:
                f = path.with_suffix(ext)
                if f.exists():
                    f.unlink()
```

### Success Criteria
- [ ] `from genotools.core import GenotypeData, setup_logging` works
- [ ] Import time < 100ms (no side effects)
- [ ] `mypy genotools/core/ --strict` passes
- [ ] Unit tests for GenotypeData format detection and conversion

---

## Phase 2: QC Migration

### Goal
Replace mutable SampleQC/VariantQC classes with pure functions composed into a pipeline.

### Issues Resolved
| ID | Issue | Resolution |
|----|-------|------------|
| A1 | Mutable state (complete) | Steps are pure functions |
| A2 | Shared execution model (QC part) | QC uses step composition pattern |
| M1 | Magic values (QC part) | qc/config.py with documented defaults |
| M2 | Duplicated logic | Extracted to step protocol |
| D1 | Type hints (QC part) | Full type hints on all steps |

### Deliverables

#### 2.1 `qc/results.py`
```python
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Optional
from ..core.genotypes import GenotypeData

@dataclass(frozen=True)
class FilterResult:
    """Result of a single QC step."""
    output: GenotypeData
    samples_removed: int
    variants_removed: int
    metrics: dict[str, Any] = field(default_factory=dict)
    log: str = ""
    pruned_samples_file: Optional[Path] = None  # Path to .outliers file

    def to_dict(self) -> dict:
        """Convert to legacy dictionary format for backward compatibility.

        Returns the standard format:
        {
            'pass': bool,
            'step': str,
            'metrics': {'outlier_count': int, ...},
            'output': {'pruned_samples': str, 'plink_out': str}
        }
        """
        return {
            'pass': True,  # If we have a result, it passed
            'step': self.metrics.get('step_name', 'unknown'),
            'metrics': {
                'outlier_count': self.samples_removed + self.variants_removed,
                **self.metrics
            },
            'output': {
                'pruned_samples': str(self.pruned_samples_file) if self.pruned_samples_file else None,
                'plink_out': str(self.output.path)
            }
        }

@dataclass(frozen=True)
class QCResult:
    """Result of full QC pipeline."""
    input: GenotypeData
    output: GenotypeData
    step_results: list[tuple[str, FilterResult]]

    @property
    def total_samples_removed(self) -> int:
        return sum(r.samples_removed for _, r in self.step_results)

    @property
    def total_variants_removed(self) -> int:
        return sum(r.variants_removed for _, r in self.step_results)

    def to_legacy_dict(self) -> dict:
        """Convert to legacy pass_fail dictionary format."""
        return {
            name: result.to_dict()
            for name, result in self.step_results
        }
```

#### 2.2 `qc/config.py`
```python
from dataclasses import dataclass

@dataclass(frozen=True)
class CallrateConfig:
    """Sample/variant call rate thresholds."""
    sample_threshold: float = 0.02  # Remove samples with >2% missing
    variant_threshold: float = 0.02  # Remove variants with >2% missing

@dataclass(frozen=True)
class SexConfig:
    """Sex check configuration."""
    female_max_f: float = 0.2   # F-statistic threshold for female
    male_min_f: float = 0.8     # F-statistic threshold for male

@dataclass(frozen=True)
class HetConfig:
    """Heterozygosity outlier detection."""
    std_devs: float = 3.0  # Remove samples beyond N std devs

@dataclass(frozen=True)
class RelatedConfig:
    """Relatedness/kinship thresholds."""
    kinship_threshold: float = 0.0884  # KING kinship for 2nd degree
    # Reference: Manichaikul et al. 2010
    # 0.354 = MZ twin/duplicate
    # 0.177 = 1st degree (parent-child, full sibling)
    # 0.0884 = 2nd degree (half-sibling, avuncular, grandparent)
    # 0.0442 = 3rd degree (first cousin)

@dataclass(frozen=True)
class HWEConfig:
    """Hardy-Weinberg equilibrium pruning."""
    p_value: float = 1e-6  # Remove variants with HWE p < threshold

@dataclass(frozen=True)
class LDConfig:
    """LD pruning parameters."""
    window_size: int = 1000    # Window size in kb
    step_size: int = 10        # Step size in variants
    r2_threshold: float = 0.02  # r² threshold

# ... additional configs for each step
```

#### 2.3 `qc/steps/callrate.py` (example step)
```python
from pathlib import Path
import logging
import pandas as pd

from ...core.genotypes import GenotypeData
from ...core.executors import run_plink2_make_pgen, get_plink2, run_command
from ..config import CallrateConfig
from ..results import FilterResult

logger = logging.getLogger(__name__)

def filter_callrate(
    data: GenotypeData,
    config: CallrateConfig,
    out_path: Path
) -> FilterResult:
    """
    Filter samples and variants by call rate.

    Removes:
    - Samples with missingness > config.sample_threshold
    - Variants with missingness > config.variant_threshold
    """
    step_name = "callrate_prune"
    logger.info(f"Filtering by callrate (sample={config.sample_threshold}, variant={config.variant_threshold})")

    # Use helper function that ensures psam-cols preservation
    result = run_plink2_make_pgen(
        input_path=data.path,
        output_path=out_path,
        input_format=data.format,
        extra_args=[
            "--mind", str(config.sample_threshold),
            "--geno", str(config.variant_threshold),
        ]
    )

    output = GenotypeData.from_path(out_path)

    # Write outlier file in correct format (tab-separated, #FID header)
    outliers_path = out_path.parent / f"{out_path.name}.outliers"
    if (irem_file := out_path.with_suffix(".mindrem.id")).exists():
        outliers = pd.read_csv(irem_file, sep=r"\s+", dtype=str)
        outliers = outliers.rename(columns={"FID": "#FID"})
        outliers.to_csv(outliers_path, sep="\t", header=True, index=False)
    else:
        outliers_path = None

    return FilterResult(
        output=output,
        samples_removed=data.sample_count - output.sample_count,
        variants_removed=data.variant_count - output.variant_count,
        metrics={
            "step_name": step_name,
            "sample_threshold": config.sample_threshold,
            "variant_threshold": config.variant_threshold,
        },
        log=result.stderr,
        pruned_samples_file=outliers_path
    )
```

#### 2.4 `qc/pipeline.py`
```python
from dataclasses import dataclass, field
from pathlib import Path
from typing import Protocol, Any
import logging

from ..core.genotypes import GenotypeData
from ..core.logging import current_step
from ..core.exceptions import QCError
from .results import FilterResult, QCResult

logger = logging.getLogger(__name__)

class QCStep(Protocol):
    """Protocol for QC step functions."""
    def __call__(
        self,
        data: GenotypeData,
        config: Any,
        out_path: Path
    ) -> FilterResult:
        ...

@dataclass
class QCPipeline:
    """Composes QC steps into a pipeline."""
    steps: list[tuple[str, QCStep, Any]]  # (name, function, config)
    warn_only: bool = False  # If True, continue on step failure (--warn flag)

    def run(self, data: GenotypeData, out_dir: Path) -> QCResult:
        """Execute all steps sequentially."""
        out_dir.mkdir(parents=True, exist_ok=True)

        current = data
        step_results: list[tuple[str, FilterResult]] = []

        for name, step_fn, config in self.steps:
            current_step.set(name)
            logger.info(f"Starting step: {name}")

            step_out = out_dir / name / "output"
            step_out.parent.mkdir(exist_ok=True)

            try:
                result = step_fn(current, config, step_out)
                step_results.append((name, result))

                logger.info(
                    f"Completed {name}: "
                    f"-{result.samples_removed} samples, "
                    f"-{result.variants_removed} variants"
                )

                current = result.output

                if current.sample_count == 0:
                    raise QCError(f"All samples removed at step: {name}")
                if current.variant_count == 0:
                    raise QCError(f"All variants removed at step: {name}")

            except QCError as e:
                if self.warn_only:
                    logger.warning(f"Step {name} failed but continuing (--warn mode): {e}")
                    # Skip this step, keep using previous data
                    continue
                else:
                    raise

        current_step.set("")

        return QCResult(
            input=data,
            output=current,
            step_results=step_results
        )
```

### Success Criteria
- [ ] Each step function works independently
- [ ] `QCPipeline.run()` produces same output as old `SampleQC`/`VariantQC`
- [ ] No mutable state in any QC code
- [ ] All thresholds documented in config.py
- [ ] `mypy genotools/qc/ --strict` passes

---

## Phase 3: Ancestry Migration

### Goal
Replace current ancestry code with ML-pattern fit/predict interface.

### Issues Resolved
| ID | Issue | Resolution |
|----|-------|------------|
| A2 | Shared execution model (complete) | Ancestry uses fit/predict pattern |
| M1 | Magic values (Ancestry part) | ancestry/config.py with documented defaults |
| D1 | Type hints (Ancestry part) | Full type hints |

### Deliverables

#### 3.1 `ancestry/config.py`
```python
from dataclasses import dataclass
from typing import Literal

# CRITICAL: umap-learn must be pinned to 0.5.3 for model compatibility
# See requirements.txt: umap-learn==0.5.3

@dataclass(frozen=True)
class PCAConfig:
    """PCA configuration."""
    n_components: int = 20
    maf: float = 0.01
    geno: float = 0.01
    hwe: float = 5e-6

@dataclass(frozen=True)
class UMAPConfig:
    """UMAP configuration.

    Note: Requires umap-learn==0.5.3 for model compatibility.
    """
    n_neighbors: int = 15
    min_dist: float = 0.1
    n_components: int = 2
    metric: str = "euclidean"

@dataclass(frozen=True)
class ClassifierConfig:
    """XGBoost classifier configuration."""
    n_estimators: int = 100
    max_depth: int = 6
    learning_rate: float = 0.1

@dataclass(frozen=True)
class AncestryConfig:
    """Full ancestry prediction configuration."""
    pca: PCAConfig = PCAConfig()
    umap: UMAPConfig = UMAPConfig()
    classifier: ClassifierConfig = ClassifierConfig()

    # Supported ancestry labels (all 10 from current implementation)
    labels: tuple[str, ...] = (
        "AFR",  # African
        "SAS",  # South Asian
        "EAS",  # East Asian
        "EUR",  # European
        "AMR",  # American/Latino
        "AJ",   # Ashkenazi Jewish
        "CAS",  # Central Asian
        "MDE",  # Middle Eastern
        "FIN",  # Finnish
        "AAC",  # African American
    )
```

#### 3.2 `ancestry/results.py`
```python
from dataclasses import dataclass
import pandas as pd

@dataclass
class AncestryPredictions:
    """Ancestry prediction results."""
    predictions: pd.DataFrame  # IID, predicted_ancestry, probabilities...

    def get_ancestry(self, sample_id: str) -> str:
        """Get predicted ancestry for a sample."""
        ...

    def filter_by_ancestry(self, ancestry: str) -> list[str]:
        """Get sample IDs for a given ancestry."""
        ...
```

#### 3.3 `ancestry/model.py`
```python
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional
import pickle
import logging

from umap import UMAP
from xgboost import XGBClassifier

from ..core.genotypes import GenotypeData
from ..core.exceptions import AncestryError
from .config import AncestryConfig
from .reference import ReferencePanel
from .results import AncestryPredictions

logger = logging.getLogger(__name__)

@dataclass
class AncestryModel:
    """Ancestry prediction model with fit/predict interface."""
    config: AncestryConfig = field(default_factory=AncestryConfig)
    _reducer: Optional[UMAP] = field(default=None, repr=False)
    _classifier: Optional[XGBClassifier] = field(default=None, repr=False)
    _is_fitted: bool = field(default=False, repr=False)

    def fit(self, reference: ReferencePanel) -> "AncestryModel":
        """Fit model on reference panel."""
        logger.info("Fitting ancestry model on reference panel")

        # 1. Run PCA on reference
        pca_result = self._run_pca(reference)

        # 2. Fit UMAP on PCA output
        self._reducer = UMAP(
            n_neighbors=self.config.umap.n_neighbors,
            min_dist=self.config.umap.min_dist,
            n_components=self.config.umap.n_components,
        )
        umap_embedding = self._reducer.fit_transform(pca_result)

        # 3. Fit XGBoost classifier
        self._classifier = XGBClassifier(
            n_estimators=self.config.classifier.n_estimators,
            max_depth=self.config.classifier.max_depth,
            learning_rate=self.config.classifier.learning_rate,
        )
        self._classifier.fit(umap_embedding, reference.labels)

        self._is_fitted = True
        logger.info("Ancestry model fitted successfully")
        return self

    def predict(self, data: GenotypeData) -> AncestryPredictions:
        """Predict ancestry for samples."""
        if not self._is_fitted:
            raise AncestryError("Model must be fitted before prediction")

        logger.info(f"Predicting ancestry for {data.sample_count} samples")
        # ... implementation
        ...

    def save(self, path: Path) -> None:
        """Serialize fitted model."""
        if not self._is_fitted:
            raise AncestryError("Cannot save unfitted model")
        with open(path, 'wb') as f:
            pickle.dump(self, f)

    @classmethod
    def load(cls, path: Path) -> "AncestryModel":
        """Load fitted model."""
        with open(path, 'rb') as f:
            model = pickle.load(f)
        if not isinstance(model, cls):
            raise AncestryError(f"Invalid model file: {path}")
        return model
```

#### 3.4 `ancestry/reference.py`
```python
from dataclasses import dataclass
from pathlib import Path
import pandas as pd

from ..core.genotypes import GenotypeData

@dataclass
class ReferencePanel:
    """Reference panel for ancestry prediction."""
    genotypes: GenotypeData
    labels: pd.Series  # Sample ID -> ancestry label

    @classmethod
    def load(cls, geno_path: Path, labels_path: Path) -> "ReferencePanel":
        """Load reference panel from files."""
        ...
```

#### 3.5 Container/Cloud Support

The current implementation supports multiple execution modes that must be preserved:

```python
from dataclasses import dataclass
from typing import Literal, Optional
from enum import Enum

class InferenceMode(Enum):
    """Ancestry inference execution modes."""
    LOCAL = "local"           # Local Python execution
    CONTAINER = "container"   # Docker container
    SINGULARITY = "singularity"  # Singularity container
    CLOUD = "cloud"           # Google Cloud AI Platform

@dataclass(frozen=True)
class InferenceConfig:
    """Configuration for ancestry inference execution."""
    mode: InferenceMode = InferenceMode.LOCAL
    container_image: Optional[str] = None
    cloud_project: Optional[str] = None
    cloud_region: str = "us-central1"
```

CLI flags to support:
- `--container` → Docker execution
- `--singularity` → Singularity execution
- `--cloud` → Google Cloud AI Platform

### Success Criteria
- [ ] `AncestryModel.fit()` works with reference panel
- [ ] `AncestryModel.predict()` produces same results as old code
- [ ] Model can be saved/loaded
- [ ] All ancestry parameters documented in config.py
- [ ] Container and cloud modes work as before
- [ ] `mypy genotools/ancestry/ --strict` passes

---

## Method Mapping (Old → New)

This table shows how current methods map to new functions.

### SampleQC Methods

| Current Method | New Function | Module |
|---------------|--------------|--------|
| `SampleQC.run_callrate_prune()` | `filter_callrate()` | `qc/steps/callrate.py` |
| `SampleQC.run_sex_prune()` | `filter_sex()` | `qc/steps/sex.py` |
| `SampleQC.run_het_prune()` | `filter_heterozygosity()` | `qc/steps/heterozygosity.py` |
| `SampleQC.run_related_prune()` | `filter_relatedness()` | `qc/steps/relatedness.py` |
| `SampleQC.run_confirming_kinship()` | `verify_kinship()` | `qc/steps/relatedness.py` |

### VariantQC Methods

| Current Method | New Function | Module |
|---------------|--------------|--------|
| `VariantQC.run_geno_prune()` | `filter_variant_missingness()` | `qc/steps/variant_missingness.py` |
| `VariantQC.run_case_control_prune()` | `filter_case_control()` | `qc/steps/case_control.py` |
| `VariantQC.run_haplotype_prune()` | `filter_haplotype()` | `qc/steps/haplotype.py` |
| `VariantQC.run_hwe_prune()` | `filter_hwe()` | `qc/steps/hwe.py` |
| `VariantQC.run_ld_prune()` | `prune_ld()` | `qc/steps/ld_prune.py` |

### Ancestry Methods

| Current Method | New Function | Module |
|---------------|--------------|--------|
| `Ancestry.predict()` | `AncestryModel.predict()` | `ancestry/model.py` |
| `Ancestry.train_umap_classifier()` | `AncestryModel.fit()` | `ancestry/model.py` |
| PCA calculation | `run_pca()` | `ancestry/reducers/pca.py` |
| UMAP transformation | `UMAPReducer.fit_transform()` | `ancestry/reducers/umap.py` |

---

## Phase 4: CLI Migration

### Goal
Replace monolithic `__main__.py` with clean CLI layer.

### Issues Resolved
| ID | Issue | Resolution |
|----|-------|------------|
| A3 | No separation of concerns | CLI → Runner → Domain layers |
| M3 | Monolithic main | Split into parser.py, runner.py, output.py |
| C1 | Boolean string parsing | argparse action='store_true' |

### Deliverables

#### 4.1 `cli/parser.py`
```python
import argparse
from pathlib import Path
from typing import Optional
from dataclasses import dataclass

@dataclass
class PipelineArgs:
    """Validated pipeline arguments."""
    geno_path: Path
    out_dir: Path

    # QC options
    run_callrate: bool = True
    run_sex: bool = True
    run_het: bool = True
    run_related: bool = True

    # Ancestry options
    run_ancestry: bool = False
    ref_panel: Optional[Path] = None

    # GWAS options
    run_gwas: bool = False
    pheno_file: Optional[Path] = None

    # Error handling
    warn_only: bool = False  # --warn flag: continue on step failure

    # Inference mode
    use_container: bool = False
    use_singularity: bool = False
    use_cloud: bool = False

def create_parser() -> argparse.ArgumentParser:
    """Create argument parser with proper types."""
    parser = argparse.ArgumentParser(
        prog="genotools",
        description="Genotype QC and ancestry prediction pipeline"
    )

    parser.add_argument("--geno", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)

    # Boolean flags - no more string parsing!
    parser.add_argument("--skip-callrate", action="store_true")
    parser.add_argument("--skip-sex", action="store_true")
    # ... other skip flags

    # Error handling - IMPORTANT: preserve --warn behavior
    parser.add_argument("--warn", action="store_true",
                        help="Continue pipeline on step failure instead of stopping")

    # Inference mode flags
    parser.add_argument("--container", action="store_true",
                        help="Run ancestry inference in Docker container")
    parser.add_argument("--singularity", action="store_true",
                        help="Run ancestry inference in Singularity container")
    parser.add_argument("--cloud", action="store_true",
                        help="Run ancestry inference on Google Cloud AI Platform")

    return parser

def parse_args(args: Optional[list[str]] = None) -> PipelineArgs:
    """Parse and validate arguments."""
    parser = create_parser()
    ns = parser.parse_args(args)
    return PipelineArgs(
        geno_path=ns.geno,
        out_dir=ns.out,
        run_callrate=not ns.skip_callrate,
        warn_only=ns.warn,
        use_container=ns.container,
        use_singularity=ns.singularity,
        use_cloud=ns.cloud,
        # ...
    )
```

#### 4.2 `cli/runner.py`
```python
from pathlib import Path
import logging

from ..core.genotypes import GenotypeData
from ..core.logging import setup_logging
from ..core.exceptions import QCError, AncestryError
from ..qc.pipeline import QCPipeline
from ..qc.steps import filter_callrate, filter_sex, ...
from ..qc.config import CallrateConfig, SexConfig, ...
from ..ancestry.model import AncestryModel
from ..ancestry.config import InferenceMode
from .parser import PipelineArgs
from .output import write_results

logger = logging.getLogger(__name__)

def run_pipeline(args: PipelineArgs) -> dict:
    """Main pipeline orchestration."""
    setup_logging(log_file=args.out_dir / "genotools.log")

    # Load data
    data = GenotypeData.from_path(args.geno_path)
    logger.info(f"Loaded {data.sample_count} samples, {data.variant_count} variants")

    results = {}
    pass_fail = {}  # Track step success/failure

    # Build and run QC pipeline
    if any([args.run_callrate, args.run_sex, args.run_het, args.run_related]):
        steps = []
        if args.run_callrate:
            steps.append(("callrate", filter_callrate, CallrateConfig()))
        if args.run_sex:
            steps.append(("sex", filter_sex, SexConfig()))
        # ... etc

        qc_pipeline = QCPipeline(steps=steps, warn_only=args.warn_only)
        try:
            qc_result = qc_pipeline.run(data, args.out_dir / "qc")
            results["qc"] = qc_result
            data = qc_result.output
            pass_fail["qc"] = {"status": True}
        except QCError as e:
            logger.error(f"QC failed: {e}")
            pass_fail["qc"] = {"status": False, "error": str(e)}
            if not args.warn_only:
                raise

    # Run ancestry if requested
    if args.run_ancestry:
        # Determine inference mode
        if args.use_container:
            mode = InferenceMode.CONTAINER
        elif args.use_singularity:
            mode = InferenceMode.SINGULARITY
        elif args.use_cloud:
            mode = InferenceMode.CLOUD
        else:
            mode = InferenceMode.LOCAL

        try:
            model = AncestryModel.load(args.ref_panel)
            ancestry_result = model.predict(data, mode=mode)
            results["ancestry"] = ancestry_result
            pass_fail["ancestry"] = {"status": True}
        except AncestryError as e:
            logger.error(f"Ancestry prediction failed: {e}")
            pass_fail["ancestry"] = {"status": False, "error": str(e)}
            if not args.warn_only:
                raise

    results["pass_fail"] = pass_fail

    # Write outputs
    write_results(results, args.out_dir)

    return results
```

#### 4.3 `cli/output.py`
```python
from pathlib import Path
import json
from dataclasses import asdict
from typing import Any

from ..qc.results import QCResult
from ..ancestry.results import AncestryPredictions

def write_results(results: dict, out_dir: Path) -> None:
    """Write pipeline results to files.

    Output JSON structure must match current format:
    {
        "ancestry_counts": {"EUR": 100, "AFR": 50, ...},
        "ancestry_labels": [{"#FID": "...", "IID": "...", "label": "EUR"}, ...],
        "QC": [
            {"step": "callrate_prune", "pruned_count": 5, "metric": "outlier_count", ...}
        ],
        "GWAS": [
            {"value": 1.02, "metric": "lambda", "ancestry": "EUR"}
        ],
        "pruned_samples": [...],
        "related_samples": [...]
    }
    """
    output = {}

    # QC results
    if "qc" in results:
        qc_result: QCResult = results["qc"]
        output["QC"] = [
            {
                "step": name,
                "pruned_count": r.samples_removed + r.variants_removed,
                "metric": "outlier_count",
                **r.metrics
            }
            for name, r in qc_result.step_results
        ]
        output["pruned_samples"] = _collect_pruned_samples(qc_result)

    # Ancestry results
    if "ancestry" in results:
        ancestry_result: AncestryPredictions = results["ancestry"]
        output["ancestry_counts"] = ancestry_result.predictions["predicted_ancestry"].value_counts().to_dict()
        output["ancestry_labels"] = ancestry_result.predictions[["#FID", "IID", "predicted_ancestry"]].rename(
            columns={"predicted_ancestry": "label"}
        ).to_dict(orient="records")

    # GWAS results
    if "gwas" in results:
        output["GWAS"] = results["gwas"]

    # Write JSON output
    with open(out_dir / "results.json", "w") as f:
        json.dump(output, f, indent=2)

def _collect_pruned_samples(qc_result: QCResult) -> list[dict[str, str]]:
    """Collect all pruned samples across QC steps."""
    # Implementation reads from each step's pruned_samples file
    ...
```

#### 4.4 Updated `__main__.py`
```python
"""GenoTools CLI entry point."""
from .cli.parser import parse_args
from .cli.runner import run_pipeline

def main():
    args = parse_args()
    run_pipeline(args)

if __name__ == "__main__":
    main()
```

### Success Criteria
- [ ] CLI works identically to old interface
- [ ] No boolean string parsing
- [ ] Clean separation: parser → runner → domain
- [ ] `python -m genotools --help` shows proper help
- [ ] Easy to add new CLI options
- [ ] `--warn` flag continues pipeline on step failure
- [ ] `--container`, `--singularity`, `--cloud` flags work
- [ ] Output JSON matches expected structure (see `cli/output.py`)
- [ ] `pass_fail` dict available for step status inspection

---

## Phase 5: GWAS & Cleanup

### Goal
Migrate GWAS module and remove deprecated code.

### Deliverables

#### 5.1 GWAS Migration
- Apply same patterns as QC (pure functions, typed configs)
- Determine if step-based or different pattern fits better

#### 5.2 Deprecation
1. Add deprecation warnings to old modules
2. Update all documentation
3. Update examples/tests to use new API

#### 5.3 Removal
1. Remove old `qc.py`, `ancestry.py`, `pipeline.py`
2. Remove `REFACTOR.md` and `REFACTOR_ISSUES.md`
3. Final cleanup pass

### Success Criteria
- [ ] All tests pass with new implementation
- [ ] No references to old modules
- [ ] Documentation updated
- [ ] Clean git history with squashed commits per phase

---

## Open Questions

1. **Backwards compatibility** - Should we maintain the old API surface with deprecation warnings?
2. **Config file format** - YAML, TOML, or Python dataclasses only?
3. **Reference panel packaging** - Bundle with package or download on demand?
4. **GWAS pattern** - Same as QC (steps), or different (more statistical focus)?

---

## Progress Tracking

| Phase | Status | Issues Resolved |
|-------|--------|-----------------|
| 0: Regression Testing | Not Started | (enables validation) |
| 1: Core Foundation | Not Started | A1, P1, R1, R2, R3, D1, D2, D3 |
| 2: QC Migration | Not Started | A1, A2, M1, M2, D1 |
| 3: Ancestry Migration | Not Started | A2, M1, D1 |
| 4: CLI Migration | Not Started | A3, M3, C1 |
| 5: GWAS & Cleanup | Not Started | - |
