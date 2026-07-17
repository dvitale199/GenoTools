# Decouple new CLI from legacy — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make the new modular pipeline (`cli/`, `core/`, `qc/`, `gwas/`, new `ancestry/`) import zero top-level legacy modules — removing the four `utils.py` imports and the `ancestry.py` importlib hack from `cli/runner.py`, and porting the three legacy `Ancestry` helper methods the new ancestry path still uses.

**Architecture:** Two phases. **6a** relocates four small helpers into `core/` (legacy-free reimplementations) and `git mv`s `ancestry.py` → `ancestry/legacy.py` so it imports normally. **6b** ports `get_raw_files` / `get_common_snps` / `split_cohort_ancestry` / `clean_up` into the new `ancestry/` package as pure functions, then rewires the `genotools-new` path to use them so it never touches legacy `Ancestry`. All ports are *faithful* (same PLINK commands, same pandas ops), proven by **differential tests** that run legacy and new on identical synthetic inputs and assert byte-identical output.

**Tech Stack:** Python 3.11, PLINK/PLINK2 (via `genotools/core/executors.py`), pandas, pytest.

## Global Constraints

- **Branch target:** `refactor/main` (this feature branch → PR into `refactor/main`; **never `main`**). Already on `refactor/decouple-from-legacy`.
- **Behavior-preserving.** The default `genotools` engine (legacy `Ancestry`) and all QC/GWAS output must be unchanged. Only `genotools-new`'s internals change, and its *output* must not (faithful port).
- **Legacy files stay.** Do **not** delete `utils.py`, `ancestry/legacy.py`, `gwas.py`, or `imputation.py` — they remain for the legacy cluster and are deleted wholesale in Phase 5/6. Leave their originals (`gt_header`, `bfiles_to_pfiles`, `vcf_to_pfiles`, `upfront_check`, `get_common_snps`) untouched.
- **New code uses `core/` only.** Every new/relocated function uses `core/executors` (`run_command`, `get_plink`, `get_plink2`) and `core/exceptions`; zero imports from `genotools.utils`.
- **`shell_do(cmd)` → `run_command(cmd, tool_name="plink2")`** (or `tool_name="plink"` for PLINK1.9). These produce equivalent argv for the unquoted commands here (`shell_do` uses `.split()`, `run_command` uses `shlex.split()`).
- **PLINK2 psam-cols flag** stays exactly `--make-pgen psam-cols=fid,parents,sex,pheno1,phenos` wherever it appears (`core.executors.PLINK2_PSAM_COLS`).
- **Deferred (do NOT implement here):** the `upfront_check` data-driven skip logic (tracked separately); flipping the default to new ancestry; deleting legacy.
- **Run tests from the repo root** (existing tests reference `tests/data/synthetic/...` relative to CWD).

**Design spec:** `docs/superpowers/specs/2026-07-17-decouple-from-legacy-design.md`.

## File Structure

**Create:**
- `genotools/core/validation.py` — `validate_input()` (pre-flight input validation + data breakdown).
- `genotools/ancestry/preprocessing.py` — `get_common_snps()`, `get_raw_files()`, `clean_up_files()`.
- `genotools/ancestry/cohort.py` — `split_cohort_by_ancestry()`.
- `genotools/ancestry/legacy.py` — via `git mv genotools/ancestry.py` (retained legacy `Ancestry`).
- `tests/unit/test_ancestry/conftest.py` — synthetic reference-panel fixtures.
- `tests/unit/test_ancestry/test_preprocessing.py` — differential tests for the preprocessing ports.
- `tests/unit/test_ancestry/test_cohort.py` — differential test for the cohort-split port.

**Modify:**
- `genotools/core/logging.py` — add `banner()`.
- `genotools/core/genotypes.py` — add `GenotypeData.from_vcf()`.
- `genotools/core/__init__.py` — export `banner`, `validate_input`.
- `genotools/cli/runner.py` — rewire `_setup_logging`, `_convert_input_format`, `_initialize_modules` (kill importlib hack), `_run_ancestry_prediction_new` / `_run_training_mode` / `_run_inference_mode`.
- `tests/unit/test_core.py` — add `TestBanner`, `TestFromVcf`, `TestValidation`.

**Delete (via git mv only):** `genotools/ancestry.py` (moves to `ancestry/legacy.py`).

---

# Phase 6a — Import decoupling (behavior-preserving)

### Task A1: `banner()` in `core/logging.py`

**Files:**
- Modify: `genotools/core/logging.py`
- Modify: `genotools/core/__init__.py`
- Modify: `genotools/cli/runner.py:298-300` (`_setup_logging`)
- Test: `tests/unit/test_core.py`

**Interfaces:**
- Produces: `banner() -> str` — returns the GenoTools ASCII header (identical string to legacy `genotools.utils.gt_header()`).

- [ ] **Step 1: Write the failing test**

Add to `tests/unit/test_core.py`:

```python
class TestBanner:
    """Tests for the ASCII banner."""

    def test_banner_returns_ascii_header(self):
        from genotools.core.logging import banner
        b = banner()
        assert isinstance(b, str)
        assert "█" in b            # box-drawing art present
        assert b.count("\n") >= 5  # multi-line

    def test_banner_matches_legacy(self):
        # Faithful copy of legacy gt_header (differential; drop when utils.py is removed)
        from genotools.core.logging import banner
        from genotools.utils import gt_header
        assert banner() == gt_header()
```

- [ ] **Step 2: Run test to verify it fails**

Run: `python -m pytest tests/unit/test_core.py::TestBanner -v`
Expected: FAIL with `ImportError: cannot import name 'banner'`.

- [ ] **Step 3: Add `banner()` to `core/logging.py`**

Append to `genotools/core/logging.py` (copy the exact string from `genotools/utils.py:30-40`):

```python
def banner() -> str:
    """Return the GenoTools ASCII banner (for the consolidated log header)."""
    return """
     ██████╗ ███████╗███╗  ██╗ █████╗ ████████╗ █████╗  █████╗ ██╗      ██████╗
    ██╔════╝ ██╔════╝████╗ ██║██╔══██╗╚══██╔══╝██╔══██╗██╔══██╗██║     ██╔════╝
    ██║  ██╗ █████╗  ██╔██╗██║██║  ██║   ██║   ██║  ██║██║  ██║██║     ╚█████╗
    ██║  ╚██╗██╔══╝  ██║╚████║██║  ██║   ██║   ██║  ██║██║  ██║██║      ╚═══██╗
    ╚██████╔╝███████╗██║ ╚███║╚█████╔╝   ██║   ╚█████╔╝╚█████╔╝███████╗██████╔╝
    ╚═════╝ ╚══════╝╚═╝  ╚══╝ ╚════╝    ╚═╝    ╚════╝  ╚════╝ ╚══════╝╚═════╝
    """
```

> Copy the string byte-for-byte from `utils.py:32-39` so `test_banner_matches_legacy` passes.

- [ ] **Step 4: Export from `core/__init__.py`**

In `genotools/core/__init__.py`, add `banner` to the `from .logging import (...)` block and to `__all__`.

- [ ] **Step 5: Rewire the runner**

In `genotools/cli/runner.py`, `_setup_logging()`, replace:

```python
        # Create new log files with header
        from ..utils import gt_header

        header = gt_header()
```

with:

```python
        # Create new log files with header
        from ..core.logging import banner

        header = banner()
```

- [ ] **Step 6: Run tests**

Run: `python -m pytest tests/unit/test_core.py::TestBanner -v`
Expected: PASS (2 tests).

- [ ] **Step 7: Commit**

```bash
git add genotools/core/logging.py genotools/core/__init__.py genotools/cli/runner.py tests/unit/test_core.py
git commit -m "refactor: add core.logging.banner, drop utils.gt_header import from runner"
```

---

### Task A2: `GenotypeData.from_vcf()` + rewire input conversion

**Files:**
- Modify: `genotools/core/genotypes.py`
- Modify: `genotools/cli/runner.py:314-325` (`_convert_input_format`, bfile + vcf branches)
- Test: `tests/unit/test_core.py`

**Interfaces:**
- Consumes: `GenotypeData.from_path`, `GenotypeData.to_pfile` (existing).
- Produces: `GenotypeData.from_vcf(vcf_path: Path | str, plink2_path: Optional[Path | str] = None) -> "GenotypeData"` — converts VCF → bfile → pfile at the VCF prefix (VCF path with `.vcf*` stripped), removes the intermediate bed/bim/fam, returns a `GenotypeData` for the pfile. Faithful to `utils.vcf_to_pfiles` (`utils.py:92-111`).

- [ ] **Step 1: Write the failing tests**

Add to `tests/unit/test_core.py` inside `TestGenotypeData` (uses the existing `test_data_path` fixture):

```python
    def test_from_vcf_creates_pfile(self, test_data_path: Path, tmp_path: Path):
        """from_vcf converts a VCF to pfile and cleans up intermediates."""
        from genotools.core import GenotypeData
        from genotools.core.executors import get_plink2
        # Make a VCF from the synthetic pfile
        vcf_prefix = tmp_path / "from_vcf_input"
        import subprocess
        subprocess.run(
            [str(get_plink2()), "--pfile", str(test_data_path),
             "--autosome", "--export", "vcf", "--out", str(vcf_prefix)],
            check=True, capture_output=True,
        )
        vcf_file = vcf_prefix.with_suffix(".vcf")
        assert vcf_file.exists()

        data = GenotypeData.from_vcf(vcf_file)

        assert data.format == "pfile"
        assert (vcf_prefix.with_suffix(".pgen")).exists()
        # intermediate bfile removed
        assert not (vcf_prefix.with_suffix(".bed")).exists()

    def test_from_vcf_raises_on_missing(self, tmp_path: Path):
        from genotools.core import GenotypeData
        with pytest.raises(FileNotFoundError):
            GenotypeData.from_vcf(tmp_path / "nope.vcf")
```

- [ ] **Step 2: Run to verify failure**

Run: `python -m pytest "tests/unit/test_core.py::TestGenotypeData::test_from_vcf_creates_pfile" -v`
Expected: FAIL with `AttributeError: type object 'GenotypeData' has no attribute 'from_vcf'`.

- [ ] **Step 3: Implement `from_vcf`**

Add this classmethod to `GenotypeData` in `genotools/core/genotypes.py` (place after `from_path`). It mirrors `utils.vcf_to_pfiles`:

```python
    @classmethod
    def from_vcf(
        cls,
        vcf_path: Path | str,
        plink2_path: Optional[Path | str] = None,
    ) -> "GenotypeData":
        """Convert a VCF to pfile format at the VCF's prefix.

        Faithful reimplementation of legacy utils.vcf_to_pfiles: VCF -> bfile
        -> pfile, removing the intermediate bed/bim/fam. Returns a GenotypeData
        for the resulting pfile.

        Raises:
            FileNotFoundError: if the VCF is missing or a conversion produces no output.
            ExternalToolError: if PLINK2 fails.
        """
        from genotools.core.executors import run_command, get_plink2

        vcf_path = Path(vcf_path)
        if not vcf_path.is_file():
            raise FileNotFoundError(f"{vcf_path} does not exist.")

        if plink2_path is None:
            plink2_path = get_plink2()

        # Prefix = VCF path with .vcf/.vcf.gz stripped (matches legacy split('.vcf'))
        prefix = Path(str(vcf_path).split(".vcf")[0])

        run_command(
            f"{plink2_path} --vcf {vcf_path} --make-bed --out {prefix}",
            tool_name="plink2",
        )
        if not prefix.with_suffix(".bed").is_file():
            raise FileNotFoundError(
                f"{prefix} bed/bim/fam files do not exist. Conversion from VCF failed"
            )

        bfile = cls.from_path(prefix)
        pfile = bfile.to_pfile(prefix, plink2_path=plink2_path)

        if prefix.with_suffix(".pgen").is_file():
            for ext in (".bed", ".bim", ".fam"):
                p = prefix.with_suffix(ext)
                if p.exists():
                    p.unlink()
        else:
            raise FileNotFoundError(
                f"{prefix} pgen/pvar/psam files do not exist. "
                f"Conversion from bed/bim/fam failed."
            )
        return pfile
```

> Note: `from_path(prefix)` after the VCF→bed step detects `bfile` (only bed/bim/fam exist at `prefix`), so `to_pfile` uses `--bfile`, matching legacy.

- [ ] **Step 4: Rewire the runner's conversion branches**

In `genotools/cli/runner.py`, `_convert_input_format()`, replace the whole `if/elif` conversion block:

```python
        if input_format == "bfile":
            from ..utils import bfiles_to_pfiles

            bfiles_to_pfiles(bfile_path=str(self.args.input.bfile))
        elif input_format == "vcf":
            from ..utils import vcf_to_pfiles

            vcf_to_pfiles(vcf_path=str(self.args.input.vcf))
```

with:

```python
        from ..core import GenotypeData

        if input_format == "bfile":
            bfile = str(self.args.input.bfile)
            GenotypeData.from_path(bfile).to_pfile(bfile)
        elif input_format == "vcf":
            GenotypeData.from_vcf(str(self.args.input.vcf))
```

- [ ] **Step 5: Run tests**

Run: `python -m pytest "tests/unit/test_core.py::TestGenotypeData" -v`
Expected: PASS (all, including the two new tests).

- [ ] **Step 6: Commit**

```bash
git add genotools/core/genotypes.py genotools/cli/runner.py tests/unit/test_core.py
git commit -m "refactor: convert bfile/vcf input via core.GenotypeData, drop utils conversion imports"
```

---

### Task A3: `core/validation.validate_input()` + rewire

**Files:**
- Create: `genotools/core/validation.py`
- Modify: `genotools/core/__init__.py`
- Modify: `genotools/cli/runner.py:327-332` (`_convert_input_format`, validation block)
- Test: `tests/unit/test_core.py`

**Interfaces:**
- Produces: `validate_input(geno_path: str | Path, out_path: str | Path, skip_fails: bool = False) -> None`. Raises `ValidationError` if `{out_path}_all_logs.log` exists and not `skip_fails`; raises `FileNotFoundError` if `{geno_path}.pgen` missing; raises `ValidationError` if the psam is missing `SEX` or `PHENO1`; otherwise prints the genetic-sex/phenotype breakdown. **Does not** implement legacy's data-driven step-skips (deferred).

- [ ] **Step 1: Write the failing tests**

Add a new class to `tests/unit/test_core.py`:

```python
class TestValidation:
    """Tests for core.validation.validate_input."""

    @pytest.fixture
    def geno(self) -> Path:
        return Path("tests/data/synthetic/genotools_test")

    def test_passes_on_valid_input(self, geno: Path, tmp_path: Path, capsys):
        from genotools.core.validation import validate_input
        validate_input(geno, tmp_path / "out", skip_fails=False)  # no raise
        out = capsys.readouterr().out
        assert "breakdown" in out.lower()

    def test_raises_on_missing_pgen(self, tmp_path: Path):
        from genotools.core.validation import validate_input
        with pytest.raises(FileNotFoundError):
            validate_input(tmp_path / "nope", tmp_path / "out", skip_fails=False)

    def test_raises_when_log_exists(self, geno: Path, tmp_path: Path):
        from genotools.core.validation import validate_input
        from genotools.core.exceptions import ValidationError
        log = tmp_path / "out_all_logs.log"
        log.write_text("prior run\n")
        with pytest.raises(ValidationError):
            validate_input(geno, tmp_path / "out", skip_fails=False)

    def test_skip_fails_ignores_existing_log(self, geno: Path, tmp_path: Path):
        from genotools.core.validation import validate_input
        (tmp_path / "out_all_logs.log").write_text("prior\n")
        validate_input(geno, tmp_path / "out", skip_fails=True)  # no raise

    def test_raises_on_missing_sex_column(self, geno: Path, tmp_path: Path):
        from genotools.core.validation import validate_input
        from genotools.core.exceptions import ValidationError
        import pandas as pd
        # Build a pfile-shaped fixture whose psam lacks SEX
        bad = tmp_path / "bad"
        bad.with_suffix(".pgen").write_bytes((geno.with_suffix(".pgen")).read_bytes())
        bad.with_suffix(".pvar").write_bytes((geno.with_suffix(".pvar")).read_bytes())
        psam = pd.read_csv(geno.with_suffix(".psam"), sep=r"\s+")
        psam.drop(columns=["SEX"]).to_csv(bad.with_suffix(".psam"), sep="\t", index=False)
        with pytest.raises(ValidationError):
            validate_input(bad, tmp_path / "out2", skip_fails=False)
```

- [ ] **Step 2: Run to verify failure**

Run: `python -m pytest tests/unit/test_core.py::TestValidation -v`
Expected: FAIL with `ModuleNotFoundError: No module named 'genotools.core.validation'`.

- [ ] **Step 3: Implement `core/validation.py`**

Create `genotools/core/validation.py`. This ports the *checks + breakdown* of `utils.upfront_check` (`utils.py:115-159`), raising `ValidationError` instead of `ValueError`/`KeyError`, and omits the deferred skip logic:

```python
"""Pre-flight validation of pipeline input (legacy-free).

Ports the validation + data-breakdown behavior of legacy utils.upfront_check.
The legacy data-driven step-skip logic is intentionally NOT ported (deferred;
see the design spec) — the current runner already discarded it.
"""

from pathlib import Path
from typing import Union

import pandas as pd

from genotools.core.exceptions import ValidationError


def validate_input(
    geno_path: Union[str, Path],
    out_path: Union[str, Path],
    skip_fails: bool = False,
) -> None:
    """Validate pipeline input and print a data breakdown.

    Raises:
        ValidationError: if the output log already exists (and not skip_fails),
            or the psam is missing SEX/PHENO1.
        FileNotFoundError: if the input pgen is missing.
    """
    geno_path = str(geno_path)
    out_path = str(out_path)

    if Path(f"{out_path}_all_logs.log").is_file() and not skip_fails:
        raise ValidationError(
            f"{out_path}_all_logs.log exists, which means the pipeline has "
            f"previously been run on this output file! Rerun with --skip_fails "
            f"to ignore this, or write output to a new file name."
        )

    if not Path(f"{geno_path}.pgen").is_file():
        raise FileNotFoundError(f"{geno_path} does not exist.")

    sam = pd.read_csv(f"{geno_path}.psam", sep=r"\s+")
    var = pd.read_csv(
        f"{geno_path}.pvar", delimiter="\t", comment="#", header=None,
        usecols=range(5), names=["#CHROM", "POS", "ID", "REF", "ALT"],
        low_memory=False,
    )

    if "SEX" not in sam.columns:
        raise ValidationError(
            f"{geno_path}.psam is missing SEX column. GenoTools requires a SEX column."
        )
    if "PHENO1" not in sam.columns:
        raise ValidationError(
            f"{geno_path}.psam is missing PHENO1 column. GenoTools requires a PHENO1 column."
        )

    sex_counts = sam["SEX"].value_counts().to_dict()
    pheno_counts = sam["PHENO1"].value_counts().to_dict()

    print("Your data has the following breakdown:")
    print("- Genetic Sex:")
    for sex in sex_counts:
        if sex == 1:
            print(f"{sex_counts[sex]} Males \n")
        if sex == 2:
            print(f"{sex_counts[sex]} Females \n")
        if sex in (0, -9):
            print(f"{sex_counts[sex]} Unknown \n")

    print("- Phenotypes:")
    for pheno in pheno_counts:
        if pheno == 2:
            print(f"{pheno_counts[pheno]} Cases \n")
        if pheno == 1:
            print(f"{pheno_counts[pheno]} Controls \n")
        if pheno in (0, -9):
            print(f"{pheno_counts[pheno]} Missing \n")
```

- [ ] **Step 4: Export from `core/__init__.py`**

Add `from .validation import validate_input` and `"validate_input"` to `__all__` in `genotools/core/__init__.py`.

- [ ] **Step 5: Rewire the runner**

In `genotools/cli/runner.py`, `_convert_input_format()`, replace:

```python
        # Run upfront validation
        if not self.args.output.skip_fails:
            from ..utils import upfront_check

            legacy_dict = self.args.to_legacy_dict()
            upfront_check(str(self.args.geno_path), legacy_dict)
```

with:

```python
        # Run upfront validation
        from ..core.validation import validate_input

        validate_input(
            self.args.geno_path,
            self.args.out_path,
            skip_fails=self.args.output.skip_fails,
        )
```

> The `if not skip_fails` guard moves *inside* `validate_input` (which still prints the breakdown and checks columns even when `skip_fails=True`, matching legacy `upfront_check`'s early section).

- [ ] **Step 6: Run tests**

Run: `python -m pytest tests/unit/test_core.py::TestValidation -v`
Expected: PASS (5 tests).

- [ ] **Step 7: Commit**

```bash
git add genotools/core/validation.py genotools/core/__init__.py genotools/cli/runner.py tests/unit/test_core.py
git commit -m "refactor: add core.validation.validate_input, drop utils.upfront_check import"
```

---

### Task A4: `git mv ancestry.py` → `ancestry/legacy.py`, kill importlib hack

**Files:**
- Move: `genotools/ancestry.py` → `genotools/ancestry/legacy.py`
- Modify: `genotools/cli/runner.py:258-272` (`_initialize_modules`)
- Test: `tests/unit/test_cli/test_runner.py` (add an import-smoke test)

**Interfaces:**
- Produces: `genotools.ancestry.legacy.Ancestry` importable normally (no importlib).

- [ ] **Step 1: Move the file with git**

```bash
git mv genotools/ancestry.py genotools/ancestry/legacy.py
```

> Its internal imports (`from genotools.utils import ...`, `from genotools.dependencies import ...`) resolve unchanged from the new location.

- [ ] **Step 2: Write the failing test**

Add to `tests/unit/test_cli/test_runner.py` (new class at end of file):

```python
class TestLegacyAncestryImport:
    def test_legacy_ancestry_importable(self):
        from genotools.ancestry.legacy import Ancestry
        assert Ancestry is not None

    def test_entry_points_import(self):
        from genotools.cli import main, main_new
        assert callable(main) and callable(main_new)
```

- [ ] **Step 3: Run to verify current state**

Run: `python -m pytest tests/unit/test_cli/test_runner.py::TestLegacyAncestryImport -v`
Expected: `test_legacy_ancestry_importable` PASSES already (git mv done); this confirms the move. Proceed to replace the hack.

- [ ] **Step 4: Replace the importlib hack in the runner**

In `genotools/cli/runner.py`, `_initialize_modules()`, replace the whole block:

```python
        # Load legacy ancestry module (ancestry.py) via importlib
        # since genotools.ancestry is the new package directory
        import importlib.util
        import sys
        from pathlib import Path as PathLib

        genotools_dir = PathLib(__file__).parent.parent
        ancestry_spec = importlib.util.spec_from_file_location(
            "genotools.legacy_ancestry", genotools_dir / "ancestry.py"
        )
        legacy_ancestry = importlib.util.module_from_spec(ancestry_spec)  # type: ignore
        sys.modules["genotools.legacy_ancestry"] = legacy_ancestry
        ancestry_spec.loader.exec_module(legacy_ancestry)  # type: ignore
        Ancestry = legacy_ancestry.Ancestry
        self._ancestry = Ancestry()
```

with:

```python
        # Legacy ancestry engine (still the default; A/B control until Phase 5/6)
        from ..ancestry.legacy import Ancestry

        self._ancestry = Ancestry()
```

- [ ] **Step 5: Run tests + import/CLI smoke**

Run:
```bash
python -m pytest tests/unit/test_cli/test_runner.py::TestLegacyAncestryImport -v
python -c "import genotools, genotools.core, genotools.qc, genotools.ancestry, genotools.gwas, genotools.cli"
python -m genotools --help >/dev/null && echo "CLI OK"
```
Expected: tests PASS; imports succeed; `CLI OK`.

- [ ] **Step 6: Verify no new-code legacy imports remain (6a gate)**

Run:
```bash
grep -rn "from ..utils import\|from genotools.utils import\|importlib" genotools/cli genotools/core genotools/qc genotools/gwas 2>/dev/null | grep "\.py:"
```
Expected: **no output** (all four `utils` imports and the importlib hack are gone).

- [ ] **Step 7: Commit**

```bash
git add genotools/ancestry/legacy.py genotools/cli/runner.py tests/unit/test_cli/test_runner.py
git commit -m "refactor: relocate legacy ancestry.py to ancestry/legacy.py, remove importlib hack"
```

---

### Task A5: Phase 6a regression gate

**Files:** none (verification only).

- [ ] **Step 1: Full unit + regression suite**

Run: `python -m pytest tests/unit tests/regression -q`
Expected: all green (prior count + the new A1–A4 tests; no failures).

- [ ] **Step 2: Parity (if `.venv-stable` present)**

Run: `python -m pytest tests/regression/test_parity.py -v`
Expected: PASS, or SKIP if `.venv-stable` absent. (Proves 6a's conversion + validation produce identical pfiles.)

- [ ] **Step 3: Commit (only if any fixup was needed)**

If everything was already green, nothing to commit. Otherwise commit fixups with `refactor: fix 6a regression`.

---

# Phase 6b — Finish the new ancestry path (faithful port)

### Task B1: Synthetic reference-panel fixtures

**Files:**
- Create: `tests/unit/test_ancestry/conftest.py`

**Interfaces:**
- Produces pytest fixtures:
  - `synth_geno_pfile -> Path` — the existing synthetic pfile prefix.
  - `synth_ref_bfile(tmp_path_factory) -> Path` — a bfile (bed/bim/fam) copy of the synthetic data, used as a reference panel.
  - `synth_ref_labels(tmp_path_factory) -> Path` — a tab-separated `FID IID label` file assigning alternating `EUR`/`AFR` labels to the ref-panel samples (no header).

- [ ] **Step 1: Write the conftest fixtures**

Create `tests/unit/test_ancestry/conftest.py`:

```python
from pathlib import Path

import pandas as pd
import pytest

from genotools.core.executors import run_command, get_plink2

SYNTH = Path("tests/data/synthetic/genotools_test")


@pytest.fixture(scope="session")
def synth_geno_pfile() -> Path:
    assert SYNTH.with_suffix(".pgen").exists(), "run from repo root"
    return SYNTH


@pytest.fixture(scope="session")
def synth_ref_bfile(tmp_path_factory) -> Path:
    """A bfile copy of the synthetic data, used as a reference panel."""
    out = tmp_path_factory.mktemp("refpanel") / "ref"
    run_command(
        f"{get_plink2()} --pfile {SYNTH} --make-bed --out {out}",
        tool_name="plink2",
    )
    return out


@pytest.fixture(scope="session")
def synth_ref_labels(synth_ref_bfile: Path, tmp_path_factory) -> Path:
    """FID/IID/label for the ref-panel samples (alternating EUR/AFR)."""
    fam = pd.read_csv(f"{synth_ref_bfile}.fam", sep=r"\s+", header=None)
    labels = pd.DataFrame({
        "FID": fam[0],
        "IID": fam[1],
        "label": ["EUR" if i % 2 == 0 else "AFR" for i in range(len(fam))],
    })
    path = tmp_path_factory.mktemp("reflabels") / "ref_labels.txt"
    labels.to_csv(path, sep="\t", header=False, index=False)
    return path
```

- [ ] **Step 2: Smoke-run the fixtures**

Add a temporary check and run it (then delete the check):

```bash
python -m pytest tests/unit/test_ancestry -q --collect-only
```
Expected: collection succeeds (fixtures import cleanly). Then verify plink2 conversion works:
```bash
python -c "from genotools.core.executors import get_plink2; print(get_plink2())"
```
Expected: prints a path.

- [ ] **Step 3: Commit**

```bash
git add tests/unit/test_ancestry/conftest.py
git commit -m "test: synthetic reference-panel fixtures for ancestry preprocessing"
```

---

### Task B2: Port `get_common_snps` → `ancestry/preprocessing.py`

**Files:**
- Create: `genotools/ancestry/preprocessing.py`
- Test: `tests/unit/test_ancestry/test_preprocessing.py`

**Interfaces:**
- Consumes: `core.executors.run_command`, `get_plink`, `get_plink2`.
- Produces: `get_common_snps(geno_path1: str, geno_path2: str, out_name: str) -> dict` returning `{"common_snps": <path>, "bed": <out_name>}`. Faithful reimplementation of `utils.get_common_snps` (`utils.py:404-459`).

- [ ] **Step 1: Write the differential failing test**

Create `tests/unit/test_ancestry/test_preprocessing.py`:

```python
from pathlib import Path

import pandas as pd


def _read_bim_snps(prefix: str) -> set:
    bim = pd.read_csv(f"{prefix}.bim", sep=r"\s+", header=None)
    return set(bim[1])


def test_get_common_snps_matches_legacy(synth_ref_bfile, tmp_path):
    """New port yields the same extracted SNP set + common_snps file as legacy."""
    from genotools.ancestry.preprocessing import get_common_snps as new_fn
    from genotools.utils import get_common_snps as legacy_fn

    # Use the same bfile as both inputs (guarantees a nonempty common set)
    g1, g2 = str(synth_ref_bfile), str(synth_ref_bfile)

    legacy_out = tmp_path / "legacy_common"
    new_out = tmp_path / "new_common"

    legacy_res = legacy_fn(g1, g2, str(legacy_out))
    new_res = new_fn(g1, g2, str(new_out))

    assert set(new_res.keys()) == set(legacy_res.keys())
    assert _read_bim_snps(str(new_out)) == _read_bim_snps(str(legacy_out))
    assert (Path(f"{new_out}.common_snps").read_text()
            == Path(f"{legacy_out}.common_snps").read_text())
```

- [ ] **Step 2: Run to verify failure**

Run: `python -m pytest tests/unit/test_ancestry/test_preprocessing.py::test_get_common_snps_matches_legacy -v`
Expected: FAIL with `ModuleNotFoundError: No module named 'genotools.ancestry.preprocessing'`.

- [ ] **Step 3: Implement the port**

Create `genotools/ancestry/preprocessing.py` with `get_common_snps`. This is `utils.py:404-459` with `plink_exec`/`plink2_exec` → resolved executables and `shell_do(cmd)` → `run_command(cmd, tool_name=...)`:

```python
"""Ancestry PLINK preprocessing (legacy-free port of the Ancestry helpers).

Faithful reimplementations of legacy Ancestry.get_raw_files / clean_up and
utils.get_common_snps, using core.executors. Genotype output is identical to
the legacy versions (proven by differential tests).
"""

import os

import numpy as np
import pandas as pd

from genotools.core.executors import run_command, get_plink, get_plink2


def get_common_snps(geno_path1: str, geno_path2: str, out_name: str) -> dict:
    """Extract SNPs common to two bfiles from geno_path1. Returns output paths."""
    print("Getting Common SNPs")
    plink = get_plink()
    plink2 = get_plink2()

    bim1 = pd.read_csv(f"{geno_path1}.bim", sep="\t", header=None, dtype={0: str}, low_memory=False)
    bim1.columns = ["chr", "rsid", "kb", "pos", "a1", "a2"]
    bim2 = pd.read_csv(f"{geno_path2}.bim", sep="\t", header=None, dtype={0: str}, low_memory=False)
    bim2.columns = ["chr", "rsid", "kb", "pos", "a1", "a2"]

    bim1["rsid"].to_csv(f"{geno_path1}.snplist", sep="\t", header=None, index=None)

    bim1["merge_id"] = bim1["chr"].astype(str) + ":" + bim1["pos"].astype(str) + ":" + bim1["a2"] + ":" + bim1["a1"]
    bim2["merge_id1"] = bim2["chr"].astype(str) + ":" + bim2["pos"].astype(str) + ":" + bim2["a2"] + ":" + bim2["a1"]
    bim2["merge_id2"] = bim2["chr"].astype(str) + ":" + bim2["pos"].astype(str) + ":" + bim2["a1"] + ":" + bim2["a2"]

    common_snps1 = bim2[["rsid", "merge_id1", "a1", "a2"]].merge(bim1, how="inner", left_on=["merge_id1"], right_on=["merge_id"])
    common_snps2 = bim2[["rsid", "merge_id2", "a1", "a2"]].merge(bim1, how="inner", left_on=["merge_id2"], right_on=["merge_id"])
    common_snps = pd.concat([common_snps1, common_snps2], axis=0)

    run_command(
        f"{plink} --bfile {geno_path1} --flip {geno_path1}.snplist --make-bed --out {geno_path1}_flip",
        tool_name="plink",
    )

    bim1_flip = pd.read_csv(f"{geno_path1}_flip.bim", sep="\t", header=None)
    bim1_flip.columns = ["chr", "rsid", "kb", "pos", "a1", "a2"]

    bim1_flip["merge_id"] = bim1_flip["chr"].astype(str) + ":" + bim1_flip["pos"].astype(str) + ":" + bim1_flip["a2"] + ":" + bim1_flip["a1"]
    common_snps1 = bim2[["rsid", "merge_id1", "a1", "a2"]].merge(bim1_flip, how="inner", left_on=["merge_id1"], right_on=["merge_id"])
    common_snps2 = bim2[["rsid", "merge_id2", "a1", "a2"]].merge(bim1_flip, how="inner", left_on=["merge_id2"], right_on=["merge_id"])

    common_snps = pd.concat([common_snps, common_snps1, common_snps2], axis=0)
    common_snps = common_snps.drop_duplicates(subset=["chr", "pos"], ignore_index=True)

    common_snps_file = f"{out_name}.common_snps"
    common_snps["rsid_y"].to_csv(f"{common_snps_file}", sep="\t", header=False, index=False)

    run_command(
        f"{plink2} --bfile {geno_path1} --extract {common_snps_file} --make-bed --out {out_name}",
        tool_name="plink2",
    )

    return {"common_snps": common_snps_file, "bed": out_name}
```

> Diff this against `utils.py:404-459` — the only changes are the two `run_command` lines (was `shell_do`) and the `get_plink()/get_plink2()` resolution. Everything pandas is identical.

- [ ] **Step 4: Run the differential test**

Run: `python -m pytest tests/unit/test_ancestry/test_preprocessing.py::test_get_common_snps_matches_legacy -v`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add genotools/ancestry/preprocessing.py tests/unit/test_ancestry/test_preprocessing.py
git commit -m "refactor: port get_common_snps into ancestry.preprocessing (differential-tested)"
```

---

### Task B3: Port `get_raw_files` + `clean_up_files` → `ancestry/preprocessing.py`

**Files:**
- Modify: `genotools/ancestry/preprocessing.py`
- Test: `tests/unit/test_ancestry/test_preprocessing.py`

**Interfaces:**
- Consumes: `get_common_snps` (Task B2), `core.executors.run_command`, `get_plink2`.
- Produces:
  - `clean_up_files(files: list) -> None` — faithful port of `Ancestry.clean_up` (`legacy.py:54-70`).
  - `get_raw_files(geno_path: str, ref_panel: str, ref_labels: str, out_path: str, train: bool, common_snps_file: str | None = None) -> dict` — faithful port of `Ancestry.get_raw_files` (`legacy.py:72-264`) returning `{"raw_ref": DataFrame, "raw_geno": DataFrame, "out_paths": dict}`. Differences from legacy, all deliberate: (1) reads `self.*` as params; (2) `shell_do` → `run_command`; (3) `concat_logs(...)` calls **dropped** (structured logging replaces raw `.log` aggregation; does not affect genotype output); (4) the **containerized branch is not ported** (new path always runs non-containerized); (5) inference mode takes `common_snps_file` directly instead of deriving it from `self.model_path`.

- [ ] **Step 1: Write the differential failing tests**

Append to `tests/unit/test_ancestry/test_preprocessing.py`:

```python
def _sorted_df(df):
    df = df.sort_values(["FID", "IID"]).reset_index(drop=True)
    return df.reindex(sorted(df.columns), axis=1)


def test_get_raw_files_train_matches_legacy(synth_ref_bfile, synth_ref_labels, synth_geno_pfile, tmp_path):
    """Train-mode raw_ref/raw_geno match legacy Ancestry.get_raw_files byte-for-byte."""
    from genotools.ancestry.preprocessing import get_raw_files as new_fn
    from genotools.ancestry.legacy import Ancestry

    ref, labels, geno = str(synth_ref_bfile), str(synth_ref_labels), str(synth_geno_pfile)

    # Legacy
    legacy_dir = tmp_path / "legacy"; legacy_dir.mkdir()
    anc = Ancestry()
    anc.geno_path = geno
    anc.ref_panel = ref
    anc.ref_labels = labels
    anc.out_path = str(legacy_dir / "out")
    anc.train = True
    anc.model_path = None
    anc.containerized = False
    legacy_res = anc.get_raw_files()

    # New
    new_dir = tmp_path / "new"; new_dir.mkdir()
    new_res = new_fn(
        geno_path=geno, ref_panel=ref, ref_labels=labels,
        out_path=str(new_dir / "out"), train=True,
    )

    pd.testing.assert_frame_equal(
        _sorted_df(legacy_res["raw_ref"]), _sorted_df(new_res["raw_ref"]),
    )
    pd.testing.assert_frame_equal(
        _sorted_df(legacy_res["raw_geno"]), _sorted_df(new_res["raw_geno"]),
    )


def test_get_raw_files_inference_matches_legacy(synth_ref_bfile, synth_ref_labels, synth_geno_pfile, tmp_path):
    """Inference-mode raw_geno matches legacy, exercising the missing-column fill loop.

    The inference geno is built MISSING chr22 while the reference panel and the
    common-SNP set retain it, so the inference-only fill loop (`np.repeat(2, ...)`
    + reorder to ref columns) actually runs — and the differential comparison
    validates it against legacy byte-for-byte.
    """
    import shutil
    from unittest.mock import patch
    import genotools.ancestry.preprocessing as prep
    from genotools.ancestry.preprocessing import get_raw_files as new_fn
    from genotools.ancestry.legacy import Ancestry
    from genotools.core.executors import run_command, get_plink2

    ref, labels, geno = str(synth_ref_bfile), str(synth_ref_labels), str(synth_geno_pfile)

    # Produce a common_snps file via a train run on the FULL geno (includes chr22)
    prep_dir = tmp_path / "prep"; prep_dir.mkdir()
    train_res = new_fn(geno_path=geno, ref_panel=ref, ref_labels=labels,
                       out_path=str(prep_dir / "out"), train=True)
    common_snps_src = train_res["out_paths"]["common_snps"]

    # Inference geno MISSING chr22 -> ref's chr22 common SNPs must be filled in
    subset_geno = tmp_path / "subset_geno"
    run_command(
        f"{get_plink2()} --pfile {geno} --not-chr 22 "
        f"--make-pgen psam-cols=fid,parents,sex,pheno1,phenos --out {subset_geno}",
        tool_name="plink2",
    )

    # Legacy inference derives the common_snps path from model_path: <dir>/model.common_snps
    ld = tmp_path / "legacy"; ld.mkdir()
    shutil.copy2(common_snps_src, ld / "model.common_snps")
    anc = Ancestry()
    anc.geno_path, anc.ref_panel, anc.ref_labels = str(subset_geno), ref, labels
    anc.out_path = str(ld / "out")
    anc.train = False
    anc.model_path = str(ld / "model.pkl")
    anc.containerized = False
    legacy_res = anc.get_raw_files()

    # New inference takes the common_snps file directly; spy proves the fill loop ran
    nd = tmp_path / "new"; nd.mkdir()
    shutil.copy2(common_snps_src, nd / "model.common_snps")
    with patch.object(prep.np, "repeat", wraps=prep.np.repeat) as spy:
        new_res = new_fn(geno_path=str(subset_geno), ref_panel=ref, ref_labels=labels,
                         out_path=str(nd / "out"), train=False,
                         common_snps_file=str(nd / "model.common_snps"))
    assert spy.called, "missing-column fill loop was not exercised (degenerate test)"

    pd.testing.assert_frame_equal(
        _sorted_df(legacy_res["raw_geno"]), _sorted_df(new_res["raw_geno"]),
    )


def test_clean_up_files_removes_extensions(tmp_path):
    from genotools.ancestry.preprocessing import clean_up_files
    prefix = tmp_path / "junk"
    for ext in ("bed", "bim", "fam", "raw"):
        (tmp_path / f"junk.{ext}").write_text("x")
    clean_up_files([str(prefix)])
    for ext in ("bed", "bim", "fam", "raw"):
        assert not (tmp_path / f"junk.{ext}").exists()
```

- [ ] **Step 2: Run to verify failure**

Run: `python -m pytest tests/unit/test_ancestry/test_preprocessing.py::test_get_raw_files_train_matches_legacy -v`
Expected: FAIL with `ImportError: cannot import name 'get_raw_files'`.

- [ ] **Step 3: Implement `clean_up_files` and `get_raw_files`**

Add to `genotools/ancestry/preprocessing.py`. `clean_up_files` (from `legacy.py:54-70`):

```python
def clean_up_files(files: list) -> None:
    """Remove intermediate PLINK files by known extensions (port of Ancestry.clean_up)."""
    extensions = ["bim", "bed", "fam", "hh", "snplist", "ref_allele", "alleles", "raw"]
    for file in files:
        for ext in extensions:
            file_ext = f"{file}.{ext}"
            if os.path.exists(file_ext):
                os.remove(file_ext)
```

`get_raw_files` — port `legacy.py:72-264` applying the five transformations from the Interfaces block. Reproduce the legacy body with these exact edits:

```python
def get_raw_files(
    geno_path: str,
    ref_panel: str,
    ref_labels: str,
    out_path: str,
    train: bool,
    common_snps_file: str | None = None,
) -> dict:
    """Process reference + genotype data into labeled raw feature matrices.

    Faithful port of legacy Ancestry.get_raw_files (train + model-inference
    branches; containerized branch intentionally omitted). In inference mode,
    pass the model's common-SNP list path via `common_snps_file`.
    """
    plink2 = get_plink2()
    out_paths = {}

    # variant prune geno before getting common snps
    geno_prune_path = f"{out_path}_variant_pruned"
    run_command(
        f"{plink2} --pfile {geno_path} --geno 0.1 --max-alleles 2 --make-bed --out {geno_prune_path}",
        tool_name="plink2",
    )
    out_paths["geno_pruned_bed"] = geno_prune_path

    if train:
        ref_common_snps = f"{out_path}_umap_linearsvc_ancestry_model"
        common_snps_files = get_common_snps(ref_panel, geno_prune_path, ref_common_snps)
        out_paths = {**out_paths, **common_snps_files}
    else:
        if common_snps_file is None or not os.path.isfile(common_snps_file):
            raise FileNotFoundError(f"{common_snps_file} file does not exist.")
        ref_common_snps = f"{os.path.dirname(out_path)}/" + os.path.basename(common_snps_file).split(".")[0]
        run_command(
            f"{plink2} --bfile {ref_panel} --extract {common_snps_file} --make-bed --out {ref_common_snps}",
            tool_name="plink2",
        )
        out_paths["common_snps"] = common_snps_file
        out_paths["bed"] = ref_common_snps

    if not os.path.exists(f"{ref_common_snps}.bed"):
        raise FileNotFoundError(f"{ref_common_snps} PLINK binaries (bed/bim/fam) do not exist.")

    # get raw version of common snps - reference panel
    run_command(
        f"{plink2} --bfile {ref_common_snps} --recode A --out {ref_common_snps}",
        tool_name="plink2",
    )
    if not os.path.exists(f"{ref_common_snps}.raw"):
        raise FileNotFoundError(f"{ref_common_snps}.raw does not exist.")

    ref_raw = pd.read_csv(f"{ref_common_snps}.raw", sep=r"\s+")
    ref_ids = ref_raw[["FID", "IID"]]
    ref_snps = ref_raw.drop(columns=["FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE"], axis=1)
    ref_snps_cols = ref_snps.columns.str.extract("(.*)_")[0]
    ref_snps.columns = ref_snps_cols
    col_names = ["FID", "IID"] + list(ref_snps_cols)
    ref_raw = pd.concat([ref_ids, ref_snps], axis=1)
    ref_raw.columns = col_names

    # read ancestry file with reference labels
    ancestry = pd.read_csv(f"{ref_labels}", sep="\t", header=None, names=["FID", "IID", "label"])
    ref_fam = pd.read_csv(f"{ref_panel}.fam", sep=r"\s+", header=None)
    ref_labeled = ref_fam.merge(ancestry, how="left", left_on=[0, 1], right_on=["FID", "IID"])

    labeled_ref_raw = ref_raw.merge(ref_labeled, how="left", on=["FID", "IID"])
    labeled_ref_raw.drop(columns=[0, 1, 2, 3, 4, 5], inplace=True)

    print()
    print()
    print("Labeled Reference Ancestry Counts:")
    print(labeled_ref_raw.label.value_counts())
    print()
    print()

    # get reference alleles from ref_common_snps
    ref_common_snps_ref_alleles = f"{ref_common_snps}.ref_allele"
    ref_common_snps_bim = pd.read_csv(f"{ref_common_snps}.bim", header=None, sep="\t")
    ref_common_snps_bim.columns = ["chr", "rsid", "kb", "pos", "a1", "a2"]
    ref_common_snps_bim[["rsid", "a1"]].to_csv(ref_common_snps_ref_alleles, sep="\t", header=False, index=False)
    out_paths["ref_alleles"] = ref_common_snps_ref_alleles

    geno_common_snps = f"{out_path}_common_snps"
    get_common_snps(geno_prune_path, ref_common_snps, geno_common_snps)

    geno_common_snps_bim = pd.read_csv(f"{geno_common_snps}.bim", sep=r"\s+", header=None)
    geno_common_snps_bim.columns = ["chr", "rsid", "kb", "pos", "a1", "a2"]

    ref_common_snps_bim["merge_id"] = ref_common_snps_bim["chr"].astype(str) + ":" + ref_common_snps_bim["pos"].astype(str)
    geno_common_snps_bim["merge_id"] = geno_common_snps_bim["chr"].astype(str) + ":" + geno_common_snps_bim["pos"].astype(str)

    merge_common_snps_bim = geno_common_snps_bim[["merge_id", "a1", "a2"]].merge(ref_common_snps_bim, how="inner", on=["merge_id"])
    merge_common_snps_bim[["chr", "rsid", "kb", "pos", "a1_x", "a2_x"]].to_csv(f"{geno_common_snps}.bim", sep="\t", header=None, index=None)

    switch = {"A": "T", "T": "A", "C": "G", "G": "C"}
    merge_common_snps_bim["a1_x_switch"] = merge_common_snps_bim["a1_x"].map(switch)
    merge_common_snps_switch = merge_common_snps_bim[
        (merge_common_snps_bim["a1_y"] != merge_common_snps_bim["a1_x"])
        & (merge_common_snps_bim["a1_y"] != merge_common_snps_bim["a1_x_switch"])
    ]
    merge_common_snps_switch[["rsid", "a2_x"]].to_csv(f"{geno_common_snps}_switch.alleles", sep="\t", header=False, index=False)

    run_command(
        f"{plink2} --bfile {geno_common_snps} --alt1-allele {geno_common_snps}_switch.alleles --recode A --out {geno_common_snps}",
        tool_name="plink2",
    )

    raw_geno = pd.read_csv(f"{geno_common_snps}.raw", sep=r"\s+")
    geno_ids = raw_geno[["FID", "IID"]]
    geno_snps = raw_geno.drop(columns=["FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE"], axis=1)
    geno_snps.columns = geno_snps.columns.str.extract("(.*)_")[0]

    # adding missing snps when not training
    missing_cols = []
    if not train:
        for col in ref_snps.columns:
            if col not in geno_snps.columns:
                missing_cols += [pd.Series(np.repeat(2, geno_snps.shape[0]), name=col)]
        if len(missing_cols) > 0:
            missing_cols = pd.concat(missing_cols, axis=1)
            geno_snps = pd.concat([geno_snps, missing_cols], axis=1)
        geno_snps = geno_snps[ref_snps.columns]

    raw_geno = pd.concat([geno_ids, geno_snps], axis=1)
    raw_geno.columns = col_names
    raw_geno["label"] = "new"

    # remove intermediate files (concat_logs dropped: structured logging replaces raw .log aggregation)
    files = [geno_prune_path, f"{geno_prune_path}_flip", f"{out_path}_common_snps_switch"]
    clean_up_files(files)

    return {
        "raw_ref": labeled_ref_raw,
        "raw_geno": raw_geno,
        "out_paths": out_paths,
    }
```

> This body is `legacy.py:148-264` with exactly these edits: `plink2_exec` → local `plink2`; `shell_do(x)` → `run_command(x, tool_name="plink2")`; `self.ref_labels/ref_panel/out_path/train` → params; `self.clean_up` → `clean_up_files`; the `listOfFiles`/`concat_logs` block removed. Diff against the original to confirm no pandas line drifted — the Step-1 differential test is the gate.

- [ ] **Step 4: Run the differential + clean_up tests**

Run: `python -m pytest tests/unit/test_ancestry/test_preprocessing.py -v`
Expected: PASS (all three tests). If `assert_frame_equal` fails, diff your `# ----` region against `legacy.py:148-262` line-by-line.

- [ ] **Step 5: Commit**

```bash
git add genotools/ancestry/preprocessing.py tests/unit/test_ancestry/test_preprocessing.py
git commit -m "refactor: port get_raw_files + clean_up into ancestry.preprocessing (differential-tested)"
```

---

### Task B4: Port `split_cohort_ancestry` → `ancestry/cohort.py`

**Files:**
- Create: `genotools/ancestry/cohort.py`
- Test: `tests/unit/test_ancestry/test_cohort.py`

**Interfaces:**
- Consumes: `core.executors.run_command`, `get_plink2`.
- Produces: `split_cohort_by_ancestry(labels_path: str, geno_path: str, out_path: str, min_samples: int, subset: list | None = None) -> dict` returning `{"labels": [...], "paths": [...], "pruned_samples": DataFrame}`. Faithful port of `Ancestry.split_cohort_ancestry` (`legacy.py:907-963`).

- [ ] **Step 1: Write the differential failing test**

Create `tests/unit/test_ancestry/test_cohort.py`:

```python
from pathlib import Path

import pandas as pd


def test_split_cohort_matches_legacy(synth_geno_pfile, tmp_path):
    from genotools.ancestry.cohort import split_cohort_by_ancestry
    from genotools.ancestry.legacy import Ancestry

    geno = str(synth_geno_pfile)
    # Build a predicted-labels file over the real sample IDs
    psam = pd.read_csv(f"{geno}.psam", sep=r"\s+")
    fid = psam.get("#FID", psam.get("FID"))
    iid = psam["IID"]
    labels = pd.DataFrame({
        "FID": fid, "IID": iid,
        "label": ["EUR" if i % 2 == 0 else "AFR" for i in range(len(iid))],
    })
    labels_path = tmp_path / "pred_labels.txt"
    labels.to_csv(labels_path, sep="\t", index=False)

    # Legacy
    ld = tmp_path / "legacy"; ld.mkdir()
    anc = Ancestry()
    anc.geno_path = geno
    anc.out_path = str(ld / "out")
    anc.subset = None
    anc.min_samples = 10
    legacy_res = anc.split_cohort_ancestry(labels_path=str(labels_path))

    # New
    nd = tmp_path / "new"; nd.mkdir()
    new_res = split_cohort_by_ancestry(
        labels_path=str(labels_path), geno_path=geno,
        out_path=str(nd / "out"), min_samples=10, subset=None,
    )

    assert sorted(new_res["labels"]) == sorted(legacy_res["labels"])
    assert (new_res["pruned_samples"].shape == legacy_res["pruned_samples"].shape)
    # Each retained ancestry produced a pgen
    for label in new_res["labels"]:
        assert Path(f"{nd / 'out'}_{label}.pgen").exists()
```

- [ ] **Step 2: Run to verify failure**

Run: `python -m pytest tests/unit/test_ancestry/test_cohort.py -v`
Expected: FAIL with `ModuleNotFoundError: No module named 'genotools.ancestry.cohort'`.

- [ ] **Step 3: Implement the port**

Create `genotools/ancestry/cohort.py` (port of `legacy.py:907-963`; `self.*` → params, `shell_do` → `run_command`, drop `concat_logs`):

```python
"""Cohort splitting by predicted ancestry (legacy-free port)."""

import pandas as pd

from genotools.core.executors import run_command, get_plink2

PLINK2_PSAM_COLS = "psam-cols=fid,parents,sex,pheno1,phenos"


def split_cohort_by_ancestry(
    labels_path: str,
    geno_path: str,
    out_path: str,
    min_samples: int,
    subset: list | None = None,
) -> dict:
    """Split a cohort into per-ancestry pfiles based on predicted labels."""
    plink2 = get_plink2()

    pred_labels = pd.read_csv(labels_path, sep="\t")
    labels_list = []
    outfiles = []

    split_labels = subset if subset else pred_labels.label.unique()
    pruned_samples = pd.DataFrame(columns=["FID", "IID", "step", "label"])

    for label in split_labels:
        if pred_labels[pred_labels.label == label].shape[0] >= min_samples:
            labels_list.append(label)
            outname = f"{out_path}_{label}"
            outfiles.append(outname)
            ancestry_group_outpath = f"{outname}.samples"
            pred_labels[pred_labels.label == label][["FID", "IID"]].to_csv(
                ancestry_group_outpath, index=False, header=False, sep="\t"
            )
            run_command(
                f"{plink2} --pfile {geno_path} --keep {ancestry_group_outpath} "
                f"--make-pgen {PLINK2_PSAM_COLS} --out {outname}",
                tool_name="plink2",
            )
        else:
            pruned_samples_label = pred_labels[pred_labels.label == label].copy()
            pruned_samples_label["step"] = "insufficient_ancestry_sample_n"
            pruned_samples_label = pruned_samples_label[["FID", "IID", "step", "label"]]
            pruned_samples = pd.concat([pruned_samples, pruned_samples_label], axis=0, ignore_index=True)

    return {"labels": labels_list, "paths": outfiles, "pruned_samples": pruned_samples}
```

> `PLINK2_PSAM_COLS` here matches `core.executors.PLINK2_PSAM_COLS`; import from core instead if preferred (`from genotools.core.executors import PLINK2_PSAM_COLS`).

- [ ] **Step 4: Run the differential test**

Run: `python -m pytest tests/unit/test_ancestry/test_cohort.py -v`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add genotools/ancestry/cohort.py tests/unit/test_ancestry/test_cohort.py
git commit -m "refactor: port split_cohort_ancestry into ancestry.cohort (differential-tested)"
```

---

### Task B5: Rewire the new ancestry path to use the ported functions

**Files:**
- Modify: `genotools/cli/runner.py` (`_run_ancestry_prediction_new`, `_run_training_mode`, `_run_inference_mode`)
- Test: `tests/unit/test_cli/test_runner.py`

**Interfaces:**
- Consumes: `ancestry.preprocessing.get_raw_files`, `ancestry.preprocessing.clean_up_files`, `ancestry.cohort.split_cohort_by_ancestry` (Tasks B2–B4).
- Produces: a `genotools-new` path that never calls `self._ancestry.*`.

- [ ] **Step 1: Write the failing test (spy that legacy is not used)**

Add to `tests/unit/test_cli/test_runner.py`:

```python
class TestNewAncestryStandalone:
    def test_new_path_uses_ported_functions_not_legacy(self, monkeypatch):
        """The new ancestry training path calls ancestry.preprocessing.get_raw_files,
        not the legacy Ancestry.get_raw_files."""
        import genotools.ancestry.preprocessing as prep
        calls = {"new_get_raw_files": 0}

        real = prep.get_raw_files
        def spy(*a, **k):
            calls["new_get_raw_files"] += 1
            return real(*a, **k)
        monkeypatch.setattr(prep, "get_raw_files", spy)

        # Guard: legacy get_raw_files must not be called from the new path
        from genotools.ancestry.legacy import Ancestry
        def boom(self):
            raise AssertionError("new path called legacy Ancestry.get_raw_files")
        monkeypatch.setattr(Ancestry, "get_raw_files", boom)

        # Static guard is the primary check here (full ancestry run needs a ref panel):
        import inspect
        from genotools.cli import runner
        src = inspect.getsource(runner.PipelineRunner._run_training_mode)
        src += inspect.getsource(runner.PipelineRunner._run_inference_mode)
        src += inspect.getsource(runner.PipelineRunner._run_ancestry_prediction_new)
        assert "self._ancestry.get_raw_files" not in src
        assert "self._ancestry.split_cohort_ancestry" not in src
        assert "self._ancestry.clean_up" not in src
```

- [ ] **Step 2: Run to verify failure**

Run: `python -m pytest tests/unit/test_cli/test_runner.py::TestNewAncestryStandalone -v`
Expected: FAIL (the `assert "self._ancestry.get_raw_files" not in src` fails — the new path still calls legacy).

- [ ] **Step 3: Rewire `_run_training_mode`**

In `genotools/cli/runner.py`, `_run_training_mode`, replace the legacy preprocessing lines:

```python
        self._ancestry.model_path = None
        self._ancestry.train = True

        # PLINK preprocessing
        raw = self._ancestry.get_raw_files()
```

with:

```python
        from ..ancestry.preprocessing import get_raw_files

        # PLINK preprocessing (ported; legacy-free)
        raw = get_raw_files(
            geno_path=str(self.state.geno_path),
            ref_panel=str(self.args.ancestry.ref_panel),
            ref_labels=str(self.args.ancestry.ref_labels),
            out_path=actual_out,
            train=True,
        )
```

- [ ] **Step 4: Rewire `_run_inference_mode`**

In `_run_inference_mode`, replace the common-SNP file round-trip + legacy call (the block writing `*.common_snps` / setting `self._ancestry.model_path` and later `raw = self._ancestry.get_raw_files()`) with a direct list write + ported call. Replace from `# Write common SNPs to temp file...` through `raw = self._ancestry.get_raw_files()` with:

```python
        # Write the model's common SNPs to a file for the ported preprocessing
        if model.common_snps is not None:
            common_snps_file = f"{actual_out}_loaded_model.common_snps"
            with open(common_snps_file, "w") as f:
                for snp in model.common_snps:
                    f.write(f"{snp}\n")
        else:
            model_dir = Path(model_path)
            snps_file = model_dir / "common_snps.txt" if model_dir.is_dir() else None
            if snps_file and snps_file.exists():
                common_snps_file = f"{actual_out}_loaded_model.common_snps"
                shutil.copy2(snps_file, common_snps_file)
            else:
                raise FileNotFoundError(f"No common_snps found in model: {model_path}")

        from ..ancestry.preprocessing import get_raw_files

        # PLINK preprocessing (ported; legacy-free)
        raw = get_raw_files(
            geno_path=str(self.state.geno_path),
            ref_panel=str(self.args.ancestry.ref_panel),
            ref_labels=str(self.args.ancestry.ref_labels),
            out_path=actual_out,
            train=False,
            common_snps_file=common_snps_file,
        )
```

- [ ] **Step 5: Rewire `_run_ancestry_prediction_new` (split + clean up)**

In `_run_ancestry_prediction_new`, replace:

```python
        # Cohort splitting (reuse legacy)
        labels_path = f"{actual_out}_umap_linearsvc_predicted_labels.txt"
        ancestry_split = self._ancestry.split_cohort_ancestry(
            labels_path=labels_path
        )
```

with:

```python
        # Cohort splitting (ported; legacy-free)
        from ..ancestry.cohort import split_cohort_by_ancestry

        labels_path = f"{actual_out}_umap_linearsvc_predicted_labels.txt"
        ancestry_split = split_cohort_by_ancestry(
            labels_path=labels_path,
            geno_path=str(self.state.geno_path),
            out_path=actual_out,
            min_samples=self.args.ancestry.min_samples,
            subset=self.args.ancestry.subset_ancestry,
        )
```

And replace the cleanup call:

```python
        self._ancestry.clean_up(files_to_clean)
```

with:

```python
        from ..ancestry.preprocessing import clean_up_files

        clean_up_files(files_to_clean)
```

> Also delete the now-dead legacy field assignments in `_run_ancestry_prediction_new` that only fed the legacy helpers (`self._ancestry.geno_path`, `.out_path`, `.final_out_path`, `.ref_panel`, `.ref_labels`, `.containerized`, `.singularity`, `.subset`, `.min_samples`, and in the mode helpers `.model_path`, `.train`). Keep any field still read by retained code; if none remain, they can all go. Verify with the Step 7 grep.

- [ ] **Step 6: Run the guard test + suite**

Run:
```bash
python -m pytest tests/unit/test_cli/test_runner.py::TestNewAncestryStandalone -v
python -m pytest tests/unit tests/regression -q
```
Expected: guard test PASS; full suite green.

- [ ] **Step 7: Verify the new path is legacy-free (6b gate)**

Run:
```bash
python - <<'PY'
import inspect
from genotools.cli import runner
src = "".join(inspect.getsource(getattr(runner.PipelineRunner, m)) for m in
    ("_run_ancestry_prediction_new", "_run_training_mode", "_run_inference_mode"))
for bad in ("self._ancestry.get_raw_files", "self._ancestry.split_cohort_ancestry", "self._ancestry.clean_up"):
    assert bad not in src, bad
print("new ancestry path is legacy-free")
PY
```
Expected: `new ancestry path is legacy-free`.

- [ ] **Step 8: Commit**

```bash
git add genotools/cli/runner.py tests/unit/test_cli/test_runner.py
git commit -m "refactor: rewire genotools-new ancestry path to ported functions (legacy-free)"
```

---

### Task B6: Final verification + tracker update

**Files:**
- Modify: `REFACTOR_HARDENING.md`

- [ ] **Step 1: Full static coupling gate**

Run:
```bash
grep -rn "from ..utils import\|from genotools.utils import\|importlib" genotools/cli genotools/core genotools/qc genotools/gwas 2>/dev/null | grep "\.py:"
grep -rn "from genotools.utils\|from ..utils" genotools/ancestry/preprocessing.py genotools/ancestry/cohort.py 2>/dev/null
```
Expected: **no output** for both (new code imports zero legacy).

- [ ] **Step 2: Import + CLI smoke**

Run:
```bash
python -c "import genotools, genotools.core, genotools.qc, genotools.ancestry, genotools.gwas, genotools.cli"
python -m genotools --help >/dev/null && echo OK
```
Expected: no errors, `OK`.

- [ ] **Step 3: Full suite + parity**

Run:
```bash
python -m pytest tests/unit tests/regression -q
python -m pytest tests/regression/test_parity.py -v   # PASS or SKIP if no .venv-stable
```
Expected: green.

- [ ] **Step 4: Update the tracker**

In `REFACTOR_HARDENING.md`, mark item #6 done (Round 5): summarize 6a (utils → core, ancestry importlib hack removed) and 6b (get_raw_files/get_common_snps/split_cohort_ancestry/clean_up ported; genotools-new legacy-free; differential-tested). Note deferred: default-flip + legacy deletion (gated on item #1), and the upfront_check skip-gap.

- [ ] **Step 5: Commit**

```bash
git add REFACTOR_HARDENING.md
git commit -m "docs: record round-5 (decouple from legacy, item #6) in tracker"
```

- [ ] **Step 6: Open PR into `refactor/main`**

```bash
git push -u origin refactor/decouple-from-legacy
gh pr create --base refactor/main --title "refactor: decouple new CLI from legacy (item #6)" --body "See docs/superpowers/specs/2026-07-17-decouple-from-legacy-design.md"
```
