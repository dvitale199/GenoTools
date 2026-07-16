# Retire Dead Abstractions Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Delete two confirmed-orphaned Python abstractions — the unused `QCPipeline` orchestration API (with `QCResult`/`QCStepProtocol`) and the unadopted `ReferencePanel` module — from the mid-cutover refactor.

**Architecture:** Deletion-only cleanup. None of the removed code is reachable at runtime, so correctness is proven by (a) zero remaining references, (b) imports still working, (c) the full unit+regression+parity suite staying green (minus the deleted self-tests). No production behavior changes.

**Tech Stack:** Python 3.11, pytest, PLINK/PLINK2 (auto-downloaded), `.venv` (new code), `.venv-stable` (pre-refactor baseline for parity).

## Global Constraints

- Branch off `refactor/main`; PR targets **`refactor/main`, never `main`**.
- Follow CLAUDE.md: no `Co-Authored-By` line in commits; keep each concern in its own commit.
- Do NOT touch `genotools/container/` — intentionally retained as a historical reference.
- Do NOT touch `ancestry/results.py`, `qc/results.py::FilterResult`, or `gwas/results.py` — all live and tested.
- Run tests with `.venv/bin/python -m pytest`.

---

### Task 1: Remove the unused QCPipeline orchestration API

**Files:**
- Delete: `genotools/qc/pipeline.py`
- Modify: `genotools/qc/results.py` (remove `QCResult` class; keep `FilterResult`)
- Modify: `genotools/qc/__init__.py` (drop `QCPipeline`/`QCResult` imports, `__all__` entries, docstring example)
- Delete: `tests/unit/test_qc/test_pipeline.py`
- Modify: `tests/unit/test_qc/test_results.py` (remove `TestQCResult`, keep `FilterResult` tests)
- Modify: `tests/regression/test_qc_steps.py` (remove `TestQCPipelineSampleQC`/`TestQCPipelineVariantQC`)

**Interfaces:**
- Consumes: nothing from other tasks.
- Produces: `genotools.qc` no longer exports `QCPipeline` or `QCResult`; still exports `FilterResult` and all `filter_*`/`prune_ld`/`verify_kinship` functions and configs.

- [ ] **Step 1: Confirm nothing outside the deletion set references these symbols**

Run:
```bash
grep -rn "QCPipeline\|QCStepProtocol" genotools/ | grep -vE "qc/pipeline.py|qc/__init__.py"
grep -rn "QCResult" genotools/ | grep -vE "qc/pipeline.py|qc/results.py|qc/__init__.py"
grep -rn "from genotools.qc.pipeline\|from genotools.qc import .*QCPipeline\|qc\.QCPipeline" genotools/ tests/ | grep -v "test_qc_steps.py\|test_pipeline.py"
```
Expected: no output (only the files being edited/deleted reference them; `test_qc_steps.py`/`test_pipeline.py` reference QCPipeline via inline imports removed in this task).

- [ ] **Step 2: Delete the pipeline module and its unit test**

Run:
```bash
git rm genotools/qc/pipeline.py tests/unit/test_qc/test_pipeline.py
```
Expected: both files staged for deletion.

- [ ] **Step 3: Remove the `QCResult` class from `genotools/qc/results.py`**

Delete the entire `QCResult` dataclass (the block beginning `@dataclass(frozen=True)` / `class QCResult:` through the end of `get_step_result`, i.e. from the blank lines after `FilterResult` to EOF). The file must end right after `FilterResult`'s `to_dict` method. Leave the imports (`from typing import Any, Optional` is still used by `FilterResult`) and `FilterResult` untouched.

- [ ] **Step 4: Update `genotools/qc/__init__.py`**

Change the import lines:
```python
from genotools.qc.pipeline import QCPipeline
from genotools.qc.results import FilterResult, QCResult
```
to:
```python
from genotools.qc.results import FilterResult
```

Remove these `__all__` entries (and the now-empty `# Pipeline` comment):
```python
    # Pipeline
    "QCPipeline",
    # Results
    "QCResult",
```
so the Results section reads:
```python
    # Results
    "FilterResult",
```

Replace the module docstring's pipeline example (the `QCPipeline(...)` / `pipeline.run(...)` block and the `QCPipeline,` line in the import example) with just the single-step usage:
```python
"""QC module for genotype quality control.

This module provides pure functions for sample and variant quality control.

Example usage:
    from genotools.qc import filter_callrate, CallrateConfig, FilterResult

    data = GenotypeData.from_path(Path("input"))
    result = filter_callrate(data, CallrateConfig(mind=0.05), Path("output"))
"""
```

- [ ] **Step 5: Remove `TestQCResult` from `tests/unit/test_qc/test_results.py`**

Change the import:
```python
from genotools.qc.results import FilterResult, QCResult
```
to:
```python
from genotools.qc.results import FilterResult
```
Delete the entire `class TestQCResult:` block (from `class TestQCResult:` up to, but not including, `class TestLegacyCompatibility:`). Keep `TestFilterResult` and `TestLegacyCompatibility`.

- [ ] **Step 6: Remove QCPipeline regression classes from `tests/regression/test_qc_steps.py`**

Delete the section from the `# QCPipeline API Regression Tests` header comment through the end of `TestQCPipelineVariantQC` (i.e. `class TestQCPipelineSampleQC:` and `class TestQCPipelineVariantQC:` and their methods), up to but not including the next section header (`TestAllSamplePipeline`). The inline `from genotools.qc import QCPipeline, ...` statements live inside those methods and are removed with them. Keep every direct step-function golden test (`TestCallrateRegression` … `TestVariantOutlierValidation`) and the CLI-driven `TestAllSamplePipeline`/`TestAllVariantPipeline`/`TestFullPipeline`.

- [ ] **Step 7: Verify no references remain and imports work**

Run:
```bash
grep -rn "QCPipeline\|QCStepProtocol" genotools/ tests/
grep -rn "QCResult" genotools/ tests/
.venv/bin/python -c "import genotools, genotools.qc, genotools.ancestry, genotools.cli, genotools.gwas; from genotools.qc import FilterResult; print('imports OK')"
```
Expected: first two greps produce no output; import prints `imports OK`.

- [ ] **Step 8: Run the QC + regression suites**

Run:
```bash
.venv/bin/python -m pytest tests/unit/test_qc tests/regression/test_qc_steps.py -q
```
Expected: all pass (fewer tests than before — the deleted `QCPipeline`/`QCResult` tests are gone, not failing).

- [ ] **Step 9: Commit**

```bash
git add -A
git commit -m "refactor: remove unused QCPipeline orchestration API

QCPipeline (qc/pipeline.py) plus QCResult and QCStepProtocol were an
orchestration API the CLI runner never used (it reimplements orchestration
in _run_single_step/_run_qc_pipeline). Their golden coverage was redundant
with the direct step-function golden tests, which remain. Removes the class,
its QCResult return type, their exports, and their self-tests. FilterResult
(returned by every step) is kept."
```

---

### Task 2: Remove the unadopted ReferencePanel module

**Files:**
- Delete: `genotools/ancestry/reference.py`
- Modify: `genotools/ancestry/__init__.py` (drop the reference import block + `__all__` entries + docstring example)
- Modify: `genotools/cli/runner.py` (drop `ReferencePanel` from the ancestry import + the dead `_new_modules` registration)

**Interfaces:**
- Consumes: nothing from Task 1.
- Produces: `genotools.ancestry` no longer exports `ReferencePanel`, `get_default_model_path`, or `validate_model_files`; the runner's `_new_modules` no longer has a `"ReferencePanel"` key (it was never read).

- [ ] **Step 1: Confirm the module's symbols are unused outside itself**

Run:
```bash
grep -rn "ReferencePanel\|get_default_model_path\|validate_model_files" genotools/ tests/ | grep -vE "ancestry/reference.py|ancestry/__init__.py|cli/runner.py"
```
Expected: at most docstring hits in `ancestry/model.py` (lines with `>>>`). No real import or call.

- [ ] **Step 2: Delete the reference module**

Run:
```bash
git rm genotools/ancestry/reference.py
```
Expected: file staged for deletion.

- [ ] **Step 3: Update `genotools/ancestry/__init__.py`**

Remove the entire reference import block:
```python
# Reference panel management
from genotools.ancestry.reference import (
    ReferencePanel,
    get_default_model_path,
    validate_model_files,
)
```
Remove these `__all__` entries:
```python
    "ReferencePanel",
```
```python
    "get_default_model_path",
    "validate_model_files",
```
If the module docstring contains a `ReferencePanel` example (`>>> from genotools.ancestry import ... ReferencePanel ...` and `>>> ref = ReferencePanel.load(...)`), remove those two lines from the docstring so it no longer references the deleted class.

- [ ] **Step 4: Update `genotools/cli/runner.py`**

Change:
```python
        from ..ancestry import AncestryModel, ReferencePanel, AncestryConfig
```
to:
```python
        from ..ancestry import AncestryModel, AncestryConfig
```
Remove the dead registration line from the `self._new_modules = {...}` dict:
```python
            "ReferencePanel": ReferencePanel,
```

- [ ] **Step 5: Verify no references remain and imports work**

Run:
```bash
grep -rn "ReferencePanel\|get_default_model_path\|validate_model_files" genotools/ tests/
.venv/bin/python -c "import genotools, genotools.ancestry, genotools.cli; from genotools.ancestry import AncestryModel, AncestryConfig; print('imports OK')"
```
Expected: grep shows at most `>>>` docstring lines in `ancestry/model.py` (harmless prose) and nothing else; import prints `imports OK`.

- [ ] **Step 6: Run the ancestry + CLI suites**

Run:
```bash
.venv/bin/python -m pytest tests/unit/test_ancestry tests/unit/test_cli -q
```
Expected: all pass.

- [ ] **Step 7: Commit**

```bash
git add -A
git commit -m "refactor: remove unadopted ReferencePanel module

ancestry/reference.py (ReferencePanel + get_default_model_path +
validate_model_files) was a public export referenced only in docstrings,
never constructed in shipped code, with zero test coverage. The new-ancestry
path uses AncestryModel directly, and the runner registered ReferencePanel in
_new_modules but never read it. Removes the module, its exports, and the dead
registration."
```

---

### Task 3: Full verification and PR

**Files:** none (verification + PR only).

**Interfaces:**
- Consumes: the two commits from Tasks 1 and 2.
- Produces: a PR into `refactor/main`.

- [ ] **Step 1: Full unit + regression suite**

Run:
```bash
.venv/bin/python -m pytest tests/unit tests/regression -q
```
Expected: all pass (parity tests included when `.venv-stable` is present; they skip cleanly otherwise). Total count is lower than the 397 baseline by exactly the deleted tests — all deletions, no failures.

- [ ] **Step 2: CLI smoke test**

Run:
```bash
.venv/bin/python -m genotools --help
```
Expected: usage text prints, exit code 0.

- [ ] **Step 3: Push and open PR into refactor/main**

```bash
git push -u origin refactor/retire-dead-abstractions
gh pr create --base refactor/main --head refactor/retire-dead-abstractions \
  --title "refactor: retire unused QCPipeline and ReferencePanel abstractions" \
  --body "<markdown summary + test plan per CLAUDE.md; note: targets refactor/main, container/ intentionally retained>"
```
Expected: PR URL printed; `gh pr view --json baseRefName` shows `refactor/main`.

- [ ] **Step 4: Confirm CI green**

Run:
```bash
gh pr checks <PR#>
```
Expected: both "Unit + regression tests" and "Old-vs-new parity" pass.

---

## Self-Review

**Spec coverage:** Spec's two removal targets → Task 1 (QCPipeline/QCResult/QCStepProtocol) and Task 2 (ReferencePanel module). Spec's "explicitly NOT touched" (container/, ancestry/results.py, FilterResult, gwas/results.py) → encoded in Global Constraints and never referenced by any task. Spec's verification (grep, import smoke, full suite, parity, CLI smoke) → Task 1 steps 7–8, Task 2 steps 5–6, Task 3 steps 1–2. ✓

**Placeholder scan:** No TBD/TODO. The only free-text placeholder is the PR body in Task 3 Step 3 (`<markdown summary...>`), which is intentionally authored at PR time. ✓

**Type consistency:** `FilterResult` kept everywhere; `QCResult`/`QCPipeline`/`QCStepProtocol`/`ReferencePanel`/`get_default_model_path`/`validate_model_files` removed consistently across source, exports, and tests. ✓
