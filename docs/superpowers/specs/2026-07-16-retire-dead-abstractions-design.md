# Design: Retire dead abstractions (refactor hardening, round 3+)

**Date:** 2026-07-16
**Branch target:** `refactor/main` (feature branch → PR into `refactor/main`; **never `main`**).
**Tracker item:** REFACTOR_HARDENING.md remaining-work #7 ("Retire dead abstractions").

## Goal

Remove two confirmed-orphaned Python abstractions from the mid-cutover refactor so
the codebase carries less unused surface before the higher-risk legacy-decouple
work (item #6). This is a **deletion-only cleanup**: none of the removed code is
reachable at runtime, so it cannot affect the in-flight real-data parity run
(the gate for `refactor/main` → `main`, owned by another dev).

## Why these two are safe to remove (evidence)

1. **`QCPipeline`** (`genotools/qc/pipeline.py`) — an orchestration API
   (`.run()`, `.sample_qc()`/`.variant_qc()` factories). The CLI runner
   (`cli/runner.py`) reimplements orchestration in
   `_run_single_step`/`_run_qc_pipeline` and **never uses `QCPipeline`**. Its
   golden regression coverage is **redundant**: `tests/regression/test_qc_steps.py`
   already exercises every pure step function directly against golden files
   (`TestCallrateRegression`, `TestSexRegression`, `TestHetRegression`,
   `TestRelatedRegression`, `TestDirectOutlierValidation`, the variant-side
   equivalents, and the CLI-driven `TestAllSample/AllVariant/FullPipeline`). The
   `QCPipeline` test classes only exercise the unused orchestration path.
   - `QCResult` (`qc/results.py`) is produced **only** by `QCPipeline.run`; it
     dies with it. `FilterResult` (returned by every step) is unrelated and stays.
   - `QCStepProtocol` (`qc/pipeline.py`) is used only within that module.

2. **`ReferencePanel`** (`genotools/ancestry/reference.py`) — a public export,
   referenced only in a docstring, **never constructed in shipped code**, **zero
   test references**. The whole module is dead: its two sibling helpers
   `get_default_model_path()` and `validate_model_files()` are likewise only
   exported and never used. The new-ancestry path uses `AncestryModel` directly.
   The runner registers `_new_modules["ReferencePanel"]` but never reads it.

**Explicitly NOT touched:**
- `genotools/container/` — although the harness (`run.py`, `Dockerfile`, bundled
  `*.pkl` models, `requirements.txt`) is orphaned and not imported anywhere, it is
  **intentionally retained** as a historical reference for rebuilding containerized
  ancestry-inference later (maintainer decision). Its `setup.py` `package_data`
  entry stays as-is.
- Verified live + tested: `ancestry/results.py` (all classes imported by
  `ancestry/model.py` + tested), `qc/results.py::FilterResult`, `gwas/results.py`.

## Plan — one PR, two commits

### Commit 1 — `refactor: remove unused QCPipeline orchestration API`
- Delete `genotools/qc/pipeline.py` (`QCPipeline`, `QCStepProtocol`).
- Remove `QCResult` from `genotools/qc/results.py` (keep `FilterResult`).
- `genotools/qc/__init__.py`: drop `QCPipeline` and `QCResult` imports + `__all__` entries.
- Tests:
  - Delete `tests/unit/test_qc/test_pipeline.py`.
  - Remove `TestQCResult` from `tests/unit/test_qc/test_results.py` (keep `FilterResult` tests).
  - Remove `TestQCPipelineSampleQC` / `TestQCPipelineVariantQC` (lines ~465–633) and
    the `QCPipeline` imports from `tests/regression/test_qc_steps.py`. Keep all
    direct step-function golden tests and the CLI-driven pipeline tests.

### Commit 2 — `refactor: remove unadopted ReferencePanel module`
- Delete `genotools/ancestry/reference.py` in full (`ReferencePanel` +
  `get_default_model_path` + `validate_model_files`).
- `genotools/ancestry/__init__.py`: remove the `from genotools.ancestry.reference import (...)`
  block + the three `__all__` entries.
- `genotools/cli/runner.py`: drop `ReferencePanel` from the `from ..ancestry import ...`
  line and the dead `"ReferencePanel": ReferencePanel` registration in `_new_modules`.

## Verification (after each commit + at the end)

1. **Static:** `grep` for every removed name — `QCPipeline`, `QCResult`,
   `QCStepProtocol`, `ReferencePanel`, `get_default_model_path`,
   `validate_model_files` — expect zero live references.
2. **Import smoke:** `python -c "import genotools, genotools.qc, genotools.ancestry, genotools.cli, genotools.gwas"`.
3. **Full suite:** `tests/unit tests/regression` green. Test count drops by the
   deleted tests (~35–40 fewer) — deletions, not failures.
4. **Parity:** `tests/regression/test_parity.py` green with `.venv-stable` present
   (proves runtime paths untouched).
5. **CLI smoke:** `python -m genotools --help`.

## Risks / mitigations

- *A "dead" symbol is actually referenced somewhere not grepped* → import smoke
  test + full suite + step-1 grep catch it before each commit lands.
- *Interference with the parity run* → none: no QC/GWAS/ancestry execution path is
  touched; changes are pure deletions of unreachable code.

## Out of scope

- Item #6 (decouple CLI from legacy `utils.py` / ancestry importlib hack) — the
  higher-value structural work; behavior-preserving but touches the startup hot
  path, so coordinate with the parity run.
- Items #9/#10 (performance) — change execution; defer until real-data parity is proven.
