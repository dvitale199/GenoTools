# Design: Flip default to new ancestry + cleave legacy (refactor hardening, round 6)

**Date:** 2026-07-20
**Branch target:** `refactor/main` (feature branch → PR/merge into `refactor/main`; **never `main`**).
**Tracker item:** REFACTOR_HARDENING.md round-6 (follows item #6, "decouple from legacy", now merged).

## Goal

Make `refactor/main` a clean, **legacy-free** candidate for real-data parity by:

1. **Flip the default** so the `genotools` command runs the new `AncestryModel`
   path (retire the `genotools-new` A/B entry point and the `use_new_ancestry`
   flag — there is only one path now).
2. **Delete** the two dead legacy modules: `genotools/ancestry/legacy.py` and the
   orphaned `genotools/gwas.py`.
3. **Golden-convert** the four differential tests that lose their legacy oracle
   when `legacy.py` is deleted.
4. **Fix the `upfront_check` skip-gap** so a real-cohort run vs `main` does not
   diverge on data-driven step skipping (currently dropped in the new pipeline).

## Validation model (recap)

Parity is `refactor/main` (new) vs `main` (git-preserved legacy baseline). The
parity harness (`tests/regression/test_parity.py`) shells out to a `.venv-stable`
built from `origin/main` and compares against `python -m genotools` on the branch.
It never imports repo legacy modules, and its 6 tests exercise **QC and GWAS only
(no ancestry)** — so flipping the ancestry default and deleting `legacy.py` cannot
break the parity gate.

## Scope

### In scope

- Flip default (runner + `cli/__init__` + `setup.py`).
- Delete `ancestry/legacy.py`, `gwas.py`, and the now-obsolete A/B script
  `tests/scripts/test_ancestry_ab.py` (it invokes `genotools-new` and imports
  nothing, but its whole purpose — legacy-vs-new A/B — is gone).
- Golden-convert 4 tests; delete/trim 1 test class.
- Port the `upfront_check` data-driven skip logic into `core/validation.py`.

### Deferred (kept in place this round)

- `genotools/utils.py` and `genotools/imputation.py` — **`imputation.py` imports
  `utils.shell_do`, so `utils.py` must stay while `imputation.py` stays.** The two
  differential tests that import `genotools.utils` (`banner`, `get_common_snps`)
  keep working because `utils.py` remains.
- `genotools/container/`.

## Detailed decisions

### 1. Flip-the-default wiring

**`genotools/cli/runner.py`:**
- `PipelineRunner.__init__`: drop the `use_new_ancestry` parameter and the
  `self._use_new_ancestry` attribute.
- `_initialize_modules`: drop the trailing
  `from ..ancestry.legacy import Ancestry` + `self._ancestry = Ancestry()` block
  (runner.py:258-261). Also drop the `self._ancestry: Any = None` line in
  `__init__`.
- Delete the legacy `_run_ancestry_prediction` method (runner.py:390-443).
- `_run_with_ancestry`: replace the `if self._use_new_ancestry: … else: …` branch
  (runner.py:336-339) with an unconditional call to
  `self._run_ancestry_prediction_new()`.
- Module-level `run_pipeline`: drop the `use_new_ancestry` parameter; call
  `PipelineRunner(args)`.

**`genotools/cli/__init__.py`:**
- Delete the `main_new()` function and remove `"main_new"` from `__all__`.

**`setup.py`:**
- Remove the `'genotools-new=genotools.cli:main_new'` console-script line.

`python -m genotools` → `genotools.cli:main` → `run_pipeline(args)` now always
uses the new `AncestryModel` path.

### 2. Deletions

- `genotools/ancestry/legacy.py` (1089 lines) — no longer imported by any
  non-test module after the flip. Verified: only importers are the 4 tests being
  converted/deleted and `tests/scripts/test_ancestry_ab.py` (also deleted).
- `genotools/gwas.py` (408 lines) — orphan module **shadowed** by the
  `genotools/gwas/` package (package wins on import). All real imports use
  `genotools.gwas.<submodule>` (resolve to the package); the two
  `from genotools.gwas import run_association` hits (`gwas/__init__.py:30`,
  `gwas/pipeline.py:229`) are inside docstrings/examples. Nothing imports
  `Assoc` from it.
- `tests/scripts/test_ancestry_ab.py` — obsolete A/B harness.

### 3. Golden conversion (A → golden, not deletion)

The four differential tests prove the ported functions equal legacy
byte-for-byte. When `legacy.py` is deleted they lose their oracle, so convert them
to **golden** tests: snapshot the (already-proven-equal-to-legacy) new output to
committed fixtures **while legacy still exists**, then compare new-vs-fixture.

Tests to convert:
- `tests/unit/test_ancestry/test_preprocessing.py`
  - `test_get_raw_files_train_matches_legacy`
  - `test_get_raw_files_inference_matches_legacy` (keep the `np.repeat` spy to
    prove the missing-column fill loop runs)
- `tests/unit/test_ancestry/test_cohort.py`
  - `test_split_cohort_matches_legacy`
  - `test_split_cohort_subset_matches_legacy`

Test to **delete** (legacy no longer importable):
- `tests/unit/test_cli/test_runner.py::TestLegacyAncestryImport` (both methods —
  `test_legacy_ancestry_importable` and `test_entry_points_import`, the latter
  references `main_new`).

**Fixture form & footprint.** DataFrames → parquet, compared with
`pd.testing.assert_frame_equal` on sorted frames; split pfiles →
`tests/regression/compare.py::compare_genotypes`.

To keep the repo lean, the golden fixtures are generated from a **reduced
reference panel restricted to chromosomes 21 + 22** (~1,110 variants) instead of
the full 40,500-variant synthetic panel. Rationale:
- The reduced panel exercises the **identical code paths** (common-SNP extraction,
  merge, the inference missing-column fill — chr22 is retained in the ref but
  excluded from the inference geno via `--not-chr 22`, so the fill loop still
  fires — cohort split, and the insufficient-N pruning branch).
- Full-panel new==legacy equivalence was already proven by the differential tests
  and is preserved in git history plus the real-data parity gate. The golden's job
  going forward is **regression** protection of the ported code, which a
  path-complete reduced panel provides at a fraction of the size.

New session-scoped fixtures live in `tests/unit/test_ancestry/conftest.py`
(reduced ref bfile + labels + golden paths). Golden data files live under
`tests/unit/test_ancestry/golden/`.

**Generation is a one-time, legacy-present step** (must run before `legacy.py` is
deleted): a generator script runs both new and legacy on the reduced inputs,
asserts equality, and writes the new output as the committed fixture — so the
fixture is trustworthy because it equaled legacy at creation time.

### 4. Skip-gap fix (faithful port of `utils.upfront_check` lines ~161-189)

Port the data-driven step-skip logic into `core/validation.py`. `validate_input`
returns a typed `ValidationDecisions`; the runner applies it. Core stays
CLI-decoupled — it takes primitive request flags, not `PipelineArgs`.

```python
@dataclass(frozen=True)
class ValidationDecisions:
    skip_sex: bool = False
    skip_case_control: bool = False
    skip_het: bool = False
    disable_filter_controls: bool = False
```

`validate_input(geno_path, out_path, skip_fails=False, *, sex_requested=False,
het_requested=False, hwe_requested=False, filter_controls=False,
case_control_requested=False) -> ValidationDecisions`

Logic (only when **not** `skip_fails`; if `skip_fails`, return all-False):
1. `skip_sex` if `sex_requested` and ((no `1` and no `2` in `sex_counts`) **or**
   (no `'23'` and no `'X'` in `chr_counts`)).
2. `disable_filter_controls` if `hwe_requested` and `filter_controls` and
   (`1 not in pheno_counts`).
3. `skip_case_control` if `case_control_requested` and (`1 not in pheno_counts`
   **or** `2 not in pheno_counts`).
4. `skip_het` if `het_requested` and `var.shape[0] < 50`.
   **QUIRK TO PRESERVE:** `var.shape[0]` is the **variant** count (read from the
   pvar), even though the legacy warning text says "samples". Replicate exactly
   for parity; do **not** "fix" it to a sample count.

`sex_counts`/`pheno_counts` come from the psam; `chr_counts` from the pvar
`#CHROM` column (via the pvar read that `validate_input` already performs but
currently ignores). Preserve legacy `warnings.warn(...)` messages.

**Runner wiring (`cli/runner.py`):**
- `_convert_input_format` captures the return into `self._validation_decisions`,
  passing the request flags from typed args:
  `sex_requested=self.args.sample_qc.run_sex`,
  `het_requested=self.args.sample_qc.run_het`,
  `hwe_requested=self.args.variant_qc.run_hwe`,
  `filter_controls=self.args.variant_qc.filter_controls`,
  `case_control_requested=self.args.variant_qc.run_case_control`.
- `__init__` initializes `self._validation_decisions = ValidationDecisions()`.
- New helper `_filter_steps_by_decisions(steps) -> steps` removes `sex`,
  `case_control`, `het` per the decisions; `run()` wraps
  `get_all_enabled_steps()` with it.
- `_run_qc_pipeline`: after `legacy_args = self.args.to_legacy_dict()`, set
  `legacy_args["filter_controls"] = False` when
  `self._validation_decisions.disable_filter_controls`.

**Regression tests** (mirror `test_core.py::TestValidation`'s column-dropped psam
fixtures — copy real pgen, write a crafted psam/pvar):
- core-level (in `test_core.py::TestValidation`): each of the 4 conditions returns
  the right decision; plus `skip_fails=True` and "flags not requested" both return
  all-False.
- runner-level (in `test_runner_regression.py`): `_filter_steps_by_decisions`
  drops the right steps; `disable_filter_controls` reaches `legacy_args`.

## Verification gates

- `.venv/bin/python -m pytest tests/unit tests/regression -q` → green.
- `.venv/bin/python -m pytest tests/regression/test_parity.py -v` → 8/8
  (`.venv-stable` present).
- Static: no `genotools-new` / `use_new_ancestry` / legacy
  `_run_ancestry_prediction` / `ancestry.legacy` / `gwas.py` references left in
  `genotools/`.
- End-to-end (background, >2 min): default `genotools … --ancestry --ref-panel
  <ref> --ref-labels <labels> --min-samples 10` trains a model (GridSearchCV +
  UMAP) and produces split pfiles + JSON + model dir — confirms the NEW model runs
  by default.

## Non-blocking follow-ups (fold in if convenient)

- Bump `setup.py` `python_requires` to `>=3.10` (PEP-604 annotations already in
  the codebase).
- Guard `str(ref_panel)` in `_run_training_mode` → clean error when `--ref-panel`
  is omitted in training mode.

## Environment / workflow

- Use `.venv/bin/python` from repo root.
- Feature branch off `refactor/main` → PR/merge to `refactor/main`. NEVER touch
  `main`.
- Conventional commits, no `Co-Authored-By`/attribution (CLAUDE.md).
- Faithful-port rule: `shell_do(cmd)` → `run_command(cmd, tool_name="plink2")`;
  `psam-cols=fid,parents,sex,pheno1,phenos` stays.
- Execute via `superpowers:subagent-driven-development`; add a round-6 section to
  `REFACTOR_HARDENING.md` when done.
