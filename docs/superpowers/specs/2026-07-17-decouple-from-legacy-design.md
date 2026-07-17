# Design: Decouple the new CLI from legacy (refactor hardening, item #6)

**Date:** 2026-07-17
**Branch target:** `refactor/main` (feature branch → PR into `refactor/main`; **never `main`**).
**Tracker item:** REFACTOR_HARDENING.md remaining-work #6 ("Decouple from legacy").

## Goal

Make the new modular pipeline (`cli/`, `core/`, `qc/`, `gwas/`, and the new
`ancestry/` package) **stand entirely on its own** — importing zero top-level
legacy modules. Today the new CLI still reaches into legacy code in two ways, both
in `cli/runner.py`:

1. Four `from ..utils import …` calls (`gt_header`, `bfiles_to_pfiles`,
   `vcf_to_pfiles`, `upfront_check`).
2. An `importlib` file-path hack that loads legacy `ancestry.py` (the `Ancestry`
   class) because the top-level `ancestry.py` *module* name-collides with the new
   `ancestry/` *package*.

Additionally, even the *new* ancestry path (`genotools-new`) is not actually
legacy-free: it delegates to the legacy `Ancestry` class for three helper methods.
This design closes both gaps so the new pipeline is self-contained, which is the
prerequisite for the eventual legacy removal (Phase 5/6).

## Background: the three subsystems are at different migration stages

| Subsystem | State | Legacy coupling |
|-----------|-------|-----------------|
| **QC** | ✅ fully migrated | No `qc.py` exists; `qc/` package is the only impl; runner calls it directly. |
| **GWAS** | ✅ fully migrated | `gwas/` package **shadows** the orphaned `gwas.py` file (unreachable dead code). |
| **Ancestry** | ⚠️ **partially** migrated | New `ancestry/` package (`AncestryModel`) exists but is **not** a complete replacement (see below). |

Ancestry is the only hard part. Two facts make it so:

- The **default** `genotools` command still runs the **legacy** `Ancestry`
  (`run_pipeline(..., use_new_ancestry=False)`).
- The **new** path (`genotools-new`, `_run_ancestry_prediction_new`) still calls
  the **legacy** `Ancestry` for `get_raw_files()`, `split_cohort_ancestry()`, and
  `clean_up()`.

## Scope

This is **one cohesive spec** delivering two parts:

- **6a — import decoupling** (behavior-preserving; safe).
- **6b — finish the new ancestry path** so `genotools-new` is 100% legacy-free
  (faithful port; changes no output).

**Explicitly OUT of scope (gated on tracker item #1, real-data parity):**

- Flipping the default `genotools` engine from legacy `Ancestry` to `AncestryModel`.
- Deleting the legacy `Ancestry` class / the `genotools-new` A/B entry point.
- Deleting the legacy cluster files (`utils.py`, `ancestry/legacy.py`, `gwas.py`,
  `imputation.py`) — that wholesale deletion is Phase 5/6.

After this work: legacy `Ancestry` stays as the **default engine** (`genotools`);
`genotools-new` becomes a fully standalone new pipeline.

**Validation is cross-branch (`main` vs `refactor/main`).** The handoff parity gate
compares `genotools` on `main` (old everything) against `genotools-new` on
`refactor/main` (new everything) — a true full old-vs-new A/B, which is exactly
what the existing parity harness does (`.venv-stable` = old baseline). Because that
comparison is on **genotype content** (split pfiles, labels, pruned samples), not
log output, it is **robust to the in-flight logging updates** being made on
`refactor/main` before handoff. This is why we do *not* flip the default now: the
final round (gated on that parity gate, owned by the other dev) flips the default
and deletes legacy.

---

## 6a — Import decoupling

Relocate the four `utils.py` helpers into `core/`, reimplemented **legacy-free**
(using `core/executors` + `core/exceptions`). The `utils.py` originals are left
untouched (still used by the legacy cluster) and are deleted wholesale in Phase 5/6.

| Legacy helper | New home | Approach |
|---|---|---|
| `gt_header()` | `core/logging.py` → `banner()` | Copy the ASCII banner string; `_setup_logging()` calls it. |
| `bfiles_to_pfiles()` (bfile input) | **already exists**: `GenotypeData.to_pfile()` | Runner bfile branch → `GenotypeData.from_path(bfile).to_pfile(bfile)`. |
| `vcf_to_pfiles()` | new `GenotypeData.from_vcf()` in `core/genotypes.py` | vcf → bed (`--vcf --make-bed`) → pfile (`to_pfile`), remove intermediate bed/bim/fam; via `run_plink2`. |
| `upfront_check()` | new `core/validation.py` → `validate_input(args)` | Typed-args validation + data breakdown print; raises `ValidationError`. |

Then move legacy `Ancestry` into the package and kill the importlib hack:

- `git mv genotools/ancestry.py genotools/ancestry/legacy.py`. This resolves the
  module/package name collision the same way `gwas/` already did — except we
  *retain* it (it is still the default `genotools` engine) rather than orphaning it.
- Its `from genotools.utils import …` / `from genotools.dependencies import …`
  keep working (legacy→legacy / shared is acceptable).
- **Not** exported from `ancestry/__init__.py` — keeps the new public API clean.
- `cli/runner.py`: replace the importlib block (lines ~258–272) with
  `from ..ancestry.legacy import Ancestry`.

### Deferred (tracked separately): the `upfront_check` skip-gap

The current new runner calls `upfront_check(...)` but **discards its return value**
(`runner.py:332`). In legacy code that return dict carried *data-driven auto-skips*
(skip sex-prune if no X chromosome; skip het if <50 samples; skip case/control if
only one class present; drop `filter_controls` if no controls). The new pipeline
has **no equivalent logic anywhere**, so those auto-skips silently don't happen —
a pre-existing latent behavior difference, independent of decoupling.

**Decision:** keep this out of scope. `validate_input` ports the checks + breakdown
only (the observable behavior of the current new path); it does **not** re-implement
the skip logic. Fixing the gap changes scientific behavior and belongs with the
parity work (item #1), not a decoupling round. Recorded here so it is not lost.

---

## 6b — Finish the new ancestry path (make `genotools-new` legacy-free)

The new path reaches into legacy `Ancestry` for exactly three methods, plus one
transitive dependency (`get_common_snps` in `utils.py`). Port surface:

| Legacy method | ~Size | Nature | Risk |
|---|---|---|---|
| `clean_up()` | 15 | deletes temp files by extension | trivial |
| `split_cohort_ancestry()` | 55 | per-ancestry `--keep --make-pgen`, collects insufficient-N pruned samples | mechanical |
| `get_raw_files()` | 190 | variant prune → common-SNP extract → allele-flip reconciliation → `.raw` export; branches on train vs inference | **high — output-critical** |
| `get_common_snps()` (`utils.py`) | 55 | PLINK1.9 `--flip` + plink2 `--extract` + allele-merge | **high — output-critical** |

### New home

- `ancestry/preprocessing.py` — `get_raw_files` + `get_common_snps` + a `clean_up`
  helper, as **pure functions taking explicit args** (geno_path, ref_panel,
  ref_labels, out_path, train, common_snps, …) instead of reading `self`.
- `ancestry/cohort.py` — `split_cohort_by_ancestry(labels_path, geno_path,
  out_path, subset, min_samples)`.
- Rewire `_run_ancestry_prediction_new` / `_run_training_mode` /
  `_run_inference_mode` to call these functions instead of `self._ancestry.*`.
  State currently set as fields on `self._ancestry` (train, model_path, subset,
  min_samples, …) becomes explicit function arguments sourced from
  `self.args.ancestry.*`.

### The principle that keeps 6b safe: faithful port

The port must issue the **same PLINK commands** and the **same pandas operations**
so that `genotools-new`'s output is **unchanged by the port**. The only new-vs-legacy
ancestry difference remains the ML model (which already existed); the port introduces
none. PLINK calls go through `core/executors.run_command` (return-code-checked — a
strict improvement on the *failure* path, identical on *success*), consistent with
the rest of the new code.

### Simplifications that fall out naturally

- **Inference round-trip collapses.** Today `_run_inference_mode` writes the model's
  common-SNP list to a `*.common_snps` file purely so legacy `get_raw_files(train=False)`
  can read it back off disk (and fakes a `model_path`). The ported function takes the
  `common_snps` list directly — no temp file, no fake path.
- **Container branch not carried over.** The new path always sets
  `containerized=False` (`runner.py:506`), so `get_raw_files`'s containerized branch
  (legacy lines ~122–134) is unreachable in the new path and is intentionally not
  ported. Recorded as a deliberate scope decision: `genotools-new` does not support
  `--container` (already true today).

### `get_common_snps` duplication (accepted)

`get_common_snps` still has legacy callers (`ancestry/legacy.py`). It is **copied**
into `ancestry/preprocessing.py` (new, `run_command`-based); the `utils.py` original
stays for the legacy cluster. Same "duplicate-then-delete-wholesale-in-Phase-5/6"
pattern as the other utils. The alternative (legacy importing from the new package)
creates backwards coupling and is rejected.

---

## Testing

Two layers: an **in-suite** proof that the code is correct/faithful (below), and
the **handoff correctness gate**, which is the cross-branch parity run — `genotools`
on `main` vs `genotools-new` on `refactor/main` — owned by the other dev (tracker
item #1). Everything here is **log-agnostic** (compares genotype content, not log
output), so it composes with the in-flight logging changes on `refactor/main`.

### 6b — before/after golden invariant (in-suite faithful-port check)

Pin `genotools-new`'s ancestry output as a golden **before** the port, then assert
**byte-identical** output **after**. Because the port is faithful, they must match
exactly. This is a within-branch before/after check — it does **not** require
`.venv-stable`, and it is independent of the cross-branch gate — and it directly
proves the port changed nothing.

- Capture on synthetic ancestry input, both **training** and **inference** modes.
- Compare: per-ancestry split pfiles (IDs + genotype content via existing
  `compare_genotypes()`), predicted-labels file, and the pruned-samples set.

### Unit tests

- `core/genotypes.py::from_vcf` — vcf → pfile round-trip; intermediate bed cleaned up;
  raises on missing/failed conversion.
- `core/validation.py::validate_input` — raises on missing `.pgen`, on
  `{out}_all_logs.log` already present (unless `skip_fails`), on missing SEX/PHENO1
  columns; prints the data breakdown.
- `core/logging.py::banner` — returns the expected header string.
- `ancestry/preprocessing.py` — `get_common_snps` and `get_raw_files` on synthetic
  ref+geno (train and inference); `clean_up` removes the right extensions.
- `ancestry/cohort.py::split_cohort_by_ancestry` — split counts, `min_samples`
  threshold moves samples to pruned, `subset` restricts groups.

### Regression / smoke

- `from genotools.ancestry.legacy import Ancestry` imports and constructs.
- Both entry points import: `genotools.cli:main`, `genotools.cli:main_new`.
- Existing parity harness (`tests/regression/test_parity.py`) green with
  `.venv-stable` present — proves 6a's conversion + validation produce identical
  pfiles, and both `genotools` and `genotools-new` still run. This is the same
  cross-branch mechanism the handoff gate uses (old baseline vs new branch); the
  ancestry A/B runs `genotools-new` on the new side.

## Verification (end of round)

1. **Static — coupling gone:** `grep` shows zero `from ..utils import` and zero
   `importlib`/`spec_from_file_location` referencing ancestry anywhere in `cli/`,
   `core/`, `qc/`, `gwas/`, or the new `ancestry/` package (excluding
   `ancestry/legacy.py`). Zero `self._ancestry.get_raw_files` /
   `split_cohort_ancestry` / `clean_up` in the *new* path.
2. **Import smoke:** `python -c "import genotools, genotools.core, genotools.qc, genotools.ancestry, genotools.gwas, genotools.cli"`.
3. **CLI smoke:** `python -m genotools --help`.
4. **Full suite:** `tests/unit tests/regression` green (current 374 + new tests).
5. **Parity:** `tests/regression/test_parity.py` green with `.venv-stable`.

## Risks / mitigations

- **A `get_raw_files`/`get_common_snps` port drifts from legacy** (changing
  ancestry output) → the before/after golden invariant catches any drift byte-for-byte,
  in both train and inference modes.
- **Missed legacy reach-back in the new path** → the static grep (verification #1)
  plus import/CLI smoke assert the new path never touches `self._ancestry`.
- **Relocating `ancestry.py` breaks a hidden importer** → import smoke + full suite;
  legacy path (`genotools`) exercised by parity harness.
- **Interference with the in-flight real-data parity run** (owned by another dev):
  6a touches the startup hot path (input conversion, validation, logging) but is
  behavior-preserving; 6b changes only the `genotools-new` path, not the default
  `genotools`. Coordinate the merge; the default engine's behavior is untouched.

## Out of scope (recap)

- Flipping the default to new ancestry; deleting legacy `Ancestry` and the
  `genotools-new` entry point — gated on item #1 (real-data parity), final round.
- The `upfront_check` skip-gap fix — tracked, belongs with the parity work.
- Phase 5/6 wholesale deletion of the legacy cluster (`utils.py`,
  `ancestry/legacy.py`, `gwas.py`, `imputation.py`).
- Items #9/#10 (performance) — separate.
