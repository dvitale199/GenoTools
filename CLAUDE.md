# GenoTools — AI Development Guide

Conventions for AI assistants working on this codebase. This file is loaded
every session, so it holds only what changes what you do. Everything else is a
pointer:

| Need | Read |
|---|---|
| CLI flags and their semantics | `docs/cli_args.md` |
| Test suite, parity harness, real-cohort runs | `TESTING.md` |
| What changed in 2.0 and why | `MIGRATION_2.0.md` |
| Refactor history, tracked remaining work | `REFACTOR.md` |
| JSON report shape | `docs/json_output_overview.md` |

---

## Architecture

2.0 replaced the `SampleQC`/`VariantQC`/`Ancestry` classes with pure functions
over frozen config dataclasses. Guidance anywhere referring to
`genotools/qc.py`, `pipeline.py`, `run_het_prune()`, `shell_do()` or
`concat_logs()` predates 2.0 and is wrong.

```
argv → parse_args() → PipelineArgs → PipelineRunner.run()
  → validate_input() → ValidationDecisions (steps the data rules out)
  → convert to pfiles → ancestry prediction (optional)
  → _run_with_ancestry (per label)  |  _run_qc_only (flat)
      → _run_qc_pipeline → _run_single_step → filter_<step>(data, config, out)
         callrate → sex → het → related → kinship_check →
         geno → case_control → haplotype → hwe → ld → assoc
  → JSON report + pfiles + {out}_all_logs.log
```

- `cli/` — `parser.py` typed args + validation, `runner.py` orchestration and
  step dispatch, `output.py` JSON report
- `qc/` — `steps/*.py` one function per step, `config.py` frozen configs,
  `results.py` `FilterResult` / `STEP_REPORT`
- `core/` — `executors.py` (`run_plink2`, `run_plink`, `run_king`,
  `run_command`), `genotypes.py` (`GenotypeData`), `logging.py` (`RunLog`,
  `step_context`), `validation.py`, `exceptions.py`
- `ancestry/` — `AncestryModel`, preprocessing, reducers, cohort split
- `gwas/` — PCA and association

**Off the 2.0 path:** `utils.py`, `imputation.py`, `container/`, `prototype/`.
Never import them from `cli`/`core`/`qc`/`gwas`/`ancestry`; that decoupling is
tested.

---

## Contracts

**Every QC step is a module-level function**, no classes or shared state:

```python
def filter_<step>(data: GenotypeData, config: <Step>Config, out_path: Path) -> FilterResult:
```

Copy the nearest existing step rather than a template in this file — they are
all the same shape, and real code can't drift out of date.

**Return a `FilterResult`** (`qc/results.py`); never hand-build its `to_dict()`
output. Its `metrics` holds **counts only** — `cli/output.py` renders each
metric as a row in a long `(step, metric, pruned_count)` table, so a mode name
or a derived bound would land under a column meaning "how many were pruned."
Descriptive context goes to the log. Register report metric names in
`STEP_REPORT` so skipped and failed steps emit matching zero rows.

**Validate thresholds in the config, not the step.** Frozen dataclasses
validating in `__post_init__` via `ThresholdConfig._validate_*`
(`core/config.py`) fail before PLINK runs, and `parse_args` surfaces it as a
clean `ERROR:` line.

**Raise, don't return, on failure** — a `GenoToolsError` subclass
(`ValidationError`, `QCError`, `ExternalToolError`, …). `_run_qc_pipeline`
decides what happens: warn-and-continue by default, fail-fast under
`--no-warn`. `main()` renders those plus `FileNotFoundError` / `ValueError` /
`TypeError` as a one-line `ERROR:` with exit 1. Anything else is a bug and
keeps its traceback.

**Log via `get_logger(__name__)` inside `with step_context(step_name):`.**
Never `print()` or `warnings.warn()`. The executors harvest PLINK's own `.log`
themselves — there is no `concat_logs()` call to make.

**Keep the pfile chain intact.** Each step reads `{geno}.pgen/.pvar/.psam` and
writes `{out}.pgen/.pvar/.psam`, or the next step has nothing to read.
`run_plink2` adds `--make-pgen psam-cols=fid,parents,sex,pheno1,phenos` by
default; pass `make_pgen=False` for report-only commands (`--het`,
`--indep-pairwise`). Outlier files are tab-separated with a `#FID` header.

GenoTools is **synchronous**. Do not introduce async.

---

## Traps

These caused real bugs. They are why this file exists.

**Two run shapes.** `_run_with_ancestry` and `_run_qc_only` are separate paths,
and anything per-step must be reachable from both. Production runs ancestry
once, then QCs each group as a separate *flat* job — so a flag read only in the
ancestry branch silently does nothing exactly where it matters (`--amr-het`,
round 11). When a setting can vary per dataset, resolve it in one shared place
(`SampleQCArgs.het_config_for()`), not at each call site.

**Never hardcode an ancestry label.** The vocabulary is user-supplied — from
the `--ref-labels` TSV via an unconstrained `LabelEncoder`, or a pickled
model's encoder. Admixture detection also invents `CAH`, which appears in no
panel. `label == "AMR"` was the second defect behind `--amr-het`.

**`str(None)` is `"None"`,** which reaches PLINK as a filename and produces
`Failed to open None.bed`. Validate optional paths up front, at the CLI
boundary — this is how `--ancestry` without `--ref-panel` used to fail
(round 12).

**Data preconditions are skips, not crashes.** A step the data rules out gets a
`<step>_skip_reason()` in `core/validation.py`, wired into `_skip_reasons()`,
so the report distinguishes "not requested" from "requested but impossible".

**KING is Linux-only.** Check `platform.system()` and warn.

**`python -m genotools` puts the cwd first on `sys.path`** — from the repo root
it imports the working tree whatever virtualenv you use. Matters when comparing
against `.venv-stable`.

**Flag policy.** Hyphenated spellings only for anything added in 2.0+.
Underscore spellings are deprecated 1.x aliases in `_DEPRECATED_SPELLINGS`; a
*removed* flag goes in `_REMOVED_FLAGS`, which gives a targeted error naming
the replacement instead of argparse's bare "unrecognized arguments".
`--container` / `--singularity` / `--cloud` follow the same pattern.

---

## Adding a QC step

Seven touch points. Miss one and the step silently never runs.

1. **Config** — frozen dataclass in `qc/config.py`, validating in `__post_init__`
2. **Function** — `qc/steps/<name>.py`, shape above
3. **Exports** — `qc/steps/__init__.py` and `qc/__init__.py` (imports + `__all__`)
4. **Report** — entry in `STEP_REPORT` (`qc/results.py`)
5. **CLI** — flag in `create_parser()`, field on `SampleQCArgs`/`VariantQCArgs`,
   `to_<step>_config()`, wiring in `parse_args()`, entry in
   `get_enabled_*_steps()`
6. **Dispatch** — register in `_initialize_modules()`, branch in
   `_run_single_step()`
7. **Tests, and `docs/cli_args.md`**

---

## Testing

```bash
.venv/bin/python -m pytest tests/unit -q        # ~25s, 616 tests
.venv/bin/python -m pytest tests/regression -q  # ~7min, needs PLINK
```

- **Revert-check any test claiming to gate a fix.** Break the production call
  site and confirm the test fails. `REFACTOR.md` item 17 records round-7 tests
  that passed against the broken code because they drove a helper rather than
  the path the bug lived in.
- **Test through the layer that owns the behaviour.** The parser has three
  worth covering — dataclass, `_parse_*` helper, end-to-end `parse_args` — and
  bugs hide in the gaps.
- **Runner regressions** go in `tests/unit/test_cli/test_runner_regression.py`.
- **Split pure logic out of PLINK-driven steps** so it can be tested against a
  hand-built frame. `select_het_outliers` is the pattern.
- **Changing the JSON report shape** means regenerating
  `tests/regression/golden/` via `tests/scripts/generate_golden.py` — a
  contract change, not a chore.

`.venv` is `pip install -e .`; `.venv-stable` is the frozen parity baseline
(currently 1.3.6). See `TESTING.md`.

---

## Writing it up

- **Commit messages:** conventional prefix (`feat:`, `fix:`, `refactor:`,
  `docs:`, `test:`), concise first line, bullets for anything non-trivial.
  **Never** add `Co-Authored-By:` or any attribution trailer.
- **PR descriptions:** markdown, `## Summary` + `## Test plan`, ready to paste
  into GitHub. **Never** include "Generated with Claude Code" or similar.
- **`REFACTOR.md`:** append a *new* round entry; never edit historical ones.
  Deferred work goes in "Remaining work" with a number.
