# Design: Logging / pipeline-visibility redesign (refactor hardening, round 7)

**Date:** 2026-07-23
**Branch:** `refactor/logging-visibility-redesign` (feature branch → PR/merge into
`refactor/main`; **never `main`**).
**Tracker item:** REFACTOR.md item **#0** — "Logging / pipeline-visibility
redesign" (do this before item #1 real-data parity and the perf items).

## Goal

Make each pipeline step's behavior **clearly visible** — on the terminal while a
run is in progress, and on disk afterwards. The current logging scored 6.5/10: a
well-built structured core (`core/logging.py`) let down by incomplete integration
(raw PLINK output dropped from disk, a dead 0-byte `cleaned_logs.log`, and three
uncoordinated user-facing channels — `print`, `logging` file-only, and
`warnings.warn`).

## Settled decisions (from the grill session)

1. **Console = curated progress stream.** One logging stream. The console shows a
   concise, step-tagged progress narrative at INFO; raw PLINK output never reaches
   the terminal. `print()` and `warnings.warn` in the pipeline path are retired
   and routed through `logging`.
2. **`{out}_cleaned_logs.log` is deleted** (file + all references). Replaced by a
   genuinely good consolidated log.
3. **Consolidated log = per-step sections, summary then raw.** `{out}_all_logs.log`
   is banner header → one section per step (`===== step =====` header, structured
   summary lines, then the verbatim harvested PLINK output for that step) → an
   end-of-run summary table.
4. **Raw source = PLINK's native `.log`.** We harvest the `{--out}.log` file PLINK
   writes automatically (authoritative: counts + warnings + command line + timing),
   not the incomplete captured `stderr`.
5. **Harvest is executor-level, every invocation.** `run_command` harvests the
   native `.log` after *every* PLINK/KING call so compound steps (GWAS/PCA,
   ancestry preprocessing) are fully captured, not just their last invocation.
6. **Logs always persist at `{out}`** regardless of `--full_output` (which now
   governs only intermediate pfiles). Raw `.log` files are harvested out of the
   temp dir before it is cleaned.
7. **Single consolidated file for everything** under `--ancestry`, sectioned by
   label; per-step raw files are `{out}_{label}_{step}.log`. (Item #9
   parallelization will revisit this; accepted.)
8. **Guard decoupled; logging set up first.** The "output already exists" re-run
   guard becomes a standalone pre-check; then logging is configured; then
   `validate_input` (now pure) logs the data breakdown + skip decisions into the
   consolidated log.
9. **Single `RunLog` writer object** owns the consolidated file handle. Structured
   records reach it via a thin custom `logging.Handler`; section headers, buffered
   raw, and the summary table go through explicit methods. One owner → deterministic
   ordering + explicit flush; unit-testable without logging globals.
10. **End-of-run summary table** (step → removed → PASS/FAIL, per ancestry group)
    to both the console and the tail of the consolidated log.
11. **Console format = concise step-tagged**: `[callrate_prune] complete: 5 samples
    removed`, with a `! ` prefix on WARNING/ERROR. No timestamp/level/logger-name.
    The file keeps the full detailed format.

### Smaller defaults (ratified)

- No new `--verbose`/`--quiet` flag this round; console is always the curated INFO
  stream.
- The consolidated file captures INFO+ structured records; the executor's
  per-command DEBUG line stays out of the file (the raw `.log` already carries the
  exact command).
- Pre-split & compound phases get their own sections (`input_preparation`,
  `ancestry_prediction`, and the PCA/GWAS sub-phases already have `step_context`),
  so their PLINK calls are attributed, not dumped into an empty-step bucket.
- Standardize on `get_logger(__name__)` (retire the bare `logging.getLogger`
  idiom in steps/runner/gwas/ancestry-model — 11 sites).
- Do **not** touch third-party warnings (sklearn/pandas) or
  `logging.captureWarnings`; only GenoTools' own `print`/`warn` in the pipeline
  path migrate.
- Out of scope: `genotools/utils.py`, `imputation.py`, `download_refs.py`,
  `dependencies.py`, `container/` (deferred/other entry points). Module-docstring
  `print(...)` examples are not touched.

## Architecture

### `core/logging.py` — new `RunLog` + wiring

```python
class RunLog:
    """Owns the consolidated {out}_all_logs.log file handle.

    Single writer for: the banner header, per-step section headers, live
    structured records (via _RunLogHandler), buffered raw PLINK blocks flushed
    at section end (also written to per-step raw files), and the final summary.
    """
    def __init__(self, path: Path) -> None: ...
    def write_banner(self) -> None: ...
    def begin_section(self, title: str) -> None: ...          # writes "===== title ====="
    def append_raw(self, command: str, text: str) -> None: ... # buffers into current section
    def end_section(self, raw_log_path: Path | None) -> None: ...# flush buffered raw → file + per-step
    def write_record(self, record: logging.LogRecord) -> None: ...# live structured line (full format)
    def write_summary(self, rows: list[SummaryRow]) -> None: ...
    def close(self) -> None: ...

class _RunLogHandler(logging.Handler):
    """Routes genotools structured records into a RunLog (INFO+)."""
    def __init__(self, runlog: RunLog): ...
    def emit(self, record): self.runlog.write_record(record)
```

- **Raw sink ContextVar.** `raw_sink: ContextVar[Optional[RunLog]]` (mirrors the
  existing `current_step`). Set when run-logging is installed; read by the executor.
  `None` → harvesting is a no-op (library/unit use).
- **`install_run_logging(out_path, level="INFO", console=True) -> RunLog`.** Builds
  the `RunLog`, writes the banner, attaches `_RunLogHandler` + a console
  `StreamHandler` (concise formatter) to the `"genotools"` logger, sets `raw_sink`,
  returns the `RunLog`. Closes/clears prior handlers first (FD-leak safe, like
  today's `setup_logging`).
- **Formatters.** File/full: reuse today's
  `"%(asctime)s [%(levelname)s] %(step)s %(name)s: %(message)s"`. Console/concise:
  a custom `Formatter` producing `"[step] message"` (from the `record.step` the
  existing `StepContextFilter` sets), prefixed `"! "` for `levelno >= WARNING`.
- **`setup_logging(...)` retained** for existing/simple callers (e.g. `test_core`
  `TestLogging`); the runtime uses `install_run_logging` instead of the plain
  `FileHandler`. Keep `banner()`, `step_context`, `current_step`,
  `StepContextFilter`, `get_logger`, `set_step`/`clear_step`.

### `core/executors.py` — universal raw harvest in `run_command`

After the subprocess returns and `cmd_result` is built, **before** the
`check`/raise branch (so a failing step's PLINK log is still captured):

```python
_harvest_raw_log(cmd_list, cmd_str)   # best-effort; never raises
```

`_harvest_raw_log`:
- Find the output prefix: token after `--out`, or after `--prefix` (KING).
- `log_path = f"{prefix}.log"`; if it exists, read it (`errors="replace"`).
- `sink = raw_sink.get()`; if `sink is not None: sink.append_raw(cmd_str, text)`.
- Wrapped in try/except → harvesting must never break execution.

Every PLINK/PLINK2 call in the codebase passes `--out <prefix>` (QC via
`run_plink2`, PCA/assoc via `run_command`, ancestry preprocessing via
`run_command`, conversions in `genotypes.py`); KING passes `--prefix`. So one hook
in `run_command` captures all invocations. Attribution to a section is by the
currently-open `RunLog` section (opened by the runner), not by parsing.

### `core/validation.py` — pure validation, logged

- **Remove** the `{out}_all_logs.log`-exists guard from `validate_input` (moves to
  the runner pre-check). `validate_input` keeps: pgen-exists check, SEX/PHENO1
  column checks, the data breakdown, and the skip decisions.
- Replace the breakdown `print(...)` (lines 76–93) with `logger.info(...)` and the
  five `warnings.warn(...)` skip messages (103–129) with `logger.warning(...)`
  (same text). `logger = get_logger(__name__)`.
- **New** `guard_output_not_exists(out_path, skip_fails)` (module-level): raises
  `ValidationError` if `{out}_all_logs.log` exists and not `skip_fails`. Called by
  the runner *before* logging setup.

### `cli/runner.py` — orchestration, sectioning, summary

- **`run()` reorder:** `guard_output_not_exists(...)` → `_setup_logging()` (build
  `RunLog`, install handlers, write banner) → `_convert_input_format()` (wrapped in
  a `input_preparation` section so conversion + `validate_input` land in the log) →
  steps.
- **`_setup_logging()`** becomes: delete any stale `{out}_all_logs.log`; **stop
  creating `cleaned_logs.log`**; `self._runlog = install_run_logging(out_path,
  console=True)`. Store `self._runlog` on the runner (and clean up / `close()` in
  the `finally`).
- **`_run_qc_pipeline`** per step: compute `raw_log_path`
  (`{out}_{step}.log` or `{out}_{label}_{step}.log`), wrap the run in
  `runlog.begin_section(section_title)` / `runlog.end_section(raw_log_path)`. Replace
  `print("Running: …")` / `print("cannot be run")` / `print("failed but
  continuing")` with `logger.info`/`logger.warning`. Section title is label-aware
  under ancestry (e.g. `EUR / callrate_prune`); the console `[step]` tag likewise.
- **`_run_single_step`** kinship "Relatedness … Linux only" `print` → `logger.warning`.
- **`_run_ancestry_prediction_new`**: wrap in `begin_section("ancestry_prediction")`
  + a `step_context("ancestry")` so the ported preprocessing logs + PLINK raw are
  attributed. `end_section({out}_ancestry.log)`.
- **End-of-run summary:** a helper assembles `SummaryRow`s from `state.pass_fail` +
  `state.step_results` (+ per-ancestry results) and calls
  `runlog.write_summary(...)`; the same rows print to the console via `logger.info`.

### `ancestry/preprocessing.py` — logged, not printed

Replace `print("Getting Common SNPs")` and the "Labeled Reference Ancestry
Counts" block (lines 18, 141–146) with `logger.info(...)`
(`logger = get_logger(__name__)`). No behavior change.

### `get_logger` standardization

Change `logger = logging.getLogger(__name__)` → `logger = get_logger(__name__)`
(and drop now-unused `import logging` where only used for that) in: the 8 QC steps,
`gwas/steps/association.py`, `gwas/steps/pca.py`, `gwas/pipeline.py`,
`cli/runner.py`, `ancestry/model.py`. Keeps every module under the `genotools.*`
namespace so records reach the `RunLog` handler.

## Resulting artifacts (after any run)

| File | Contents |
|------|----------|
| `{out}_all_logs.log` | banner → per-step sections (header, structured summary, raw PLINK) → summary table |
| `{out}_{step}.log` (`{out}_{label}_{step}.log`) | isolated verbatim PLINK output for that step (all its invocations) |
| `{out}_cleaned_logs.log` | **gone** |

Console: concise step-tagged progress + final summary table. No raw PLINK.

## Verification gates

- `.venv/bin/python -m pytest tests/unit tests/regression -q` → green.
- `tests/regression/test_logging.py` (rewritten) asserts: consolidated log has
  section headers + structured `[callrate_prune]`/`[geno_prune]` markers +
  "filtering complete" + **real PLINK content** (e.g. a variant/sample count phrase
  from the harvested `.log`); per-step raw files exist and contain PLINK output; **no**
  `cleaned_logs.log`; a summary table line is present; console (captured
  stdout/stderr) shows a concise `[step]` line and no full-format timestamp noise.
- Unit: `RunLog` section/raw/summary ordering (no disk-race, deterministic bytes);
  `_harvest_raw_log` parses `--out`/`--prefix` and no-ops without a sink;
  `validate_input` emits `logger.warning` (via `caplog`) instead of `warnings.warn`;
  `guard_output_not_exists` raises/skips correctly.
- Parity unaffected: `tests/regression/test_parity.py` compares genotype content +
  IDs + lambda, **not** logs, and shells out to `.venv-stable`. Still 8/8 with
  `.venv-stable` present.
- End-to-end (background): a `--callrate --geno` run and an `--ancestry` run produce
  the expected files; eyeball the consolidated log reads cleanly top-to-bottom.
- Static: no `cleaned_logs` references left in `genotools/` (except deferred
  `utils.py`); no `warnings.warn` in `core/validation.py`; no pipeline-path
  `print(` outside module-docstring examples.

## Environment / workflow

- Use `.venv/bin/python` from repo root.
- Feature branch off `refactor/main` → PR/merge to `refactor/main`. NEVER touch `main`.
- Conventional commits, no `Co-Authored-By`/attribution (CLAUDE.md).
- Faithful-port rule preserved: `psam-cols=fid,parents,sex,pheno1,phenos` stays;
  no genotype-affecting changes.
- Add a round-7 section to `REFACTOR.md` when done; mark item #0 resolved.
