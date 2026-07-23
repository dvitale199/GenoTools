# Plan: Logging / pipeline-visibility redesign (round 7)

Test-first, task-by-task. Each task: write/adjust tests → implement → run the
task's tests green → move on. Whole-diff review + full suite + e2e at the end,
then PR into `refactor/main`. Spec:
`docs/superpowers/specs/2026-07-23-logging-visibility-redesign-design.md`.

Run tests with `.venv/bin/python -m pytest`. Branch:
`refactor/logging-visibility-redesign`.

---

## Task 1 — `RunLog` writer + custom handler (`core/logging.py`)

**Tests** (`tests/unit/test_core.py::TestRunLog`, new):
- `write_banner` then `begin_section("callrate_prune")` writes the `===== callrate_prune =====` header.
- `write_record` (INFO record with a `.step`) appends a full-format structured line inside the open section.
- `append_raw` buffers; nothing hits disk until `end_section`; `end_section(raw_path)` appends the buffered raw to the consolidated file **after** the structured lines AND writes the same raw to `raw_path`.
- Ordering is deterministic: banner < header < structured < raw (assert byte order in the file).
- `end_section(None)` writes raw into the consolidated only (no per-step file).
- `write_summary([...])` appends a recognizable summary block.
- `close()` is idempotent; re-`install` doesn't leak FDs (reuse existing pattern).

**Implement:** `RunLog`, `_RunLogHandler`, `raw_sink` ContextVar,
`install_run_logging(out_path, level, console) -> RunLog`, concise console
`Formatter` (with `! ` WARNING+ prefix). Keep `setup_logging` and all existing
exports intact. Update `core/__init__` / `core/logging` exports as needed.

---

## Task 2 — Universal raw harvest in `run_command` (`core/executors.py`)

**Tests** (`tests/unit/test_executors.py` or `test_core.py`):
- `_harvest_raw_log`: given a cmd list with `--out /p/x` and an existing `/p/x.log`, calls `sink.append_raw` with the file text; with `--prefix` (KING) likewise; with no `--out`/no `.log`, no-op; with `raw_sink` unset, no-op; never raises if the file read fails.
- Integration: with a `RunLog` installed as the sink and a real `run_plink2` call (real PLINK2, skip if absent), the harvested text contains a known PLINK phrase.

**Implement:** `_harvest_raw_log(cmd_list, cmd_str)` called in `run_command`
after building `cmd_result`, before the `check` raise, in try/except. Reads
`raw_sink.get()`.

---

## Task 3 — Pure `validate_input` + decoupled guard (`core/validation.py`)

**Tests** (`tests/unit/test_core.py::TestValidation`):
- `validate_input` no longer raises on an existing `{out}_all_logs.log` (guard moved); still raises on missing pgen / missing SEX / missing PHENO1.
- Skip decisions now emitted via logging: use `caplog` to assert a WARNING record (not `warnings.warn`) for each of the 4 skip conditions; decisions dataclass unchanged.
- Breakdown emitted via `logger.info` (assert a "breakdown"/"Genetic Sex" INFO record via `caplog`).
- New `guard_output_not_exists(out, skip_fails)`: raises `ValidationError` when the log exists and not `skip_fails`; returns quietly otherwise. (Adjust the two existing `test_raises_when_log_exists` / `test_skip_fails_ignores_existing_log` to target the guard.)

**Implement:** move the guard out; swap `print`→`logger.info`,
`warnings.warn`→`logger.warning`; `logger = get_logger(__name__)`; drop `import
warnings`.

---

## Task 4 — Runner reorder, sectioning, channel migration (`cli/runner.py`)

**Tests** (`tests/unit/test_cli/test_runner_regression.py`):
- `run()` calls `guard_output_not_exists` before logging setup (assert via ordering spy / monkeypatch, or that a pre-existing log raises before any handler is installed).
- `_setup_logging` no longer creates `{out}_cleaned_logs.log` (assert it does not exist after a run) and installs a `RunLog` (assert `self._runlog` is set / banner present).
- No pipeline-path `print(` remains in `runner.py` (static assert via source scan in the test, or capture stdout and assert the old "Running: …" line is absent from stdout while present as a log record).
- Existing `--warn` / final-step-failure / GWAS-arg-threading regression tests still pass unchanged.

**Implement:** reorder `run()`; wrap `_convert_input_format` in an
`input_preparation` section; per-step `begin_section`/`end_section` with
`raw_log_path`; ancestry section + `step_context("ancestry")`; migrate the 4
`print`s; `logger = get_logger(__name__)`; store/close `self._runlog`.

---

## Task 5 — End-of-run summary table (`cli/runner.py`)

**Tests** (`test_runner_regression.py`):
- After a 2-step QC run (hermetic or real-PLINK), `runlog.write_summary` receives rows with step name + removed count + PASS; per-ancestry runs produce per-label rows.
- The summary appears in the consolidated log tail and as console/log records.

**Implement:** `SummaryRow` (or a plain tuple/dict), a
`_emit_run_summary()` helper reading `state.pass_fail`/`step_results` (+ ancestry
results), called once at the end of `run()`.

---

## Task 6 — `print`→log in ancestry preprocessing + `get_logger` sweep

**Tests:** light — `caplog` asserts `get_common_snps`/`get_raw_files` emit an
INFO record where they used to print (real-PLINK differential/golden tests already
cover behavior; assert no stdout leak). Static: a test scanning the touched
modules asserts they bind `get_logger(__name__)`.

**Implement:** `ancestry/preprocessing.py` prints → `logger.info`; standardize
`logging.getLogger(__name__)` → `get_logger(__name__)` across the 11 sites in the
spec. Pure-refactor; no behavior change.

---

## Task 7 — Rewrite `tests/regression/test_logging.py` (lock the new behavior)

**Tests** (rewrite; real CLI, skip without PLINK2/data):
- Consolidated `{out}_all_logs.log`: banner + `===== ` section headers + `[callrate_prune]`/`[geno_prune]` markers + "filtering complete" + **harvested PLINK content** (a known count/warning phrase) + a summary-table marker.
- Per-step raw files `{out}_callrate.log` / `{out}_geno.log` exist and contain PLINK output; they are NOT the consolidated file.
- `{out}_cleaned_logs.log` does **not** exist.
- Logs persist without `--full_output` (default run).
- Console: captured stdout/stderr contains a concise `[callrate_prune]` line and does not dump raw PLINK.
- (If cheap) an `--ancestry` smoke asserts `{out}_{label}_{step}.log` naming.

---

## Task 8 — Wrap-up

- Full suite: `.venv/bin/python -m pytest tests/unit tests/regression -q` green.
- Parity (if `.venv-stable` present): `test_parity.py` 8/8.
- Static gates (spec): no `cleaned_logs` in `genotools/` except `utils.py`; no
  `warnings.warn` in `core/validation.py`; no pipeline-path `print(`.
- E2e (background): `--callrate --geno` run + `--ancestry` run; eyeball the
  consolidated log top-to-bottom.
- Whole-diff self-review (or `code-reviewer` subagent).
- Add a **Round 7** section to `REFACTOR.md`; flip item #0 to resolved; update the
  "Start here" pointer to item #1 (real-data parity).
- PR into `refactor/main` (markdown body, no attribution).

---

## Notes / faithful-port guardrails

- No genotype-affecting changes; `psam-cols=fid,parents,sex,pheno1,phenos` stays.
- Harvest is best-effort and must never break a run (try/except around file read +
  sink call).
- `raw_sink` unset ⇒ harvesting is a silent no-op, so library/unit callers and the
  existing step unit tests are unaffected.
- Keep `FilterResult.log` / `GWASResult.log` fields as-is (still carry stderr);
  they're orthogonal to the disk harvest.
