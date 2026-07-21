# Flip default to new ancestry + cleave legacy — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make `refactor/main` legacy-free: flip the default `genotools` command to the new `AncestryModel`, delete the dead legacy `ancestry/legacy.py` and orphan `gwas.py`, convert the 4 differential tests to golden, and port the `upfront_check` data-driven step-skip logic into `core/validation.py`.

**Architecture:** The new modular pipeline (`cli/`, `core/`, `qc/`, `gwas/`, `ancestry/`) already runs standalone. This round removes the last legacy escape hatches (the `use_new_ancestry` A/B switch and the two dead legacy modules) and restores one behavior that was dropped in the rewrite (data-driven step skipping), so a real-cohort run matches `main`.

**Tech Stack:** Python 3.8+, pandas, PLINK2 (via `genotools.core.executors.run_command`), pytest, pyarrow (parquet golden fixtures).

## Global Constraints

- Work on a feature branch off `refactor/main`; PR/merge into `refactor/main`. **NEVER touch `main`.**
- Conventional commits; **no** `Co-Authored-By`/attribution lines (CLAUDE.md).
- Use `.venv/bin/python` from repo root `/Users/vitaled2/Desktop/projects/GenoTools`.
- Faithful-port rule: `shell_do(cmd)` → `run_command(cmd, tool_name="plink2")`; the `--make-pgen` flag stays `psam-cols=fid,parents,sex,pheno1,phenos`.
- Preserve legacy `warnings.warn(...)` message text when porting skip logic.
- **QUIRK to preserve:** the het-skip threshold uses the **variant** count (`var.shape[0]` from the pvar), even though the legacy warning says "samples". Do not "fix" it.

---

### Task 1: Port data-driven step-skip logic into `core/validation.py`

**Files:**
- Modify: `genotools/core/validation.py`
- Test: `tests/unit/test_core.py` (class `TestValidation`)

**Interfaces:**
- Consumes: nothing new.
- Produces:
  - `ValidationDecisions` dataclass (frozen) with bool fields `skip_sex`, `skip_case_control`, `skip_het`, `disable_filter_controls` (all default `False`).
  - `validate_input(geno_path, out_path, skip_fails=False, *, sex_requested=False, het_requested=False, hwe_requested=False, filter_controls=False, case_control_requested=False) -> ValidationDecisions` — same raising behavior as today, now returning decisions. Existing positional callers `validate_input(geno, out, skip_fails=...)` keep working (new params are keyword-only with defaults) and get an all-`False` `ValidationDecisions`.

- [ ] **Step 1: Write failing tests** in `tests/unit/test_core.py`, appended inside `class TestValidation`.

```python
    # --- data-driven step-skip decisions (ported from legacy upfront_check) ---

    def _pfile_with(self, geno: Path, dst: Path, *, psam=None, pvar_rows=None):
        """Copy geno's pgen; write a crafted psam (DataFrame) and/or truncate pvar."""
        import pandas as pd
        dst.with_suffix(".pgen").write_bytes(geno.with_suffix(".pgen").read_bytes())
        if psam is None:
            dst.with_suffix(".psam").write_bytes(geno.with_suffix(".psam").read_bytes())
        else:
            psam.to_csv(dst.with_suffix(".psam"), sep="\t", index=False)
        src_pvar = geno.with_suffix(".pvar").read_text().splitlines()
        if pvar_rows is None:
            dst.with_suffix(".pvar").write_text("\n".join(src_pvar) + "\n")
        else:
            header = [ln for ln in src_pvar if ln.startswith("#")]
            body = [ln for ln in src_pvar if not ln.startswith("#")][:pvar_rows]
            dst.with_suffix(".pvar").write_text("\n".join(header + body) + "\n")

    def test_skip_sex_when_no_sex_data(self, geno: Path, tmp_path: Path):
        import pandas as pd
        from genotools.core.validation import validate_input
        psam = pd.read_csv(geno.with_suffix(".psam"), sep=r"\s+")
        psam["SEX"] = 0
        bad = tmp_path / "nosex"
        self._pfile_with(geno, bad, psam=psam)
        d = validate_input(bad, tmp_path / "o1", sex_requested=True)
        assert d.skip_sex is True

    def test_skip_sex_when_no_x_chrom(self, geno: Path, tmp_path: Path):
        from genotools.core.validation import validate_input
        # Keep only chr-1 variants (no X) by truncating to the first block.
        bad = tmp_path / "nox"
        self._pfile_with(geno, bad, pvar_rows=100)
        d = validate_input(bad, tmp_path / "o2", sex_requested=True)
        assert d.skip_sex is True

    def test_disable_filter_controls_when_no_controls(self, geno: Path, tmp_path: Path):
        import pandas as pd
        from genotools.core.validation import validate_input
        psam = pd.read_csv(geno.with_suffix(".psam"), sep=r"\s+")
        psam["PHENO1"] = 2  # all cases, no controls
        bad = tmp_path / "nocontrol"
        self._pfile_with(geno, bad, psam=psam)
        d = validate_input(bad, tmp_path / "o3", hwe_requested=True, filter_controls=True)
        assert d.disable_filter_controls is True

    def test_skip_case_control_when_only_cases(self, geno: Path, tmp_path: Path):
        import pandas as pd
        from genotools.core.validation import validate_input
        psam = pd.read_csv(geno.with_suffix(".psam"), sep=r"\s+")
        psam["PHENO1"] = 2
        bad = tmp_path / "onlycases"
        self._pfile_with(geno, bad, psam=psam)
        d = validate_input(bad, tmp_path / "o4", case_control_requested=True)
        assert d.skip_case_control is True

    def test_skip_het_when_few_variants(self, geno: Path, tmp_path: Path):
        from genotools.core.validation import validate_input
        bad = tmp_path / "fewvar"
        self._pfile_with(geno, bad, pvar_rows=40)  # < 50 variants
        d = validate_input(bad, tmp_path / "o5", het_requested=True)
        assert d.skip_het is True

    def test_no_decisions_when_skip_fails(self, geno: Path, tmp_path: Path):
        import pandas as pd
        from genotools.core.validation import validate_input
        psam = pd.read_csv(geno.with_suffix(".psam"), sep=r"\s+")
        psam["SEX"] = 0
        bad = tmp_path / "sf"
        self._pfile_with(geno, bad, psam=psam)
        d = validate_input(bad, tmp_path / "o6", skip_fails=True, sex_requested=True)
        assert d.skip_sex is False

    def test_no_decisions_when_not_requested(self, geno: Path, tmp_path: Path):
        from genotools.core.validation import validate_input
        d = validate_input(geno, tmp_path / "o7")  # nothing requested
        assert (d.skip_sex, d.skip_case_control, d.skip_het, d.disable_filter_controls) \
            == (False, False, False, False)
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `.venv/bin/python -m pytest tests/unit/test_core.py::TestValidation -q`
Expected: FAIL (validate_input returns `None`, has no keyword params / `ValidationDecisions` undefined).

- [ ] **Step 3: Implement** `ValidationDecisions` + skip logic in `genotools/core/validation.py`.

Add `import warnings` and `from dataclasses import dataclass` at the top. Add the dataclass above `validate_input`:

```python
@dataclass(frozen=True)
class ValidationDecisions:
    """Data-driven step-skip decisions derived from the input breakdown."""
    skip_sex: bool = False
    skip_case_control: bool = False
    skip_het: bool = False
    disable_filter_controls: bool = False
```

Change `validate_input`'s signature and append the decision logic after the existing breakdown prints (replace the current `-> None` and the tail of the function):

```python
def validate_input(
    geno_path: Union[str, Path],
    out_path: Union[str, Path],
    skip_fails: bool = False,
    *,
    sex_requested: bool = False,
    het_requested: bool = False,
    hwe_requested: bool = False,
    filter_controls: bool = False,
    case_control_requested: bool = False,
) -> ValidationDecisions:
    ...  # (all existing validation + breakdown-print code unchanged)

    chr_counts = var["#CHROM"].value_counts().to_dict()

    if skip_fails:
        return ValidationDecisions()

    skip_sex = False
    if sex_requested:
        if (1 not in sex_counts) and (2 not in sex_counts):
            warnings.warn(
                "You tried calling sex prune but no sample sex data is available. "
                "Skipping...", stacklevel=2)
            skip_sex = True
        elif ("23" not in chr_counts) and ("X" not in chr_counts):
            warnings.warn(
                "You tried calling sex prune but no X chromosome data is "
                "available. Skipping...", stacklevel=2)
            skip_sex = True

    disable_filter_controls = False
    if hwe_requested and filter_controls and (1 not in pheno_counts):
        warnings.warn(
            "You tried calling hwe prune with controls filtered but no controls "
            "are available. Skipping...", stacklevel=2)
        disable_filter_controls = True

    skip_case_control = False
    if case_control_requested and ((1 not in pheno_counts) or (2 not in pheno_counts)):
        warnings.warn(
            "You tried calling case-control prune but only cases or controls are "
            "available, not both. Skipping...", stacklevel=2)
        skip_case_control = True

    skip_het = False
    if het_requested and (var.shape[0] < 50):
        warnings.warn(
            "You tried calling het prune with less than 50 samples. Skipping...",
            stacklevel=2)
        skip_het = True

    return ValidationDecisions(
        skip_sex=skip_sex,
        skip_case_control=skip_case_control,
        skip_het=skip_het,
        disable_filter_controls=disable_filter_controls,
    )
```

Note: `chr_counts` uses the already-read `var` frame's `#CHROM` column (previously unused). The `pheno_counts`/`sex_counts` dicts already exist from the breakdown code.

- [ ] **Step 4: Run tests to verify they pass**

Run: `.venv/bin/python -m pytest tests/unit/test_core.py::TestValidation -q`
Expected: PASS (all, including the pre-existing raise tests).

- [ ] **Step 5: Commit**

```bash
git add genotools/core/validation.py tests/unit/test_core.py
git commit -m "feat: port data-driven step-skip decisions into validate_input"
```

---

### Task 2: Wire skip decisions into the runner

**Files:**
- Modify: `genotools/cli/runner.py`
- Test: `tests/unit/test_cli/test_runner_regression.py`

**Interfaces:**
- Consumes: `ValidationDecisions`, updated `validate_input` (Task 1).
- Produces: `PipelineRunner._validation_decisions: ValidationDecisions`; `PipelineRunner._filter_steps_by_decisions(steps: List[str]) -> List[str]`.

- [ ] **Step 1: Write failing tests** in `tests/unit/test_cli/test_runner_regression.py` (append at end).

```python
from genotools.core.validation import ValidationDecisions


class TestValidationDecisionsApplied:
    """Round-6: data-driven step-skip decisions from validate_input must be
    applied by the runner (skipped steps dropped; filter_controls forced off)."""

    def test_filter_steps_drops_decided_steps(self, tmp_path: Path) -> None:
        runner, _, _ = _make_runner(tmp_path, warn_only=False)
        runner._validation_decisions = ValidationDecisions(
            skip_sex=True, skip_case_control=True, skip_het=True
        )
        kept = runner._filter_steps_by_decisions(
            ["callrate", "sex", "het", "case_control", "hwe"]
        )
        assert kept == ["callrate", "hwe"]

    def test_filter_steps_noop_by_default(self, tmp_path: Path) -> None:
        runner, _, _ = _make_runner(tmp_path, warn_only=False)
        steps = ["callrate", "sex", "het"]
        assert runner._filter_steps_by_decisions(steps) == steps

    def test_disable_filter_controls_reaches_legacy_args(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        runner, geno, out = _make_runner(tmp_path, warn_only=False)
        # Args request filter_controls=True; the decision must force it False,
        # so the assertion is non-trivial (default is already False).
        runner.args.variant_qc.filter_controls = True
        runner._validation_decisions = ValidationDecisions(disable_filter_controls=True)
        captured: Dict[str, Any] = {}

        def fake_single_step(step, step_input, step_output, legacy_args):
            captured["filter_controls"] = legacy_args["filter_controls"]
            _touch_pfiles(Path(step_output))
            return {"pass": True, "step": step, "metrics": {}, "output": {}}

        monkeypatch.setattr(runner, "_run_single_step", fake_single_step)
        runner._run_qc_pipeline(steps=["hwe"], geno_path=str(geno), out_path=str(out))
        assert captured["filter_controls"] is False
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `.venv/bin/python -m pytest tests/unit/test_cli/test_runner_regression.py::TestValidationDecisionsApplied -q`
Expected: FAIL (`_validation_decisions`/`_filter_steps_by_decisions` do not exist).

- [ ] **Step 3: Implement runner wiring** in `genotools/cli/runner.py`.

(a) In `__init__`, after `self.args = args`, add (import `ValidationDecisions` at top: `from ..core.validation import ValidationDecisions`):

```python
        self._validation_decisions = ValidationDecisions()
```

(b) In `_convert_input_format`, replace the `validate_input(...)` call with:

```python
        from ..core.validation import validate_input

        self._validation_decisions = validate_input(
            self.args.geno_path,
            self.args.out_path,
            skip_fails=self.args.output.skip_fails,
            sex_requested=self.args.sample_qc.run_sex,
            het_requested=self.args.sample_qc.run_het,
            hwe_requested=self.args.variant_qc.run_hwe,
            filter_controls=self.args.variant_qc.filter_controls,
            case_control_requested=self.args.variant_qc.run_case_control,
        )
```

(c) Add the helper method (place near `_run_qc_only`):

```python
    def _filter_steps_by_decisions(self, steps: List[str]) -> List[str]:
        """Drop steps the input breakdown says to skip (ported upfront_check)."""
        d = self._validation_decisions
        result = list(steps)
        if d.skip_sex and "sex" in result:
            result.remove("sex")
        if d.skip_case_control and "case_control" in result:
            result.remove("case_control")
        if d.skip_het and "het" in result:
            result.remove("het")
        return result
```

(d) In `run()`, change `steps = self.args.get_all_enabled_steps()` to:

```python
        steps = self._filter_steps_by_decisions(self.args.get_all_enabled_steps())
```

(e) In `_run_qc_pipeline`, right after `legacy_args = self.args.to_legacy_dict()` (and the existing `het_values` override), add:

```python
        if self._validation_decisions.disable_filter_controls:
            legacy_args["filter_controls"] = False
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `.venv/bin/python -m pytest tests/unit/test_cli/test_runner_regression.py -q`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add genotools/cli/runner.py tests/unit/test_cli/test_runner_regression.py
git commit -m "feat: apply data-driven step-skip decisions in the pipeline runner"
```

---

### Task 3: Golden-convert the preprocessing differential tests

**Files:**
- Modify: `tests/unit/test_ancestry/conftest.py`
- Modify: `tests/unit/test_ancestry/test_preprocessing.py`
- Create: `tests/unit/test_ancestry/golden/ref21_22.{bed,bim,fam}`, `ref21_22_labels.txt`, `geno21_22.{pgen,pvar,psam}`, `raw_ref_train.parquet`, `raw_geno_train.parquet`, `raw_geno_inference.parquet`
- Scratchpad (not committed): a one-time generator script.

**Interfaces:**
- Consumes: `genotools.ancestry.preprocessing.get_raw_files`, `genotools.ancestry.legacy.Ancestry` (generator only, while legacy still exists), `genotools.core.executors.{run_command,get_plink2}`.
- Produces: committed reduced inputs + golden parquet fixtures; golden tests with **zero** legacy imports in committed code.

> **Why reduced (chr21+22):** keeps committed fixtures tiny while exercising the identical code paths (common-SNP extraction, merge, the inference missing-column fill loop — chr22 stays in the ref but is excluded from the inference geno, so the fill still fires). Full-panel new==legacy equivalence was already proven differentially and is preserved by the parity gate.

- [ ] **Step 1: Add reduced-input fixtures** to `tests/unit/test_ancestry/conftest.py`.

```python
GOLDEN = Path(__file__).parent / "golden"


@pytest.fixture(scope="session")
def ref21_22_bfile() -> Path:
    """Committed chr21+22 reference bfile (small; retains chr22)."""
    p = GOLDEN / "ref21_22"
    assert p.with_suffix(".bed").exists(), "run the golden generator first"
    return p


@pytest.fixture(scope="session")
def ref21_22_labels() -> Path:
    return GOLDEN / "ref21_22_labels.txt"


@pytest.fixture(scope="session")
def geno21_22_pfile() -> Path:
    """Committed chr21+22 geno pfile (all 500 samples, small)."""
    p = GOLDEN / "geno21_22"
    assert p.with_suffix(".pgen").exists(), "run the golden generator first"
    return p
```

- [ ] **Step 2: Write the consuming golden tests** — replace the three legacy-oracle tests in `tests/unit/test_ancestry/test_preprocessing.py` (`test_get_raw_files_train_matches_legacy`, `test_get_raw_files_inference_matches_legacy`) with golden versions. Keep `test_get_common_snps_matches_legacy` (imports `genotools.utils`, which stays) and `test_clean_up_files_removes_extensions` untouched. Keep the `_sorted_df` helper.

```python
from pathlib import Path  # (already imported)
from .conftest import GOLDEN


def test_get_raw_files_train_golden(ref21_22_bfile, ref21_22_labels, geno21_22_pfile, tmp_path):
    """Train-mode raw_ref/raw_geno match the committed golden (== legacy at capture)."""
    from genotools.ancestry.preprocessing import get_raw_files
    new_dir = tmp_path / "new"; new_dir.mkdir()
    res = get_raw_files(
        geno_path=str(geno21_22_pfile), ref_panel=str(ref21_22_bfile),
        ref_labels=str(ref21_22_labels), out_path=str(new_dir / "out"), train=True,
    )
    pd.testing.assert_frame_equal(
        _sorted_df(res["raw_ref"]), pd.read_parquet(GOLDEN / "raw_ref_train.parquet"))
    pd.testing.assert_frame_equal(
        _sorted_df(res["raw_geno"]), pd.read_parquet(GOLDEN / "raw_geno_train.parquet"))


def test_get_raw_files_inference_golden(ref21_22_bfile, ref21_22_labels, geno21_22_pfile, tmp_path):
    """Inference-mode raw_geno matches golden; the np.repeat spy proves the
    missing-column fill loop ran (chr22 excluded from the inference geno)."""
    import shutil
    from unittest.mock import patch
    import genotools.ancestry.preprocessing as prep
    from genotools.ancestry.preprocessing import get_raw_files
    from genotools.core.executors import run_command, get_plink2

    # common_snps from a train run on the full (chr21+22) reduced geno
    prep_dir = tmp_path / "prep"; prep_dir.mkdir()
    train_res = get_raw_files(
        geno_path=str(geno21_22_pfile), ref_panel=str(ref21_22_bfile),
        ref_labels=str(ref21_22_labels), out_path=str(prep_dir / "out"), train=True)
    common_snps_src = train_res["out_paths"]["common_snps"]

    # inference geno MISSING chr22 -> ref chr22 common SNPs must be filled
    subset_geno = tmp_path / "subset_geno"
    run_command(
        f"{get_plink2()} --pfile {geno21_22_pfile} --not-chr 22 "
        f"--make-pgen psam-cols=fid,parents,sex,pheno1,phenos --out {subset_geno}",
        tool_name="plink2")

    nd = tmp_path / "new"; nd.mkdir()
    shutil.copy2(common_snps_src, nd / "model.common_snps")
    with patch.object(prep.np, "repeat", wraps=prep.np.repeat) as spy:
        res = get_raw_files(
            geno_path=str(subset_geno), ref_panel=str(ref21_22_bfile),
            ref_labels=str(ref21_22_labels), out_path=str(nd / "out"), train=False,
            common_snps_file=str(nd / "model.common_snps"))
    assert spy.called, "missing-column fill loop was not exercised (degenerate test)"
    pd.testing.assert_frame_equal(
        _sorted_df(res["raw_geno"]), pd.read_parquet(GOLDEN / "raw_geno_inference.parquet"))
```

- [ ] **Step 3: Run the tests to confirm they fail** (golden files absent).

Run: `.venv/bin/python -m pytest tests/unit/test_ancestry/test_preprocessing.py -q`
Expected: FAIL (missing `golden/ref21_22.bed` → fixture assert, or missing parquet).

- [ ] **Step 4: Generate + verify + commit the golden fixtures** with a scratchpad script (run while `legacy.py` still exists). Write to the scratchpad dir, run it, and confirm every `assert_frame_equal(new, legacy)` passes.

```python
# scratchpad/gen_prep_golden.py  — run once: .venv/bin/python scratchpad/gen_prep_golden.py
import shutil
from pathlib import Path
import pandas as pd
from genotools.core.executors import run_command, get_plink2
from genotools.ancestry.preprocessing import get_raw_files
from genotools.ancestry.legacy import Ancestry

REPO = Path("/Users/vitaled2/Desktop/projects/GenoTools")
SYNTH = REPO / "tests/data/synthetic/genotools_test"
GOLDEN = REPO / "tests/unit/test_ancestry/golden"
GOLDEN.mkdir(parents=True, exist_ok=True)
WORK = REPO / "scratchpad/prepwork"; WORK.mkdir(parents=True, exist_ok=True)

def sorted_df(df):
    df = df.sort_values(["FID", "IID"]).reset_index(drop=True)
    return df.reindex(sorted(df.columns), axis=1)

# 1. Build committed reduced inputs (deterministic from SYNTH)
run_command(f"{get_plink2()} --pfile {SYNTH} --chr 21,22 --make-bed --out {GOLDEN/'ref21_22'}", tool_name="plink2")
run_command(f"{get_plink2()} --pfile {SYNTH} --chr 21,22 --make-pgen psam-cols=fid,parents,sex,pheno1,phenos --out {GOLDEN/'geno21_22'}", tool_name="plink2")
fam = pd.read_csv(f"{GOLDEN/'ref21_22'}.fam", sep=r"\s+", header=None)
labels = pd.DataFrame({"FID": fam[0], "IID": fam[1],
                       "label": ["EUR" if i % 2 == 0 else "AFR" for i in range(len(fam))]})
labels.to_csv(GOLDEN / "ref21_22_labels.txt", sep="\t", header=False, index=False)

ref, lbl, geno = str(GOLDEN/"ref21_22"), str(GOLDEN/"ref21_22_labels.txt"), str(GOLDEN/"geno21_22")

# 2. TRAIN: new vs legacy, assert equal, write golden
new = get_raw_files(geno_path=geno, ref_panel=ref, ref_labels=lbl, out_path=str(WORK/"new_out"), train=True)
anc = Ancestry(); anc.geno_path, anc.ref_panel, anc.ref_labels = geno, ref, lbl
anc.out_path = str(WORK/"leg_out"); anc.train = True; anc.model_path = None; anc.containerized = False
leg = anc.get_raw_files()
for key, fname in (("raw_ref", "raw_ref_train.parquet"), ("raw_geno", "raw_geno_train.parquet")):
    g = sorted_df(new[key])
    pd.testing.assert_frame_equal(g, sorted_df(leg[key]))   # new == legacy at capture
    g.to_parquet(GOLDEN / fname)
    pd.testing.assert_frame_equal(g, pd.read_parquet(GOLDEN / fname))  # parquet round-trips exactly

# 3. INFERENCE: build subset geno (no chr22), new vs legacy, write golden
common_src = new["out_paths"]["common_snps"]
run_command(f"{get_plink2()} --pfile {geno} --not-chr 22 --make-pgen psam-cols=fid,parents,sex,pheno1,phenos --out {WORK/'subset'}", tool_name="plink2")
li = WORK/"leg_inf"; li.mkdir(exist_ok=True)
shutil.copy2(common_src, li/"model.common_snps")
anc2 = Ancestry(); anc2.geno_path, anc2.ref_panel, anc2.ref_labels = str(WORK/"subset"), ref, lbl
anc2.out_path = str(li/"out"); anc2.train = False; anc2.model_path = str(li/"model.pkl"); anc2.containerized = False
leg_inf = anc2.get_raw_files()
ni = WORK/"new_inf"; ni.mkdir(exist_ok=True)
shutil.copy2(common_src, ni/"model.common_snps")
new_inf = get_raw_files(geno_path=str(WORK/"subset"), ref_panel=ref, ref_labels=lbl,
                        out_path=str(ni/"out"), train=False, common_snps_file=str(ni/"model.common_snps"))
g = sorted_df(new_inf["raw_geno"])
pd.testing.assert_frame_equal(g, sorted_df(leg_inf["raw_geno"]))
g.to_parquet(GOLDEN / "raw_geno_inference.parquet")
pd.testing.assert_frame_equal(g, pd.read_parquet(GOLDEN / "raw_geno_inference.parquet"))
print("PREP GOLDEN OK")
```

Run: `.venv/bin/python scratchpad/gen_prep_golden.py`
Expected: prints `PREP GOLDEN OK` (every `assert_frame_equal` passed → fixtures equal legacy).

- [ ] **Step 5: Run the golden tests to verify they pass, then commit**

```bash
.venv/bin/python -m pytest tests/unit/test_ancestry/test_preprocessing.py -q   # PASS
ls -lh tests/unit/test_ancestry/golden/                                        # confirm fixtures are small (< ~1 MB total)
git add tests/unit/test_ancestry/conftest.py tests/unit/test_ancestry/test_preprocessing.py tests/unit/test_ancestry/golden
git commit -m "test: golden-convert ancestry preprocessing differential tests"
```

---

### Task 4: Golden-convert the cohort-split differential tests

**Files:**
- Modify: `tests/unit/test_ancestry/test_cohort.py`
- Create: `tests/unit/test_ancestry/golden/split_EUR.{pgen,pvar,psam}`, `split_AFR.{pgen,pvar,psam}`, `cohort_pruned_samples.parquet`
- Scratchpad (not committed): a one-time generator script.

**Interfaces:**
- Consumes: `genotools.ancestry.cohort.split_cohort_by_ancestry`, `tests.regression.compare.compare_genotypes`, `genotools.ancestry.legacy.Ancestry` (generator only), the `geno21_22_pfile` fixture from Task 3.
- Produces: golden split pfiles + pruned-samples fixture; golden tests with zero legacy imports.

> The split is driven by a **deterministic, programmatically built** labels frame (a tiny "EAS" group below `min_samples` exercises the insufficient-N prune; EUR/AFR retained). Reduced geno (chr21+22) keeps the golden split pfiles tiny.

- [ ] **Step 1: Write the consuming golden tests** — replace both legacy-oracle tests in `tests/unit/test_ancestry/test_cohort.py`. Use the shared deterministic label builder so the generator and test agree exactly.

```python
from pathlib import Path
import pandas as pd
from tests.regression.compare import compare_genotypes
from .conftest import GOLDEN


def _write_labels(geno: Path, out: Path):
    """Deterministic 3-group labels: tiny EAS (<min_samples) + EUR/AFR."""
    psam = pd.read_csv(f"{geno}.psam", sep=r"\s+")
    fid = psam.get("#FID", psam.get("FID")); iid = psam["IID"]
    label = ["EAS" if i < 3 else ("EUR" if i % 2 == 0 else "AFR") for i in range(len(iid))]
    df = pd.DataFrame({"FID": fid, "IID": iid, "label": label})
    df.to_csv(out, sep="\t", index=False)


def _sorted(d):
    return d.sort_values(["FID", "IID"]).reset_index(drop=True)


def test_split_cohort_golden(geno21_22_pfile, tmp_path):
    """Retained + insufficient-N pruned branches match the committed golden."""
    from genotools.ancestry.cohort import split_cohort_by_ancestry
    labels_path = tmp_path / "pred_labels.txt"
    _write_labels(geno21_22_pfile, labels_path)
    nd = tmp_path / "new"; nd.mkdir()
    res = split_cohort_by_ancestry(
        labels_path=str(labels_path), geno_path=str(geno21_22_pfile),
        out_path=str(nd / "out"), min_samples=10, subset=None)

    assert sorted(res["labels"]) == ["AFR", "EUR"]
    assert "EAS" not in res["labels"]
    assert set(res["pruned_samples"]["label"]) == {"EAS"}
    pd.testing.assert_frame_equal(
        _sorted(res["pruned_samples"]),
        pd.read_parquet(GOLDEN / "cohort_pruned_samples.parquet"))
    for lbl in ("EUR", "AFR"):
        cmp = compare_genotypes(GOLDEN / f"split_{lbl}", Path(f"{nd/'out'}_{lbl}"))
        assert cmp.equal, cmp.message


def test_split_cohort_subset_golden(geno21_22_pfile, tmp_path):
    """subset=['EUR'] restricts which ancestries are split."""
    from genotools.ancestry.cohort import split_cohort_by_ancestry
    labels_path = tmp_path / "pred_labels.txt"
    psam = pd.read_csv(f"{geno21_22_pfile}.psam", sep=r"\s+")
    fid = psam.get("#FID", psam.get("FID")); iid = psam["IID"]
    df = pd.DataFrame({"FID": fid, "IID": iid,
                       "label": ["EUR" if i % 2 == 0 else "AFR" for i in range(len(iid))]})
    df.to_csv(labels_path, sep="\t", index=False)
    nd = tmp_path / "new"; nd.mkdir()
    res = split_cohort_by_ancestry(
        labels_path=str(labels_path), geno_path=str(geno21_22_pfile),
        out_path=str(nd / "out"), min_samples=10, subset=["EUR"])
    assert res["labels"] == ["EUR"]
```

- [ ] **Step 2: Run the tests to confirm they fail** (golden files absent).

Run: `.venv/bin/python -m pytest tests/unit/test_ancestry/test_cohort.py -q`
Expected: FAIL (missing `golden/split_EUR.*` / `cohort_pruned_samples.parquet`).

- [ ] **Step 3: Generate + verify + commit the golden fixtures** with a scratchpad script (legacy still present). Verify new==legacy for labels, pruned samples, and **genotype content of each split** via `compare_genotypes`.

```python
# scratchpad/gen_cohort_golden.py — run: .venv/bin/python scratchpad/gen_cohort_golden.py
from pathlib import Path
import shutil
import pandas as pd
from genotools.ancestry.cohort import split_cohort_by_ancestry
from genotools.ancestry.legacy import Ancestry
import sys; sys.path.insert(0, str(Path("/Users/vitaled2/Desktop/projects/GenoTools")))
from tests.regression.compare import compare_genotypes

REPO = Path("/Users/vitaled2/Desktop/projects/GenoTools")
GOLDEN = REPO / "tests/unit/test_ancestry/golden"
GENO = GOLDEN / "geno21_22"          # committed in Task 3
WORK = REPO / "scratchpad/cohortwork"; WORK.mkdir(parents=True, exist_ok=True)

psam = pd.read_csv(f"{GENO}.psam", sep=r"\s+")
fid = psam.get("#FID", psam.get("FID")); iid = psam["IID"]
label = ["EAS" if i < 3 else ("EUR" if i % 2 == 0 else "AFR") for i in range(len(iid))]
labels_path = WORK / "pred_labels.txt"
pd.DataFrame({"FID": fid, "IID": iid, "label": label}).to_csv(labels_path, sep="\t", index=False)

# new
nd = WORK / "new"; nd.mkdir(exist_ok=True)
new = split_cohort_by_ancestry(labels_path=str(labels_path), geno_path=str(GENO),
                               out_path=str(nd / "out"), min_samples=10, subset=None)
# legacy
ld = WORK / "leg"; ld.mkdir(exist_ok=True)
anc = Ancestry(); anc.geno_path = str(GENO); anc.out_path = str(ld / "out"); anc.subset = None; anc.min_samples = 10
leg = anc.split_cohort_ancestry(labels_path=str(labels_path))

def s(d): return d.sort_values(["FID", "IID"]).reset_index(drop=True)
assert sorted(new["labels"]) == sorted(leg["labels"]) == ["AFR", "EUR"]
pd.testing.assert_frame_equal(s(new["pruned_samples"]), s(leg["pruned_samples"]))
for lbl in ("EUR", "AFR"):
    cmp = compare_genotypes(Path(f"{ld/'out'}_{lbl}"), Path(f"{nd/'out'}_{lbl}"))
    assert cmp.equal, cmp.message                       # new split == legacy split
    for ext in (".pgen", ".pvar", ".psam"):
        shutil.copy2(f"{nd/'out'}_{lbl}{ext}", GOLDEN / f"split_{lbl}{ext}")
s(new["pruned_samples"]).to_parquet(GOLDEN / "cohort_pruned_samples.parquet")
print("COHORT GOLDEN OK")
```

Run: `.venv/bin/python scratchpad/gen_cohort_golden.py`
Expected: prints `COHORT GOLDEN OK`.

- [ ] **Step 4: Run the golden tests to verify they pass, then commit**

```bash
.venv/bin/python -m pytest tests/unit/test_ancestry/test_cohort.py -q   # PASS
git add tests/unit/test_ancestry/test_cohort.py tests/unit/test_ancestry/golden
git commit -m "test: golden-convert ancestry cohort-split differential tests"
```

---

### Task 5: Flip the default to the new ancestry path

**Files:**
- Modify: `genotools/cli/runner.py`
- Modify: `genotools/cli/__init__.py`
- Modify: `setup.py`
- Modify: `tests/unit/test_cli/test_runner.py` (delete `TestLegacyAncestryImport`)

**Interfaces:**
- Consumes: existing `_run_ancestry_prediction_new`.
- Produces: `run_pipeline(args)` and `PipelineRunner(args)` with no `use_new_ancestry` parameter; `genotools.cli` exports only `main` (no `main_new`).

- [ ] **Step 1: Update the entry-point / import tests first (red)** in `tests/unit/test_cli/test_runner.py`.

Delete the entire `TestLegacyAncestryImport` class (both `test_legacy_ancestry_importable` and `test_entry_points_import`). Add a replacement guard that pins the flip:

```python
class TestDefaultEntryPoint:
    def test_only_main_entry_point(self):
        import genotools.cli as cli
        assert callable(cli.main)
        assert not hasattr(cli, "main_new")

    def test_run_pipeline_has_no_ab_flag(self):
        import inspect
        from genotools.cli import run_pipeline
        assert "use_new_ancestry" not in inspect.signature(run_pipeline).parameters
```

Run: `.venv/bin/python -m pytest tests/unit/test_cli/test_runner.py::TestDefaultEntryPoint -q`
Expected: FAIL (`main_new` still present; `use_new_ancestry` still in signature).

- [ ] **Step 2: Edit `genotools/cli/runner.py`.**

- `PipelineRunner.__init__`: change signature to `def __init__(self, args: PipelineArgs) -> None:`; delete the `use_new_ancestry` docstring arg, `self._use_new_ancestry = use_new_ancestry`, and `self._ancestry: Any = None`.
- `_initialize_modules`: delete the trailing block:

```python
        # Legacy ancestry engine (still the default; A/B control until Phase 5/6)
        from ..ancestry.legacy import Ancestry

        self._ancestry = Ancestry()
```

- Delete the whole legacy `_run_ancestry_prediction` method (the one that configures `self._ancestry` and calls `self._ancestry.run_ancestry()`).
- `_run_with_ancestry`: replace

```python
        if self._use_new_ancestry:
            ancestry_result = self._run_ancestry_prediction_new()
        else:
            ancestry_result = self._run_ancestry_prediction()
```

with

```python
        ancestry_result = self._run_ancestry_prediction_new()
```

- Module-level `run_pipeline`: change to

```python
def run_pipeline(args: PipelineArgs) -> PipelineOutput:
    """Run the GenoTools pipeline.

    This is the main entry point for running the pipeline programmatically.

    Args:
        args: Validated pipeline arguments.

    Returns:
        PipelineOutput containing all results.
    """
    runner = PipelineRunner(args)
    return runner.run()
```

- [ ] **Step 3: Edit `genotools/cli/__init__.py`.** Delete the `main_new()` function and remove `"main_new"` from `__all__` (leave `"main"`).

- [ ] **Step 4: Edit `setup.py`.** Remove the line `'genotools-new=genotools.cli:main_new',` from `console_scripts`.

- [ ] **Step 5: Run the flip tests + the standalone guard.**

Run: `.venv/bin/python -m pytest tests/unit/test_cli/test_runner.py -q`
Expected: PASS (including `TestNewAncestryStandalone`, which still holds).

- [ ] **Step 6: Commit**

```bash
git add genotools/cli/runner.py genotools/cli/__init__.py setup.py tests/unit/test_cli/test_runner.py
git commit -m "refactor: make new AncestryModel the default; drop genotools-new A/B path"
```

---

### Task 6: Delete legacy modules, verify end-to-end, document

**Files:**
- Delete: `genotools/ancestry/legacy.py`, `genotools/gwas.py`, `tests/scripts/test_ancestry_ab.py`
- Modify: `REFACTOR_HARDENING.md` (add round-6 section)

**Interfaces:**
- Consumes: nothing (removals + docs).
- Produces: a legacy-free branch.

- [ ] **Step 1: Delete the dead modules.**

```bash
git rm genotools/ancestry/legacy.py genotools/gwas.py tests/scripts/test_ancestry_ab.py
```

- [ ] **Step 2: Static-check that no live references remain.**

Run:
```bash
grep -rn "genotools-new\|use_new_ancestry\|ancestry\.legacy\|_run_ancestry_prediction\b\|main_new" genotools/ | grep -v "_run_ancestry_prediction_new"
grep -rn "genotools/gwas\.py\|from genotools import gwas\b" genotools/ tests/
```
Expected: no output (the only `gwas` hits are the `genotools.gwas.<submodule>` package imports, which are legitimate).

- [ ] **Step 3: Reinstall (console-scripts changed) and run the full unit + regression suite.**

```bash
.venv/bin/pip install -e . -q
.venv/bin/python -m pytest tests/unit tests/regression -q
```
Expected: green (no collection errors from the deleted modules; `import genotools.ancestry.legacy` appears in no committed test).

- [ ] **Step 4: Run the parity gate.**

Run: `.venv/bin/python -m pytest tests/regression/test_parity.py -v`
Expected: 8/8 pass (`.venv-stable` present).

- [ ] **Step 5: End-to-end default-ancestry smoke run (background, >2 min).** Build a ref panel + labels the way `tests/unit/test_ancestry/conftest.py` does (full synthetic → bfile + alternating labels), then run the **default** command and confirm the NEW model trained.

```bash
# in scratchpad (resolve plink2 via genotools' own resolver — it is cached under ~/.genotools):
PLINK2=$(.venv/bin/python -c "from genotools.core.executors import get_plink2; print(get_plink2())")
$PLINK2 --pfile tests/data/synthetic/genotools_test --make-bed --out scratchpad/e2e_ref
# build labels (FID IID label alternating EUR/AFR) into scratchpad/e2e_labels.txt, then:
.venv/bin/python -m genotools --pfile tests/data/synthetic/genotools_test \
    --out scratchpad/e2e_out --ancestry \
    --ref-panel scratchpad/e2e_ref --ref-labels scratchpad/e2e_labels.txt \
    --min-samples 10 --full-output
```
Expected: exit 0; `scratchpad/e2e_out*.json` present; `scratchpad/e2e_out_ancestry_model/` model dir created (proves training via the new `AncestryModel`); split pfiles per predicted ancestry present. **Run in the background** (trains GridSearchCV + UMAP).

- [ ] **Step 6: Add a round-6 section to `REFACTOR_HARDENING.md`** summarizing: default flipped to new `AncestryModel`; `genotools-new`/`use_new_ancestry` retired; `legacy.py`/`gwas.py`/`test_ancestry_ab.py` deleted; 4 differential tests golden-converted (reduced chr21+22 fixtures); `upfront_check` skip-gap ported into `core/validation.py` (with the preserved variant-count quirk). Note remaining work (#1 real-data parity, #9 parallelize, #10 PCA cache) and the still-deferred `utils.py`/`imputation.py`/`container/`.

- [ ] **Step 7: Commit**

```bash
git add -A
git commit -m "refactor: delete legacy ancestry/gwas modules; document round 6"
```

---

## Self-Review notes

- **Spec coverage:** Task 1+2 = skip-gap fix; Task 3+4 = golden conversion; Task 5 = flip default; Task 6 = deletions + verification + docs. All four spec goals covered.
- **Ordering dependency:** golden generation (Tasks 3–4) and the flip (Task 5) must precede legacy deletion (Task 6). Generators import `genotools.ancestry.legacy`, which exists until Task 6.
- **Legacy in committed code:** after Task 6, no committed test imports `ancestry.legacy` (generators live only in scratchpad). The `TestNewAncestryStandalone` static guard and `genotools.utils`-based tests (`banner`, `get_common_snps`) remain valid — `utils.py` is deferred, not deleted.
- **Type consistency:** `ValidationDecisions` field names (`skip_sex`, `skip_case_control`, `skip_het`, `disable_filter_controls`) are used identically in `validate_input` (Task 1), `_filter_steps_by_decisions` and `_run_qc_pipeline` (Task 2), and the runner tests.
