# Migrating to GenoTools 2.0

GenoTools 2.0 is the refactor release. QC and ancestry **results** are unchanged
— see [Verification](#verification) — but the **CLI spelling** and parts of the
**JSON report schema** changed, and two long-standing flag bugs were fixed in
ways that alter behavior.

```bash
pip install --upgrade the_real_genotools     # 2.0.0
pip install 'the_real_genotools<2.0'         # stay on 1.3.6 while you migrate
```

QC pruning behavior and defaults are unchanged. The changes that need attention
are the flag spellings (old ones still work, with a warning), one renamed JSON
column, and slightly shifted GWAS p-values.

---

## Command-line flags

Every flag renamed from `underscore_style` to `hyphen-style`. **The old
spellings still work** and emit a deprecation warning naming the replacement,
so existing scripts keep running. They will be removed in a future release.

| 1.x | 2.0 |
|---|---|
| `--all_sample` | `--all-sample` |
| `--all_variant` | `--all-variant` |
| `--case_control` | `--case-control` |
| `--covar_names` | `--covar-names` |
| `--duplicated_cutoff` | `--duplicated-cutoff` |
| `--filter_controls` | `--filter-controls` |
| `--full_output` | `--full-output` |
| `--kinship_check` | `--kinship-check` |
| `--maf_lambdas` | `--maf-lambdas` |
| `--min_samples` | `--min-samples` |
| `--prune_related` | `--prune-related` |
| `--ref_labels` | `--ref-labels` |
| `--ref_panel` | `--ref-panel` |
| `--related_cutoff` | `--related-cutoff` |
| `--skip_fails` | `--skip-fails` |
| `--subset_ancestry` | `--subset-ancestry` |

To find them in your own scripts:

```bash
grep -rnE '\-\-(all|case|covar|duplicated|filter|full|kinship|maf|min|prune|ref|related|skip|subset)_' .
```

### Boolean flags no longer accept a value

In 1.x every boolean flag was declared `type=str, nargs='?', const='True'`, so
`--all_sample False` parsed and quietly meant "off". That was never the intent —
these are presence flags, and the pattern was copied from the threshold flags
(`--callrate` and friends) where an optional value *is* wanted. 2.0 declares
them `action="store_true"`, so a value is an error:

```
$ genotools ... --all_sample False
error: --all_sample takes no value (got 'False'). In 1.x every boolean flag
accepted True/False; they are now presence-only. pass the flag on its own to
enable it, or omit it to disable it.
```

To migrate: **pass the flag alone to enable, omit it to disable.** Affects
`--all_sample`, `--all_variant`, `--ancestry`,
`--filter_controls`, `--full_output`, `--gwas`, `--kinship_check`,
`--maf_lambdas`, `--prune_related`, `--related`, `--skip_fails`, `--warn`,
`--prune_duplicated`. (`--container` and `--singularity` are also presence-only
now, but 2.0 rejects them outright — see below.)

Two of those defaulted to *on* in 1.x, so "disable" needs an explicit flag
rather than omission:

| 1.x | 2.0 |
|---|---|
| `--warn` / omitted (continue past errors) | default; pass `--no-warn` to stop on the first error |
| `--prune_duplicated` / omitted (prune duplicates) | default; pass `--no-prune-duplicated` to disable |

`--warn` and `--prune_duplicated` are still accepted on their own as no-ops,
since what they requested is now the default.

### `--amr-het` was removed, replaced by `--het-ancestry`

Not a rename — a different flag with different reach. Both spellings
(`--amr-het`, `--amr_het`) are rejected with a message naming the replacement:

```
$ genotools ... --all-sample --amr-het
ERROR: --amr-het was removed in GenoTools 2.0.1. Use '--het-ancestry AMR sd'
instead, which also works in a flat run (--amr-het silently did nothing
without --ancestry). See MIGRATION_2.0.md.
```

| 1.x / 2.0.0 | 2.0.1 |
|---|---|
| `--ancestry --all-sample --amr_het` | `--ancestry --all-sample --het-ancestry AMR sd` |

Three reasons it changed rather than being carried forward:

- **It was silently inert outside `--ancestry`.** The flag was read in exactly
  one place, inside the ancestry branch; a flat run never consulted it. The
  per-ancestry production workflow runs ancestry once and then QCs each group
  as a separate flat job — precisely the path where it did nothing, with no
  warning and a normal-looking JSON. `--het-ancestry` works in both run shapes,
  and errors rather than being ignored when it cannot apply.
- **It hardcoded one reference panel's label.** `label == "AMR"` was the only
  user-facing feature in the codebase that assumed a particular panel's naming.
  The label vocabulary is user-supplied — from `--ref-labels`, or from a
  pickled model's encoder — so the flag was unusable on any panel spelling the
  group differently, and could never reach `CAH`, the synthetic admixed label
  that admixture detection invents.
- **The multiplier was invisible.** `--amr-het` was described as "auto-detect",
  but only the location and scale were derived from the data; the `3` was
  hardcoded. `sd [N]` makes it a knob.

**The replacement is not bit-for-bit identical in principle, though it was in
practice here.** `--amr-het` thresholded the derived heterozygosity *rate*;
`sd` thresholds `F`, so that both spellings of `--het` bound the same
statistic. The two are near-perfectly anti-correlated within a group, so the
`mean ± 3σ` rule picks the same samples either way: measured on the GP2 r12 10k
subset (9,771 samples post-callrate) the two rules select **identical** sets —
99 samples each, zero disagreements cohort-wide and zero within AMR's 340 —
despite bounds on quite different scales (`F: [-0.081, 0.160]` against
`rate: [0.258, 0.332]`). Nothing guarantees that on other data, so treat a
borderline sample moving as possible rather than expected.

`--het` itself gained the same spec grammar (`--het sd [N]` alongside
`LOWER UPPER`), and the 1.x `--het -1 -1` sentinel still works while warning
that `--het sd` is the spelling now. See
[docs/cli_args.md](docs/cli_args.md) for the full grammar.

### New flags

`--het-ancestry`, `--quiet`, `--debug`, `--no-warn`, `--no-prune-duplicated`.
`--het` accepts a new `sd [N]` form.

---

## Behavior changes

### Relatedness pruning is unchanged

Worth stating explicitly because it is easy to misread the flag rename as a
behavior change. It is not. Under `--all-sample`, both versions **report**
related pairs and prune only duplicates:

| | 1.x | 2.0 |
|---|---|---|
| `prune_related` default | `False` | `False` |
| `prune_duplicated` default | `True` | `True` |

1.x declared both as `type=str` (`'False'`/`'True'`) and then mapped them onto
real booleans in `__main__.py`, so the documented defaults were the effective
ones. Related pairs land in `{out}_{ancestry}.related` and the
`related_samples` JSON block either way; they are labeled with a degree
(`unrel`/`second_deg`/`first_deg`/`duplicate`), not removed. Pass
`--prune-related` to remove them, in either version.

### Steps that could not run are now reported

In 1.x a step that was requested but impossible was handled two different ways,
neither of them visible:

- **Skipped** (the data ruled it out — sex prune with no X chromosome, het on
  too few samples): the step was silently dropped. No row in `QC`, no entry in
  `pass_fail`. Indistinguishable from a step you never asked for.
- **Failed** (the step ran and errored): a `QC` row appeared with counts zeroed
  and `pass: false`, but the reason was never written to the report.

2.0 reports both, with a reason. See
[JSON report schema](#json-report-schema) for the new `outcome` and `reason`
fields. Nothing about *which* steps run changed — only whether you can tell
what happened to them.

Three checks are now re-decided **per dataset** rather than inherited from the
cohort, because `--ancestry` splits the cohort *after* the cohort-level check and
a group's data can differ sharply from the whole:

| step | why a group can differ | 1.x | 2.0 |
|---|---|---|---|
| `het` | group falls under PLINK's 50-sample LD floor | failed inside PLINK | skipped |
| `sex` | group has no recorded sample sex | failed inside PLINK | skipped |
| `case_control` | group holds only cases, or only controls | raised in the step | skipped |

In 1.x these were reported as failures for a group even though the cohort-level
version of the same finding was reported as a skip — one decision with two
behaviors depending on where it was noticed. Only sample-derived checks are
re-decided: the X-chromosome half of the sex check cannot change, since the split
keeps samples and every group inherits the cohort's pvar. All three read the psam
as it stands *before* the QC chain runs, so a precondition broken by an earlier
prune (callrate removing the last control, say) still fails in the step.

`--skip-fails` suppresses the cohort-level decisions but not these, so a cohort
that 1.x would have run-and-failed is now skipped. Samples and variants carried
forward are unaffected either way — the step did not run in either version.

One skip is new, because the guard behind it never worked. 1.x meant to skip
het below 50 samples but tested the **variant** count instead of the sample
count (`utils.py:185`), and a real dataset never has fewer than 50 variants, so
the guard never fired. Het then ran and died inside PLINK, which refuses to
estimate LD from fewer than 50 samples. 2.0 tests the sample count, and — since
`--ancestry` splits the cohort *after* that check — re-tests it per ancestry
group, where small groups actually occur.

If you run `--ancestry --het` on a cohort with a small ancestry group, that
group's het step changes from a failure to a skip:

```
1.x:  het_prune  FIN  outlier_count  0  pass=false
2.0:  het_prune  FIN  outlier_count  0  pass=false  outcome="skipped"  reason="12 samples is fewer than..."
```

Samples and variants carried forward are identical — het did not run in either
version. Only the reporting changed. Anything reading `pass` still works; read
`outcome` to tell a skip from a failure.

### `--container`, `--singularity` and `--cloud` are rejected

2.0 does not run ancestry prediction remotely. Passing any of the three is an
error rather than a silent no-op:

```
$ genotools ... --ancestry --container
ERROR: --container is not supported in GenoTools 2.0. 1.x ran prediction in a
Docker image built around a 1.x model, which 2.0 cannot load. Drop the flag to
predict locally, or pin 'genotools<2.0' to keep the 1.x container.
```

**`--container` / `--singularity`** worked in 1.x. They wrote the projected PCs
to `genotools/container/`, ran `mkoretsky1/genotools_ancestry:python3.11`, and
read predicted labels back. That image's `run.py` unpickles
`GP2_merge_release6_NOVEMBER_..._umap_linearsvc_ancestry_model.pkl` — a 1.x
model, which 2.0's `AncestryModel` cannot load. Restoring the flags therefore
needs a rebuilt and republished image carrying a 2.0-format model; until that
exists, the flags fail rather than quietly predicting locally.

**`--cloud` never did anything, in any version.** It is not a lost 1.x feature:
1.3.6 has no `--cloud` flag and no cloud code path at all. 2.0 added the flag
name without an implementation behind it.

To migrate: drop the flag — prediction runs locally, which is what 2.0 would
have done anyway. If you need the 1.x container, pin `genotools<2.0`.

### 1.x ancestry models cannot be loaded

`--model` takes a model directory, or a single `.pkl` written by 2.0. A model
from 1.x holds an `sklearn.pipeline.Pipeline` rather than an `AncestryModel` and
is rejected:

```
Invalid model file: expected AncestryModel, got <class 'sklearn.pipeline.Pipeline'>.
This looks like a GenoTools 1.x model, which 2.0 cannot load. Pass a model
directory written by 2.0, or retrain by dropping --model and passing
--ref-panel/--ref-labels.
```

Retrain once against your reference panel and reuse the resulting directory.
(1.x's `--model` also required a sibling `.common_snps` file; 2.0 keeps
`common_snps.txt` inside the model directory instead.)

### A model records the libraries it was fitted under

`metadata.json` now carries a `versions` block — umap-learn, scikit-learn,
xgboost, numpy, pandas, scipy, GenoTools and Python — captured at fit time.
Loading a model compares them against the environment and warns on any
difference:

```
Model version drift: umap-learn 0.5.3 -> 0.5.7. This model was fitted under the
recorded versions, and the embedding can differ under different ones, so
ancestry calls may not match what this model was validated on. Reinstall the
recorded versions, or retrain, to reproduce them.
```

This is a **warning, never an error** — the load always succeeds. It exists
because the failure it describes is otherwise silent: a model fitted under one
umap and loaded under another unpickles cleanly and embeds differently, so the
run finishes with no error and different ancestry calls.

A model trained before this block existed loads with a *provenance unknown*
warning instead, since "cannot tell" and "no drift" are different answers.
Retrain to record it.

### GWAS p-values shift slightly

PCA now prunes high-LD and MHC regions that 1.x left in, so association
p-values differ marginally across the board. This is intentional (ratified as
"decision B" in `REFACTOR.md`). Genomic-inflation lambda is unchanged within
0.05 and the tested-variant set is identical; do not expect bit-identical
p-values against 1.x output.

---

## JSON report schema

Most of the report is unchanged. The differences:

| Key | Change |
|---|---|
| `QC[].pruned_count` | unchanged (a pre-release refactor renamed it to `count`; reverted) |
| `QC[].step` | unchanged — a pre-release refactor shortened it to the flag name (`callrate`) on runs without `--ancestry`; reverted to the reported name (`callrate_prune`) |
| `QC[].outcome`, `QC[].reason` | **new** — `pass` / `fail` / `skipped`, and why. See [Steps that could not run are now reported](#steps-that-could-not-run-are-now-reported) |
| `pass_fail[].outcome`, `pass_fail[].reason` | **new** — same two fields on the per-step status block. `status` still a boolean |
| `pruned_samples[].label` | unchanged — still carries the sample's ancestry group |
| `total_umap`, `ref_umap`, `new_samples_umap` | **columns renamed** `"0".."24"` → `"UMAP1".."UMAP25"` |
| `common_snps` | **new** — count of variants shared between input and reference panel |
| `projected_pcs` | same columns and values; column *order* differs (`label` moved earlier) |

The UMAP rename is the one breaking change here. If you read those columns by
name:

```python
# 1.x
umap = pd.DataFrame(report["total_umap"])[["0", "1"]]

# 2.0
umap = pd.DataFrame(report["total_umap"])[["UMAP1", "UMAP2"]]
```

`label` and `dataset` columns in those blocks are unchanged. Numeric positional
keys were replaced because `"0"` is ambiguous as a JSON key and sorts
unintuitively; the values themselves are identical.

---

## Verification

2.0 was validated against 1.3.6 on a 10,000-sample subset of GP2 release 12,
run as `--ancestry --all-sample --all-variant`:

- **All 11 ancestry groups produced byte-identical genotypes and sample/variant
  IDs** (verified with `plink2 --pgen-diff` plus allele-coding comparison)
- Per-sample ancestry labels identical for all 10,000 samples; identical
  per-label counts, test accuracy (0.985037), and confusion matrix
- QC metrics, pruned-sample IDs, related pairs, and per-step pass/fail identical
- Both versions reproduce the released full-cohort labels at the same 98.68%,
  with identical disagreement patterns

Reproduce with `tests/scripts/compare_ancestry_run.py`.

Known issue carried over from 1.x: `het` pruning fails on very small ancestry
groups (observed on a 12-sample FIN group in both versions). It is reported as a
failed step and the pipeline continues.
