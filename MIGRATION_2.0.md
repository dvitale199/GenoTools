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
| `--amr_het` | `--amr-het` |
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
grep -rnE '\-\-(all|amr|case|covar|duplicated|filter|full|kinship|maf|min|prune|ref|related|skip|subset)_' .
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
`--all_sample`, `--all_variant`, `--amr_het`, `--ancestry`, `--container`,
`--filter_controls`, `--full_output`, `--gwas`, `--kinship_check`,
`--maf_lambdas`, `--prune_related`, `--related`, `--singularity`,
`--skip_fails`, `--warn`, `--prune_duplicated`.

Two of those defaulted to *on* in 1.x, so "disable" needs an explicit flag
rather than omission:

| 1.x | 2.0 |
|---|---|
| `--warn` / omitted (continue past errors) | default; pass `--no-warn` to stop on the first error |
| `--prune_duplicated` / omitted (prune duplicates) | default; pass `--no-prune-duplicated` to disable |

`--warn` and `--prune_duplicated` are still accepted on their own as no-ops,
since what they requested is now the default.

### New flags

`--quiet`, `--debug`, `--no-warn`, `--no-prune-duplicated`, `--cloud`.

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
