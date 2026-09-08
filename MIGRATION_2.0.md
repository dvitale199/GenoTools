# Migrating to GenoTools 2.x

**2.1.0 is the first published 2.x release.** 2.0.0 and 2.0.1 were development
versions and never reached PyPI, so if you are upgrading you are coming from
**1.3.6**, and this whole document applies to you. Everything below describes
the difference between 1.3.6 and 2.1.0.

2.x is the refactor release. QC **results** are unchanged — see
[Verification](#verification) — but the **CLI spelling** and parts of the
**JSON report schema** changed, and several long-standing bugs were fixed in
ways that alter behavior.

The QC path is the safe part: same pruning, same defaults, same counts. **The
ancestry path is where the real changes are.** Four independent fixes land here
— the training race, the dosage-2 fill, palindromic SNP matching, and the umap
unpin — and three of them can move ancestry calls. If you predict ancestry,
read [Behavior changes](#behavior-changes) before upgrading; if you only run
QC, the flag spellings are most of what matters to you.

```bash
pip install --upgrade the_real_genotools     # 2.1.0
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
ERROR: --amr-het was removed in GenoTools 2.0. Use '--het-ancestry AMR sd'
instead, which also works in a flat run (--amr-het silently did nothing
without --ancestry). See MIGRATION_2.0.md.
```

| 1.x | 2.x |
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

`--het-ancestry`, `--quiet`, `--debug`, `--no-warn`, `--no-prune-duplicated`,
`--ancestry-min-fit-accuracy`, `--ancestry-fit-fallbacks`.
`--het` accepts a new `sd [N]` form.

Ancestry prediction adds `--no-admixture-detection`,
`--ancestry-missing-fill`, `--ancestry-max-missing-snps`, `--ancestry-plots`
and `--ancestry-self-test`. See
[Absent model SNPs are no longer filled with dosage 2](#absent-model-snps-are-no-longer-filled-with-dosage-2)
and `docs/cli_args.md`.

---

## Behavior changes

### Ancestry training is now deterministic, and picks different hyperparameters

`ClassifierConfig.learning_rate` was declared, documented and validated but
never passed to `XGBClassifier`, in 1.x and in 2.0 alike. XGBoost therefore
used its own gblinear default of 0.5, which sits at the edge of numerical
divergence — and gblinear's default `updater="shotgun"` is Hogwild, so with the
thread count unset a race decided per run which side of that edge a fit landed
on. A diverged fit reaches `|intercept| ~1e15`, saturates the softmax, and
predicts one label for every sample while reporting an accuracy exactly equal
to that label's prevalence in the reference panel. It was pickled and used
anyway.

Measured over 20 identical repeats: on dense long-read WGS panel PCs, 6/20 and
19/20 collapses on two draws of the *same* configuration. On GP2 array PCs,
0/20 — array cohorts sit inside the stable region, which is why this went
unnoticed. Having *more* overlapping variants than usual is what caused the
failure.

What changes for you:

- **Retraining is now reproducible.** With `n_jobs=1` and `learning_rate=0.1`
  wired through, repeated fits are bit-identical.
- **The search will select different hyperparameters.** Where fits collapsed
  at random, each grid point was scored partly by luck, so `best_params_` was
  whichever candidate drew the luckiest folds. With the race gone, model
  selection is meaningful — and its outcome may differ from what a 1.x run
  chose. Expect equal or better quality, not identical labels.
- **`n_estimators` rose from 100 to 200.** `learning_rate=0.1` takes five times
  smaller steps, and 100 rounds no longer reach the same place: measured on GP2
  panel PCs, held-out balanced accuracy was 0.9570 at (0.5, 100), 0.9150 at
  (0.1, 100), and 0.9560 at (0.1, 200).
- **A collapsed fit now fails the run instead of being saved.** Tune with
  `--ancestry-min-fit-accuracy` and `--ancestry-fit-fallbacks`; see
  [docs/cli_args.md](docs/cli_args.md).
- **The report gained an `ancestry_fit` block** carrying the selected
  hyperparameters, balanced accuracies, the fitted model's numerical health,
  cheap k-NN and nearest-centroid baselines, and the attempt table.

**Audit models you already have.** The defect is present in every 1.x and
pre-2.1.0 model, so a saved model may be collapsed:

```bash
python tests/scripts/check_model_health.py <model dir or 1.x .pkl>
```

It reports `|coef|`, `|intercept|` and the distinct-class count, and exits
non-zero if any model checked has diverged. A healthy GP2 model measures
`|coef| 1.03 / |intercept| 4.24`; a collapsed one measured `5.32 / 2.0e15`.

It also says which side of the fix a model was trained on: every pre-2.1.0
model pickled `learning_rate=None`, because the field was never passed. A
converged pre-fix model is still evidence about luck rather than about the
process, so retrain when you next can.

### Retraining moves about 1.3% of ancestry labels — and that is retraining, not the fix

If you retrain an ancestry model after upgrading and compare its calls against
labels a 1.x model produced, expect roughly **1.3% of samples to change group**.
Measured on the full GP2 release 12 (129,831 samples) against the released 1.x
labels: **1,695 samples moved, 1.306%**.

**This is not caused by the determinism fix, and it is not evidence that the
new labels are better or worse.** Three independent measurements say so:

- Holding the training data fixed and changing only the code, old-vs-new moves
  **0.08%** — not 1.3%.
- All three determinism arms (as-shipped, fixed, and the coordinate-descent
  alternative) label the cohort **identically**.
- **1.3.6 pays the same ~1.3%** when retrained against those released labels.
  The old code has the property too.

What actually moves the labels is **retraining itself**. Two compounding causes,
both measured:

- **The hyperparameter grid's top is a plateau.** 48 of 216 candidates score
  within one fold-to-fold standard deviation of the winner, and five tie to ten
  decimal places. At full scale the search won by 0.000667 over the candidate a
  10,000-sample run picked, while its own fold-to-fold std was 0.014178 — **21×
  the gap**. Selection among near-equals is decided by noise. The substantive
  parameters did reproduce (`n_components=25`, `n_neighbors=5`,
  `lambda=0.001`); only the UMAP `a`/`b` shape pair moved.
- **The common-SNP list depends on the cohort.** `--geno 0.1` and the
  lowest-missingness tie-break both read per-variant missingness, so a different
  cohort yields a slightly different SNP list. One variant differed between a
  10,000-sample run and the full release — enough to reshuffle an argmax on a
  plateau that flat.

**Neither labeling is demonstrably more correct.** A model-free nearest-centroid
arbiter cannot rank them: 31.8% of the moved samples match the released label,
33.0% the new one, 35.2% neither — because at cohort scale that baseline agrees
with *either* labeling only ~87% of the time. Do not report the upgrade as
having improved ancestry calls; there is no measurement supporting that.

**What to do about it.** If you need labels stable across releases, do not
retrain — keep predicting with your existing model, which is unaffected. If you
do retrain, treat it as a new model with its own labels rather than an update to
the old ones, and expect ~1% of a large cohort to move at every retrain,
including retrains that change nothing but the cohort.

### `ClassifierConfig` lost two fields

`max_depth` and `reg_lambda` were never passed to `XGBClassifier` and are gone.
`max_depth` does nothing under the gblinear booster by its own docstring, and
the grid's `xgb__lambda` travels through `**kwargs` and never touched
`reg_lambda`, so reading it invited tuning a value that changed nothing. Only
code constructing `ClassifierConfig` directly is affected; no CLI flag exposed
either.

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

Two ways forward. For the NeuroBooster array, download the 2.x model:

```bash
genotools-download --model nba_gp2_r12
```

`nba_gp2_r12` is trained on GP2 release 12 (43,173 SNPs, 10 ancestry labels)
and is the `genotools-download` default. It was fitted under the library
versions 2.1.0 itself requires, so it loads without a drift warning; when those
libraries later move, the warning it then emits is the real thing and worth
reading. The 1.x models (`nba_v1`, `nba_v2`,
`neurochip_v1`) remain available for anyone still on 1.x, and asking for one
now prints a warning saying it will not load in 2.x.

Otherwise retrain once against your reference panel and reuse the resulting
directory. (1.x's `--model` also required a sibling `.common_snps` file; 2.x
keeps `common_snps.txt` inside the model directory instead.)

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

**GenoTools' own version is held apart from that check.** It moves on every
release, including one that changes nothing numerical, so a GenoTools-only
difference is reported as provenance at `INFO` rather than as drift:

```
Model provenance: fitted under GenoTools 2.0.1, running 2.1.0. The libraries
that determine the embedding are unchanged, so this is not library drift. Check
that release's changelog if it changed ancestry behaviour.
```

It is not treated as harmless — a release *can* change ancestry behaviour, and
2.x's SNP tie-break and absent-SNP fill both did — but the changelog is the
authority on whether a given one did, and a warning that fires for every model
after every release is one users learn to skip. When a library moved too, the
drift warning fires as normal and names the GenoTools move alongside it.

A model trained before this block existed loads with a *provenance unknown*
warning instead, since "cannot tell" and "no drift" are different answers.
Retrain to record it.

### umap_learn is no longer pinned — ancestry calls shift by ~1.2%

**This changes ancestry calls. Retrain your models and revalidate.**

1.x and early 2.0 pinned `umap_learn==0.5.3`. That pin is now removed
(`umap_learn>=0.5.5`), and the `setuptools` runtime dependency with it.

**Why it had to go.** umap-learn below 0.5.5 runs `import pkg_resources` at
import time, and setuptools deleted `pkg_resources` in 82.0.0. A fresh install
that resolves a current setuptools therefore cannot import umap at all, so
`--ancestry` fails before any of your data is touched. The pin had stopped
protecting reproducibility and started preventing installation.

**What it costs.** Measured on a 10,000-sample GP2 subset across 11 ancestry
groups, comparing umap-learn 0.5.3 / numpy 2.3.5 / pandas 2.3.3 against
umap-learn 0.5.12 / numpy 2.4.6 / pandas 3.0.5:

| Comparison | Calls changed |
|---|---|
| Same model, both environments (inference drift only) | 129 / 10,000 — 1.29% |
| Each environment trains its own model (the upgrade path) | 122 / 10,000 — 1.22% |

Model quality is unchanged: test balanced accuracy 0.9850 before, 0.9838 after.

**Retraining does not avoid this.** 113 of the moved samples are the same in
both comparisons, and the largest single shift — 33 samples from AFR to AAC —
is the *identical 33 samples* either way. The drift is a systematic property of
the newer library stack at population boundaries, not an artifact of a stale
fit, so retraining under the new stack reproduces most of it.

Where the calls move (retrain comparison, groups of 5+):

```
AFR -> AAC  33      MDE -> AFR  11      EUR -> AMR   7
EUR -> AFR  17      AJ  -> EUR   9      AJ  -> MDE   6
EUR -> MDE  13      CAS -> SAS   9
```

These are adjacent and admixed groups, which is where a slightly different
embedding would be expected to tip samples across a boundary — but 1.2% is a
real change to your results, not rounding.

**Reproducing prior results.** Nothing already produced is altered by
upgrading. A model trained by 2.1.0 or later ships a `requirements.txt` beside
`pipeline.pkl` recording the exact environment it was fitted under; recreate it
with `pip install -r <model_dir>/requirements.txt` to get the original calls
back. For models predating that, `requirements-lock.txt` in the repo root
pins a validated environment. Development installs stay unpinned so CI keeps
catching upstream breakage early.

### Ancestry SNP matching excludes palindromes and prefers better-called probes

**This can change ancestry calls, by very little.** Measured at 1 call in
10,000 on a 10,000-sample GP2 subset; model accuracy and every QC count were
unchanged.

Two corrections to how a cohort is matched against the reference panel:

- **Palindromic (A/T, C/G) sites are excluded.** Their alleles survive strand
  complement unchanged, so nothing downstream can tell which strand they came
  from, and such a site was previously accepted at whatever orientation it
  arrived in. If you build panels with the recipe in
  `docs/prep_reference_panel.md` they were already excluded and nothing changes
  for you — the GP2 panel has 0 of 209,517. If you use a panel that kept them,
  this is a correctness fix and your calls may move more than the figure above.

- **A position offering several probes now resolves to the best-called one.**
  Cohorts routinely carry the same site under multiple probe IDs
  (`rs301801`, `IlmnSeq_rs301801`, `seq_rs301801`). Previously whichever the
  merge emitted first won, which was arbitrary and sensitive to the order of
  variants in your input file; now the one with the lowest missingness wins.

**Existing models are not invalidated.** Only the cohort side of the match
changed; the reference-panel side is byte-identical, so a saved model's
common-SNP list, its fitted classifier and its accuracy are all unaffected. You
do not need to retrain, and predictions from an existing model shift only by
the amount above.

### Absent model SNPs are no longer filled with dosage 2

**This can change ancestry calls a lot, but only for a cohort that was already
being predicted badly.**

When predicting with a saved model, the cohort's feature matrix has to have
exactly the model's SNPs in the model's order, so any SNP the cohort does not
carry has to be invented. 1.x — and 2.0 up to this point — filled every such
SNP with dosage **2** for every sample. Dosage 2 is not a neutral value: it is
one end of the range, applied identically to everybody, so a cohort missing a
large share of the model's SNP list acquires a large *shared* offset in PC
space, lands somewhere no reference sample sits, and can end up nearer the
global centroid than to any ancestry centroid — which is the rule that labels a
sample `CAH`. A high-accuracy model predicting `CAH` for an entire cohort is
this, not a model failure.

The fill is now a missing value, which `PCAReducer`'s reference-fitted imputer
replaces with the panel's own mean for that SNP. That is already how a
genuinely missing call at a SNP the cohort *does* carry is treated, so the
change also makes the matrix internally consistent. A filled SNP now
contributes nothing to the projection instead of pulling it.

Three things follow:

- **A cohort that matched the model's SNPs well barely moves.** The filled
  columns were a small part of the matrix, and they now sit at the panel mean
  instead of at 2.
- **A cohort that matched badly may move a great deal** — that is the fix, not
  a regression. The old labels for such a cohort were describing the fill.
- **Above 50% filled, prediction is refused** rather than reported, naming the
  counts and the likely cause. Raise `--ancestry-max-missing-snps` to predict
  anyway; every run logs the overlap and writes the filled IDs to
  `{out}_filled_snps.txt` either way.

`--ancestry-missing-fill constant` restores the old behaviour exactly, int
dosages included, for reproducing an older run.

**Existing models are not invalidated** — the model is untouched; this is
entirely about the matrix handed to it at prediction time.

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
| `ancestry_diagnostics` | **new** — what the prediction path measured about itself: `snp_overlap` (how much of the model's SNP list the cohort carried, and how much was filled), `allele_frequency` (panel/cohort concordance over matched sites, with a swap-signature count), `pc_drift` (per-PC cohort displacement in reference SDs), `admixture` (the `CAH` count with the classifier's labels from before the override), `self_test` (with `--ancestry-self-test`), `plots`, and `warnings` |
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

2.x was validated against 1.3.6 twice: a full differential comparison on a
10,000-sample subset, and a production-scale run on the whole release.

**Differential parity, 10,000-sample subset of GP2 release 12**, run as
`--ancestry --all-sample --all-variant`:

- **All 11 ancestry groups produced byte-identical genotypes and sample/variant
  IDs** (verified with `plink2 --pgen-diff` plus allele-coding comparison)
- Per-sample ancestry labels identical for all 10,000 samples; identical
  per-label counts, test accuracy (0.985037), and confusion matrix
- QC metrics, pruned-sample IDs, related pairs, and per-step pass/fail identical
- Both versions reproduce the released full-cohort labels at the same 98.68%,
  with identical disagreement patterns

Reproduce with `tests/scripts/compare_ancestry_run.py`.

**Production scale, the full GP2 release 12** — 129,831 samples, 4h50m, run
with the production configuration:

- **All 88 QC steps passed across all 11 ancestry groups**
- The trained model is healthy and converged (max abs coefficient 0.866,
  intercept 0.676, all 10 labels predicted), with none of the collapse
  signature the fix removes
- **Bit-identical across 20 repeated fits** — the determinism claim holds at
  full scale, not just on the small panels it was found on
- QC reproduces the release almost exactly: callrate 2,906 vs 2,906 with
  perfect per-sample agreement, sex 1,910 vs 1,909. The one real difference,
  heterozygosity, is a *threshold* difference and not a code difference — the
  release gave AMR `sd`-derived bounds where the comparison run used the base
  fixed window, and rebuilding that group reproduces the release's counts
  exactly under the release's own setting
- Labels moved 1.306% against the released 1.x labels, for the reasons in
  [Retraining moves about 1.3% of ancestry labels](#retraining-moves-about-13-of-ancestry-labels--and-that-is-retraining-not-the-fix)

**Memory note.** Predicting a cohort this size currently needs a high-memory
machine: the preprocessing step materializes the cohort as a dense 8-byte
matrix and peaks near 187 GiB at 129,831 samples. This is not new in 2.x — 1.x
has the identical code path — but plan hardware accordingly for full-release
prediction. Tracked as REFACTOR.md item 40.

Known issue carried over from 1.x: `het` pruning fails on very small ancestry
groups (observed on a 12-sample FIN group in both versions). It is reported as a
failed step and the pipeline continues.
