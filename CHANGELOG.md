# Changelog

Notable changes per release. Releases before 2.1.0 are recorded in the git
history and the GitHub releases page.

## 2.1.0

**The first published 2.x release.** 2.0.0 and 2.0.1 were development versions
and never reached PyPI, so this release supersedes **1.3.6** and carries the
whole 2.x change set. Read [MIGRATION_2.0.md](MIGRATION_2.0.md) before
upgrading.

QC results are unchanged from 1.3.6 and verified so. The CLI spelling, the JSON
report schema, and the ancestry prediction path all changed.

### Upgrade in one paragraph

Flags moved from `underscore_style` to `hyphen-style`; the old spellings still
work but warn. If you only run QC, that is most of what affects you. If you
predict ancestry, three separate fixes can move calls — read the behavior
changes below. If you retrain a model, expect ~1.3% of labels to move for
reasons that are not the fixes; see the note under Known behavior.

### Fixed

- **Ancestry training was nondeterministic and could silently collapse.**
  `learning_rate` was declared, validated and documented but never passed to
  XGBoost, which fell back to a default at the edge of numerical divergence;
  with gblinear's Hogwild updater and an unset thread count, a race decided per
  run whether a fit diverged. A collapsed fit predicted one label for every
  sample while reporting an accuracy equal to that label's prevalence, and was
  saved and used anyway. Measured at 6/20 and 19/20 collapses on dense WGS
  panel PCs. Array cohorts never collapsed but still produced a **different
  model on every run**. Fits are now bit-identical across repeats, verified at
  129,831 samples, and a collapsed fit fails the run instead of being saved.
  `n_estimators` rose 100 → 200 to match the corrected learning rate.
- **Absent model SNPs were filled with dosage 2.** A SNP the cohort does not
  carry has to be invented to match the model's fixed feature width; filling it
  with dosage 2 applied a large shared offset to every sample and could push a
  whole cohort to the `CAH` admixed label. The fill is now a missing value that
  the reference-fitted imputer replaces with the panel mean.
  `--ancestry-missing-fill constant` restores the old behavior exactly.
- **Ancestry SNP matching used an order-dependent tie-break.** The second
  `get_common_snps` call passed its arguments reversed, so which probe won at a
  duplicated position depended on row order. It decided 1,078 of 43,173 common
  SNPs. Palindromic sites are now excluded as well.
- `--amr-het` was silently inert without `--ancestry` — which is exactly the
  path production uses. Replaced by `--het-ancestry` (see Removed).
- `--ancestry` without `--ref-panel`/`--ref-labels` failed inside PLINK with
  `Failed to open None.bed`; it now fails at the CLI boundary with a clear
  message.
- `--container`, `--singularity` and `--cloud` were accepted and silently did
  nothing; they now exit with an error.
- Binary-GWAS summary crash, and `--warn` step-failure handling.

### Changed

- **All flags are hyphenated.** Underscore spellings still work and warn.
- **`umap_learn` is no longer pinned.** The 0.5.3 pin had become an install
  blocker — setuptools removed `pkg_resources` in 82.0.0. Unpinning shifts
  ~1.2% of ancestry calls, and retraining does not avoid it.
- **GWAS p-values shift slightly.** PCA now excludes high-LD and MHC regions
  that 1.x left in. Genomic-inflation lambda is unchanged within 0.05 and the
  tested-variant set is identical.
- **1.x ancestry models cannot be loaded.** Retrain, or stay on 1.3.6.
- Logging redesigned: one consolidated run log with per-step sections and
  inlined PLINK output, plus per-step raw logs that persist regardless of
  `--full-output`.
- Internals: the `SampleQC`/`VariantQC`/`Ancestry` classes were replaced by
  pure functions over frozen config dataclasses. Anything importing those
  classes directly will need updating; the CLI is unaffected.

### Added

- `--het-ancestry`, and a `sd [N]` form for `--het`, so heterozygosity bounds
  can be derived per ancestry group in both flat and per-ancestry runs.
- Ancestry diagnostics: `--ancestry-plots`, `--ancestry-self-test`,
  `--ancestry-missing-fill`, `--ancestry-max-missing-snps`,
  `--no-admixture-detection`, `--ancestry-min-fit-accuracy`,
  `--ancestry-fit-fallbacks`.
- `--quiet`, `--debug`, `--no-warn`, `--no-prune-duplicated`.
- JSON report gains `outcome`/`reason` on every step (distinguishing "not
  requested" from "requested but impossible"), a `parameters` section recording
  the settings a run actually used, an `ancestry_diagnostics` block (SNP
  overlap and fill, panel/cohort allele-frequency concordance, PC drift,
  admixture decisions), an `ancestry_fit` block (selected hyperparameters,
  accuracies, model health, baselines), `common_snps`, and a `software` section
  naming each external tool's resolved path and version.
- `tests/scripts/check_model_health.py` — reports whether a saved model, 1.x or
  2.x, is collapsed.

### Removed

- `--amr-het` / `--amr_het`. Use `--het-ancestry AMR sd`, which works in a flat
  run and on any label vocabulary. Rejected with a message naming the
  replacement.
- `ClassifierConfig.max_depth` and `reg_lambda` — neither was ever passed to
  the classifier.

### Known behavior

- **Retraining moves ~1.3% of ancestry labels, and this is retraining, not the
  fixes above.** Holding training data fixed, old-vs-new moves 0.08%; all
  determinism arms label identically; and 1.3.6 pays the same ~1.3% against the
  same released labels. The causes are a flat hyperparameter plateau (48 of 216
  candidates within one fold-std) and a cohort-dependent common-SNP list.
  Neither labeling is demonstrably more correct. Predicting with an existing
  model is unaffected. See MIGRATION_2.0.md.
- **Full-release ancestry prediction needs a high-memory machine.**
  Preprocessing materializes the cohort as a dense 8-byte matrix and peaks near
  187 GiB at 129,831 samples. Not new in 2.x — 1.x shares the code path.
- `het` pruning fails on very small ancestry groups (observed on a 12-sample
  group). Reported as a failed step; the pipeline continues. Carried over from
  1.x.
