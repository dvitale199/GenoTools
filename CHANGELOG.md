# Changelog

Notable changes per release. Releases before 2.1.0 are recorded in the git
history and the GitHub releases page.

## Unreleased — next release is 2.2.0

A **minor** version, not a patch: the change below alters what a default run
leaves on disk. Nothing about your results, final outputs or logs changes, but
a default that deletes files users previously got to keep is more than a patch.

### Changed

- **Intermediate pfiles are now deleted as the run progresses, at defaults.**

  Each QC step reads a genotype file and writes a new one — callrate reads your
  input and writes `…_callrate`, sex reads that and writes `…_sex`, and so on —
  so a ten-step run creates ten copies of your data. At GP2 release scale each
  copy is tens of GiB.

  Cleanup used to be skipped whenever `warn_only` was set, and warn-and-continue
  is the **default**, so at defaults every copy survived for the whole run. The
  only way to reclaim the disk was `--no-warn`, which also turns
  warn-and-continue into fail-fast — so you could not have both. That is a bad
  trade on a long run, where warn-and-continue exists precisely so an eight-hour
  job does not die on one recoverable step.

  The two flags no longer interfere:

  | | before | now |
  |---|---|---|
  | *(default)* | keeps every copy | **deletes each copy once the next is written** |
  | `--full-output` | keeps every copy | keeps every copy |
  | `--no-warn` | deletes + fail-fast | fail-fast only |
  | `--no-warn --full-output` | keeps + fail-fast | keeps + fail-fast |

  `--full-output` is now the only switch governing retention; `--no-warn` is
  purely failure policy. Peak disk during a run drops from roughly one working
  copy per step to one overall — on the round-19 full-GP2 run, ~98 GiB written
  against a ~50 GiB peak.

  **If you relied on the old default to inspect intermediates, pass
  `--full-output`** — it reproduces the previous behavior exactly. A **failed**
  step's input is still kept either way, so a broken run can still be examined.

## 2.1.1

A packaging and diagnostics patch. No QC, ancestry or report behavior changes.

### Fixed

- **`python_requires` was `>=3.8`, three minor versions below what the stack
  needs.** `pandas` and `scikit-learn` both declare `requires_python >=3.11`,
  so the floor admitted interpreters the dependencies do not support. On
  3.8-3.10 pip does not report an unsupported interpreter: it backtracks,
  resolving years-old pandas and sklearn that satisfy their own floors, and
  installs a combination nothing tests. The floor is now `>=3.11`,
  matching what CI tests. `README.md` advertised 3.8/3.9/3.10 badges and none
  for 3.11; it now shows 3.11+.

### Added

- **`genotools --version`.** Previously the only ways to ask which GenoTools was
  installed were `pip show the_real_genotools` and `importlib.metadata`. The
  obvious `genotools.__version__` is misleading: from any directory holding a
  `genotools/` package it reports the source tree rather than the installed
  distribution, and 1.x never defined the attribute, so on a real 1.3.6 install
  it raises `AttributeError` instead of answering. The flag reads the installed
  distribution's metadata, and falls back to the package version — labeled as
  coming from a source tree — when nothing is installed.

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
- **1.x ancestry models cannot be loaded.** Download `nba_gp2_r12` (the new
  default) for the NeuroBooster array, retrain against your own reference
  panel, or stay on 1.3.6.
- Logging redesigned: one consolidated run log with per-step sections and
  inlined PLINK output, plus per-step raw logs that persist regardless of
  `--full-output`.
- **A model's recorded GenoTools version no longer counts as library drift.**
  It changes on every release, so comparing it meant every distributed model
  warned that its ancestry calls might have moved after every release —
  including the model GenoTools ships, on first use. A GenoTools-only
  difference is now reported as provenance at `INFO` and points at the
  changelog; drift in the libraries that determine the embedding warns exactly
  as before, and names the GenoTools move alongside it.
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
- **`nba_gp2_r12`, a 2.x-format ancestry model** for the NeuroBooster array,
  trained on GP2 release 12 (43,173 SNPs, 10 labels) through the fixed
  deterministic path, and fitted under the library versions this release
  requires — so it loads without a drift warning rather than one on first use. Now the `genotools-download` default, since the three 1.x
  models it used to serve cannot be loaded by 2.x. Asking for one of those now
  warns instead of failing later at load time, and an unknown name lists what
  is available rather than raising `KeyError`. `--model default` /
  `--ref default` now work, as the help text has always claimed. A cached
  archive that fails checksum validation is now re-downloaded instead of
  failing forever: `download_data_from_gcs` returned early whenever the
  destination existed, so re-publishing an artifact stranded anyone holding the
  previous copy, and the error never named the file to delete.

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
