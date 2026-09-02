# GenoTools Command Line Arguments Documentation

## Overview

Detailed reference for the `genotools` command-line arguments, covering input,
output, quality control (QC), ancestry prediction, and GWAS.

Flags use **hyphens**, not underscores (`--full-output`, not `--full_output`).
The 1.x underscore spellings are still accepted for backward compatibility and
emit a deprecation warning; they will be removed in a future release. See
[MIGRATION_2.0.md](../MIGRATION_2.0.md) for the full mapping.

`genotools --help` is the authoritative list; this document covers the same flags
in more detail and is kept in sync with `genotools/cli/parser.py`.

**Threshold flags are optional-argument flags.** `--callrate`, `--geno`,
`--case-control`, `--haplotype`, `--hwe`, and `--pca` each run their step with a
default when passed bare, and accept a value to override it: `--callrate` uses
0.02, `--callrate 0.05` uses 0.05. Omitting the flag entirely skips the step.

---

### File I/O Arguments

Exactly one input format is required. If more than one is given, `--pfile` wins
over `--bfile`, which wins over `--vcf`.

- **`--bfile`**
  - *Type*: `path`
  - *Default*: `None`
  - *Description*: PLINK 1.9 format genotype file path, without the
    `.bed`/`.bim`/`.fam` extension. Converted to PLINK 2 pfiles at pipeline start.

- **`--pfile`**
  - *Type*: `path`
  - *Default*: `None`
  - *Description*: PLINK 2 format genotype file path, without the
    `.pgen`/`.pvar`/`.psam` extension. The pipeline's native format.

- **`--vcf`**
  - *Type*: `path`
  - *Default*: `None`
  - *Description*: VCF format genotype file path (`.vcf` or `.vcf.gz`). Converted
    to PLINK 2 pfiles at pipeline start.

- **`--out`** *(required)*
  - *Type*: `path`
  - *Default*: `None`
  - *Description*: Prefix for output files, including the path.

---

### Output and Logging Arguments

- **`--full-output`**
  - *Type*: `flag`
  - *Default*: `False`
  - *Description*: Keep every intermediate pfile instead of only the final
    output. Does **not** affect logs — the consolidated and per-step logs are
    always written.

- **`--skip-fails`**
  - *Type*: `flag`
  - *Default*: `False`
  - *Description*: Bypass the up-front checks. Two effects: (1) re-running over an
    already-used output prefix is allowed — the previous run's consolidated log is
    rotated to `{out}_all_logs.log.N` rather than overwritten (per-step logs are
    replaced); and (2) the data-driven step auto-skips described below are
    **disabled**, so a step runs even when the input can't support it.

- **`--no-warn`**
  - *Type*: `flag`
  - *Default*: `False` (warn-and-continue is the default behavior)
  - *Description*: Stop the pipeline at the first failing step. Without this
    flag, a failing step is recorded as failed and the pipeline continues from
    the last step that passed.

- **`--quiet`**
  - *Type*: `flag`
  - *Default*: `False`
  - *Description*: Suppress the console progress stream. Log files are still
    written in full: the consolidated `{out}_all_logs.log` and the per-step
    `{out}_{step}.log` raw PLINK logs. Useful for batch/cluster jobs.

- **`--debug`**
  - *Type*: `flag`
  - *Default*: `False`
  - *Description*: Log at `DEBUG` level instead of `INFO`, on both the console
    and the consolidated log, and re-raise errors so the full Python traceback is
    shown (without it, expected failures print a one-line `ERROR:` message).
    Combine with `--quiet` for verbose file-only logs.

---

### Data-Driven Step Auto-Skips

A requested step is **skipped automatically** when the input can't support it.
`--skip-fails` disables these checks, forcing the step to run anyway.

| Step | Skipped when |
|------|--------------|
| `--sex` | No sample sex data (no males and no females recorded), or no X chromosome in the input |
| `--het` | The input has fewer than 50 samples (PLINK's LD floor) |
| `--case-control` | Cases or controls are missing (both are required) |
| `--filter-controls` | No controls present — the HWE filter runs on everyone instead |

Each skip is logged as a `!` warning on the console and recorded in the
consolidated run log.

---

### Sample-Level QC Arguments

- **`--callrate [THRESHOLD]`**
  - *Type*: `float`
  - *Default*: `0.02` when the flag is passed bare; step skipped if omitted
  - *Description*: Maximum per-sample missingness (PLINK `--mind`). Samples above
    the threshold are pruned.

- **`--sex [FEMALE_MAX_F MALE_MIN_F]`**
  - *Type*: `float` (zero or two values)
  - *Default*: `0.25 0.75` when the flag is passed bare; step skipped if omitted
  - *Description*: Sex-check F-statistic cutoffs. Samples whose genetic sex is
    discordant with the recorded sex are pruned. Requires 0 or exactly 2 values.
    Subject to auto-skip (see above).

- **`--het [LOWER UPPER | sd [N]]`**
  - *Type*: `str` (see the spec grammar below)
  - *Default*: `-0.15 0.15` when the flag is passed bare; step skipped if omitted
  - *Description*: Heterozygosity bounds on PLINK's F coefficient, applied to
    every dataset. Samples outside the bounds are pruned. Subject to auto-skip
    (see above).

- **`--het-ancestry LABEL [LOWER UPPER | sd [N]]`**
  - *Type*: `str`, repeatable (one ancestry group per use)
  - *Default*: none; `--het`'s value applies to every group
  - *Description*: Override `--het` for one ancestry group. Requires
    `--ancestry` — a run without ancestry prediction has no labels to match, so
    passing it in a flat run is an error. Using it implies the het step runs,
    even without `--het` or `--all-sample`. A label that no group in the run
    carries produces a warning naming the groups that were predicted.

#### The heterozygosity spec grammar

Both flags take the same four-row spec, so they cannot mean different things:

| Tokens | Meaning |
|---|---|
| *(none)* | fixed `-0.15 0.15` (the defaults) |
| `LOWER UPPER` | fixed bounds on PLINK's `F` |
| `sd` | bounds at `mean ± 3σ` of this dataset's `F` |
| `sd N` | bounds at `mean ± N·σ` |

`sd` bounds `F`, the same statistic the fixed bounds do — only the location and
scale come from the data. A per-group override always beats the base `--het`.

| Command | Effect |
|---|---|
| `--het` | every group: fixed `-0.15 0.15` |
| `--het -0.2 0.2` | every group: fixed `-0.2 0.2` |
| `--het sd` | every group / the whole input: `mean ± 3σ` |
| `--het sd 2.5` | every group: `mean ± 2.5σ` |
| `--het -0.2 0.2 --het-ancestry AMR sd` | AMR adaptive, everything else `-0.2 0.2` |
| `--het sd 2 --het-ancestry AMR sd 1.5` | 2σ everywhere, AMR tightened to 1.5σ |
| `--het-ancestry AMR sd --het-ancestry CAH -0.3 0.3` | adaptive and custom ranges, mixed per group |

#### When to reach for `sd`

Fixed bounds are one absolute cut applied to groups whose dispersion differs by
an order of magnitude. Measured on GP2 r12 (10k subset, 9,759 samples across 10
groups), `sd(F)` is 0.0669 for AMR against 0.0088 for EUR, so `±0.15` is a 2.24σ
cut for AMR and a 17σ cut for EUR. In practice `--het -0.15 0.15` is an AMR-only
filter: 31 of the 35 samples it prunes cohort-wide are AMR.

`sd` expresses the cut in multiples of a group's own spread instead. Three
things follow that are worth knowing before using it:

- **`sd` and fixed are not interchangeable defaults.** On the same data a global
  `sd 3` prunes 137 samples where fixed prunes 35, and EUR alone goes from 2 to
  85. That is why the per-group override exists — `sd` is not the global default.
- **`sd N` is a dispersion multiplier, not an expected-flag-rate dial.** The
  "3σ ≈ 0.3%" intuition assumes normality; heterozygosity within an admixed
  group is often skewed and multimodal by admixture proportion, which is the
  premise of the feature.
- **`sd` is not robust, deliberately.** Mean and σ are computed over every
  sample, outliers included, so a real subpopulation inflates σ, widens the
  bounds and spares itself. For AMR at `sd 3` the bounds work out to
  `-0.177 .. +0.225` — wider than the observed range on both sides, so het
  pruning is effectively a no-op for that group. That is the intended outcome
  (see below), but it means het pruning is **not** a contamination screen for a
  widely-dispersed group. On small groups, σ estimated near the 50-sample floor
  makes the boundary itself wobble in a way a fixed bound does not.

#### Why AMR, and why `sd 3` prunes nothing there

The long-standing rationale for special-casing AMR was that "heterozygosity is
higher in admixed populations." On GP2 r12 that is not what the data shows: AMR
has the *lowest* mean het rate and the *highest* mean F of all ten groups, and
30 of the 31 samples fixed bounds prune are on the heterozygote-*deficit* side
(`F >= +0.15`), 28 of them from a single contributing cohort.

AMR is a pooled label spanning heterogeneous admixture proportions; pooling
those allele frequencies produces an apparent heterozygote deficit (Wahlund),
pushing a structured subgroup past the fixed **upper** bound. Those samples sit
at only ~2.5σ of AMR's own spread. So `sd 3` pruning zero AMR samples is the
correct outcome rather than under-pruning — `sd 2` would prune 21 of them, which
is to say it would discard the very samples the feature exists to keep. The
per-group multiplier exists so base and group *can* disagree, not because AMR
should be tightened.

- **`--related`**
  - *Type*: `flag`
  - *Default*: `False`
  - *Description*: Run relatedness analysis (KING). By default it reports and
    prunes duplicates but retains related samples — see `--prune-related` and
    `--no-prune-duplicated`.

- **`--related-cutoff`**
  - *Type*: `float`
  - *Default*: `0.0884`
  - *Description*: Kinship coefficient above which a pair counts as related.

- **`--duplicated-cutoff`**
  - *Type*: `float`
  - *Default*: `0.354`
  - *Description*: Kinship coefficient above which a pair counts as duplicated.

- **`--prune-related`**
  - *Type*: `flag`
  - *Default*: `False`
  - *Description*: Also prune related (not just duplicated) samples.

- **`--no-prune-duplicated`**
  - *Type*: `flag`
  - *Default*: `False` (duplicates are pruned by default)
  - *Description*: Report duplicated samples without pruning them.

- **`--kinship-check`**
  - *Type*: `flag`
  - *Default*: `False`
  - *Description*: Confirm reported family relationships against genotype-derived
    kinship. **Linux only** (requires KING).

- **`--all-sample`**
  - *Type*: `flag`
  - *Default*: `False`
  - *Description*: Run callrate, sex, het, and related with default thresholds.
    Note that the callrate threshold under this flag is **0.05**, not the 0.02
    used by a bare `--callrate`. Does not include `--kinship-check`.
    Passing `--het` alongside it still applies — `--all-sample` turns the step
    on, `--het` supplies the bounds — so `--all-sample --het sd` runs the full
    sample battery with derived heterozygosity bounds.

---

### Variant-Level QC Arguments

- **`--geno [THRESHOLD]`**
  - *Type*: `float`
  - *Default*: `0.05` when the flag is passed bare; step skipped if omitted
  - *Description*: Maximum per-variant missingness. Variants above the threshold
    are pruned.

- **`--case-control [P_THRESHOLD]`**
  - *Type*: `float`
  - *Default*: `1e-4` when the flag is passed bare; step skipped if omitted
  - *Description*: P-value threshold for the case/control missingness
    differential test. Subject to auto-skip (see above).

- **`--haplotype [P_THRESHOLD]`**
  - *Type*: `float`
  - *Default*: `1e-4` when the flag is passed bare; step skipped if omitted
  - *Description*: P-value threshold for the missingness-by-haplotype test.

- **`--hwe [P_THRESHOLD]`**
  - *Type*: `float`
  - *Default*: `1e-4` when the flag is passed bare; step skipped if omitted
  - *Description*: Hardy-Weinberg equilibrium p-value threshold.

- **`--filter-controls`**
  - *Type*: `flag`
  - *Default*: `False`
  - *Description*: Apply the HWE filter using controls only. Subject to
    auto-disable (see above).

- **`--ld [WINDOW_SIZE STEP_SIZE R2_THRESHOLD]`**
  - *Type*: `float` (zero or three values)
  - *Default*: `50 5 0.5` when the flag is passed bare; step skipped if omitted
  - *Description*: Linkage-disequilibrium pruning parameters. Requires 0 or
    exactly 3 values.

- **`--all-variant`**
  - *Type*: `flag`
  - *Default*: `False`
  - *Description*: Run every variant-level QC step with default thresholds
    (geno, case-control, haplotype, hwe) — **excluding** LD pruning, which must
    be requested explicitly with `--ld`.

---

### Ancestry Arguments

- **`--ancestry`**
  - *Type*: `flag`
  - *Default*: `False`
  - *Description*: Run ancestry prediction and split the cohort by predicted
    ancestry, then run the requested QC steps once per ancestry group. Requires
    either `--ref-panel` + `--ref-labels` (train a new model) or `--model` (use a
    pre-trained one).

- **`--ref-panel`**
  - *Type*: `path`
  - *Default*: `None`
  - *Description*: Reference panel genotype file path (PLINK format, without
    extension). **Required whenever `--ancestry` is used** — training builds
    the reference PCA from it, and inference (`--model`) re-derives that PCA by
    subsetting the panel to the model's common SNPs.

- **`--ref-labels`**
  - *Type*: `path`
  - *Default*: `None`
  - *Description*: Reference panel ancestry labels file: tab-separated
    `FID IID label`, no header. **Required whenever `--ancestry` is used**,
    inference included.

- **`--model`**
  - *Type*: `path`
  - *Default*: `None`
  - *Description*: Pre-trained ancestry model directory (or a legacy `.pkl`
    file). Runs inference instead of training. Does **not** remove the need for
    `--ref-panel` and `--ref-labels`.

- **`--container`**, **`--singularity`**, **`--cloud`**
  - *Type*: `flag`
  - *Default*: `False`
  - *Description*: **Not supported — rejected with an error.** Ancestry
    prediction always runs locally, in process. The container images were built
    around a 1.x model format that 2.0 cannot load, and cloud prediction was
    never implemented in any version. See
    [MIGRATION_2.0.md](../MIGRATION_2.0.md).

- **`--subset-ancestry [LABEL ...]`**
  - *Type*: `str` (zero or more labels)
  - *Default*: `None` (all predicted groups)
  - *Description*: Restrict downstream QC to the named ancestry labels, e.g.
    `--subset-ancestry EUR AFR`.

- **`--min-samples`**
  - *Type*: `int`
  - *Default*: `0`
  - *Description*: Minimum samples an ancestry group needs to be carried
    forward. Groups below the threshold get no split pfiles at all; their samples
    are recorded in the pruned-samples output with the step
    `insufficient_ancestry_sample_n`.

#### Prediction-path diagnostics

Training accuracy is measured on reference samples that went through the
*training* preprocessing branch. A cohort goes through the *inference* branch,
which aligns it to the model's SNP list — and can fail there in ways training
never sees. The classic symptom is a model with high test accuracy predicting a
single label (often `CAH`) for an entire cohort, because the cohort matched
almost none of the model's SNPs and the rest were invented.

Every ancestry run now measures that path and reports it in the log and in the
JSON report's `ancestry_diagnostics` section: how much of the model's SNP list
the cohort carried, whether matched sites agree with the panel on allele
frequency, where the cohort landed in reference PC space, and what the `CAH`
override did to the classifier's labels. The flags below tune or extend it.
Each one requires `--ancestry`; passing one without it is an error rather than a
silent no-op.

- **`--ancestry-missing-fill`**
  - *Type*: `str` (`ref-mean` or `constant`)
  - *Default*: `ref-mean`
  - *Description*: How to fill a model SNP the cohort does not carry.
    `ref-mean` writes a missing value, which the reference-fitted imputer fills
    with the panel's own mean for that SNP, so the filled SNP contributes
    nothing to the projection. `constant` writes dosage 2 — one end of the
    range, applied identically to every sample — which is what 1.x and 2.0
    before this did, and is the mechanism behind whole cohorts collapsing onto
    one label. Kept only for reproducing an older run.

- **`--ancestry-max-missing-snps`**
  - *Type*: `float` (fraction in `[0, 1]`)
  - *Default*: `0.5`
  - *Description*: Refuse to predict when more than this share of the model's
    SNPs is absent from the cohort. Above roughly half, a predicted label
    describes the fill rather than the sample. The error names the counts and
    the usual causes (chromosome naming, genome build); raise the limit to
    predict anyway. Below the limit, exceeding 5% logs a warning.

- **`--no-admixture-detection`**
  - *Type*: `flag`
  - *Default*: `False` (detection on)
  - *Description*: Report the classifier's own labels instead of relabeling as
    `CAH` the samples that sit nearer the global centroid than any ancestry
    centroid. The first thing to try on an all-`CAH` result: it separates a
    broken classifier from a broken override. The per-sample
    `{out}_decisions.txt` file records both labels either way.

- **`--ancestry-plots`**
  - *Type*: `flag`
  - *Default*: `False`
  - *Description*: Write two diagnostic PNGs: `{out}_pca_diagnostic.png`, the
    reference panel coloured by label with the cohort overlaid across the three
    leading PC pairs, and `{out}_centroid_distance.png`, the margin behind each
    `CAH` decision against the line where the decision flips.

- **`--ancestry-self-test`**
  - *Type*: `flag`
  - *Default*: `False`
  - *Description*: Predict the reference panel's *own* ancestries through the
    prediction path, as a positive control. Accuracy is optimistic by
    construction (most of these samples were fitted on), so the load-bearing
    number is the panel's `CAH` rate: the panel defines the ancestry centroids,
    so a panel that comes back largely `CAH` indicts the override or the
    projection, not the cohort. A panel that predicts itself correctly here
    leaves only the cohort's own preprocessing to blame.

**Files an ancestry run writes for debugging** (alongside the report, at the
final output prefix — they survive a run without `--full-output`):

| File | Contents |
|---|---|
| `{out}_decisions.txt` | Per sample: the classifier's label, the label after the `CAH` override, the nearest centroid, and one distance column per centroid |
| `{out}_filled_snps.txt` | IDs of the model SNPs the cohort did not carry |
| `{out}_pca_diagnostic.png` | With `--ancestry-plots` |
| `{out}_centroid_distance.png` | With `--ancestry-plots` |

---

### GWAS and PCA Arguments

- **`--pca [N_PCS]`**
  - *Type*: `int`
  - *Default*: `10` when the flag is passed bare; PCA skipped if omitted
  - *Description*: Run PCA with N principal components. High-LD and MHC regions
    are excluded before pruning.

- **`--build`**
  - *Type*: `str` (`hg19` or `hg38`)
  - *Default*: `hg38`
  - *Description*: Genome build, which selects the PCA exclusion-region file.

- **`--gwas`**
  - *Type*: `flag`
  - *Default*: `False`
  - *Description*: Run the association test. Uses the PCA eigenvectors as
    covariates unless `--covars` is given.

- **`--covars`**
  - *Type*: `path`
  - *Default*: `None`
  - *Description*: External covariates file. When supplied, it replaces the PCA
    eigenvectors as the covariate set.

- **`--covar-names`**
  - *Type*: `str`
  - *Default*: `None`
  - *Description*: Comma-separated covariate column names to use from the
    external covariates file.

- **`--maf-lambdas`**
  - *Type*: `flag`
  - *Default*: `False`
  - *Description*: Apply a MAF prune before computing genomic-inflation lambdas.

---

### Examples

```bash
# All sample-level QC with defaults
genotools --pfile data/mydata --out results/output --all-sample

# Ancestry prediction only
genotools --pfile data/mydata --out results/output --ancestry \
          --ref-panel refs/panel --ref-labels refs/labels.txt

# Full pipeline, split by ancestry
genotools --pfile data/mydata --out results/output --ancestry \
          --ref-panel refs/panel --ref-labels refs/labels.txt \
          --all-sample --all-variant

# Custom thresholds, PCA + GWAS, quiet console
genotools --pfile data/mydata --out results/output \
          --callrate 0.05 --geno 0.01 --hwe 1e-6 \
          --pca 20 --gwas --quiet

# Debugging a prediction that came back one label for every sample
genotools --pfile data/mydata --out results/debug --ancestry \
          --model models/my_model --ref-panel refs/panel \
          --ref-labels refs/labels.txt \
          --ancestry-self-test --ancestry-plots --no-admixture-detection
```
