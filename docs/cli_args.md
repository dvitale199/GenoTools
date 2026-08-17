# GenoTools Command Line Arguments Documentation

## Overview

Detailed reference for the `genotools` command-line arguments, covering input,
output, quality control (QC), ancestry prediction, and GWAS.

Flags use **hyphens**, not underscores (`--full-output`, not `--full_output`).
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
| `--het` | The input has fewer than 50 variants |
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

- **`--het [LOWER UPPER]`**
  - *Type*: `float` (zero or two values)
  - *Default*: `-0.15 0.15` when the flag is passed bare; step skipped if omitted
  - *Description*: Heterozygosity rate (F coefficient) bounds. Samples outside
    the range are pruned. Requires 0 or exactly 2 values. Subject to auto-skip
    (see above).

- **`--amr-het`**
  - *Type*: `flag`
  - *Default*: `False`
  - *Description*: Use auto-detected heterozygosity bounds for the AMR ancestry
    group instead of the fixed `--het` values. Only meaningful with `--ancestry`.

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
    extension), used to train the ancestry model.

- **`--ref-labels`**
  - *Type*: `path`
  - *Default*: `None`
  - *Description*: Reference panel ancestry labels file: tab-separated
    `FID IID label`, no header.

- **`--model`**
  - *Type*: `path`
  - *Default*: `None`
  - *Description*: Pre-trained ancestry model directory (or a legacy `.pkl`
    file). Runs inference instead of training.

- **`--container`**
  - *Type*: `flag`
  - *Default*: `False`
  - *Description*: Run ancestry prediction inside a Docker container. Mutually
    exclusive with `--model`.

- **`--singularity`**
  - *Type*: `flag`
  - *Default*: `False`
  - *Description*: Run ancestry prediction inside a Singularity container.

- **`--cloud`**
  - *Type*: `flag`
  - *Default*: `False`
  - *Description*: Run ancestry prediction on Google Cloud AI Platform.

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
```
