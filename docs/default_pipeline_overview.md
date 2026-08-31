# Default Pipeline in GenoTools

## Overview
GenoTools' default pipeline commences with ancestry prediction, followed by a series of structured steps for quality control at both the sample and variant levels. This systematic methodology ensures the processing of genotype data is thorough, safeguarding data integrity and relevance for advanced analysis.

## Pipeline Flow

### Ancestry Prediction Pipeline

The `AncestryModel` class (`genotools/ancestry/model.py`) provides a comprehensive pipeline for ancestry predictions in genetic data. Key features include:

- **Data Preparation**: Processes and prunes genotype data, aligns it with reference panels, and prepares labeled datasets for analysis.
- **Principal Component Analysis**: Performs PCA to transform the genetic data for downstream analysis, crucial for capturing genetic variance.
- **Machine Learning Pipeline**: Utilizes a combination of UMAP and XGBoost for sophisticated ancestry prediction, including training, parameter tuning, and evaluation.
- **Local Prediction**: Prediction runs locally, in process. `--container`,
  `--singularity` and `--cloud` were removed in 2.0 and are rejected with an
  error — the container images were built around a 1.x model format that 2.0
  cannot load, and cloud prediction was never implemented in any version. See
  [MIGRATION_2.0.md](../MIGRATION_2.0.md).
- **Admixture Analysis**: Adjusts labels for samples with complex genetic backgrounds using clustering techniques.
- **Cohort Splitting**: Efficiently splits data based on predicted ancestries, facilitating focused analysis on specific genetic groups.

Designed for robustness and versatility, `AncestryModel` is a key component of GenoTools, streamlining complex tasks in ancestry estimation and genetic data analysis.


### QC Steps
Following ancestry analysis, the pipeline undergoes these sample-level and variant-level QC procedures:

### Sample-Level QC Steps
- **Callrate**: Default threshold: 0.05
- **Sex Check**: Default cutoffs: [0.25, 0.75]
- **Relatedness Check**: Enabled by default
- **Heterozygosity (Het)**: Default range: `[-0.15, 0.15]` on PLINK's F
  coefficient. Bounds can instead be derived from each dataset's own spread
  with `--het sd [N]` (`mean ± N·σ` of F, N defaults to 3), and overridden for
  a single ancestry group with `--het-ancestry LABEL ...`. Under `--ancestry`
  the base applies to every group unless a group has an override. Fixed bounds
  cut groups of very different dispersion at very different strictness, which
  is what the `sd` form exists to address — see
  [cli_args.md](cli_args.md#the-heterozygosity-spec-grammar) for the grammar
  and for when to prefer each.

### Variant-Level QC Steps
- **Case-Control Check**: Default threshold: 1e-4
- **Haplotype Check**: Default threshold: 1e-4
- **Hardy-Weinberg Equilibrium (HWE)**: Default threshold: 1e-4
- **Genotype Missingness (Geno)**: Default threshold: 0.05
