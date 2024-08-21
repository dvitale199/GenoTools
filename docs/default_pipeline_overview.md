# Default Pipeline in GenoTools

## Overview
GenoTools' default pipeline commences with ancestry prediction, followed by a series of structured steps for quality control at both the sample and variant levels. This systematic methodology ensures the processing of genotype data is thorough, safeguarding data integrity and relevance for advanced analysis.

## Pipeline Flow

### Ancestry Prediction Pipeline

The `Ancestry` class in GenoTools provides a comprehensive pipeline for ancestry predictions in genetic data. Key features include:

- **Data Preparation**: Processes and prunes genotype data, aligns it with reference panels, and prepares labeled datasets for analysis.
- **Principal Component Analysis**: Performs PCA to transform the genetic data for downstream analysis, crucial for capturing genetic variance.
- **Machine Learning Pipeline**: Utilizes a combination of UMAP and XGBoost for sophisticated ancestry prediction, including training, parameter tuning, and evaluation.
- **Flexible Prediction Environment**: Supports predictions both in local, containerized, and cloud environments, catering to diverse computational needs.
- **Admixture Analysis**: Adjusts labels for samples with complex genetic backgrounds using clustering techniques.
- **Cohort Splitting**: Efficiently splits data based on predicted ancestries, facilitating focused analysis on specific genetic groups.

Designed for robustness and versatility, the `Ancestry` class is a key component of GenoTools, streamlining complex tasks in ancestry estimation and genetic data analysis.


### QC Steps
Following ancestry analysis, the pipeline undergoes these sample-level and variant-level QC procedures:

### Sample-Level QC Steps
- **Callrate**: Default threshold: 0.05
- **Sex Check**: Default cutoffs: [0.25, 0.75]
- **Relatedness Check**: Enabled by default
- **Heterozygosity Rate (Het)**: Default range: [-0.15, 0.15]

### Variant-Level QC Steps
- **Case-Control Check**: Default threshold: 1e-4
- **Haplotype Check**: Default threshold: 1e-4
- **Hardy-Weinberg Equilibrium (HWE)**: Default threshold: 1e-4
- **Genotype Missingness (Geno)**: Default threshold: 0.05
