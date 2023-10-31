# SampleQC Class Functions Guide

## Introduction

This guide details the methods available in the `SampleQC` class for quality control on genotype data.

## Table of Contents
- [Constructor](#constructor)
- [run_callrate_prune](#run_callrate_prune)
- [run_sex_prune](#run_sex_prune)
- [run_het_prune](#run_het_prune)
- [run_related_prune](#run_related_prune)

---

## run_callrate_prune

\`\`\`python
run_callrate_prune(self, mind=0.02)
\`\`\`

### Description

Executes call rate pruning on genotype data using PLINK. Removes individuals based on missing genotype rate.

### Parameters

- **mind**: Threshold for permissible missing genotype rate. Defaults to 0.02 (2%).

### Returns

- **dict**: Contains
    * 'pass': Successful completion flag.
    * 'step': Label for this procedure ('callrate_prune').
    * 'metrics': Metrics such as 'outlier_count'.
    * 'output': Paths to the generated output files.

---

## run_sex_prune

\`\`\`python
run_sex_prune(self, check_sex=[0.25,0.75])
\`\`\`

### Description

Executes sex-based pruning on genotype data using PLINK.

### Parameters

- **check_sex**: Lower and upper bounds for sex discrepancy checks. Defaults to [0.25, 0.75].

### Returns

- **dict**: Contains
    * 'pass': Successful completion flag.
    * 'step': Label for this procedure ('sex_prune').
    * 'metrics': Metrics such as 'outlier_count'.
    * 'output': Paths to the generated output files.

---

## run_het_prune

\`\`\`python
run_het_prune(self, het_threshold=0.1)
\`\`\`

### Description

Executes heterozygosity-based pruning on genotype data using PLINK.

### Parameters

- **het_threshold**: Threshold for acceptable heterozygosity rate. Defaults to 0.1 (10%).

### Returns

- **dict**: Contains
    * 'pass': Successful completion flag.
    * 'step': Label for this procedure ('het_prune').
    * 'metrics': Metrics such as 'outlier_count'.
    * 'output': Paths to the generated output files.

---

## run_related_prune

\`\`\`python
run_related_prune(self, pi_hat=0.125)
\`\`\`

### Description

Executes relatedness pruning on genotype data using PLINK.

### Parameters

- **pi_hat**: Threshold for removing related individuals based on PI_HAT value. Defaults to 0.125.

### Returns

- **dict**: Contains
    * 'pass': Successful completion flag.
    * 'step': Label for this procedure ('related_prune').
    * 'metrics': Metrics such as 'outlier_count'.
    * 'output': Paths to the generated output files.


# VariantQC Class Functions Guide

## Introduction

This guide details the methods available in the `VariantQC` class for quality control on variant data.

## Table of Contents
- [run_geno_prune](#run_geno_prune)
- [run_case_control_prune](#run_case_control_prune)
- [run_haplotype_prune](#run_haplotype_prune)
- [run_hwe_prune](#run_hwe_prune)
- [run_ld_prune](#run_ld_prune)

---

## run_geno_prune

\`\`\`python
run_geno_prune(self, geno_threshold=0.05)
\`\`\`

### Description

Prunes SNPs based on missing genotype data.

### Parameters

- **geno_threshold**: Maximum allowable proportion of missing genotyping for a SNP. Defaults to 0.05 (5%).

### Returns

- **dict**: Contains
    * 'pass': Successful completion flag.
    * 'step': Label for this procedure ('geno_prune').
    * 'metrics': Metrics such as 'geno_removed_count'.
    * 'output': Paths to the generated output files.

---

## run_case_control_prune

\`\`\`python
run_case_control_prune(self, p_threshold=1e-4)
\`\`\`

### Description

Prunes SNPs based on missing genotype data, stratified by case-control status.

### Parameters

- **p_threshold**: P-value threshold for testing the difference in missing rates between cases and controls. Defaults to 1e-4.

### Returns

- **dict**: Contains
    * 'pass': Successful completion flag.
    * 'step': Label for this procedure ('case_control_missingness_prune').
    * 'metrics': Metrics such as 'mis_removed_count'.
    * 'output': Paths to the generated output files.

## run_haplotype_prune

\`\`\`python
run_haplotype_prune(self, p_threshold=1e-4)
\`\`\`

### Description

Prunes SNPs based on missing haplotype data.

### Parameters

- **p_threshold**: P-value threshold for missingness by haplotype. Defaults to 1e-4.

### Returns

- **dict**: Contains
    * 'pass': Successful completion flag.
    * 'step': Label for this procedure ('haplotype_prune').
    * 'metrics': Metrics such as 'haplotype_removed_count'.
    * 'output': Paths to the generated output files.

## run_hwe_prune

\`\`\`python
run_hwe_prune(self, hwe_threshold=1e-4, filter_controls=False)
\`\`\`

### Description

Prunes SNPs based on Hardy-Weinberg equilibrium (HWE) p-values.

### Parameters

- **hwe_threshold**: P-value threshold for SNP filtering based on HWE. Defaults to 1e-4.
- **filter_controls**: If set to True, perform HWE test only on control samples. Defaults to False.

### Returns

- **dict**: Contains
    * 'pass': Successful completion flag.
    * 'step': Label for this procedure ('hwe_prune').
    * 'metrics': Metrics such as 'hwe_removed_count'.
    * 'output': Paths to the generated output files.

## run_ld_prune

\`\`\`python
run_ld_prune(self, window_size=50, step_size=5, r2_threshold=0.5)
\`\`\`

### Description

Prunes SNPs based on Linkage Disequilibrium.

### Parameters

- **window_size (int, optional)**: Size of the sliding window over the genome in variant count. Defaults to 50.
- **step_size (int, optional)**: Step of the sliding window over the genome. Defaults to 5.
- **r2_threshold (float, optional)**: \( r^2 \) threshold to include a variant. Defaults to 0.5.

### Returns

- **dict**: Contains
    * 'pass': Boolean indicating the successful completion of the process.
    * 'step': The label for this procedure ('ld_prune').
    * 'metrics': Metrics associated with the pruning, such as 'ld_removed_count'.
    * 'output': Dictionary containing paths to the generated output files.














# Assoc Class Functions Guide

## Introduction

This guide details the methods available in the `Assoc` class for running various genomic analyses.

## Table of Contents
- [Constructor](#constructor)
- [write_exclusion_file](#write_exclusion_file)
- [run_pca_pruning](#run_pca_pruning)
- [run_plink_pca](#run_plink_pca)
- [calculate_inflation](#calculate_inflation)
- [run_gwas](#run_gwas)
- [run_prs](#run_prs)
- [run_association](#run_association)

---

## Constructor

\`\`\`python
def __init__(self, geno_path=None, out_path=None, pca=10, build='hg38', gwas=True, pheno_name='PHENO1', covar_path=None, covar_names=None)
\`\`\`

### Parameters

- **geno_path**: The path to the genotype file.
- **out_path**: The output directory.
- **pca**: Number of principal components to calculate. Defaults to 10.
- **build**: Genome build, either 'hg19' or 'hg38'.
- **gwas**: Flag to run GWAS analysis.
- **pheno_name**: Phenotype column name.
- **covar_path**: Path to the covariate file.
- **covar_names**: Names of the covariates.

---

## write_exclusion_file

\`\`\`python
def write_exclusion_file(self)
\`\`\`

### Description

Writes genomic regions to be excluded based on the genome build (`hg19` or `hg38`).

### Returns

- **exclusion_file**: Path to the generated exclusion file.

---

## run_pca_pruning

\`\`\`python
def run_pca_pruning(self, maf=0.01, geno=0.01, hwe=5e-6, indep_pairwise=[1000, 10, 0.02])
\`\`\`

### Parameters

- **maf**: Minor allele frequency.
- **geno**: Genotype missing rate.
- **hwe**: Hardy-Weinberg Equilibrium p-value.
- **indep_pairwise**: Parameters for SNP pruning.

---

## run_plink_pca

\`\`\`python
def run_plink_pca(self)
\`\`\`

### Description

Runs PCA on pruned data.

---

## calculate_inflation

\`\`\`python
def calculate_inflation(self, pval_array, normalize=False, ncases=None, ncontrols=None)
\`\`\`

### Parameters

- **pval_array**: Array of p-values.
- **normalize**: Flag to normalize the inflation factor.
- **ncases**: Number of cases.
- **ncontrols**: Number of controls.

---

## run_gwas

\`\`\`python
def run_gwas(self, covars=True)
\`\`\`

### Parameters

- **covars**: Whether to use covariates in the analysis.

---

## run_prs

\`\`\`python
def run_prs(self)
\`\`\`

### Description

Method coming soon.

---

## run_association

\`\`\`python
def run_association(self)
\`\`\`

### Description

Wrapper function to run PCA and GWAS analysis.

---
