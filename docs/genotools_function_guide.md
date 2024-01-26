# GenoTools Functions Guide

## Introduction

This guide details the methods available in the `GenoTools` package for quality control, ancestry prediction, and analysis of genotype data.

## Table of Contents
- [SampleQC](#SampleQC)
    - [run_callrate_prune](#run_callrate_prune)
    - [run_sex_prune](#run_sex_prune)
    - [run_het_prune](#run_het_prune)
    - [run_related_prune](#run_related_prune)
- [VariantQC](#VariantQC)
    - [run_geno_prune](#run_geno_prune)
    - [run_case_control_prune](#run_case_control_prune)
    - [run_haplotype_prune](#run_haplotype_prune)
    - [run_hwe_prune](#run_hwe_prune)
    - [run_ld_prune](#run_ld_prune)
- [Assoc](#Assoc)
    - [write_exclusion_file](#write_exclusion_file)
    - [run_pca_pruning](#run_pca_pruning)
    - [run_plink_pca](#run_plink_pca)
    - [calculate_inflation](#calculate_inflation)
    - [run_gwas](#run_gwas)
    - [run_prs](#run_prs)
    - [run_association](#run_association)
- [Ancestry](#Ancestry)
    - [get_raw_files](#get_raw_files)
    - [munge_training_data](#munge_training_data)
    - [calculate_pcs](#calculate_pcs)
    - [transform](#transform)
    - [train_umap_classifier](#train_umap_classifier)
    - [load_umap_classifier](#load_umap_classifier)
    - [predict_ancestry_from_pcs](#predict_ancestry_from_pcs)
    - [get_containerized_predictions](#get_containerized_predictions)
    - [get_cloud_predictions](#get_cloud_predictions)
    - [predict_admixed_samples](#predict_admixed_samples)
    - [umap_transform_with_fitted](#umap_transform_with_fitted)
    - [split_cohort_ancestry](#split_cohort_ancestry)
    - [run_ancestry](#run_ancestry)


---

# Sample Quality Control

## SampleQC

```python
__init__(self, geno_path=None, out_path=None)
```

### Description

Initializes the `SampleQC` class with paths to genotype data and output path.

### Parameters
- **geno_path**: Path to the genotype data.
- **out_path**: Path for output data.

---

## run_callrate_prune

```python
run_callrate_prune(self, mind=0.02)
```

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

```python
run_sex_prune(self, check_sex=[0.25,0.75])
```

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

```python
run_het_prune(self, het_threshold=0.1)
```

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

```python
run_related_prune(self, related_cutoff=0.0884, duplicated_cutoff=0.354, prune_related=True, prune_duplicated=True)
```

### Description

Execute pruning based on relatedness and duplication checks on genotype data using PLINK and KING.

### Parameters

- **related_cutoff** (float, optional): Threshold for relatedness check. Default value is 0.0884.
- **duplicated_cutoff** (float, optional): Threshold for duplication check. Default value is 0.354.
- **prune_related** (bool, optional): Whether to prune related samples. Default is True.
- **prune_duplicated** (bool, optional): Whether to prune duplicated samples. Default is True.

### Returns

- **dict**: A structured dictionary containing:
  * 'pass': Boolean indicating the successful completion of the process.
  * 'step': The label for this procedure ('related_prune').
  * 'metrics': Metrics associated with the pruning.
  * 'output': Dictionary containing paths to the generated output files.

---

# Variant Quality Control

## VariantQC

```python
__init__(self, geno_path=None, out_path=None)
```

### Description

Initializes the `VariantQC` class with paths to genotype data and output path.

### Parameters
- **geno_path**: Path to the genotype data.
- **out_path**: Path for output data.

---

## run_geno_prune

```python
run_geno_prune(self, geno_threshold=0.05)
```

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

```python
run_case_control_prune(self, p_threshold=1e-4)
```

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

```python
run_haplotype_prune(self, p_threshold=1e-4)
```

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

```python
run_hwe_prune(self, hwe_threshold=1e-4, filter_controls=False)
```

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

```python
run_ld_prune(self, window_size=50, step_size=5, r2_threshold=0.5)
```

### Description

Prunes SNPs based on Linkage Disequilibrium.

### Parameters

- **window_size**: Size of the sliding window over the genome in variant count. Defaults to 50.
- **step_size**: Step of the sliding window over the genome. Defaults to 5.
- **r2_threshold**: \(r^2\) threshold to include a variant. Defaults to 0.5.

### Returns

- **dict**: Contains
    * 'pass': Boolean indicating the successful completion of the process.
    * 'step': The label for this procedure ('ld_prune').
    * 'metrics': Metrics associated with the pruning, such as 'ld_removed_count'.
    * 'output': Dictionary containing paths to the generated output files.


# Associations

---

## Assoc

```python
__init__(self, geno_path=None, out_path=None, pca=10, build='hg38', gwas=True, pheno_name='PHENO1', covar_path=None, covar_names=None)
```
### Description

Initializes the `Assoc` class with paths to genotype data and output path, pca, build, gwas, pheno_name, covar_path, and covar_names.

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

```python
write_exclusion_file(self)
```

### Description

Writes genomic regions to be excluded based on the genome build (`hg19` or `hg38`).

### Returns

- **exclusion_file**: Path to the generated exclusion file.

---

## run_pca_pruning

```python
run_pca_pruning(self, maf=0.01, geno=0.01, hwe=5e-6, indep_pairwise=[1000, 10, 0.02])
```

### Parameters

- **maf**: Minor allele frequency.
- **geno**: Genotype missing rate.
- **hwe**: Hardy-Weinberg Equilibrium p-value.
- **indep_pairwise**: Parameters for SNP pruning.

---

## run_plink_pca

```python
run_plink_pca(self)
```

### Description

Runs PCA on pruned data.

---

## calculate_inflation

```python
calculate_inflation(self, pval_array, normalize=False, ncases=None, ncontrols=None)
```

### Parameters

- **pval_array**: Array of p-values.
- **normalize**: Flag to normalize the inflation factor.
- **ncases**: Number of cases.
- **ncontrols**: Number of controls.

---

## run_gwas

```python
def run_gwas(self, covars=True)
```

### Parameters

- **covars**: Whether to use covariates in the analysis.

---

## run_prs

```python
def run_prs(self)
```
### Description

Method coming soon.

---

## run_association

```python
def run_association(self)
```

### Description

Wrapper function to run PCA and GWAS analysis.

---


# Ancestry Predictions

---

## Ancestry

```python
__init__(self, geno_path=None, ref_panel=None, ref_labels=None, out_path=None, model_path=None, containerized=False, singularity=False, subset=None)
```

### Description

Initializes the `Ancestry` class with paths to genotype data, reference panels, and other options.

### Parameters

- **geno_path**: Path to the genotype data.
- **ref_panel**: Path to the reference panel.
- **ref_labels**: Path to the reference labels.
- **out_path**: Path for output data.
- **model_path**: Path to the trained model.
- **containerized**: Boolean indicating if the computation is containerized.
- **singularity**: Boolean indicating if Singularity is used.
- **subset**: Any subsets to be considered.

---

## get_raw_files

```python
get_raw_files(self)
```

### Description

Processes reference and genotype data for prediction, including variant pruning, extracting common SNPs, obtaining raw versions of common SNPs, and preparing data with labels.

### Returns

- **dict**: Contains
    * 'raw_ref': A DataFrame with labeled raw reference data.
    * 'raw_geno': A DataFrame with labeled raw genotype data.
    * 'out_paths': A dictionary with output file paths.

---

## munge_training_data

```python
munge_training_data(self, labeled_ref_raw)
```

### Description

Preprocesses labeled raw data for prediction by performing train/test split, label encoding, and data formatting.

### Parameters

- **labeled_ref_raw**: A DataFrame with labeled raw reference data and PCA information.

### Returns

- **dict**: Contains
    * 'X_train': Features of the training set.
    * 'X_test': Features of the test set.
    * 'y_train': Labels of the training set.
    * 'y_test': Labels of the test set.
    * 'train_ids': IDs of samples in the training set.
    * 'test_ids': IDs of samples in the test set.
    * 'label_encoder': A fitted LabelEncoder for label transformation.
    * 'X_all': All features (including both training and test sets).
    * 'y_all': All labels (including both training and test sets).

## calculate_pcs

```python
calculate_pcs(self, X_train, X_test, y_train, y_test, train_ids, test_ids, raw_geno, label_encoder)
```

### Description

Calculates principal components (PCs) for training and testing datasets, transforms the data, and projects new samples onto the reference panel. The function also saves the PCA data to files.

### Parameters

- **X_train** (DataFrame): Features of the training set.
- **X_test** (DataFrame): Features of the test set.
- **y_train** (Series): Labels of the training set.
- **y_test** (Series): Labels of the test set.
- **train_ids** (DataFrame): IDs of samples in the training set.
- **test_ids** (DataFrame): IDs of samples in the test set.
- **raw_geno** (DataFrame): Raw genotype data.
- **label_encoder** (LabelEncoder): A fitted LabelEncoder for label transformation.

### Returns

- **dict**: Contains
    * 'X_train': Transformed features of the training set.
    * 'X_test': Transformed features of the test set.
    * 'labeled_train_pca': Labeled PCA data for the training set.
    * 'labeled_ref_pca': Labeled PCA data for the reference panel.
    * 'new_samples_projected': Projected PCA data for new samples.
    * 'out_paths': Dictionary containing output file paths.


## transform

```python
transform(self, data, mean, sd, pca, col_names, fit=False)
```

### Description

Applies flashPCA-style scaling and PCA transformation to the input data.

### Parameters

- **data** (DataFrame): Input data to be transformed.
- **mean** (Series): Mean values for flashPCA-style scaling.
- **sd** (Series): Standard deviation values for flashPCA-style scaling.
- **pca** (PCA): PCA model for transforming data.
- **col_names** (list): List of column names for the transformed data.
- **fit** (bool, optional): If True, fit-transform the data with PCA. If False, only transform. Default is False.

### Returns

- **DataFrame**: Transformed data with named columns.

---

## train_umap_classifier

```python
train_umap_classifier(self, X_train, X_test, y_train, y_test, label_encoder)
```

### Description

Train a UMAP to linear XGBoost classifier pipeline.

### Parameters

- **X_train** (DataFrame): Training features.
- **X_test** (DataFrame): Testing features.
- **y_train** (Series): Training labels.
- **y_test** (Series): Testing labels.
- **label_encoder**: LabelEncoder object for encoding and decoding labels.

### Returns

- **dict**: Dictionary containing classifier, label encoder, parameters, confusion matrix, fitted grid, train accuracy, test accuracy, and model path.

---

## load_umap_classifier

```python
load_umap_classifier(self, X_test, y_test)
```

### Description

Load a trained UMAP classifier from a pickle file and evaluate its performance on the test set.

### Parameters

- **X_test** (DataFrame): Testing features.
- **y_test** (Series): Testing labels.

### Returns

- **dict**: Dictionary containing classifier, confusion matrix, test accuracy, and model parameters.

---

## predict_ancestry_from_pcs

```python
predict_ancestry_from_pcs(self, projected, pipe_clf, label_encoder)
```

### Description

Predict ancestry labels for new samples based on their projected principal components.

### Parameters

- **projected** (DataFrame): Testing features.
- **pipe_clf** (Series): Testing labels.
- **label_encoder** LabelEncoder object for encoding and decoding labels.

### Returns

- **dict**: Dictionary containing predicted labels, output data, and metrics.

---

## get_containerized_predictions

```python
get_containerized_predictions(self, X_test, y_test, projected, label_encoder)
```

### Description

Get predictions using a containerized environment for UMAP and XGBoost classifier.

### Parameters

- **X_test** (DataFrame): Testing features.
- **y_test** (Series): Testing labels.
- **projected** (DataFrame): Testing features.
- **label_encoder** LabelEncoder object for encoding and decoding labels.

### Returns

- **tuple**: Two dictionaries containing trained classifier results and prediction results.

---

## get_cloud_predictions

```python
get_cloud_predictions(self, X_test, y_test, projected, label_encoder, train_pca)
```

### Description

Get predictions using a cloud environment for UMAP and XGBoost classifier.

### Parameters

- **X_test** (DataFrame): Test data.
- **y_test** (Series): True labels for the test data.
- **projected** (DataFrame): Projected principal components of new samples.
- **label_encoder**: Label encoder used for encoding ancestry labels.
- **train_pca**: Labeled PCs for training data.

### Returns

- **tuple**: Two dictionaries containing trained classifier results and prediction results.

---

## predict_admixed_samples

```python
predict_admixed_samples(self, projected, train_pca)
```

### Description

Change labels of samples with complex admixture, calculated based off training PCs.

### Parameters

- **projected** (DataFrame): Projected principal components of new samples.
- **train_pca**: Labeled PCs for training data.

### Returns

- **DataFrame**: Projected principal components of new samples with updated labels.

---

## umap_transform_with_fitted

```python
umap_transform_with_fitted(self, ref_pca, X_new, y_pred, params=None)
```

### Description

Transform data using a fitted UMAP components.

### Parameters

- **ref_pca** (DataFrame): Reference PCA data with labels.
- **X_new** (DataFrame): New samples to be transformed.
- **y_pred** (DataFrame): Predicted labels for new samples.
- **params** (dict, optional): UMAP parameters. Defaults to None.

### Returns

- **dict**:  dict: Dictionary containing UMAP-transformed data.

---

## split_cohort_ancestry

```python
split_cohort_ancestry(self, labels_path)
```

### Description

Split a cohort based on predicted ancestries.

### Parameters
- **labels_path** (str): Path to the file containing predicted labels.
- **subset** (list, optional): List of ancestries to continue analysis for. Defaults to False.

### Returns
- **dict**: Dictionary containing labels and paths for each split ancestry.

---

## run_ancestry

```python
run_ancestry(self)
```

### Description

Run the ancestry prediction pipeline.

###  Returns
- **dict**: Dictionary containing data, metrics, and output information.
