# GenoTools JSON Output Documentation

## Overview
This documentation provides a detailed description of the keys in JSON output for the `GenoTools` package, a comprehensive tool for genetic data analysis. This output file provides users with important output and summary statistics for processes like quality control (QC), ancestry analysis, and genome-wide association studies (GWAS). Please note that output associated with steps that are not completed in a `GenoTools` run will not be output to the JSON file.

---

### JSON Dictionary Keys

- **`input_samples`**   
  - *Process*: `All`  
  - *Description*: Input .psam file.

- **`ancestry_counts`**  
  - *Process*: `Ancestry`  
  - *Description*: Counts of the number of samples in each predicted ancestry group.

- **`ancestry_labels`**  
  - *Process*: `Ancestry`  
  - *Description*: DataFrame containing IDs and predicted label for each sample.

- **`confusion_matrix`**  
  - *Process*: `Ancestry`  
  - *Description*: Reference panel test set confusion matrix.
 
- **`test_accuracy`**  
  - *Process*: `Ancestry`  
  - *Description*: Reference panel test set accuracy.
 
- **`ref_pcs`**  
  - *Process*: `Ancestry`  
  - *Description*: DataFrame containing the reference panel prinicpal components.

- **`projected_pca`**  
  - *Process*: `Ancestry`  
  - *Description*: DataFrame containing the input samples prinicpal components.

- **`total_umap`**  
  - *Process*: `Ancestry`  
  - *Description*: DataFrame containing the UMAP representation for the reference panel and input samples.

- **`ref_umap`**  
  - *Process*: `Ancestry`  
  - *Description*: DataFrame containing the UMAP representation for the reference panel.

- **`new_samples_umap`**  
  - *Process*: `Ancestry`  
  - *Description*: DataFrame containing the UMAP representation for the input samples.

- **`QC`**  
  - *Process*: `Quality Control (QC)`  
  - *Description*: DataFrame containing summary statistics from completed QC steps.

- **`GWAS`**  
  - *Process*: `Genome-wide Assoication Study (GWAS)`  
  - *Description*: DataFrame containing lambda and lambda1000 values from completed association analysis/analyses.

- **`pruned_samples`**  
  - *Process*: `Ancestry/Quality Control (QC)`  
  - *Description*: DataFrame containing samples pruned due to insufficient ancestry count or sample-level QC.

- **`related_samples`**  
  - *Process*: `Quality Control (QC)`  
  - *Description*: DataFrame containing relatedness information output by PLINK/KING from completed relatedness pruning.
