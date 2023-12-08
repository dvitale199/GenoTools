# GenoTools
## Acknowledgements
This tool relies on PLINK, a whole genome association analysis toolset, for various genetic data processing functionalities. We gratefully acknowledge the developers of PLINK for their foundational contributions to the field of genetic analysis. More about PLINK can be found at [their website](https://www.cog-genomics.org/plink/2.0/).


## Getting Started

GenoTools is a suite of automated genotype data processing steps written in Python. The core pipeline was built for Quality Control and Ancestry estimation of data in the Global Parkinson's Genetics Program (GP2)

Setup just requires:
```
git clone https://github.com/dvitale199/GenoTools
cd GenoTools
pip install .
```

Modify the paths in the following command to run the standard GP2 pipeline:
```
genotools \
  --geno_path /path/to/genotypes/for/qc \
  --out_path /path/to/qc/output \
  --ancestry \
  --ref_panel /path/to/reference/panel \
  --ref_labels /path/to/reference/ancestry/labels \
  --container \
  --all_sample \
  --all_variant
```
Note: add the ```--singularity``` flag to run ancestry predictions on HPC

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


### Sample-Level QC Steps
Following ancestry analysis, the pipeline undergoes these sample-level QC procedures:

### Sample-Level QC Steps
- **Callrate**: Default threshold: 0.05
- **Sex Check**: Default cutoffs: [0.25, 0.75]
- **Relatedness Check**: Enabled by default
- **Heterozygosity Rate (Het)**: Default range: [-0.25, 0.25]

### Variant-Level QC Steps
- **Case-Control Check**: Default threshold: 1e-4
- **Haplotype Check**: Default threshold: 1e-4
- **Hardy-Weinberg Equilibrium (HWE)**: Default threshold: 1e-4
- **Genotype Missingness (Geno)**: Default threshold: 0.05
- **Linkage Disequilibrium (LD)**: Default setting: None (custom parameters may be provided)



# GenoTools Command Line Arguments Documentation

## Overview
This documentation provides detailed descriptions of the command-line arguments for the `GenoTools` package, a comprehensive tool for genetic data analysis. These arguments enable users to specify various inputs, outputs, and settings for processes like quality control (QC), ancestry analysis, and genome-wide association studies (GWAS).

---

### File I/O Arguments

- **`--bfile`**  
  - *Type*: `str`  
  - *Default*: `None`  
  - *Description*: Path to PLINK 1.9 format genotype file (before the *.bed/bim/fam). Used for specifying genotype data input.

- **`--pfile`**  
  - *Type*: `str`  
  - *Default*: `None`  
  - *Description*: Path to PLINK 2 format genotype file (before the *.pgen/pvar/psam). Another format for genotype data input.

- **`--vcf`**  
  - *Type*: `str`  
  - *Default*: `None`  
  - *Description*: Path to VCF format genotype file. An alternative format for genotype data input.

- **`--out`** *(required)*  
  - *Type*: `str`  
  - *Default*: `None`  
  - *Description*: Prefix for output files, including the path.

- **`--full_output`**  
  - *Type*: `str`  
  - *Default*: `True`  
  - *Description*: If set to true, outputs all results. Otherwise, outputs a subset.

- **`--skip_fails`**  
  - *Type*: `str`  
  - *Default*: `False`  
  - *Description*: If true, skips checks for input file errors.

- **`--warn`**  
  - *Type*: `str`  
  - *Default*: `False`  
  - *Description*: If true, provides warnings on errors but continues the pipeline.

---

### Ancestry Arguments

- **`--ancestry`**  
  - *Type*: `str`  
  - *Default*: `False`  
  - *Description*: Splits the analysis by ancestry if set to true.

- **`--ref_panel`**  
  - *Type*: `str`  
  - *Default*: `None`  
  - *Description*: Path to the PLINK format reference genotype file for ancestry analysis.

- **`--ref_labels`**  
  - *Type*: `str`  
  - *Default*: `None`  
  - *Description*: File with tab-separated PLINK-style IDs and ancestry labels, without a header.

- **`--model`**  
  - *Type*: `str`  
  - *Default*: `None`  
  - *Description*: Path to the pickle file with a trained ancestry model.

- **`--container`**  
  - *Type*: `str`  
  - *Default*: `False`  
  - *Description*: Enables running predictions in a container.

- **`--singularity`**  
  - *Type*: `str`  
  - *Default*: `False`  
  - *Description*: If true, runs containerized predictions using Singularity.

- **`--cloud`**  
  - *Type*: `str`  
  - *Default*: `False`  
  - *Description*: Enables running predictions in Google Cloud.

- **`--cloud_model`**  
  - *Type*: `str`  
  - *Default*: `NeuroBooster`  
  - *Description*: Specifies the model for Google Cloud predictions.

- **`--subset_ancestry`**  
  - *Type*: `None`  
  - *Description*: Subsets the analysis for specified ancestries.

---

### Sample-Level QC Arguments

- **`--callrate`**  
  - *Type*: `float`  
  - *Default*: `None`  
  - *Description*: Minimum call rate threshold for sample-level QC.

- **`--sex`**  
  - *Type*: `None`  
  - *Description*: Specifies sex prune cutoffs.

- **`--related`**  
  - *Type*: `str`  
  - *Default*: `False`  
  - *Description*: Enables relatedness pruning if set to true.

- **`--related_cutoff`**  
  - *Type*: `float`  
  - *Default*: `0.0884`  
  - *Description*: Cutoff for relatedness pruning.

- **`--duplicated_cutoff`**  
  - *Type*: `float`  
  - *Default*: `0.354`  
  - *Description*: Cutoff for duplicated samples pruning.

- **`--prune_related`**  
  - *Type*: `str`  
  - *Default*: `False`  
  - *Description*: If true, prunes related samples.

- **`--prune_duplicated`**  
  - *Type*: `str`  
  - *Default*: `True`  
  - *Description*: If true, prunes duplicated samples.

- **`--het`**  
  - *Type*: `None`  
  - *Description*: Specifies heterozygosity prune cutoffs.

- **`--all_sample`**  
  - *Type*: `str`  
  - *Default*: `False`  
  - *Description*: If true, runs all sample-level QC.

---

### Variant-Level QC Arguments

- **`--geno`**  
  - *Type*: `float`  
  - *Default*: `None`  
  - *Description*: Minimum missingness threshold for variant-level QC.

- **`--case_control`**  
  - *Type*: `float`  
  - *Default*: `None`  
  - *Description*: Threshold for case-control prune.

- **`--haplotype`**  
  - *Type*: `float`  
  - *Default*: `None`  
  - *Description*: Threshold for haplotype prune.

- **`--hwe`**  
  - *Type*: `float`  
  - *Default*: `None`  
  - *Description*: Threshold for Hardy-Weinberg Equilibrium pruning.

- **`--filter_controls`**  
  - *Type*: `str`  
  - *Default*: `False`  
  - *Description*: If true, filters controls during HWE pruning.

- **`--ld`**  
  - *Type*: `None`  
  - *Description*: Specifies parameters for linkage disequilibrium pruning.

- **`--all_variant`**  
  - *Type*: `str`  
  - *Default*: `False`  
  - *Description*: If true, runs all variant-level QC.

---

### GWAS and PCA Arguments

- **`--pca`**  
  - *Type*: `int`  
  - *Default*: `None`  
  - *Description*: Number of principal components for PCA analysis.

- **`--build`**  
  - *Type*: `str`  
  - *Default*: `hg38`  
  - *Description*: Genome build used for PCA.

- **`--gwas`**  
  - *Type*: `str`  
  - *Default*: `False`  
  - *Description*: Enables running GWAS if set to true.

- **`--covars`**  
  - *Type*: `str`  
  - *Default*: `None`  
  - *Description*: Path to external covariates file.

- **`--covar_names`**  
  - *Type*: `str`  
  - *Default*: `None`  
  - *Description*: Names of covariates to use from the external file.

---

This documentation provides a comprehensive overview of the available command-line arguments for `GenoTools`. It is designed to assist users in effectively utilizing the tool for their genetic data analysis needs.
