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

### Ancestry Prediction
The pipeline begins with **ancestry prediction**, a critical step for discerning the genetic background of samples, which greatly influences subsequent analytical processes.

### Sample-Level QC Steps
Following ancestry analysis, the pipeline undergoes these sample-level QC procedures:

- **Callrate**: Implements checks and filters based on the call rate threshold, crucial for ensuring data completeness.
- **Sex Check**: Conducts verification and correction of sex information within the genotype data, a vital step for precise downstream analysis.
- **Relatedness Check**: Identifies and manages related individuals to avoid biases in the analysis.
- **Heterozygosity Rate (Het)**: Evaluates and filters samples based on heterozygosity rates, crucial for detecting potential data outliers or inconsistencies.

### Variant-Level QC Steps
After completing sample-level QC, the pipeline transitions to variant-level QC, encompassing:

- **Case-Control Check**: Filters variants on the basis of case-control discrepancies, pivotal for association studies.
- **Haplotype Check**: Assesses haplotype stability, essential for preserving genomic integrity.
- **Hardy-Weinberg Equilibrium (HWE)**: Confirms allele frequencies are in equilibrium, a fundamental element in population genetics.
- **Genotype Missingness (Geno)**: Excludes variants with substantial missingness, key to maintaining data quality.
- **Linkage Disequilibrium (LD)**: Reviews and prunes variants in linkage disequilibrium, vital for minimizing redundancy and bias.



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
