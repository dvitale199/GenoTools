# GenoTools

## Getting Started

GenoTools is a suite of automated genotype data processing steps written in Python. The core pipeline was built for Quality Control and Ancestry estimation of data in the Global Parkinson's Genetics Program (GP2)

The core pipeline can be called as:

`python3 run_qc_pipeline.py --geno <genotype file to be QC'd (plink format)> --ref <genotype file of reference panel (plink format)> --ref_labels <labels for reference panel ancestry> --out <path and prefix for output>`

### options
`--geno`: Path to genotype to be processed in Plink (.bed/.bim/.fam) format. Include everything up to .bed/.bim/.fam. MUST HAVE 
PHENOTYPES OR SEVERAL STEPS WILL FAIL

`--ref`: Path to reference panel genotype in Plink (.bed/.bim/.fam) format. Include everything up to .bed/.bim/.fam. For GP2, we use a combination of 1kGenomes + Ashkenazi Jewish Reference panel

`--ref_labels`: Path to a tab-separated (Plink-style) file containing ancestry labels for the reference panel with the following columns: `FID  IID  label` with NO HEADER. 

`--out`: Path and prefix to output QC'd and ancestry-predicted genotype in Plink (.bed/.bim/.fam) format. Include everything up to .bed/.bim/.fam. For now, this just outputs a log file but will include function to move output files soon.

## Core Pipeline Overview

The core pipeline is broken down into 3 main pieces:
1. Sample-level Quality Control
2. Ancestry Estimation
3. Variant-level Quality Control

The quality control steps have been developed in large part by: Cornelis Blauwendraat, Mike Nalls, Hirotaka Iwaki, Sara Bandres-Ciga, Mary Makarious, Ruth Chia, Frank Grenn, Hampton Leonard, Monica Diez-Fairen, Jeff Kim of the Laboratory of Neurogenetics and Center for Alzheimer's and Related Dementias at the National Institute on Aging, NIH and has been adapted into an automated Python package by Dan Vitale.

### Sample-level Quality Control
1. `callrate_prune(geno_path, out_path, mind=0.02)`
2. `sex_prune(geno_path, out_path, check_sex=[0.25,0.75])`


## Function Reference
### QC.qc

Python Function | Parameters | Returns |
------------ | ------------- | ------------- |
`callrate_prune(geno_path, out_path, mind=0.02)` <br /><br /> Call Rate Pruning | `geno_path` *str*:  Path to the Plink genotypes (everything before .bed/.bim/.fam). <br /><br /> `out_path` *str*: Path to the output Plink genotypes (everything before .bed/.bim/.fam). <br /><br />`mind` *float*: excludes with more than 2% missing genotypes by default. This is much more stringent than Plink default missingness threshold of 10% | *dict* <br /><br />`'step'` *str*:  Name of step in pipeline <br /><br /> `'metrics'` *dict*: metrics from step with keys:  <br />`'outlier_count'` <br /><br />`'output'` *dict*: paths to output files with keys:  <br />`'outliers_path'`<br />`'plink_out'`<br />`'phenos_path'`
`sex_prune(geno_path, out_path, check_sex=[0.25,0.75])` <br /><br /> Sex Pruning. Done in 2 steps:<br /> 1. Plink `--check-sex` on whole genotype <br /> 2. Plink `--check-sex` on X chromosome (`--chr 23 --from-bp 2699520 --to-bp 154931043`) | `geno_path` *str*:  Path to the Plink genotypes (everything before .bed/.bim/.fam). <br /><br /> `out_path` *str*: Path to the output Plink genotypes (everything before .bed/.bim/.fam). <br /><br /> `check_sex` *list*: two values indicating F threshold. A male call is made if F is more than 0.75; a femle call is made if F is less than 0.25, which is less stringent than Plink default of 0.8 and 0.2, respectively. | *dict* <br /><br />`'step'` *str*:  Name of step in pipeline <br /><br /> `'metrics'` *dict*: metrics from step with keys:  <br />`'outlier_count'` <br /><br />`'output'` *dict*: paths to output files with keys:  <br />`'sex_fails'`<br />`'plink_out'`


# More Coming Soon!!!!
