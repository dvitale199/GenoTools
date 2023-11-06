# GenoTools

# WARNING: MAIN BRANCH IS UNDER DEVELOPMENT. PLEASE PULL FROM 'WORKING' FOR THE LATEST VERSION OF THE FUNCTIONAL PIPELINES

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
genotools --geno_path /path/to/genotypes/for/qc --out_path /path/to/qc/output --ancestry --ref_panel /path/to/reference/panel --ref_labels /path/to/reference/ancestry/labels --container --all_sample --all_variant
```
Note: add the ```--singularity``` flag to run ancestry predictions on HPC
