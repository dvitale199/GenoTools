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
COMING SOON

