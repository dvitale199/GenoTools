# GenoTools
## Published in G3: [https://www.biorxiv.org/content/10.1101/2024.03.26.586362v1.full.pdf](https://doi.org/10.1093/g3journal/jkae268)
[![DOI](https://zenodo.org/badge/337965715.svg)](https://zenodo.org/doi/10.5281/zenodo.10443257)
[![PyPI version](https://badge.fury.io/py/the-real-genotools.svg)](https://badge.fury.io/py/the-real-genotools)
[![PRs Welcome](https://img.shields.io/badge/PRs-welcome-brightgreen.svg?style=flat-square)](http://makeapullrequest.com)
![GitHub License](https://img.shields.io/github/license/dvitale199/GenoTools)
![Python](https://img.shields.io/badge/python-3.8-blue.svg)
![Python](https://img.shields.io/badge/python-3.9-blue.svg)
![Python](https://img.shields.io/badge/python-3.10-blue.svg)


## Documentation
You can find the full documentation with the following links:
- [GenoTools Command Line Arguments](https://github.com/dvitale199/GenoTools/blob/main/docs/cli_args.md)
- [Default Pipeline Overview](https://github.com/dvitale199/GenoTools/blob/main/docs/default_pipeline_overview.md)
- [Package Function Guide (for developers)](https://github.com/dvitale199/GenoTools/blob/main/docs/genotools_function_guide.md)
- [JSON output guide](https://github.com/dvitale199/GenoTools/blob/main/docs/json_output_overview.md)

## Getting Started

GenoTools is a suite of automated genotype data processing steps written in Python. The core pipeline was built for Quality Control and Ancestry estimation of data in the Global Parkinson's Genetics Program (GP2)

To download the most current version from pip:
```
pip install the-real-genotools
```
Alternatively, if you'd like to download from github:
```
git clone https://github.com/dvitale199/GenoTools.git
cd GenoTools
pip install .
```
you can pull the most current references by running:
```
genotools-download
```
By default, the reference panel will be downloaded to ~/.genotools/ref. but can be download to a location of choice with `--destination`.

To download specific references/models, you can run the download with the following options:
```
genotools-download --ref 1kg_30x_hgdp_ashk_ref_panel --model nba_v1 --destination /path/to/download_directory/
```

Currently, `1kg_30x_hgdp_ashk_ref_panel` is the only available reference panel. Available models are `nba_v1` for the NeuroBooster array and `neurochip_v1` for the NeuroChip Array and both are in GRCh38. If using a different array, we would suggest training a new model by running the standard command below. Please ensure the reference panel and your genotypes are in the same build. If you're using our reference panel, your genotypes must be in GRCh38.

Modify the paths in the following command to run the standard GP2 pipeline:
```
genotools \
  --pfile /path/to/genotypes/for/qc \
  --out /path/to/qc/output \
  --ancestry \
  --ref_panel /path/to/reference/panel \
  --ref_labels /path/to/reference/ancestry/labels \
  --all_sample \
  --all_variant
```
This will find common snps between your genotype data and the reference panel, run PCA, UMAP-transform PCs, and train a new XGBoost classifier specific to your data/ref panel.

if you'd like to run the pipeline using an existing model, you can do that like so (take note of the `--model` option):
```
genotools \
  --pfile /path/to/genotypes/for/qc \
  --out /path/to/qc/output \
  --ancestry \
  --ref_panel /path/to/reference/panel \
  --ref_labels /path/to/reference/ancestry/labels \
  --all_sample \
  --all_variant
  --model /path/to/nba_v1/model
```

if you'd like to run the pipeline using the default nba_v1 model in a Docker container, you can do that like so:
```
genotools \
  --pfile /path/to/genotypes/for/qc \
  --out /path/to/qc/output \
  --ancestry \
  --ref_panel /path/to/reference/panel \
  --ref_labels /path/to/reference/ancestry/labels \
  --container \
  --all_sample \
  --all_variant
```
Note: add the ```--singularity``` flag to run containerized ancestry predictions on HPC

genotools accept `--pfile`, `--bfile`, or `--vcf`. Any bfile or vcf will be converted to a pfile before running any steps. 

Please consult the docs links listed at the top of the README for the full argument guide, function guide, Default pipeline overview, and guide for navigating the output JSON.

## Acknowledgements
GenoTools was developed as the core genotype and wgs processing pipeline for the Global Parkinson's Genetics Program (GP2) at the Center for Alzheimer's and Related Dementias (CARD) at the National Institutes of Health.

This tool relies on PLINK, a whole genome association analysis toolset, for various genetic data processing functionalities. We gratefully acknowledge the developers of PLINK for their foundational contributions to the field of genetics. More about PLINK can be found at [their website](https://www.cog-genomics.org/plink/2.0/).



