# GenoTools

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


## Documentation
- [GenoTools Command Line Arguments](docs/cli_args.md)
- [Default Pipeline Overview](docs/default_pipeline_overview.md)
- [Package Function Guide (for developers)](docs/genotools_function_guide.md)

## Acknowledgements
GenoTools was developed as the core genotype and wgs processing pipeline for the Global Parkinson's Genetics Program (GP2) at the Center for Alzheimer's and Related Dementias (CARD) at the National Institutes of Health.

This tool relies on PLINK, a whole genome association analysis toolset, for various genetic data processing functionalities. We gratefully acknowledge the developers of PLINK for their foundational contributions to the field of genetics. More about PLINK can be found at [their website](https://www.cog-genomics.org/plink/2.0/).



