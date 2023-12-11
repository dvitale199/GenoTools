# GenoTools

## Getting Started

GenoTools is a suite of automated genotype data processing steps written in Python. The core pipeline was built for Quality Control and Ancestry estimation of data in the Global Parkinson's Genetics Program (GP2)

Setup just requires:
```
git clone https://github.com/dvitale199/GenoTools
cd GenoTools
pip install .
```
you can pull references by running:
```
genotools-download
```
To download specific references/models, you can run the download with the following options:
```
genotools-download --ref 1kg_30x_hgdp_ashk_ref_panel --model nba_v1
```
Currently, `1kg_30x_hgdp_ashk_ref_panel` is the only available reference panel. Available models are `nba_v1` for the NeuroBooster array and `neurochip_v1` for the NeuroChip Array. If using a different array, we would suggest training a new model by running the standard command below.

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

if you'd like to run the pipeline using an existing model, you can do that like so (take note of the `--model` option):
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
  --model nba_v1
```

This will find common snps between your genotype data and the reference panel, run PCA, UMAP-transform PCs, and train a new XGBoost classifier specific to your data/ref panel.

## Documentation
- [GenoTools Command Line Arguments](docs/cli_args.md)
- [Default Pipeline Overview](docs/default_pipeline_overview.md)
- [Package Function Guide (for developers)](docs/genotools_function_guide.md)

## Acknowledgements
GenoTools was developed as the core genotype and wgs processing pipeline for the Global Parkinson's Genetics Program (GP2) at the Center for Alzheimer's and Related Dementias (CARD) at the National Institutes of Health.

This tool relies on PLINK, a whole genome association analysis toolset, for various genetic data processing functionalities. We gratefully acknowledge the developers of PLINK for their foundational contributions to the field of genetics. More about PLINK can be found at [their website](https://www.cog-genomics.org/plink/2.0/).



