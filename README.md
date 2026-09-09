# GenoTools
## Published in G3: [https://www.biorxiv.org/content/10.1101/2024.03.26.586362v1.full.pdf](https://doi.org/10.1093/g3journal/jkae268)
[![DOI](https://zenodo.org/badge/337965715.svg)](https://zenodo.org/doi/10.5281/zenodo.10443257)
[![PyPI version](https://badge.fury.io/py/the-real-genotools.svg)](https://badge.fury.io/py/the-real-genotools)
[![PRs Welcome](https://img.shields.io/badge/PRs-welcome-brightgreen.svg?style=flat-square)](http://makeapullrequest.com)
![GitHub License](https://img.shields.io/github/license/dvitale199/GenoTools)
![Python](https://img.shields.io/badge/python-3.8-blue.svg)
![Python](https://img.shields.io/badge/python-3.9-blue.svg)
![Python](https://img.shields.io/badge/python-3.10-blue.svg)


> [!IMPORTANT]
> **2.1.0 is the first 2.x release and supersedes 1.3.6.** Upgrading requires
> `--upgrade`; a plain `pip install` leaves an existing install on 1.3.6.
>
> ```bash
> pip install --upgrade the_real_genotools
> pip show the_real_genotools          # confirm 2.1.0
> ```
>
> This is a breaking upgrade: flags moved from `underscore_style` to
> `hyphen-style` (old spellings still work but warn), the JSON report gained
> fields, and several ancestry fixes can move predicted labels. **1.x ancestry
> models cannot be loaded** — see below. Read
> [MIGRATION_2.0.md](MIGRATION_2.0.md) before upgrading;
> [CHANGELOG.md](CHANGELOG.md) has the summary.

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

Currently, `1kg_30x_hgdp_ashk_ref_panel` is the only available reference panel. Available models, all in GRCh38:

| Model | Array | Format |
|---|---|---|
| `nba_gp2_r12` (default) | NeuroBooster, trained on GP2 release 12 | **2.x** |
| `nba_v1`, `nba_v2` | NeuroBooster | 1.x only |
| `neurochip_v1` | NeuroChip | 1.x only |

The two formats are mutually incompatible: GenoTools 2.x loads only `nba_gp2_r12`, and 1.x loads only the other three. `genotools-download` warns if you ask for a 1.x model. If using a different array, we would suggest training a new model by running the standard command below. Please ensure the reference panel and your genotypes are in the same build. If you're using our reference panel, your genotypes must be in GRCh38.

> **2.x cannot load 1.x ancestry models.** `nba_v1`, `nba_v2` and
> `neurochip_v1` were trained by 1.x and are rejected with an explanatory
> error. Use `nba_gp2_r12`, the 2.x-format NeuroBooster model and now the
> `genotools-download` default, or train your own with
> `--ref-panel`/`--ref-labels` as in the standard command below.
>
> **If you retrain a model, expect ~1.3% of ancestry labels to move.** That is
> a property of retraining rather than of the fixes, and neither labeling is
> demonstrably more correct — see
> [MIGRATION_2.0.md](MIGRATION_2.0.md#retraining-moves-about-13-of-ancestry-labels--and-that-is-retraining-not-the-fix).

Modify the paths in the following command to run the standard GP2 pipeline:
```
genotools \
  --pfile /path/to/genotypes/for/qc \
  --out /path/to/qc/output \
  --ancestry \
  --ref-panel /path/to/reference/panel \
  --ref-labels /path/to/reference/ancestry/labels \
  --all-sample \
  --all-variant
```
This will find common snps between your genotype data and the reference panel, run PCA, UMAP-transform PCs, and train a new XGBoost classifier specific to your data/ref panel.

if you'd like to run the pipeline using an existing model, you can do that like so (take note of the `--model` option):
```
genotools \
  --pfile /path/to/genotypes/for/qc \
  --out /path/to/qc/output \
  --ancestry \
  --ref-panel /path/to/reference/panel \
  --ref-labels /path/to/reference/ancestry/labels \
  --all-sample \
  --all-variant \
  --model /path/to/nba_v1/model
```

Note: `--container`, `--singularity` and `--cloud` are **not supported in 2.0** and exit with an error if passed. Ancestry prediction runs locally, in process. 1.x's containerized prediction relied on a Docker image built around a 1.x model that 2.0 cannot load; `--cloud` was never implemented in any version. See [MIGRATION_2.0.md](MIGRATION_2.0.md).

genotools accepts `--pfile`, `--bfile`, or `--vcf`. Any bfile or vcf will be converted to a pfile before running any steps. 

Note: multiallelic pfiles will be converted to biallelic format by excluding multiallelic variants before running '--ancestry' steps. If you would prefer to not remove multiallelic snps, please pre-split the SNPs using bcftools prior to running genotools.  

Please consult the docs links listed at the top of the README for the full argument guide, function guide, Default pipeline overview, and guide for navigating the output JSON.

## Acknowledgements
GenoTools was developed as the core genotype and wgs processing pipeline for the Global Parkinson's Genetics Program (GP2) at the Center for Alzheimer's and Related Dementias (CARD) at the National Institutes of Health.

This tool relies on PLINK, a whole genome association analysis toolset, for various genetic data processing functionalities. We gratefully acknowledge the developers of PLINK for their foundational contributions to the field of genetics. More about PLINK can be found at [their website](https://www.cog-genomics.org/plink/2.0/).



