# GenoTools

## Getting Started

GenoTools is a suite of automated genotype data processing steps written in Python. The core pipeline was built for Quality Control and Ancestry estimation of data in the Global Parkinson's Genetics Program (GP2)

Setup just requires:
```
git clone https://github.com/dvitale199/GenoTools
cd GenoTools
pip install .
```

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
`het_prune(geno_path, out_path)` <br /><br /> Heterozygosity Pruning | `geno_path` *str*:  Path to the Plink genotypes (everything before .bed/.bim/.fam). <br /><br /> `out_path` *str*: Path to the output Plink genotypes (everything before .bed/.bim/.fam). <br /><br /> | *dict* <br /><br />`'step'` *str*:  Name of step in pipeline <br /><br /> `'metrics'` *dict*: metrics from step with keys:  <br />`'outlier_count'` <br /><br />`'output'` *dict*: paths to output files with keys:  <br />`'het_outliers'`<br />`'plink_out'`
`related_prune(geno_path, out_path, related_grm_cutoff=0.125, duplicated_grm_cutoff=0.95)` <br /><br /> Relatedness Pruning. Done using GCTA --grm | `geno_path` *str*:  Path to the Plink genotypes (everything before .bed/.bim/.fam). <br /><br /> `out_path` *str*: Path to the output Plink genotypes (everything before .bed/.bim/.fam). <br /><br /> `related_grm_cutoff` *float*: GRM cutoff for related samples <br /><br /> `duplicated_grm_cutoff` *float*: GRM cutoff for duplicated samples <br /><br />| *dict* <br /><br />`'step'` *str*:  Name of step in pipeline <br /><br /> `'metrics'` *dict*: metrics from step with keys:  <br />`'related_count'`, `'duplicated_count'` <br /><br />`'output'` *dict*: paths to output files with keys:  <br />`'relateds'`<br />`'plink_out'`
`variant_prune(geno_path, out_path)` <br /><br /> Variant Pruning. Missingness by: 1. case/control, 2. haplotype. Filtering controls for Hardy-Weinberg Equilibrium (remove hwe_p > 1e-4 | `geno_path` *str*:  Path to the Plink genotypes (everything before .bed/.bim/.fam). <br /><br /> `out_path` *str*: Path to the output Plink genotypes (everything before .bed/.bim/.fam). <br /><br /> | *dict* <br /><br />`'step'` *str*:  Name of step in pipeline <br /><br /> `'metrics'` *dict*: metrics from step with keys:  <br />`'geno_removed_count'`, `'mis_removed_count'`,`'haplotype_removed_count'`,`'hwe_removed_count'`,`'total_removed_count'` <br /><br />`'output'` *dict*: paths to output files with keys:  <br />`'plink_out'`
`avg_miss_rates(geno_path, out_path)` <br /><br /> Calculate average missingness rates (sample-level and variant-level. | `geno_path` *str*:  Path to the Plink genotypes (everything before .bed/.bim/.fam). <br /><br /> `out_path` *str*: Path to the output Plink genotypes (everything before .bed/.bim/.fam). <br /><br /> | *dict* <br /><br />`'step'` *str*:  Name of step in pipeline <br /><br /> `'metrics'` *dict*: metrics from step with keys:  <br />`'avg_lmiss'`, `'avg_imiss'` <br /><br />

### QC.utils

Python Function | Parameters | Returns |
------------ | ------------- | ------------- |
`shell_do(command, log=False, return_log=False)` <br /><br /> Run shell commands from Python  | `command` *str*:  Command to be run in shell. <br /><br /> `log` *str*: Default=False. If True, print stdout <br /><br /> `return_log` *str*: Default=False. if True, return stdout | `stdout` <br /><br /> *datatype dependent on input command*
`merge_genos(geno_path1, geno_path2, out_name)` <br /><br /> Merge 2 Plink Genotypes <br /><br /> *NEEDS TO BE FIXED TO RETURN OUTPUT FILE PATHS AND IMPORTANT METRICS. WILL CHANGE `out_name` to `out_path`* | `geno_path1` *str*:  Path to the Plink genotypes (everything before .bed/.bim/.fam). `geno_path2` *str*:  Path to the Plink genotypes (everything before .bed/.bim/.fam). | `None` <br /><br /> *Currently does not return anything but outputs Plink files to out_name*
`ld_prune(geno_path, out_name, window_size=1000, step_size=50, rsq_thresh=0.05)` <br /><br /> Prune for Linkage Disequilibrium. Produces a subset of markers that are in LD with eachother based on r^2 threshold provided  <br /><br /> *NEEDS TO BE FIXED TO RETURN OUTPUT FILE PATHS AND IMPORTANT METRICS. WILL CHANGE `out_name` to `out_path`* | `geno_path` *str*:  Path to the Plink genotypes (everything before .bed/.bim/.fam). `out_name` *str*:  Path to the output Plink genotypes (everything before .bed/.bim/.fam). <br /><br />`window_size` *int*:  Default=1000. window size in kilo-bases <br /><br />`step_size` *int*:  Default=50. step size in kilo-bases <br /><br />`rsq_threshold` *float*:  Default=o.05. r^2 threshold for inclusion| `None` <br /><br /> *Currently does not return anything but outputs Plink files to out_name*

### GWAS.gwas

Python Function | Parameters | Returns |
------------ | ------------- | ------------- |
`plink_pca(geno_path, out_path, n_pcs=10)` <br /><br /> Principal Component Analysis | `geno_path` *str*:  Path to the Plink genotypes (everything before .bed/.bim/.fam). <br /><br /> `out_path` *str*: Path to the output Plink pca files (everything before .eigenval/.eigenvec). <br /><br /> `n_pcs` *int*: Number of PCs calculated. 10 is the Plink default. | *dict* <br /><br /> `'step'` *str*: Name of step in pipeline <br /><br /> `'metrics'` *dict*: metrics from step with keys: <br />`'n_pcs'` <br /><br /> `'output'` *dict*: paths to output files with keys: <br />`'plink_out'` |
`assoc(geno_path, covar_path, out_path, model)` <br /><br /> Association Analysis | `geno_path` *str*: Path to the Plink genotypes (everything before .bed/.bim/.fam). <br /><br /> `covar_path` *str*: Path to tab-separated covariate file with `#FID IID` in the first two columns (including file extension). <br /><br /> `out_path` *str*: Path to the output Plink assocation files (everything before .PHENO1.glm.logistic/linear). <br /><br /> `model` *str*: Type of association to be run. `logistic` or `linear`. | *dict* <br /><br /> `'step'` *str*: Name of step in pipeline <br /><br /> `'metrics'` *dict*: metrics from step with keys: <br />`'hits'` <br /><br /> `'output'` *dict*: paths to output files with keys: <br />`'pheno_counts'`<br />`'hits'`<br />`'hits_info'`<br />`'plink_out'` |
`prs(geno_path, out_path, assoc, clump_p1=1e-3, clump_r2=0.50, clump_kb=250)` <br /><br /> PRS Analysis with LD Clumping | `geno_path` *str*: Path to the Plink genotypes (everything before .bed/.bim/.fam). <br /><br /> `out_path` *str*: Path to the output Plink clump and PRS files (everything before .\*). <br /><br /> `assoc` *str*: Path to the association file (/path/to/file/name.PHENO1.glm.logistic/linear). <br /><br /> `clump_p1` *float*: Clumping p-value threshold. <br /><br /> `clump_r2` *float*: Clumping r^2 threshold. <br /><br /> `clump_kb` *int*: Clumping kb threshold. | *dict* <br /><br /> `'step'` *str*: Name of step in pipeline <br /><br /> `'metrics'` *dict*: metrics from step with keys: <br />`'clump_pval'`, `'clump_r2'`, `'clump_kb'`, `'num_clumps'` <br /><br /> `'output'` *dict*: paths to output files with keys: <br />`'SNP_weights'`<br />`'clump_SNPs'`<br />`'SNP_pvals'`<br />`'ranges'`<br />`'assoc'`<br />`'plink_out'` |
`calculate_inflation(pval_array, normalize=False, ncases=False, ncontorls=False)` <br /><br /> Lambda/Genomic Inflation Calculation | `pval_array` *numpy array*: P-values from GWAS summary statistics. <br /><br /> `normalize` *bool*: Normalize to 1000 cases and 1000 controls. Recommended if there is a large discrepancy between cases and controls. <br /><br /> `ncases` *int*: Number of cases. Required if normalize=True. <br /><br /> `ncontrols` *int*: Number of controls. Required if normalize=True. | *dict* <br /><br /> `'step'` *str*: Name of step in pipeline <br /><br /> `'metrics'` *dict*: metrics from step with keys: <br /> `'inflation'` |
`munge(geno_path, out_path, assoc, ref_panel)` <br /><br /> Munge Summary Statistics | `geno_path` *str* Path to Plink genotypes (everything before .bed/.bim/.fam). <br /><br /> `out_path` *str*: Path to output Plink frequency report files (everything before .afreq). <br /><br /> `assoc` *str* Path to the association file output by PRS (/path/to/file/name.assoc). <br /><br /> `ref_panel` *str* Path to the Plink-format reference panel used to identify non-rsID SNPs (everything before .bed/.bim/.fam) | *dict* <br /><br /> `'step'` *str*: Name of step in pipeline <br /><br /> `'metrics'` *dict* metrics from step with keys: <br />`'num_snps'` <br /><br /> `'data'` *dict* dataframes from step with keys: <br />`'ma_format_df'`<br />`'coordinates'` <br /><br /> `'output'` *dict*: paths to output files with keys: <br />`'plink_out'` |

### GWAS.utils

Python Function | Parameters | Returns |
------------ | ------------- | ------------- |
`zscore_pval_conversion(zscores=None, pvals=None, stats=None)` <br /><br /> Convert Between Z-score and P-values | `zscores` *numpy array*: Z-scores to be converted to P-values. <br /><br /> `pvals` *numpy array*: P-values to be converted to Z-scores. <br /><br /> `stats` *numpy array*: Summary stats that P-values are based off of (required if converting from P-values to Z-scores). | *numpy array* <br /><br /> Array of either Z-scores or P-values depending on direction of conversion |


## Genotype Calling via Illumina Gencall CLI

Genotypes can be called from .idats in parallel as follows:

`iaap-cli gencall {bpm} {cluster_file} {ped_dir} -f {idat} -p -t 8`

`iaap-cli` can be found in directory: `executables/iaap-cli-linux-x64-1.1.0-sha.80d7e5b3d9c1fdfc2e99b472a90652fd3848bbc7/iaap-cli/`

`bpm` is the Illumina manifest file (.bpm)

`cluster_file` is the clusterfile included with the genotypes (.etg)

`ped_dir` is the directory which the .ped file will be output to 

`-f {idat_dir}` is the directory that contains .idat files. In our case, each directory contains all of the idats for a single chip.

`-p` means "output .ped"

`-t` is number of threads

More information about how the iaap-cli works is in GenoTools/executables/iaap-cli-linux-x64-1.1.0-sha.80d7e5b3d9c1fdfc2e99b472a90652fd3848bbc7 directory.

The way this is used for the GP2 pipeline can be found in the GenoTools/GP2_data_processing/GP2_shulman_processing.ipynb notebook, in which I launch a swarm job in NIH's Biowulf slurm system like so:

```
with open(f'{swarm_scripts_dir}/idat_to_ped.swarm', 'w') as f:
    for code in manifest.SentrixBarcode_A.unique():
        idat_to_ped_cmd = f'\
{iaap} gencall \
{bpm} \
{cluster_file} \
{ped_dir}/ \
-f {idat_dir} \
-p \
-t 8'
        f.write(f'{idat_to_ped_cmd}\n')
f.close()
```

In this example, I write one command per sample (individual .idat) and launch a swarm job

# More Coming Soon!!!!
