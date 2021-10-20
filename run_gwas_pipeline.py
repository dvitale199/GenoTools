import pandas as pd
import argparse
import os

# local imports
from GWAS.gwas import *
from GWAS.utils import *

parser = argparse.ArgumentParser(description='Arguments for GWAS and PRS (data in Plink .bim/.bam/.fam format)')
parser.add_argument('--geno', type=str, default='nope', help='Genotype: (string file path). Path to QC\'d and imputed PLINK format genotype file, everything before the *.bed/bim/fam [defualt: nope].')
parser.add_argument('--cov', type=str, default='nope', help='tab-seperated plink-style covariate file with FID IID as first two columns.')
parser.add_argument('--model', type=str, default='nope', help='Association being run: logistic/linear')
parser.add_argument('--ref_panel', type=str, default='nope', help='Path for PLINK format reference panel file, everything before the *.bed/bim/fam')
parser.add_argument('--out', type=str, default='nope', help='Prefix for output (indcluding path)')

args = parser.parse_args()

geno_path = args.geno
covar_path = args.cov
model = args.model
ref_panel = args.ref_panel
out_path = args.out

# generate pcs
pca = plink_pca(geno_path, out_path)

# if no covariate file provided, set covar path to PCA result
if covar_path == 'nope':
    covar_path = f'{out_path}.eigenvec'

# otherwise, merge pcs and covariates
else:
    # read covariates - need to have header row
    covars = pd.read_csv(covar_path, 
                         sep='\s+',
                         dtype={'#FID':str, 'IID':str})
    
    # read pcs
    pcs = pd.read_csv(f'{out_path}.eigenvec', 
                      sep='\s+',
                      dtype={'#FID':str,'IID':str})
    
    # merge 
    covar_merged = pcs.merge(covars, how='inner', on=['#FID','IID'])
    covar_merged.to_csv(f'{out_path}.cov', sep='\t', header=True, index=False)

    # set new covar path
    covar_path = f'{out_path}.cov'
    
# run association
glm = assoc(geno_path, covar_path, out_path, model.lower())

# grab full association file
assoc_file = f"{glm['output']['plink_out']}.PHENO1.glm.{model}"

# if full association file doesn't exist then try to grab hybrid file
if not os.path.isfile(assoc_file):
    assoc_file = f'{assoc_file}.hybrid'

# calculate inflation
assoc_file_df = pd.read_csv(assoc_file, sep='\s+')
inflation = calculate_inflation(assoc_file_df['P'])

# run prs
score = prs(geno_path, out_path, assoc_file)

# grab clean association file for munging
clean_assoc_file = f"{score['output']['assoc']}"

# munge summary stats
stats = munge(geno_path, out_path, clean_assoc_file, ref_panel)

# creating metrics dataframe
metrics_df = pd.DataFrame()

steps = [pca, glm, inflation, score, stats]

for item in steps:
    step = item['step']
    pf = item['pass']

    for metric, value in item['metrics'].items():
        tmp_metrics_df = pd.DataFrame({'step':[step], 'metric':[metric], 'value': [value], 'pass':[pf]})
        metrics_df = metrics_df.append(tmp_metrics_df)

metrics_df.reset_index(drop=True, inplace=True)

# reading score reports to output to .h5 file
s1 = pd.read_csv(f'{out_path}.PRS.s1.sscore', sep='\s+')
s2 = pd.read_csv(f'{out_path}.PRS.s2.sscore', sep='\s+')
s3 = pd.read_csv(f'{out_path}.PRS.s3.sscore', sep='\s+')

# calculate z-scores
for scores in [s1, s2, s3]:
    scores_mean = scores['SCORE1_AVG'].mean()
    scores_std = scores['SCORE1_AVG'].std()
    scores['Z'] = (scores['SCORE1_AVG']-scores_mean)/scores_std
    scores['P'] = zscore_pval_conversion(scores['Z'])


# merging ma format data and coordinates to output to .h5 file
ma_df = stats['data']['ma_format_df']
coords_df = stats['data']['coordinates']

# build output hdf
metrics_outfile = f'{out_path}.GWAS.metrics.h5'

metrics_df.to_hdf(metrics_outfile, key='GWAS', mode='w')
ma_df.to_hdf(metrics_outfile, key='ma_output')
coords_df.to_hdf(metrics_outfile, key='coordinates')
s1.to_hdf(metrics_outfile, key='score_1')
s2.to_hdf(metrics_outfile, key='score_2')
s3.to_hdf(metrics_outfile, key='score_3')