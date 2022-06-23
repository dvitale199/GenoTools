import pandas as pd
import os
import argparse
import numpy as np

parser = argparse.ArgumentParser(description='Arguments for snp metric processing')
parser.add_argument('--infile', type=str, default='nope', help='SNP metrics file generated from cnv pipeline [default: nope].')
parser.add_argument('--outfile', type=str, default='nope', help='output path prefix to split per-plate file to per-chrom, per-sample info')

args = parser.parse_args()
infile = args.infile
outfile = args.outfile

df = pd.read_csv(infile,
                 dtype={
                   'chromosome':str,
                   'position':int,
                   'snpID':str,
                   'Sample_ID':str,
                   'Allele1':str,
                   'Allele2':str,
                   'BAlleleFreq':float,
                   'LogRRatio':float,
                   'R':float,
                   'Theta':float,
                   'GenTrain_Score':str,
                   'GType':str}
                )

# work around for 'ASSAY_TYPE=0' string in GenTrain_Score column
df.loc[:,'GenTrain_Score'] = pd.to_numeric(df.GenTrain_Score, errors='coerce')

for iid in df.Sample_ID.unique():
    for chrom in sorted(df.chromosome.unique()):
        
        outfile_name = f'{outfile}_{iid}_chr{chrom}.csv'
        out_df = df.loc[(df.chromosome==chrom) & (df.Sample_ID==iid)]
        out_df.to_csv(outfile_name, header=True, index=False)
