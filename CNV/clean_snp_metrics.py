import pandas as pd
import os
import argparse
import numpy as np

parser = argparse.ArgumentParser(description='Arguments for snp metric processing')
parser.add_argument('--infile', type=str, default='nope', help='SNP metrics file generated from cnv pipeline [default: nope].')
parser.add_argument('--outpath', type=str, default='nope', help='output path prefix to split per-plate file to per-chrom, per-sample info')

args = parser.parse_args()
infile = args.infile
outpath = args.outpath

df = pd.read_csv(infile, sep='\t')


cols = ['CHROM','POS','ID','REF','ALT','sampleid','BAF','LRR']

for iid in df.sampleid.unique():
    for chrom in df.CHROM.unique():
        
        outfile = f'{outpath}_{iid}_chr{chrom}.txt' 
        out_df = df.loc[(df.CHROM==chrom) & (df.sampleid==iid)]
        out_df[['CHROM','ID','POS','BAF','LRR']].to_csv(outfile, sep='\t', header=True, index=False)