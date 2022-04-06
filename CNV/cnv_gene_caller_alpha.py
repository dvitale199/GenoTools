# -*- coding: utf-8 -*-
"""# Intro + imports + args.

Simple summary.

Grab a CSV that contains the right data you need exported from GenomeStudio ro similar.

Call the CNVs on a gene by gene level.
"""

# Imports here.

import numpy as np
from numpy.core.numeric import NaN
import pandas as pd
import argparse

# Supress copy warning.

pd.options.mode.chained_assignment = None

# Args incoming.

parser = argparse.ArgumentParser(description='Arguments for volcano plotting.')    
parser.add_argument('--infile', type=str, default='input.csv', help='CSV file with data to analyze for one sample. Header is [chromosome,position,snpID,snpName,Sample_ID,Allele1,Allele2,BAlleleFreq,LogRRatio]')
parser.add_argument('--outfile', type=str, default='./test_sample', help='Path and prefix to the sample you are analyzing usually, this will automatically have *.csv appended to it.')
parser.add_argument('--intervals', type=str, default="intervals.csv", help='Gene or other feature intervals to analyze. Header is [NAME,CHR,START,STOP], one line per itnervals.Autosomes only.')
parser.add_argument('--min_variants', type=int, default=10, help='Minimum number of variants to run the CNV algorithm on per gene.')
parser.add_argument('--kb_window', type=int, default=100, help='Kilobase window around each interval, a value of 100 would mean +/- 100kb.')
args = parser.parse_args()

# Options here.

infile = args.infile
outfile = args.outfile
intervals = args.intervals
min_variants = args.min_variants
kb_window = args.kb_window

# Load in the data.

sample_df = pd.read_csv(infile, engine='c')

temp_interval_df = pd.read_csv(intervals, engine='c')
temp_interval_df.drop_duplicates(subset = ["NAME"], inplace=True, keep='first')
intervals_df = temp_interval_df[temp_interval_df.CHR.apply(lambda x: x.isnumeric())] # This deletes non-numeric CHRs.

"""# Now reduce just to the intervals of interest and summarize each interval."""

# Break down L2R and BAF per gene.

print("Remember, we are only calling CNVs for genes with more than " + str(min_variants) + " variants.")

results = []

interval_list = intervals_df['NAME'].unique()

for INTERVAL in interval_list:
  interval_CHR = int(intervals_df.loc[intervals_df['NAME'] == INTERVAL, 'CHR'].item())
  interval_START_gene = intervals_df.loc[intervals_df['NAME'] == INTERVAL, 'START'].item()
  interval_STOP_gene = intervals_df.loc[intervals_df['NAME'] == INTERVAL, 'STOP'].item()
  interval_START = interval_START_gene - (kb_window*1000)
  interval_STOP = interval_STOP_gene + (kb_window*1000)
  temp_df = sample_df[(sample_df['chromosome'] == interval_CHR) & (sample_df['position'] >= interval_START) & (sample_df['position'] <= interval_STOP)]
  print("Working on interval " + INTERVAL + " on CHR " + str(interval_CHR) + " from " + str(interval_START) + " to " + str(interval_STOP) + " containing " + str(temp_df.shape[0]) + " variants within an window of +/- " + str(kb_window) + "kb.")
  if temp_df.shape[0] < min_variants:
    print("This interval does not meet the minimum variant count requirement.")
    results.append((INTERVAL, temp_df.shape[0], NaN, NaN, NaN, interval_START, interval_START_gene, interval_STOP_gene, interval_STOP))
  else:
    temp_df['BAF_insertion'] = np.where( (temp_df['BAlleleFreq'].between(0.65, 0.85, inclusive=False)) | (temp_df['BAlleleFreq'].between(0.15, 0.35, inclusive=False)), 1, 0)
    temp_df['L2R_deletion'] = np.where( temp_df['LogRRatio'] < -0.2, 1, 0)
    temp_df['L2R_insertion'] = np.where( temp_df['LogRRatio'] > 0.2, 1, 0)
    PERCENT_BAF_INSERTION = temp_df['BAF_insertion'].mean()
    PERCENT_L2R_DELETION = temp_df['L2R_deletion'].mean()
    PERCENT_L2R_INSERTION = temp_df['L2R_insertion'].mean()
    results.append((INTERVAL, temp_df.shape[0], PERCENT_BAF_INSERTION, PERCENT_L2R_DELETION, PERCENT_L2R_INSERTION, interval_START, interval_START_gene, interval_STOP_gene, interval_STOP))

output = pd.DataFrame(results, columns=('INTERVAL', 'NUM_VARIANTS', 'PERCENT_BAF_INSERTION', 'PERCENT_L2R_DELETION','PERCENT_L2R_DUPLICATION','START_PLUS_WINDOW','START','STOP','STOP_PLUS_WINDOW'))

outpath = outfile.replace('.csv','')
output.to_csv(outpath + "-CNVs.csv", index=False)
pd.options.display.max_columns = 10
print("A summary of your results for this sample is below.")
print("Thanks for calling CNVs from genotypes with us!")
desc = output.describe().T
desc['count'] = desc['count'].astype(int)
print(desc)
