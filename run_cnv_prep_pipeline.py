import argparse
import pandas as pd

parser = argparse.ArgumentParser(description='Arguments for Running CNV Pipeline.')    
parser.add_argument('--idat_path', type=str, default='Nope.', help='Path to idats. terrible method. will change soon')
parser.add_argument('--sample', type=str, default="Nope.", help='sampleid in sentrixbarcode_sentrixposition format')
parser.add_argument('--bim', type=str, default="Nope.", help='path to bim file containing snps of interest')
parser.add_argument('--out_path', type=str, default='Nope.', help='Path to output snp metrics per sample report.')

args = parser.parse_args()

idat_path = args.idat_path
sample = args.sample
code = sample.split('_')[0]
out_path = args.out_path
bim_path = args.bim
bim = pd.read_csv(bim_path, sep='\s+', header=None, names=['chr','id','pos','bp','a1','a2'], usecols=['id'])

chroms = [str(i) for i in range(1,23)]
total_sample_df = pd.DataFrame()

for chrom in chroms:
    
    mfile = f'{idat_path}/{code}/snp_metrics_{sample}_chr{chrom}.csv'
    metrics = pd.read_csv(mfile)
    metrics_tmp = metrics.loc[metrics.snpID.isin(bim.id)]
    total_sample_df = pd.concat([total_sample_df,metrics_tmp], ignore_index=True)

total_sample_df.to_csv(out_path, index=False, header=True)