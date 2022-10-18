import pandas as pd
import argparse
import numpy as np


parser = argparse.ArgumentParser(description='Arguments for Running CNV Pipeline.')    
# parser.add_argument('--metrics_in', type=str, default='Nope.', help='Path to directory containing metrics files. Assumes metrics are stored by SentrixBarcode_A i.e. /path/to/123456789, which contains all metrics under that barcode such as 123456789_R01C01_chr1, etc. all split by chromosome')
parser.add_argument('--dosagefile', type=str, default='Nope.', help='dosages in. sample ids in first column, each subsequent column is a gene. values are dosages')
parser.add_argument('--key', type=str, default='Nope.', help='key containing sample ids and phenotypes')
parser.add_argument('--out_path', type=str, default='Nope.', help='Path to prefix for output files.')
parser.add_argument('--pheno_out', type=str, default='Nope.', help='Path to phenotype outfile.')
args = parser.parse_args()

# dosagefile = f'{cnv_path}/CNV_{label}_chr{chrom}_{cnv_type}.csv'
# dosagefile_out = f'{cnv_path}/CNV_{label}_chr{chrom}_{cnv_type}_gp2ids.csv'
# pheno_out = f'{cnv_path}/GP2_round2_{label}_chr{chrom}_{cnv_type}.pheno'
dosagefile = args.dosagefile
key_path = args.key
key = pd.read_csv(key_path)
dosagefile_out = args.out_path
pheno_out = args.pheno_out



dosage = pd.read_csv(dosagefile)
dosage_merge = dosage.merge(key[['GP2sampleID','IID']], left_on='sampleid', right_on='IID')
dosage_out = dosage_merge.drop(columns=['sampleid','IID']).set_index('GP2sampleID').reset_index().rename(columns={'GP2sampleID':'sampleid'})

dosage_pheno = dosage_out.merge(key.loc[:,['GP2sampleID','pheno']], left_on='sampleid', right_on='GP2sampleID', how='left')
dosage_pheno2 = dosage_pheno.drop(columns=["GP2sampleID"])
dosage_pheno_out = dosage_pheno2.loc[dosage_pheno2.pheno != -9]
dosage_pheno_out.loc[:,'pheno'] = np.where(dosage_pheno_out.pheno == 1, 0, 1)
dosage_pheno_out[['sampleid','pheno']].to_csv(pheno_out, sep='\t', header=True, index=False)

dosage_final = dosage_out.loc[dosage_out.sampleid.isin(dosage_pheno_out.sampleid)]
dosage_final.to_csv(dosagefile_out, header=True, index=False)