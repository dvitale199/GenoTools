from CNV.cnv import create_cnv_dosage_matrices
import argparse
import pandas as pd

parser = argparse.ArgumentParser(description='Arguments for Running CNV Pipeline.')    
# parser.add_argument('--metrics_in', type=str, default='Nope.', help='Path to directory containing metrics files. Assumes metrics are stored by SentrixBarcode_A i.e. /path/to/123456789, which contains all metrics under that barcode such as 123456789_R01C01_chr1, etc. all split by chromosome')
parser.add_argument('--files', type=str, default='Nope.', help='df with one column with sample ids and no header')
parser.add_argument('--chrom', type=str, default='Nope.', help='chromosome string.')
parser.add_argument('--label', type=str, default='Nope.', help='ancestry label')
parser.add_argument('--cnv_type', type=str, default='Nope.', help='type of cnv.')
parser.add_argument('--out_path', type=str, default='Nope.', help='Path to prefix for output files.')
args = parser.parse_args()

# metrics_in = args.metrics_in
# chrom = args.chrom
files = args.files 
files_df = pd.read_csv(files, header=None)
files_list = list(files_df.loc[:,0])
chrom = args.chrom
label = args.label
cnv_type = args.cnv_type
# samples = args.samples
# samples_df = pd.read_csv(samples)
# samples_list = list(samples_df.iloc[:,0])
out_path = args.out_path

dosage_ = pd.DataFrame()
for sample_chrom_file in files_list:


    sample = sample_chrom_file.split('/')[-1].replace(f'CNV_{label}_','').replace(f'_chr{chrom}.csv','')
    cnvs = pd.read_csv(sample_chrom_file, usecols=['INTERVAL', cnv_type, 'NUM_VARIANTS'])
    cnvs.loc[:,'sampleid'] = sample
    cnvs_final = cnvs.loc[cnvs.NUM_VARIANTS>=10]
    dosage = cnvs_final.loc[:,[cnv_type,'INTERVAL','sampleid']]
    dosage_pivot = dosage.pivot(index='sampleid', columns='INTERVAL', values=cnv_type)

    dosage_ = pd.concat([dosage_, dosage_pivot])
dosage_.to_csv(out_path)

# dosages = create_cnv_dosage_matrices(in_path=metrics_in, samples_list=samples_list, out_path=out_path)
# dosages = create_cnv_dosage_matrices(in_path=metrics_in, samples_list=samples_list, chromosome=chrom, out_path=out_path)

# baf_df = dosages['baf_df']
# del_df = dosages['l2r_del_df']
# dup_df = dosages['l2r_dup_df']

# baf_df_merged = baf_df.merge(samples, how='left', on='sampleid')
# baf_path = dosages['baf_path']
# del_path = dosages['l2r_del_path']
# dup_path = dosages['l2r_dup_path']
