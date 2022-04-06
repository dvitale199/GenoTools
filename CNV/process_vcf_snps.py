import pandas as pd
import argparse
import numpy as np

parser = argparse.ArgumentParser(description='Arguments for vcf processing (non-gz vcf file)')
parser.add_argument('--vcf', type=str, default='nope', help='Genotype: (string file path). Path to vcf format genotype file, NOT ZIPPED [default: nope].')
# parser.add_argument('--gene_ref', type=str, default='nope', help='plink-format gene reference file: https://www.cog-genomics.org/plink/1.9/resources')
parser.add_argument('--outfile', type=str, default='nope', help='output path to csv file containing per-snp, per-sample info')

args = parser.parse_args()
vcf = args.vcf
# gene_ref = args.gene_ref
out_path = args.outfile

def get_vcf_names(vcf_path):
    with open(vcf_path, "r") as ifile:
        for line in ifile:
            if line.startswith("#CHROM"):
                vcf_names = [x.strip('\n') for x in line.split('\t')]
                break
    ifile.close()
    return vcf_names

# gene_list = pd.read_csv(gene_ref, sep='\s+', header=None, names=['chr','start','end','symbol'], dtype={'chr':str,'start':int,'end':int})

# out_colnames = ['CHROM','POS','ID','REF','ALT','sampleid','BAF','LRR']
out_colnames = ['chromosome', 'position', 'snpID', 'Sample_ID', 'Allele1', 'Allele2', 'BAlleleFreq', 'LogRRatio', 'R', 'THETA']

variant_metrics_out_df = pd.DataFrame(columns=out_colnames)
variant_metrics_out_df.to_csv(out_path, header=True, index=False)


names = get_vcf_names(vcf)        
vcf = pd.read_csv(vcf, comment='#', chunksize=10000, delim_whitespace=True, header=None, names=names, dtype={'#CHROM':str})
IIDs = [x for x in names if x not in ['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT']]

for chunk in vcf:

    chunk.rename(columns={'#CHROM':'CHROM'}, inplace=True)
    chunk_melt = chunk.melt(id_vars=['CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO'], value_vars=IIDs, value_name='metrics')
    chunk_melt[['GT','GQ','IGC','BAF','LRR','NORMX','NORMY','R','THETA','X','Y']] = chunk_melt.metrics.str.split(':', expand=True)
    chunk_melt.drop(columns=['QUAL','FILTER','INFO','GT','GQ','IGC','NORMX','NORMY','X','Y','metrics'], inplace=True)
    chunk_melt.rename(columns={'variable':'sampleid'}, inplace=True)
#     print(chunk_melt)
    chunk_melt.loc[:,'CHROM'] = chunk_melt['CHROM'].astype(str).str.replace('chr','')
    chunk_final = chunk_melt.loc[:,['CHROM','POS','ID','sampleid','REF','ALT','BAF','LRR', 'R', 'THETA']]
    chunk_final.columns = ['chromosome', 'position', 'snpID', 'Sample_ID', 'Allele1', 'Allele2', 'BAlleleFreq', 'LogRRatio', 'R', 'Theta']
    
    chunk_final.to_csv(out_path, header=False, index=False, mode='a')

    
    
#     chunk_melt.loc[:,'gene'] = np.nan
    
#     # now append gene name
#     for i, gene in enumerate(gene_list.symbol):

#         chrom = gene_list.loc[i,'chr']
#         start = gene_list.loc[i,'start'] - 100000
#         end = gene_list.loc[i, 'end'] + 100000
#         chunk_melt.loc[((chunk_melt.CHROM==chrom) & (chunk_melt.POS>=start) & (chunk_melt.POS<=end)), 'gene'] = gene
    
    
    
