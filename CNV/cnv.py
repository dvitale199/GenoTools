import os
import subprocess
import numpy as np
from numpy.core.numeric import NaN
import pandas as pd
import glob
import shutil

# Supress copy warning.

pd.options.mode.chained_assignment = None

from QC.utils import shell_do



def get_vcf_names(vcf_path):
    with open(vcf_path, "r") as ifile:
        for line in ifile:
            if line.startswith("#CHROM"):
                vcf_names = [x.strip('\n') for x in line.split('\t')]
                break
    ifile.close()
    return vcf_names


def process_vcf_snps(vcf, out_path):

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

        
def clean_snp_metrics(metrics_in, out_path):
    '''splits snp metrics files by chromosome and individual'''
    
    df = pd.read_csv(metrics_in,
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
                         'Theta':float
                     })

    for iid in df.Sample_ID.unique():
        for chrom in sorted(df.chromosome.unique()):

            outfile_name = f'{out_path}_{iid}_chr{chrom}.csv'
            out_df = df.loc[(df.chromosome==chrom) & (df.Sample_ID==iid)]
            out_df.to_csv(outfile_name, header=True, index=False)
        
        


def idat_snp_metrics(idat_path, bpm, bpm_csv, egt, ref_fasta, out_path, iaap):
    '''
    current structure of idat storage is such that a directory of each SentrixBarcode_A with all idats for that barcode in it
    for ex.
    1112223334
        --> 1112223334_R01C01_Red.idat
        --> 1112223334_R01C01_Grn.idat
        --> 1112223334_R01C02_Red.idat
        --> 1112223334_R01C02_Grn.idat
        etc.
        
    '''
    out_tmp = f'{out_path}/tmp'
    os.makedirs(out_tmp, exist_ok=True)
    barcode = idat_path.split('/')[-1].split('_')[0]
    barcode_out_path = f'{out_path}/{barcode}'
    os.makedirs(barcode_out_path, exist_ok=True)
    
    
    idat_to_gtc_cmd = f'\
{iaap} gencall \
{bpm} \
{egt} \
{barcode_out_path} \
-f {idat_path} \
-g \
-t 8'

    # export path to plugins temporarily for biowulf. will figure this out later
    gtc2vcf_cmd = f'\
export BCFTOOLS_PLUGINS="/data/vitaled2/bin"; \
bcftools +gtc2vcf \
--no-version -Ob \
--bpm {bpm} \
--csv {bpm_csv} \
--egt {egt} \
--gtcs {barcode_out_path} \
--fasta-ref {ref_fasta} | \
bcftools norm --no-version -Oz -c w -f {ref_fasta} > {barcode_out_path}/{barcode}.vcf.gz'
# use --extra to output .tsv of other info from gtc

    # sort vcf
    sort_cmd = f'\
bcftools \
sort {barcode_out_path}/{barcode}.vcf.gz \
-T {out_tmp}/ \
-Oz -o {barcode_out_path}/{barcode}_sorted.vcf.gz'

    # split indels and snps in vcf
    ext_snps_cmd = f'\
vcftools --gzvcf \
{barcode_out_path}/{barcode}_sorted.vcf.gz \
--remove-indels \
--recode \
--recode-INFO-all \
--out {barcode_out_path}/{barcode}_sorted_snps'


#     split indels and snps in vcf
# can bring this in later if needed. for now, only snps
#     keep_indels_cmd = f'\
# vcftools --gzvcf \
# {barcode_out_path}/{barcode}_sorted.vcf.gz \
# --keep-only-indels \
# --recode \
# --recode-INFO-all \
# --out {barcode_out_path}/{barcode}_sorted_snps'

    cmds = [idat_to_gtc_cmd, gtc2vcf_cmd, sort_cmd, ext_snps_cmd]
    for cmd in cmds:
        if cmd == gtc2vcf_cmd:
            subprocess.call(cmd, shell=True)
        else:
            shell_do(cmd)
            
    # get snp info from each vcf
    vcf_in = f'{barcode_out_path}/{barcode}_sorted_snps.recode.vcf'
    snp_metrics_out = f'{barcode_out_path}/snp_metrics_{barcode}.csv'
    process_vcf_snps(vcf=vcf_in, out_path=snp_metrics_out)
    
    metrics_in = f'{barcode_out_path}/snp_metrics_{barcode}.csv'
    clean_metrics_out = f'{barcode_out_path}/snp_metrics'
    clean_snp_metrics(metrics_in=metrics_in, out_path=clean_metrics_out)
 
    # output snp metrics paths
    chroms = [str(i) for i in range(1,23)] + ['X','Y']
    samples = [s.split('/')[-1].replace('_Red.idat','') for s in glob.glob(f'{idat_path}/*_Red.idat')]
    
    outfiles = []
    for sample in samples:
        for chrom in chroms:
            outfile = f'{barcode_out_path}/snp_metrics_{sample}_chr{chrom}.csv'
            outfiles.append(outfile)
            
    out_df = pd.DataFrame({'path':outfiles})
    out_df.to_csv(f'{barcode_out_path}/snp_metrics_{barcode}.log', header=False, index=False)
            
    return outfiles


def call_cnvs(snp_metrics_file, out_path, intervals_file, min_variants=10, kb_window=100):

    # Load in the data.
    sample_df = pd.read_csv(snp_metrics_file, engine='c')

    temp_interval_df = pd.read_csv(intervals_file, engine='c')
    temp_interval_df.drop_duplicates(subset = ["NAME"], inplace=True, keep='first')
    intervals_df = temp_interval_df.copy()
#     intervals_df = temp_interval_df[temp_interval_df.CHR.apply(lambda x: x.isnumeric())] # This deletes non-numeric CHRs.

    """# Now reduce just to the intervals of interest and summarize each interval."""

    # Break down L2R and BAF per gene.

#     print(f"Remember, we are only calling CNVs for genes with more than {str(min_variants)} variants.")

    results = []

    interval_list = intervals_df['NAME'].unique()

    for INTERVAL in interval_list:
      interval_CHR = int(intervals_df.loc[intervals_df['NAME'] == INTERVAL, 'CHR'].item())
      interval_START_gene = intervals_df.loc[intervals_df['NAME'] == INTERVAL, 'START'].item()
      interval_STOP_gene = intervals_df.loc[intervals_df['NAME'] == INTERVAL, 'STOP'].item()
      interval_START = interval_START_gene - (kb_window*1000)
      interval_STOP = interval_STOP_gene + (kb_window*1000)
      temp_df = sample_df[(sample_df['chromosome'] == interval_CHR) & (sample_df['position'] >= interval_START) & (sample_df['position'] <= interval_STOP)]
#       print(f"Working on interval {INTERVAL} on CHR {str(interval_CHR)} from {str(interval_START)} to {str(interval_STOP)} containing {str(temp_df.shape[0])} variants within an window of +/- {str(kb_window)} kb.")
      if temp_df.shape[0] < min_variants:
#         print("This interval does not meet the minimum variant count requirement.")
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
    output.to_csv(out_path, index=False)
    pd.options.display.max_columns = 10
#     print("A summary of your results for this sample is below.")
#     print("Thanks for calling CNVs from genotypes with us!")
#     desc = output.describe().T
#     desc['count'] = desc['count'].astype(int)
#     print(desc)

   