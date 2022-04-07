import os
import subprocess
import numpy as np
from numpy.core.numeric import NaN
import pandas as pd

# Supress copy warning.

pd.options.mode.chained_assignment = None

from QC.utils import shell_do

def idat_snp_metrics(idat_path, bpm, bpm_csv, egt, ref_fasta, out_path, iaap=iaap):
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
    
    # --extra {vcf_path}/gp2_snps_{code}_metadata.tsv | \

    # sort vcf
    sort_cmd = f'\
bcftools \
sort {barcode_out_path}/{barcode}.vcf.gz \
-T {out_tmp}/ \
-Oz -o {barcode_out_path}/{barcode}_sorted.vcf.gz'
# cd {barcode_out_path} && \

    # split indels and snps in vcf
    ext_snps_cmd = f'\
vcftools --gzvcf \
{barcode_out_path}/{barcode}_sorted.vcf.gz \
--remove-indels \
--recode \
--recode-INFO-all \
--out {barcode_out_path}/{barcode}_sorted_snps'
# cd {barcode_out_path} && \


    # split indels and snps in vcf
#     keep_indels_cmd = f'\
# cd {barcode_out_path}; \
# vcftools --gzvcf \
# {barcode}_sorted.vcf.gz \
# --keep-only-indels \
# --recode \
# --recode-INFO-all \
# --out {code}_indels'

    # get snp info from each vcf
    get_logr_baf = f'\
python3 process_vcf_snps.py \
--vcf {barcode_out_path}/{barcode}_sorted_snps.recode.vcf \
--outfile {barcode_out_path}/snp_metrics_{barcode}.csv'


    # get snp info from each vcf
    clean_snp_metrics = f'\
python3 clean_snp_metrics.py \
--infile {barcode_out_path}/snp_metrics_{barcode}.csv \
--outfile {barcode_out_path}/snp_metrics'

    cmds = [idat_to_gtc_cmd, gtc2vcf_cmd, sort_cmd, ext_snps_cmd, get_logr_baf, clean_snp_metrics]
#     cmds = [get_logr_baf, clean_snp_metrics]
#     cmds = [get_logr_baf]
#     cmds = [clean_snp_metrics]
    for cmd in cmds:
        if cmd == gtc2vcf_cmd:
            subprocess.call(cmd, shell=True)
        else:
            shell_do(cmd)
            
 
    # output snp metrics paths
    chroms = [str(i) for i in range(1,23)] + ['X','Y']
    samples = [s.split('/')[-1].replace('_Red.idat','') for s in glob.glob(f'{idat_path}/*_Red.idat')]
    
    outfiles = []
    for sample in samples:
        for chrom in chroms:
            outfile = f'{barcode_out_path}/snp_metrics_{sample}_chr{chrom}.csv'
            outfiles.append(outfile)
            
    return outfiles


# idat_path = '/data/vitaled2/cnv_test/206046180074'
# test_out = '/data/vitaled2/cnv_test'
# # snp_metrics_files = idat_snp_metrics(idat_path, bpm, bpm_csv, egt, ref_fasta, test_out)

# with open(f'{swarm_scripts_dir}/call_cnvs.swarm', 'w') as f:
    
#     for mfile in snp_metrics_files:
#         out_prefix = mfile.replace('.csv','').replace('snp_metrics_', '')
# #         chrom = out_prefix.split('_')[-1]
# #         gene_chrom_list = f'/data/CARD/PD/GP2/ref_panel/glist_hg38_{chrom}.csv'

#         call_cnvs_cmd = f'\
# python3 cnv_gene_caller_alpha.py \
# --infile {mfile} \
# --outfile {out_prefix} \
# --intervals {gene_list}'

#         f.write(f'{call_cnvs_cmd}\n')
# f.close()





def call_cnvs(snp_metrics_file, out_path, intervals_file, min_variants=10, kb_window=100):

    # Load in the data.
    sample_df = pd.read_csv(snp_metrics_file, engine='c')

    temp_interval_df = pd.read_csv(intervals_file, engine='c')
    temp_interval_df.drop_duplicates(subset = ["NAME"], inplace=True, keep='first')
    intervals_df = temp_interval_df[temp_interval_df.CHR.apply(lambda x: x.isnumeric())] # This deletes non-numeric CHRs.

    """# Now reduce just to the intervals of interest and summarize each interval."""

    # Break down L2R and BAF per gene.

    print(f"Remember, we are only calling CNVs for genes with more than {str(min_variants)} variants.")

    results = []

    interval_list = intervals_df['NAME'].unique()

    for INTERVAL in interval_list:
      interval_CHR = int(intervals_df.loc[intervals_df['NAME'] == INTERVAL, 'CHR'].item())
      interval_START_gene = intervals_df.loc[intervals_df['NAME'] == INTERVAL, 'START'].item()
      interval_STOP_gene = intervals_df.loc[intervals_df['NAME'] == INTERVAL, 'STOP'].item()
      interval_START = interval_START_gene - (kb_window*1000)
      interval_STOP = interval_STOP_gene + (kb_window*1000)
      temp_df = sample_df[(sample_df['chromosome'] == interval_CHR) & (sample_df['position'] >= interval_START) & (sample_df['position'] <= interval_STOP)]
      print(f"Working on interval {INTERVAL} on CHR {str(interval_CHR)} from {str(interval_START)} to {str(interval_STOP)} containing {str(temp_df.shape[0])} variants within an window of +/- {str(kb_window)} kb.")
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

#     outpath = outfile.replace('.csv','')
#     output.to_csv(outpath + "-CNVs.csv", index=False)
    output.to_csv(out_path, index=False)
    pd.options.display.max_columns = 10
    print("A summary of your results for this sample is below.")
    print("Thanks for calling CNVs from genotypes with us!")
    desc = output.describe().T
    desc['count'] = desc['count'].astype(int)
    print(desc)

   