import os
import subprocess
import numpy as np
from numpy.core.numeric import NaN
import pandas as pd
import glob
import shutil
from sklearn.preprocessing import MinMaxScaler
import statsmodels.api as sm

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
    out_colnames = ['chromosome', 'position', 'snpID', 'Sample_ID', 'Ref', 'Alt','ALLELE_A','ALLELE_B', 'BAlleleFreq', 'LogRRatio', 'R', 'Theta', 'GenTrain_Score', 'GType']

    variant_metrics_out_df = pd.DataFrame(columns=out_colnames)
    variant_metrics_out_df.to_csv(out_path, header=True, index=False)


    names = get_vcf_names(vcf)        
    vcf = pd.read_csv(vcf, comment='#', chunksize=10000, delim_whitespace=True, header=None, names=names, dtype={'#CHROM':str})
    IIDs = [x for x in names if x not in ['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT']]

    for chunk in vcf:
        chunk.rename(columns={'#CHROM':'CHROM'}, inplace=True)
        chunk_melt = chunk.melt(id_vars=['CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO'], value_vars=IIDs, value_name='metrics')
        chunk_melt[['GT','GQ','IGC','BAF','LRR','NORMX','NORMY','R','THETA','X','Y']] = chunk_melt.metrics.str.split(':', expand=True)
        chunk_melt.loc[:,'GenTrain_Score'] = chunk_melt.INFO.str.split(';',expand=True).iloc[:,10].str.replace('GenTrain_Score=','')
        chunk_melt.loc[:,'ALLELE_A'] = chunk_melt.INFO.str.split(';',expand=True).iloc[:,1].str.replace('ALLELE_A=','')
        chunk_melt.loc[:,'ALLELE_B'] = chunk_melt.INFO.str.split(';',expand=True).iloc[:,2].str.replace('ALLELE_B=','')
        chunk_melt.drop(columns=['QUAL','FILTER','INFO','GQ','IGC','NORMX','NORMY','X','Y','metrics'], inplace=True)
        chunk_melt.rename(columns={'variable':'sampleid'}, inplace=True)
        chunk_melt.loc[:,'CHROM'] = chunk_melt['CHROM'].astype(str).str.replace('chr','')
        chunk_final = chunk_melt.loc[:,['CHROM','POS','ID','sampleid','REF','ALT','GT','ALLELE_A','ALLELE_B','BAF','LRR', 'R', 'THETA', 'GenTrain_Score']]
        gtype_map = {'0/0':'AA', '0/1':'AB', '1/1':'BB', './.':'NC'}
        
        chunk_final.loc[:,'GType'] = chunk_final['GT'].map(gtype_map)
        chunk_final.drop(columns=['GT'], inplace=True)
        chunk_final.columns = ['chromosome', 'position', 'snpID', 'Sample_ID', 'Ref', 'Alt','ALLELE_A','ALLELE_B', 'BAlleleFreq', 'LogRRatio', 'R', 'Theta', 'GenTrain_Score', 'GType']
        

#         chunk.rename(columns={'#CHROM':'CHROM'}, inplace=True)
#         chunk_melt = chunk.melt(id_vars=['CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO'], value_vars=IIDs, value_name='metrics')
#         chunk_melt[['GT','GQ','IGC','BAF','LRR','NORMX','NORMY','R','THETA','X','Y']] = chunk_melt.metrics.str.split(':', expand=True)
#         chunk_melt.drop(columns=['QUAL','FILTER','INFO','GQ','IGC','NORMX','NORMY','X','Y','metrics'], inplace=True)
#         chunk_melt.rename(columns={'variable':'sampleid'}, inplace=True)
#     #     print(chunk_melt)
#         chunk_melt.loc[:,'CHROM'] = chunk_melt['CHROM'].astype(str).str.replace('chr','')
#         chunk_final = chunk_melt.loc[:,['CHROM','POS','ID','sampleid','REF','ALT','GT','BAF','LRR', 'R', 'THETA']]
#         chunk_final.columns = ['chromosome', 'position', 'snpID', 'Sample_ID', 'Allele1', 'Allele2', 'GT', 'BAlleleFreq', 'LogRRatio', 'R', 'Theta']

        chunk_final.to_csv(out_path, header=False, index=False, mode='a')

        
def calculate_maf(gtype_df):
    
    gtypes_map = {
        'AA': 0,
        'AB': 1,
        'BA': 1,
        'BB': 2,
        'NC': np.nan
        }

    gtypes = gtype_df.pivot(index='snpID', columns='Sample_ID', values='GT').replace(gtypes_map)
    
    # count only called genotypes
    N = gtypes.shape[1]-gtypes.isna().sum(axis=1)
    freq = pd.DataFrame({'freq': gtypes.sum(axis=1)/(2*N)})
    freq.loc[:,'maf'] = np.where(freq < 0.5, freq, 1-freq)
    maf_out = freq.drop(columns=['freq']).reset_index()

    return maf_out


# fix clean_snp_metrics()
def clean_snp_metrics(metrics_in, out_path):
    '''splits snp metrics files by chromosome and individual'''
    
    snp_metrics = pd.read_csv(metrics_in,
                     dtype={
                         'chromosome':str,
                         'position':int,
                         'snpID':str,
                         'Sample_ID':str,
                         'Ref':str,
                         'Alt':str,
                         'ALLELE_A':int,
                         'ALLELE_B':int,
                         'BAlleleFreq':float,
                         'LogRRatio':float,
                         'R':float,
                         'Theta':float,
                         'GenTrain_Score':float,
                         'GType':str
                     })
    
    
    alt_split = snp_metrics.loc[:,'Alt'].str.split(',', expand=True)
    snp_metrics.loc[:,'Alt1'], snp_metrics.loc[:,'Alt2'] = alt_split.loc[:,0], alt_split.loc[:,1]

    snp_metrics.loc[(snp_metrics['GType']=='AA') & (snp_metrics['ALLELE_A']==1), 'GT'] = 'BB'
    snp_metrics.loc[(snp_metrics['GType']=='AA') & (snp_metrics['ALLELE_A']==0), 'GT'] = 'AA'
    snp_metrics.loc[(snp_metrics['GType']=='BB') & (snp_metrics['ALLELE_B']==1), 'GT'] = 'BB'
    snp_metrics.loc[(snp_metrics['GType']=='BB') & (snp_metrics['ALLELE_B']==0), 'GT'] = 'AA'

    snp_metrics.loc[(snp_metrics['GType']=='AB'), 'GT'] = 'AB'
    snp_metrics.loc[(snp_metrics['GType']=='NC'), 'GT'] = 'NC'
    snp_metrics.loc[:,'GT'] = snp_metrics.loc[:,'GT'].fillna('NC')

    # drop snps where gentrain score, theta, and r isna
#     snp_metrics = snp_metrics.loc[(~snp_metrics['GenTrain_Score'].isna()) & (~snp_metrics['Theta'].isna()) & (~snp_metrics['R'].isna())]

    snp_metrics.loc[snp_metrics['ALLELE_A']==0, 'a1'] = snp_metrics.loc[snp_metrics['ALLELE_A']==0,'Ref']
    snp_metrics.loc[snp_metrics['ALLELE_A']==1, 'a1'] = snp_metrics.loc[snp_metrics['ALLELE_A']==1,'Alt1']
    snp_metrics.loc[snp_metrics['ALLELE_B']==0, 'a2'] = snp_metrics.loc[snp_metrics['ALLELE_B']==0,'Ref']
    snp_metrics.loc[snp_metrics['ALLELE_B']==1, 'a2'] = snp_metrics.loc[snp_metrics['ALLELE_B']==1,'Alt1']
    snp_metrics.loc[snp_metrics['ALLELE_B']==2, 'a2'] = snp_metrics.loc[snp_metrics['ALLELE_B']==2,'Alt2']
    
    # calculate maf for full 
    maf_df = calculate_maf(snp_metrics)
    snp_metrics_full = snp_metrics.merge(maf_df, how='left', on='snpID')
        
    
    # output metrics file per sample, per chrom
    for iid in snp_metrics_full.Sample_ID.unique():
        for chrom in sorted(snp_metrics_full.chromosome.unique()):

            outfile_name = f'{out_path}_{iid}_chr{chrom}.csv'
            out_df = snp_metrics_full.loc[(snp_metrics_full.chromosome==chrom) & (snp_metrics_full.Sample_ID==iid)]
            out_df.to_csv(outfile_name, header=True, index=False)
            


def idat_snp_metrics(idat_path, bpm, bpm_csv, egt, ref_fasta, out_path, iaap, bcftools_plugins_path="/data/vitaled2/bin"):
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
export BCFTOOLS_PLUGINS={bcftools_plugins_path}; \
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


    """# Now reduce just to the intervals of interest and summarize each interval."""

    # Break down L2R and BAF per gene.

    results = []

    interval_list = intervals_df['NAME'].unique()

    for INTERVAL in interval_list:
      interval_CHR = intervals_df.loc[intervals_df['NAME'] == INTERVAL, 'CHR'].item()
      interval_START_gene = intervals_df.loc[intervals_df['NAME'] == INTERVAL, 'START'].item()
      interval_STOP_gene = intervals_df.loc[intervals_df['NAME'] == INTERVAL, 'STOP'].item()
      interval_START = interval_START_gene - (kb_window*1000)
      interval_STOP = interval_STOP_gene + (kb_window*1000)
      temp_df = sample_df[(sample_df['chromosome'] == interval_CHR) & (sample_df['position'] >= interval_START) & (sample_df['position'] <= interval_STOP)]

      if temp_df.shape[0] < min_variants:

        results.append((INTERVAL, temp_df.shape[0], NaN, NaN, NaN, interval_START, interval_START_gene, interval_STOP_gene, interval_STOP))
      else:
        temp_df['BAF_insertion'] = np.where( (temp_df['BAlleleFreq'].between(0.65, 0.85, inclusive='neither')) | (temp_df['BAlleleFreq'].between(0.15, 0.35, inclusive='neither')), 1, 0)
        temp_df['L2R_deletion'] = np.where( temp_df['LogRRatio'] < -0.2, 1, 0)
        temp_df['L2R_insertion'] = np.where( temp_df['LogRRatio'] > 0.2, 1, 0)
        PERCENT_BAF_INSERTION = temp_df['BAF_insertion'].mean()
        PERCENT_L2R_DELETION = temp_df['L2R_deletion'].mean()
        PERCENT_L2R_INSERTION = temp_df['L2R_insertion'].mean()
        results.append((INTERVAL, temp_df.shape[0], PERCENT_BAF_INSERTION, PERCENT_L2R_DELETION, PERCENT_L2R_INSERTION, interval_START, interval_START_gene, interval_STOP_gene, interval_STOP))

    output = pd.DataFrame(results, columns=('INTERVAL', 'NUM_VARIANTS', 'PERCENT_BAF_INSERTION', 'PERCENT_L2R_DELETION','PERCENT_L2R_DUPLICATION','START_PLUS_WINDOW','START','STOP','STOP_PLUS_WINDOW'))
    output.to_csv(out_path, index=False)
    

def create_cnv_dosage_matrices(in_path, samples_list, out_path):
    
#     chrom = str(chromosome)
    
    baf_out = f'{out_path}_BAF.csv'
    l2r_del_out = f'{out_path}_L2R_DEL.csv'
    l2r_dup_out = f'{out_path}_L2R_DUP.csv'
    
    baf_ = pd.DataFrame()
    l2r_del_ = pd.DataFrame()
    l2r_dup_ = pd.DataFrame()

    
    for sample in samples_list:
        code = sample.split('_')[0]

        cnv_file = f'{in_path}/CNV_{sample}.csv'
        cnvs = pd.read_csv(cnv_file)
        cnvs.loc[:,'sampleid'] = sample
        cnvs_final = cnvs.loc[cnvs.NUM_VARIANTS>=10]
        
        baf = cnvs_final.loc[:,['PERCENT_BAF_INSERTION','INTERVAL','sampleid']]
        baf_pivot = baf.pivot(index='sampleid',columns='INTERVAL',values='PERCENT_BAF_INSERTION')
        
        l2r_del = cnvs_final.loc[:,['PERCENT_L2R_DELETION', 'INTERVAL', 'sampleid']]
        l2r_del_pivot = l2r_del.pivot(index='sampleid', columns='INTERVAL', values='PERCENT_L2R_DELETION')
        
        l2r_dup = cnvs_final.loc[:,['PERCENT_L2R_DUPLICATION', 'INTERVAL', 'sampleid']]
        l2r_dup_pivot = l2r_dup.pivot(index='sampleid',columns='INTERVAL',values='PERCENT_L2R_DUPLICATION')
        
        baf_ = pd.concat([baf_, baf_pivot])
        l2r_del_ = pd.concat([l2r_del_, l2r_del_pivot])
        l2r_dup_ = pd.concat([l2r_dup_, l2r_dup_pivot])
        
#         baf_ = baf_.append(baf_pivot)
#         l2r_del_ = l2r_del_.append(l2r_del_pivot)
#         l2r_dup_ =l2r_dup_.append(l2r_dup_pivot)

    baf_.columns = [x.replace('-','_') for x in baf_.columns]
    l2r_del_.columns = [x.replace('-','_') for x in l2r_del_.columns]
    l2r_dup_.columns = [x.replace('-','_') for x in l2r_dup_.columns]
    baf_.columns = [x.replace('.','_') for x in baf_.columns]
    l2r_del_.columns = [x.replace('.','_') for x in l2r_del_.columns]
    l2r_dup_.columns = [x.replace('.','_') for x in l2r_dup_.columns]
    
    baf_.to_csv(baf_out, index=True, header=True)
    l2r_del_.to_csv(l2r_del_out, index=True, header=True)
    l2r_dup_.to_csv(l2r_dup_out, index=True, header=True)
    
    out_dict = {
        'baf_df': baf_,
        'l2r_del_df': l2r_del_,
        'l2r_dup_df': l2r_dup_,
        'baf_path': baf_out,
        'l2r_del_path': l2r_del_out,
        'l2r_dup_path': l2r_dup_out
                 }
    
    return out_dict
 
    
def CNV_WAS(cnv_dosage_file, pheno, covar, out_path):
    scaler = MinMaxScaler()
    dosage_df = pd.read_csv(cnv_dosage_file)
    pheno_df = pd.read_csv(pheno, sep='\t')
    covar_df = pd.read_csv(covar, sep='\t')


    if covar_df.age_of_onset.isna().all():
        covar_df.drop(columns=['age_of_onset'], inplace=True)
    else:
        covar_df.loc[:,'age_of_onset'] = scaler.fit_transform(covar_df[['age_of_onset']])

    if covar_df.age.isna().all():
        covar_df.drop(columns=['age'], inplace=True)
    else:
        covar_df.loc[:,'age'] = scaler.fit_transform(covar_df[['age']])

    if covar_df.sex.isna().all():
        covar_df.drop(columns=['sex'], inplace=True)

    covar_df.drop(columns=['FID'], inplace=True)
    covar_df.rename(columns={'GP2sampleID':'sampleid'}, inplace=True)

    data_df = dosage_df.merge(covar_df, on='sampleid', how='left').merge(pheno_df, on='sampleid', how='left').set_index('sampleid')

    rm_pred = [f'PC{i}' for i in range(1,21)] + ['sex','age_of_onset','age','pheno']

    pred_list = [x for x in data_df.columns if x not in rm_pred]
    covars_list = [x for x in data_df.columns if x not in pred_list + [f'PC{i}' for i in range(11,21)] + ['pheno']]

    results = []
    fails = []

    for pred in range(len(pred_list)):
        pred_name = pred_list[pred]
        formula = "pheno ~ " + pred_name + " + " + ' + '.join(covars_list)

        fitted = sm.formula.glm(formula=formula, family=sm.families.Binomial(), data=data_df).fit()
        beta_coef  = fitted.params.loc[pred_name]
        beta_se  = fitted.bse.loc[pred_name]
        p_val = fitted.pvalues.loc[pred_name]

        results.append((pred_name, beta_coef, beta_se, p_val))


    output = pd.DataFrame(results, columns=('PREDICTOR', 'BETA_COEF', 'BETA_SE','P_VAL'))
    output.to_csv(out_path, sep='\t', header=True, index=False)
    
    


    
    
    
    


######## UNDER DEVELOPMENT #########

# create dosage matrices
# release2_cohorts = ['BCM','UMD','SYNAPS-KZ','MDGAP-QSBB','CORIELL']

# release2_key = key.loc[key.study.isin(release2_cohorts)]
# release2_key[['FID','GP2sampleID']].to_csv(f'{cnv_path}/release2.samples', sep='\t', header=False, index=False)
# release2_key[['GP2sampleID','IID']].to_csv(f'{cnv_path}/release2_sample_id_key.csv')
# release2_covars = release2_key.loc[:,['FID', 'GP2sampleID','sex_for_qc', 'age', 'age_of_onset']]
# # release2_key[['sex_for_qc', 'age', 'age_of_onset']].to_csv(f')


# labels = ['AAC','AFR','AJ','EAS','EUR','FIN','SAS','AMR','AMR_KZ']
# # labels = ['AAC']

# with open(f'{swarm_scripts_dir}/cnv_dosages.swarm', 'w') as f:
#     for label in labels:
#         geno = f'{cnv_path}/GP2_round2_{label}'

#         cmd1 = f'\
#     plink \
#     --bfile {geno} \
#     --keep {cnv_path}/release2.samples \
#     --make-bed \
#     --out {geno}_release2'

#         cmd2 = f'plink \
#     --bfile {geno}_release2 \
#     --pca \
#     --out {geno}_release2'

#         cmds = [cmd1, cmd2]

#         for cmd in cmds:
#             shell_do(cmd)

#         pcs = pd.read_csv(f'{geno}_release2.eigenvec', sep='\s+')
#         pc_num = pcs.iloc[:, 2:].shape[1]
#         pc_names = ['FID','GP2sampleID'] + [f'PC{i}' for i in range(1, pc_num+1)]
#         pcs.columns = pc_names

#         cov = pcs.merge(release2_covars, on=['FID','GP2sampleID'], how='left')
#         cov.age.fillna(cov.age.mean(), inplace=True)
#         cov.age_of_onset.fillna(cov.age_of_onset.mean(), inplace=True)
#         cov.sex_for_qc.fillna(cov.sex_for_qc.median(), inplace=True)
#         cov.rename(columns={'sex_for_qc':'sex'})
#         cov.to_csv(f'{geno}_release2.cov', sep='\t', header=True, index=False)
        
#         samples = cov.merge(release2_key[['GP2sampleID','IID']], on='GP2sampleID', how='left')
#         samples['IID'].to_csv(f'{geno}_release2_barcode.samples', header=False, index=False)
        
#         for chrom in chroms:
#             dosage_cmd = f'\
# python /data/vitaled2/GenoTools/run_cnv_dosage_pipeline.py \
# --metrics_in {idat_path} \
# --chrom {chrom} \
# --samples {geno}_release2_barcode.samples \
# --out_path {geno}'
#             f.write(f'{dosage_cmd}\n')
# f.close()















# cnv_dir = f'{basedir}/cnvs'
# for chrom in chroms:
#     create_cnv_dosage_matrices(in_path=ibx_idat_dir, samples_list=samples, chromosome=chrom, out_path=cnv_dir)
    
#     baf = f'{cnv_dir}/CNV_chr{chrom}_BAF.csv'
#     lrr_del = f'{cnv_dir}/CNV_chr{chrom}_LRR_DEL.csv'
#     lrr_dup = f'{cnv_dir}/CNV_chr{chrom}_LRR_DUP.csv'
    
    
    
