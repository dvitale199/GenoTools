import pandas as pd
import numpy as np
import subprocess
import os
import sys
import random
from scipy.stats import ncx2

from QC.utils import shell_do


def plink_pca(geno_path, out_path, n_pcs=10):

    # what step are we running?
    step = 'plink_pca'
    print()
    print(f'RUNNING: {step}')
    print()

    # run pca
    pca_cmd = f' \
    plink2 \
    --bfile {geno_path} \
    --pca {n_pcs} \
    --out {out_path}'

    shell_do(pca_cmd)

    # check if .eigenvec is created
    if os.path.isfile(f'{out_path}.eigenvec'):
        process_complete = True

        outfiles_dict = {
            'plink_out': out_path
        }

        metrics_dict = {
            'pcs': n_pcs
        }
    
    else:
        print(f'PCA failed!')
        print(f'Check {out_path}.log for more information')

        process_complete = False

        outfiles_dict = {
            'plink_out': out_path
        }

        metrics_dict = {
            'pcs': 0
        }

    out_dict = {
        'pass': process_complete,
        'step': step,
        'metrics': metrics_dict,
        'output': outfiles_dict
    }

    return out_dict


def assoc(geno_path, covar_path, out_path, model):

    # what step are we running?
    step = 'assoc'
    print()
    print(f'RUNNING: {step}')
    print()

    pheno_counts_out = f'{out_path}.pheno.counts'
    hits_info_out = f'{out_path}.hits.info'
    hits_id_out = f'{out_path}.hits'

    # read in .fam file
    fam = pd.read_csv(f'{geno_path}.fam', sep='\s+', header=None)
    
    # get present phenotypes
    phenos = fam[5].unique()

    # get phenotype counts for debugging
    pheno_counts = fam[5].value_counts()
    pheno_counts.to_csv(pheno_counts_out, sep='\t', header=False)

    # check if multiple phenotypes present
    if len(phenos)== 1:
        # check if phenotypes exist
        if phenos[0] == -9:
            print(f'Association failed!')
            print(f'No phenotypes present in the {geno_path}.fam file')
        # phenotypes exist but only one present
        else:
            print(f'Association failed!')
            print(f'Only one phenotype present in the {geno_path}.fam file')
        
        print(f'Phenotype Counts: {pheno_counts.to_dict()}')
        
        process_complete = False

        outfiles_dict = {
            'pheno_counts': pheno_counts_out,
            'hits': 'Association Failed!',
            'hits_info': 'Association Failed!',
            'plink_out': out_path
        }

        metrics_dict = {
            'hits': 0
        }
    
    else:
        # check if binary phenotypes are correct
        binary_phenos = [-9,1,2]
        if model == 'logistic' and (not all(x in binary_phenos for x in phenos)):
            print(f'Association failed!')
            print(f'Binary phenotypes coded incorrectly in the {geno_path}.fam file')
            print(f'Phenotype Counts: {pheno_counts.to_dict()}')

            process_complete = False

            outfiles_dict = {
                'pheno_counts': pheno_counts_out,
                'hits': 'Association Failed!',
                'hits_info': 'Association Failed!',
                'plink_out': out_path
            }

            metrics_dict = {
                'hits': 0
            }

        else:
            # run association
            assoc_cmd = f'\
            plink2 \
            --bfile {geno_path} \
            --covar {covar_path} \
            --{model} \
            --allow-no-sex \
            --adjust cols=+qq \
            --covar-variance-standardize \
            --out {out_path}'

            shell_do(assoc_cmd)

            # output file names 
            glm_file = f'{out_path}.PHENO1.glm.{model}.adjusted'
            glm_hybrid_file = f'{out_path}.PHENO1.glm.{model}.hybrid.adjusted'

            # check if assoc output file exists
            if os.path.isfile(glm_file) or os.path.isfile(glm_hybrid_file):
                # read adjusted (or hybrid) file
                try:
                    adj = pd.read_csv(glm_file, sep='\s+')
                except:
                    adj = pd.read_csv(glm_hybrid_file, sep='\s+')
                # getting hits and associated information
                hits_info = adj[adj.BONF <= 0.05]
                hits_info.to_csv(hits_info_out, sep='\t', header=True, index=False)
                hits_id = adj.loc[adj.BONF <= 0.05, 'ID']
                hits_id.to_csv(hits_id_out, sep='\t', header=False, index=False)
                hits_count = hits_id.shape[0]

                process_complete = True

                outfiles_dict = {
                    'pheno_counts': pheno_counts_out,
                    'hits': hits_id_out,
                    'hits_info': hits_info_out, 
                    'plink_out': out_path
                }

                metrics_dict = {
                    'hits': hits_count
                }
            
            else:
                print(f'Association failed!')
                print(f'Check {out_path}.log for more information')

                process_complete = False

                outfiles_dict = {
                    'pheno_counts': pheno_counts_out,
                    'hits': 'Association Failed!',
                    'hits_info': 'Association Failed!',
                    'plink_out': out_path
                }

                metrics_dict = {
                    'hits': 0
                }
    
    out_dict = {
        'pass': process_complete,
        'step': step,
        'metrics': metrics_dict,
        'output': outfiles_dict
    }

    return out_dict


def prs(geno_path, out_path, assoc, clump_p1=1e-3, clump_r2=0.50, clump_kb=250):
    # clump and prs as one function keeps things simple for user
    # if user wants to LD clump for another purpose the clump output is availiable

    # what step are we running?
    step = 'prs'
    print()
    print(f'RUNNING: {step}')
    print()

    assoc_file = f'{out_path}.assoc'
    weights = f'{out_path}_weights.tab'
    clump_snps = f'{out_path}_clump.snps'
    snp_pvals = f'{out_path}_pvals.snps'
    range = f'{out_path}.range'

    # run clump - not yet supported by plink2
    clump_cmd = f'\
    plink \
    --bfile {geno_path} \
    --clump-p1 {clump_p1} \
    --clump-r2 {clump_r2} \
    --clump-kb {clump_kb} \
    --clump {assoc} \
    --clump-snp-field ID \
    --clump-field P \
    --allow-no-sex \
    --out {out_path}_clump'

    shell_do(clump_cmd)

    # check if .clumped is created
    if os.path.isfile(f'{out_path}_clump.clumped'):
        # get weights
        glm = pd.read_csv(assoc, sep='\s+')
        glm_hits = glm.loc[(glm.TEST=='ADD')].copy()
        glm_hits = glm_hits.dropna()

        # convert OR if logistic
        if 'logistic' in assoc:
            glm_hits.loc[:,'BETA'] = np.log10(glm_hits.loc[:,'OR'])
        
        glm_hits.loc[:,['ID','A1','BETA']].to_csv(weights, sep='\t', header=True, index=False)
        
        # get clump SNPs
        clump = pd.read_csv(f'{out_path}_clump.clumped', sep='\s+')
        num_clumps = clump.shape[0]
        clump[['SNP']].to_csv(clump_snps, sep='\t', header=False, index=False)

        # get SNP pvals
        glm_hits.loc[:,['ID','P']].to_csv(snp_pvals, sep='\t', header=False, index=False)

        # get assoc file
        glm_hits.to_csv(assoc_file, sep='\t', header=True, index=False)

        # writing ranges
        lines = ["s1 0 0.001","s2 0 0.05","s3 0 0.1"]
        for line in lines:
            os.system(f'echo {line} >> {range}')

        # run PRS
        prs_cmd = f'\
        plink2 \
        --bfile {geno_path} \
        --score {weights} 1 2 3 \
        --q-score-range {range} {snp_pvals} \
        --extract {clump_snps} \
        --allow-no-sex \
        --out {out_path}.PRS \
        --memory 50000'

        shell_do(prs_cmd)

        s1 = f'{out_path}.PRS.s1.sscore'
        s2 = f'{out_path}.PRS.s2.sscore'
        s3 = f'{out_path}.PRS.s3.sscore'

        # check if three .sscore files are created
        if os.path.isfile(s1) and os.path.isfile(s2) and os.path.isfile(s3):
            process_complete = True
        
        else:
            print(f'PRS failed!')
            print(f'Check {out_path}.PRS.log for more information')

            process_complete = False

        outfiles_dict = {
            'SNP_weights': weights,
            'clump_SNPs': clump_snps,
            'SNP_pvals': snp_pvals,
            'ranges': range,
            'assoc': assoc_file,
            'plink_out': out_path
        }

        metrics_dict = {
            'clump_pval': clump_p1,
            'clump_r2': clump_r2,
            'clump_kb': clump_kb,
            'num_clumps': num_clumps
        }

    else:
        print(f'PRS failed!')
        print(f'Check {out_path}_clump.log for more information')

        process_complete = False

        outfiles_dict = {
            'SNP_weights': 'PRS Failed!',
            'clump_SNPs': 'PRS Failed!',
            'SNP_pvals': 'PRS Failed!',
            'ranges': 'PRS Failed!',
            'assoc': 'PRS Failed!',
            'plink_out': out_path
        }

        metrics_dict = {
                'clump_pval': clump_p1,
                'clump_r2': clump_r2,
                'clump_kb': clump_kb,
                'num_clumps': 0
        }
    
    out_dict = {
        'pass': process_complete,
        'step': step,
        'metrics': metrics_dict,
        'output': outfiles_dict
    }

    return out_dict


def calculate_inflation(pval_array, normalize=False, ncases=None, ncontrols=None):

    # what step are we running?
    step = 'lambda_calculation'
    print()
    print(f'RUNNING: {step}')
    print()

    # need cases and controls if normalizing
    if normalize and (not ncases or not ncontrols):
        print(f'Inflation Calculation failed!')
        print(f'If normalizing, please add ncases and ncontrols')

        process_complete = False

        inflation = 0
    
    else:
        # calculate inflation
        num = ncx2.ppf(1-pval_array, 1, nc=0)
        denom = ncx2.ppf(0.5, 1, nc = 0)
        inflation = np.nanmedian(num)/denom

        # normalize when necessary
        if normalize:
            inflation1000 = 1 + (inflation -1) * (1/ncases+ 1/ncontrols)/(1/1000 + 1/1000)
            inflation = inflation1000
        
        process_complete = True

    metrics_dict = {
        'inflation': inflation
    }

    out_dict = {
        'pass': process_complete,
        'step': step,
        'metrics': metrics_dict
    }

    return out_dict


def munge(geno_path, out_path, assoc, ref_panel):

    # what step are we running?
    step = 'munge'
    print()
    print(f'RUNNING: {step}')
    print()

    # get plink frequency report
    freq_cmd = f'\
    plink2 \
    --bfile {geno_path} \
    --freq \
    --out {out_path}'

    shell_do(freq_cmd)

    # making sure .afreq is created
    if os.path.isfile(f'{out_path}.afreq'):
        # read frequnecy report
        freq_df = pd.read_csv(f'{out_path}.afreq', sep='\s+')

        # read association file
        assoc_df = pd.read_csv(assoc, sep='\s+')

        # merge on ID
        merge_df = assoc_df.merge(freq_df, how='inner', on='ID')

        # creating A2 column
        # when A1 (effect allele) is the ALT allele, A2 is the REF
        # when A1 (effect allele) is the REF allele, A2 is the ALT
        merge_df = merge_df.assign(A2=np.where(merge_df.A1 == merge_df.ALT_x, merge_df.REF_x, merge_df.ALT_x))

        # creating freq column - assigning freq to A1 (effect allele)
        # when A1 (effect allele) is the ALT allele, freq is ALT_FREQS
        # when A1 (effect allele) is the REF allele, freq is 1-ALT_FREQS
        merge_df = merge_df.assign(freq=np.where(merge_df.A1 == merge_df.ALT_x, merge_df.ALT_FREQS, 1-merge_df.ALT_FREQS))

        # fix MAF - make sure A1 (effect allele) is the minor allele
        # if freq > 0.5, flip alleles and sign of B
        maf_df = merge_df.copy()
        maf_df.A1, maf_df.A2 = np.where(maf_df.freq > 0.5, (maf_df.A2, maf_df.A1), (maf_df.A1, maf_df.A2))
        maf_df.BETA = np.where(maf_df.freq > 0.5, -maf_df.BETA, maf_df.BETA)
        maf_df.freq = np.where(maf_df.freq > 0.5, 1-maf_df.freq, maf_df.freq)

        # rename columns and isolate needed ones
        maf_df = maf_df.rename(columns={'#CHROM_x':'chr','POS':'bp','ID':'SNP',
                                        'OBS_CT_x':'n','P':'p','BETA':'b','LOG(OR)_SE':'se'})
        maf_df = maf_df.loc[:,['SNP','chr','bp','A1','A2','freq','b','se','p','n']]

        # read in reference panel .bim and create ID
        bim1k = pd.read_csv(f'{ref_panel}.bim', sep='\s+', header=None)
        bim1k.columns = ['chr','rsid','loc','bp','a1','a2']
        bim1k['ID'] = bim1k['chr'].astype('str') + ':' + bim1k['bp'].astype('str')

        # isolate ID, rsid for merging
        ref_rsids = bim1k.loc[:,['ID','rsid']]

        # seperate maf_df into rsid snps and non-rsid snps
        rsid = maf_df[maf_df['SNP'].str.contains('rs')]
        non_rsid = maf_df[~maf_df['SNP'].str.contains('rs')]

        # create ID for non-rs SNPs
        non_rsid_cp = non_rsid.copy()
        non_rsid_cp['ID'] = non_rsid_cp['chr'].astype('str') + ':' + non_rsid_cp['bp'].astype('str')

        # merge with reference
        add_rsid = non_rsid_cp.merge(ref_rsids, how='inner', on='ID')
        add_rsid = add_rsid.loc[:,['rsid','chr','bp','A1','A2','freq','b','se','p','n']]
        add_rsid = add_rsid.rename(columns={'rsid':'SNP'})

        # append to rsid
        new_rsid_df = rsid.append(add_rsid)
        new_rsid_df = new_rsid_df.reset_index(drop=True)

        # sort by chr and bp
        new_rsid_df.sort_values(['chr','bp'], ascending=[True,True], inplace=True)

        # getting correctly formatted output and coordinates
        ma_out = new_rsid_df.loc[:, ['SNP','A1','A2','freq','b','se','p','n']]
        coords = new_rsid_df.loc[:, ['SNP','chr','bp']]

        # get num_snps
        num_snps = ma_out.shape[0]

        process_complete = True

        data_dict = {
            'ma_format_df': ma_out,
            'coordinates': coords
        }

        outfiles_dict = {
            'plink_out': out_path
        }

        metrics_dict = {
            'num_snps': num_snps
        }

    else:
        print(f'Munge Failed!')
        print(f'Check {out_path}.log for more information')

        process_complete = False

        data_dict = {
            'ma_format_df': 'Munge Failed!',
            'coordinates': 'Munge Failed!'
        }

        outfiles_dict = {
            'plink_out': out_path
        }

        metrics_dict = {
            'num_snps': 0
        }
    
    out_dict = {
        'pass': process_complete,
        'step': step,
        'data': data_dict,
        'metrics': metrics_dict,
        'output': outfiles_dict
    }

    return out_dict