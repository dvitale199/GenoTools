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
    pca_cmd = f'plink2 --bfile {geno_path} --pca {n_pcs} --out {out_path}'
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

    hits_info_out = f'{out_path}.hits.info'
    hits_id_out = f'{out_path}.hits'

    # read in .fam file
    fam = pd.read_csv(f'{geno_path}.fam', sep='\s+', header=None)
    
    # get present phenotypes
    phenos = fam[5].unique()

    # get phenotype counts for debugging
    pheno_counts = fam[5].value_counts().to_dict()

    # check if multiple phenotypes present
    if(len(phenos) == 1):
        # check if phenotypes exist
        if(phenos[0] == -9):
            print(f'Association failed!')
            print(f'No phenotypes present in the {geno_path}.fam file')
        # phenotypes exist but only one present
        else:
            print(f'Association failed!')
            print(f'Only one phenotype present in the {geno_path}.fam file')
        
        process_complete = False

        outfiles_dict = {
            'hits': 'Association Failed!',
            'hits_info': 'Association Failed!',
            'plink_out': out_path
        }

        metrics_dict = {
            'phenotype_counts': pheno_counts,
            'hits': 0
        }
    
    else:
        # check if binary phenotypes are correct
        binary_phenos = [-9,1,2]
        if(model == 'logistic' and (not all(x in binary_phenos for x in phenos))):
            print(f'Association failed!')
            print(f'Binary phenotypes coded incorrectly in the {geno_path}.fam file')

            process_complete = False

            outfiles_dict = {
                'hits': 'Association Failed!',
                'hits_info': 'Association Failed!',
                'plink_out': out_path
            }

            metrics_dict = {
                'phenotype_counts': pheno_counts,
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
                    'hits': hits_id_out,
                    'hits_info': hits_info_out, 
                    'plink_out': out_path
                }

                metrics_dict = {
                    'phenotype_counts': pheno_counts,
                    'hits': hits_count
                }
            
            else:
                print(f'Association failed!')
                print(f'Check {out_path}.log for more information')

                process_complete = False

                outfiles_dict = {
                    'hits': 'Association Failed!',
                    'hits_info': 'Association Failed!',
                    'plink_out': out_path
                }

                metrics_dict = {
                    'phenotype_counts': pheno_counts,
                    'hits': 0
                }
    
    out_dict = {
        'pass': process_complete,
        'step': step,
        'metrics': metrics_dict,
        'output': outfiles_dict
    }

    return out_dict