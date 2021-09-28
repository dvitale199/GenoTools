import pandas as pd
import numpy as np
import subprocess
import os
import sys
import random
from scipy.stats import ncx2

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