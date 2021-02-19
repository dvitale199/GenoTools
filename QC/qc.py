import subprocess
import argparse
import pandas as pd
import os
import glob
import shutil
import sys

# local imports
from QC.utils import shell_do, rm_tmps


def call_rate_pruning(geno_path, out_path, mind=0.05):

    # what step are we running?
    step = "call_rate_prune"
    print()
    print(f"RUNNING: {step}")
    print()
    
    outliers_out = f'{out_path}.outliers'
    phenos_out = f'{geno_path}.phenos'


    fam = pd.read_csv(f'{geno_path}.fam', sep='\s+', header=None)
    fam[[0,1,5]].to_csv(phenos_out, sep='\t', header=False, index=False)
    
    plink_cmd1 = f"plink --bfile {geno_path} --mind {mind} --make-bed --out {out_path}"

    shell_do(plink_cmd1)
    
    shutil.move(f'{out_path}.irem', outliers_out)
    
    outlier_count = sum(1 for line in open(f'{outliers_out}'))

    outfiles_dict = {
        'outliers_path': f'{outliers_out}',
        'plink_out': f'{out_path}',
        'phenos_path': f'{phenos_out}'
    }

    out_dict = {
        'step': step,
        'outlier_count': outlier_count,
        'output': outfiles_dict
    }

    return out_dict


def sex_check(geno_path, out_path, check_sex=[0.25,0.75]):
    
    # what step are we running?
    step = "sex_check"
    print()
    print(f"RUNNING: {step}")
    print()

    # create filenames
    sex_tmp1 = f"{out_path}_tmp1"
    sex_tmp2 = f"{out_path}_tmp2"
    sex_fails = f"{out_path}.fails"

    # check sex 2 methods
    plink_cmd1 = f"plink --bfile {geno_path} --check-sex 0.25 0.75 --maf 0.05 --out {sex_tmp1}"
    plink_cmd2 = f"plink --bfile {geno_path} --chr 23 --from-bp 2699520 --to-bp 154931043 --maf 0.05 --geno 0.05 --hwe 1E-5 --check-sex  0.25 0.75 --out {sex_tmp2}"

    cmds = [plink_cmd1, plink_cmd2]
    for cmd in cmds:
        shell_do(cmd)

    # grab fails from .sexcheck files
    sex1 = pd.read_csv('/data/vitaled2/test_data/mcgill/MCGILL_all_call_rate_het_sex_tmp1.sexcheck', sep='\s+')
    sex_fail1 = sex1[sex1.STATUS=='PROBLEM']

    sex2 = pd.read_csv('/data/vitaled2/test_data/mcgill/MCGILL_all_call_rate_het_sex_tmp2.sexcheck', sep='\s+')
    sex_fail2 = sex2[sex2.STATUS=='PROBLEM']

    # combine and output
    sex_fail_df = sex_fail1.append(sex_fail2)
    sex_fail_ids = sex_fail_df.loc[:,['FID','IID']].drop_duplicates(subset=['FID','IID'])
    sex_fail_count = sex_fail_ids.shape[0]
    sex_fail_ids.to_csv(sex_fails, sep='\t', index=False)

    # remove sex fail samples from geno
    plink_cmd3 = f"plink --bfile {geno_path} --remove {sex_fails} --make-bed --out {out_path}"
    
    shell_do(plink_cmd3)

    # remove tmp files
    tmps = [sex_tmp1, sex_tmp2]
    rm_tmps(tmps)

    # log outputs
    outfiles_dict = {
        'sex_fails': sex_fails,
        'plink_out': out_path
    }

    out_dict = {
        'step': step,
        'outlier_count': sex_fail_count,
        'output': outfiles_dict
    }

    return out_dict

    
def het_pruning(geno_path, out_path):
    
    # what step are we running?
    step = "het_prune"
    print()
    print(f"RUNNING: {step}")
    print()

    het_tmp = f"{out_path}_tmp"
    het_tmp2 = f"{out_path}_tmp2"
    het_tmp3 = f"{out_path}_tmp3"
    outliers_out = f"{out_path}.outliers"
    
    plink_cmd1 = f"plink --bfile {geno_path} --geno 0.01 --maf 0.05 --indep-pairwise 50 5 0.5 --out {het_tmp}"
    plink_cmd2 = f"plink --bfile {geno_path} --extract {het_tmp}.prune.in --make-bed --out {het_tmp2}"
    plink_cmd3 = f"plink --bfile {het_tmp2} --het --out {het_tmp3}"

    cmds1 = [plink_cmd1, plink_cmd2, plink_cmd3]

    for cmd in cmds1:
        shell_do(cmd)

    het = pd.read_csv(f'{het_tmp3}.het', sep='\s+')
    het_outliers = het[((het.F <= -0.25) | (het.F >= 0.25))]
    outlier_count = het_outliers.shape[0]
    het_outliers.to_csv(f'{outliers_out}', sep='\t', index=False)
    
    plink_cmd4 = f"plink --bfile {geno_path} --remove {outliers_out} --make-bed --out {out_path}"

    shell_do(plink_cmd4)

    # remove tmp files
    tmps = [het_tmp, het_tmp2, het_tmp3]
    rm_tmps(tmps)

    outfiles_dict = {
        'het_outliers': outliers_out,
        'plink_out': out_path
    }

    out_dict = {
        'step': step,
        'outlier_count': outlier_count,
        'output': outfiles_dict
    }

    return out_dict



    

