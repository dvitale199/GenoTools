import subprocess
import argparse
import pandas as pd
import os
import glob
import shutil
import sys

# local imports
from QC.utils import shell_do, rm_tmps


################ Sample pruning methods ####################
def callrate_prune(geno_path, out_path, mind=0.02):

    # what step are we running?
    step = "callrate_prune"
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


def sex_prune(geno_path, out_path, check_sex=[0.25,0.75]):
    
    # what step are we running?
    step = "sex_prune"
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
    sex1 = pd.read_csv(f'{sex_tmp1}.sexcheck', sep='\s+')
    sex_fail1 = sex1[sex1.STATUS=='PROBLEM']

    sex2 = pd.read_csv(f'{sex_tmp2}.sexcheck', sep='\s+')
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


def het_prune(geno_path, out_path):
    
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


def related_prune(geno_path, out_path, grm_cutoff=0.125):
    
    # what step are we running?
    step = "related_prune"
    print()
    print(f"RUNNING: {step}")
    print()

    # make temp filenames
    related_out = f"{out_path}.related" 
    grm1 = f"{out_path}_grm1"
    grm2 = f"{out_path}_grm2"
    grm3 = f"{out_path}_grm3"

    # calculate grm and select relatedness <= grm_cutoff
    gcta_cmd1 = f"gcta --bfile {geno_path} --autosome --maf 0.05 --make-grm --out {grm1}" 
    gcta_cmd2 = f"gcta --grm {grm1} --grm-cutoff {grm_cutoff} --make-grm --out {grm2}"
    plink_cmd1 = f"plink --bfile {geno_path} --keep {grm2}.grm.id --make-bed --out {out_path}"
    # see if any samples are duplicated (grm cutoff >= 0.95)
    gcta_cmd3 = f"gcta --grm {grm1} --grm-cutoff 0.95 --make-grm --out {grm3}"

    cmds = [gcta_cmd1,gcta_cmd2,plink_cmd1, gcta_cmd3]
    for cmd in cmds:
        shell_do(cmd)
    
    # get sample counts
    total_count = sum(1 for line in open(f'{geno_path}.fam'))
    unrelated_count = sum(1 for line in open(f'{grm2}.grm.id'))
    related_count = total_count - unrelated_count
    duplicated_count = sum(1 for line in open(f'{grm3}.grm.id'))

    # get related sample ids
    fam = pd.read_csv(f'{geno_path}.fam', sep='\s+', header=None, usecols=[0,1], names=['FID','IID'])
    unrelated_ids = pd.read_csv(f'{grm2}.grm.id', sep='\t', header=None, names=['FID','IID'])
    merged = fam.merge(unrelated_ids.drop_duplicates(), on=['FID','IID'], how='left', indicator=True)
    related = merged[merged['_merge'] == 'left_only'].drop(columns=['_merge'])
    related['status'] = 'related'

    # get duplicated sample ids
    duplicated = pd.read_csv(f'{grm3}.grm.id', sep='\t', header=None, names=['FID','IID'])
    duplicated['status'] = 'duplicated'

    # append duplicated sample ids to related sample ids, drop_duplicates(keep='last) because all duplicated would also be considered related
    grm_pruned = related.append(duplicated)
    grm_pruned.drop_duplicates(subset=['FID', 'IID'], keep='last', inplace=True)
    grm_pruned.to_csv(related_out, sep='\t', header=True, index=False)
    
    # remove temp files
    tmps = [grm1, grm2]
    rm_tmps(tmps)

    outfiles_dict = {
        'relateds': related_out,
        'plink_out': out_path
    }

    out_dict = {
        'step': step,
        'related_count': related_count,
        'duplicated_count': duplicated_count,
        'output': outfiles_dict
    }

    return out_dict

    

################ Variant pruning methods ####################
##### SPLIT THESE INTO SEPARATE FUNCTIONS!!!!

# def variant_pruning(self, geno_path):

#     out_path = self.out_path
    
#     step = "VARIANT-LEVEL PRUNING"
#     print(step)

#     # variant missingness
#     bash1 = "plink --bfile " + geno_path + " --make-bed --out " + geno_path + "_geno --geno 0.05"

#     #missingness by case control (--test-missing), using P > 1E-4
#     bash2 = "plink --bfile " + geno_path + "_geno --test-missing --out " + out_path + "missing_snps" 
#     bash3 = "awk '{if ($5 <= 0.0001) print $2 }' " + out_path + "missing_snps.missing > " + out_path + "missing_snps_1E4.txt"
#     bash4 = "plink --bfile " + geno_path + "_geno --exclude " + out_path + "missing_snps_1E4.txt --make-bed --out " + geno_path + "_geno_missingsnp"

#     #missingness by haplotype (--test-mishap), using P > 1E-4
#     bash5 = "plink --bfile " + geno_path + "_geno_missingsnp --test-mishap --out " + out_path + "missing_hap" 
#     bash6 = "awk '{if ($8 <= 0.0001) print $9 }' " + out_path + "missing_hap.missing.hap > " + out_path + "missing_haps_1E4.txt"
#     bash7 = "cat " + out_path + "missing_haps_1E4.txt | tr '|' '\n' > " + out_path + "missing_haps_1E4_final.txt"
#     bash8 = "plink --bfile " + geno_path + "_geno_missingsnp --exclude " + out_path + "missing_haps_1E4_final.txt --make-bed --out " +  geno_path + "_geno_missingsnp_missinghap"


#     ###### THIS DOES NOT WORK WITHOUT PHENOTYPES!!!!!!!!
#     #HWE from controls only using P > 1E-4
#     bash9 = "plink --bfile " + geno_path + "_geno_missingsnp_missinghap --filter-controls --hwe 1E-4 --write-snplist --out " + out_path + "hwe"
#     bash10 = "plink --bfile " + geno_path + "_geno_missingsnp_missinghap --extract " + out_path + "hwe.snplist --make-bed --out " + geno_path + "_geno_missingsnp_missinghap_hwe"
#     bash11 = "#### moved " + geno_path + "_geno_missingsnp_missinghap_hwe to " + geno_path + "_variant #####"
#     cmds = [bash1, bash2, bash3, bash4, bash5, bash6, bash7, bash8, bash9, bash10, bash11]

#     self.run_cmds(cmds, step)

#     exts = [".bed",".bim",".fam",".log"]
#     for ext in exts:
        
#         shutil.move(geno_path + "_geno_missingsnp_missinghap_hwe" + ext, geno_path + "_variant" + ext)