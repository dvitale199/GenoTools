import subprocess
import argparse
import pandas as pd
import os
import glob
import shutil
import sys

# local imports
from QC.utils import shell_do, rm_tmps, count_file_lines


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
    
    if os.path.isfile(f'{out_path}.irem'):
        shutil.move(f'{out_path}.irem', outliers_out)

        outlier_count = sum(1 for line in open(f'{outliers_out}'))
        
    else:
        outlier_count = 0
        
    process_complete = True
    
    outfiles_dict = {
        'pruned_samples': f'{outliers_out}',
        'plink_out': f'{out_path}',
        'phenos_path': f'{phenos_out}'
    }

    metrics_dict = {
        'outlier_count': outlier_count
    }

    out_dict = {
        'pass': process_complete,
        'step': step,
        'metrics': metrics_dict,
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
#     tmps = [sex_tmp1, sex_tmp2]
#     rm_tmps(tmps)
    
    process_complete = True
    
    # log outputs
    outfiles_dict = {
        'pruned_samples': sex_fails,
        'plink_out': out_path
    }

    metrics_dict = {
        'outlier_count': sex_fail_count
    }

    out_dict = {
        'pass': process_complete,
        'step': step,
        'metrics': metrics_dict,
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
    
    # check if het_tmp3 is created. if not, skip this.
    # NOTE: there may be legitimate reasons for this, for example, if only one sample in genotype file (sometimes happens in ancestry split method)
    hetpath = f'{het_tmp3}.het'
    if os.path.isfile(hetpath):
        het = pd.read_csv(hetpath, sep='\s+')
        het_outliers = het[((het.F <= -0.25) | (het.F >= 0.25))]
        outlier_count = het_outliers.shape[0]
        het_outliers.to_csv(f'{outliers_out}', sep='\t', index=False)

        plink_cmd4 = f"plink --bfile {geno_path} --remove {outliers_out} --make-bed --out {out_path}"

        shell_do(plink_cmd4)

#         # remove tmp files
#         tmps = [het_tmp, het_tmp2, het_tmp3]
#         rm_tmps(tmps)

        outfiles_dict = {
            'pruned_samples': outliers_out,
            'plink_out': out_path
        }

        metrics_dict = {
            'outlier_count': outlier_count
        }
        
        process_complete = True
    
    else:
        print(f'Heterozygosity pruning failed!')
        print(f'Check {het_tmp}.log, {het_tmp2}.log, or {het_tmp3}.log for more information')
        
        outfiles_dict = {
            'pruned_samples': 'Heterozygosity Pruning Failed!',
            'plink_out': [het_tmp, het_tmp2, het_tmp3]
        }

        metrics_dict = {
            'outlier_count': 0
        }
    
        process_complete = False
    
    out_dict = {
        'pass': process_complete,
        'step': step,
        'metrics': metrics_dict,
        'output': outfiles_dict
    }

    return out_dict


def related_prune(geno_path, out_path, related_grm_cutoff=0.125, duplicated_grm_cutoff=0.95):
    
    # what step are we running?
    step = "related_prune"
    print()
    print(f"RUNNING: {step}")
    print()

    # make temp filenames
    related_out = f"{out_path}.related" 
    grm1 = f"{out_path}_total_grm_tmp"
    grm2 = f"{out_path}_unrelated_grm_tmp"
    grm3 = f"{out_path}_duplicated_grm_tmp"

    # calculate grm and select relatedness <= grm_cutoff
    gcta_cmd1 = f"gcta --bfile {geno_path} --autosome --maf 0.05 --make-grm --out {grm1}" 
    gcta_cmd2 = f"gcta --grm {grm1} --grm-cutoff {related_grm_cutoff} --make-grm --out {grm2}"
    plink_cmd1 = f"plink --bfile {geno_path} --keep {grm2}.grm.id --make-bed --out {out_path}"
    # see if any samples are duplicated (grm cutoff >= 0.95)
    gcta_cmd3 = f"gcta --grm {grm1} --grm-cutoff {duplicated_grm_cutoff} --make-grm --out {grm3}"

    cmds = [gcta_cmd1, gcta_cmd2, plink_cmd1, gcta_cmd3]
    for cmd in cmds:
        shell_do(cmd)
    
    
    # get sample counts
    total_count = sum(1 for line in open(f'{geno_path}.fam'))
    unrelated_count = sum(1 for line in open(f'{grm2}.grm.id'))
    nonduplicated_count = sum(1 for line in open(f'{grm3}.grm.id'))
    duplicated_count = total_count - nonduplicated_count
    related_count = total_count - unrelated_count - duplicated_count

    # get related sample ids
    fam = pd.read_csv(f'{geno_path}.fam', sep='\s+', header=None, usecols=[0,1], names=['FID','IID'])
    unrelated_ids = pd.read_csv(f'{grm2}.grm.id', sep='\t', header=None, names=['FID','IID'])
    merged_rel = fam.merge(unrelated_ids.drop_duplicates(), on=['FID','IID'], how='left', indicator=True)
    related = merged_rel[merged_rel['_merge'] == 'left_only'].drop(columns=['_merge'])
    related['status'] = 'related'

    # get duplicated sample ids
    nonduplicated = pd.read_csv(f'{grm3}.grm.id', sep='\t', header=None, names=['FID','IID'])
    merge_dup = fam.merge(nonduplicated.drop_duplicates(), on=['FID','IID'], how='left', indicator=True)
    duplicated =  merge_dup[merge_dup['_merge'] == 'left_only'].drop(columns=['_merge'])
    duplicated['status'] = 'duplicated'

    # append duplicated sample ids to related sample ids, drop_duplicates(keep='last) because all duplicated would also be considered related
    grm_pruned = related.append(duplicated)
    grm_pruned.drop_duplicates(subset=['FID', 'IID'], keep='last', inplace=True)
    grm_pruned.to_csv(related_out, sep='\t', header=True, index=False)
    
#     # remove temp files
#     tmps = [grm1, grm2, grm3]
#     rm_tmps(tmps)
    
    process_complete = True
    
    outfiles_dict = {
        'pruned_samples': related_out,
        'plink_out': out_path
    }

    metrics_dict = {
        'related_count': related_count,
        'duplicated_count': duplicated_count
    }

    out_dict = {
        'pass': process_complete,
        'step': step,
        'metrics': metrics_dict,
        'output': outfiles_dict
    }

    return out_dict

    

################ Variant pruning methods ####################
# this will eventually be broken up into multiple functions!!!!
def variant_prune(geno_path, out_path):
    
    step = "variant_prune"
    print()
    print(f"RUNNING: {step}")
    print()

    # make tmp names
    # geno
    geno_tmp1 = f'{out_path}_geno_tmp1'
    # missingness by case-control
    mis_tmp1 = f'{out_path}_mis_tmp1'
    mis_tmp2 = f'{out_path}_mis_tmp2'
    # missingness by haplotype
    hap_tmp1 = f'{out_path}_hap_tmp1'
    hap_tmp2 = f'{out_path}_hap_tmp2'
    # HWE
    hwe_tmp1 = f'{out_path}_hwe_tmp1'

    # get initial snp count
    initial_snp_count = count_file_lines(f'{geno_path}.bim')
    
    # variant missingness
    plink_cmd1 = f"plink --bfile {geno_path} --geno 0.05 --make-bed --out {geno_tmp1}"
    shell_do(plink_cmd1)
    
    # geno pruned count
    geno_snp_count = count_file_lines(f'{geno_tmp1}.bim')
    geno_rm_count = initial_snp_count - geno_snp_count
    
    
    fam = pd.read_csv(f'{geno_path}.fam', sep='\s+', header=None, usecols=[5], names=['case'])
    # check if this contains both cases and controls
    if  all(x in fam['case'].unique() for x in [1, 2]):
        #missingness by case control (--test-missing), using P > 1E-4
        plink_cmd2 = f"plink --bfile {geno_tmp1} --test-missing --out {mis_tmp1}"
        shell_do(plink_cmd2)

        mis1 = f'{mis_tmp1}.missing'
        if os.path.isfile(mis1):

            mis = pd.read_csv(mis1, sep='\s+')
            exclude = mis[mis.P <= 0.0001].loc[:,'SNP']
            exclude.to_csv(f'{mis_tmp1}.exclude', sep='\t', header=False, index=False)

            plink_cmd3 = f"plink --bfile {geno_tmp1} --exclude {mis_tmp1}.exclude --make-bed --out {mis_tmp2}"
            shell_do(plink_cmd3)

            # mis purned count
            mis_snp_count = count_file_lines(f'{mis_tmp2}.bim')
            mis_rm_count = geno_snp_count - mis_snp_count

        else:
            print()
            print(f'Case/Control Missingness pruning failed!')
            print(f'Check {mis_tmp1}.log for more information')
            print()
            mis_snp_count = count_file_lines(f'{geno_tmp1}.bim')
            mis_rm_count = 0

            # if this step fails because case/control status unknown, set mis_tmp1=geno_tmp1 from previous step to continue pipeline
            mis_tmp2 = geno_tmp1

        # missingness by haplotype (--test-mishap), using P > 1E-4
        plink_cmd4 = f"plink --bfile {mis_tmp2} --test-mishap --out {hap_tmp1}"
        shell_do(plink_cmd4)

        # read .missing.hap file and grab flanking snps for P <= 0.0001. write flanking snps to file to exclude w bash
        mis_hap = pd.read_csv(f'{hap_tmp1}.missing.hap', sep='\s+')
        mis_hap_snps = list(mis_hap[mis_hap.P <= 0.0001].loc[:,'FLANKING'].str.split('|'))
        snp_ls_df = pd.DataFrame({'snp':[rsid for ls in mis_hap_snps for rsid in ls]})
        snp_ls_df['snp'].to_csv(f'{hap_tmp1}.exclude',sep='\t', header=False, index=False)

        plink_cmd5 = f"plink --bfile {mis_tmp2} --exclude {hap_tmp1}.exclude --make-bed --out {hap_tmp2}"

        # HWE from controls only using P > 1E-4
        plink_cmd6 = f"plink --bfile {hap_tmp2} --filter-controls --hwe 1E-4 --write-snplist --out {hwe_tmp1}"
        plink_cmd7 = f"plink --bfile {hap_tmp2} --extract {hwe_tmp1}.snplist --make-bed --out {out_path}"

        cmds3 = [plink_cmd5, plink_cmd6, plink_cmd7]
        for cmd in cmds3:
            shell_do(cmd)

        # hap pruned count
        hap_snp_count = count_file_lines(f'{hap_tmp2}.bim')
        hap_rm_count = mis_snp_count - hap_snp_count
        # hwe pruned count
        final_snp_count = count_file_lines(f'{out_path}.bim')
        hwe_rm_count = hap_snp_count - final_snp_count
        # total pruned count
        total_rm_count = initial_snp_count - final_snp_count

        # remove temp files
    #     tmps = [geno_tmp1, mis_tmp1, mis_tmp2, hap_tmp1, hap_tmp2, hwe_tmp1]
    #     rm_tmps(tmps)
        
        process_complete = True
        
        outfiles_dict = {
            'plink_out': out_path
        }

        metrics_dict = {
            'geno_removed_count': geno_rm_count,
            'mis_removed_count': mis_rm_count,
            'haplotype_removed_count': hap_rm_count,
            'hwe_removed_count': hwe_rm_count,
            'total_removed_count': total_rm_count
        }
        
    else:
        print()
        print(f'Case/control pruning failed! May be missing Controls!')
        print(f'Check number of Cases/Controls')
        print()
        
        process_complete = False
        
        outfiles_dict = {
            'plink_out': 'FAILED! Check Case/Control Status'
        }

        metrics_dict = {
            'geno_removed_count': 0,
            'mis_removed_count': 0,
            'haplotype_removed_count': 0,
            'hwe_removed_count': 0,
            'total_removed_count': 0
        }

    out_dict = {
        'pass': process_complete,
        'step': step,
        'metrics': metrics_dict,
        'output': outfiles_dict
    }

    return out_dict


def avg_miss_rates(geno_path, out_path):
    plink_miss_cmd = f'\
plink \
--bfile {geno_path} \
--missing \
--out {out_path}'

    shell_do(plink_miss_cmd)

    # get average call rate
    lmiss = pd.read_csv(f'{out_path}.lmiss', sep='\s+')
    imiss = pd.read_csv(f'{out_path}.imiss', sep='\s+')
    avg_lmiss = lmiss.F_MISS.mean()
    avg_imiss = imiss.F_MISS.mean()
    print(f'Average Missing Call Rate (lmiss): {avg_lmiss}')
    print(f'Average Missing Genotyping Rate (imiss): {avg_imiss}')

    metrics = {
        'avg_lmiss': avg_lmiss,
        'avg_imiss': avg_imiss
    }

    return metrics