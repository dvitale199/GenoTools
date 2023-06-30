import subprocess
import argparse
import pandas as pd
import numpy as np
import os
import glob
import shutil
import sys

# local imports
from QC.utils import shell_do, rm_tmps, count_file_lines, concat_logs

from utils.dependencies import check_plink, check_plink2

plink_exec = check_plink()
plink2_exec = check_plink2()


################ Sample pruning methods ####################
def callrate_prune(geno_path, out_path, mind=0.02):
    
    # what step are we running?
    step = "callrate_prune"
    print()
    print(f"RUNNING: {step}")
    print()
    
    outliers_out = f'{out_path}.outliers'

    fam = pd.read_csv(f'{geno_path}.fam', sep='\s+', header=None)
    
    plink_cmd1 = f"{plink2_exec} --bfile {geno_path} --mind {mind} --make-bed --out {out_path}"

    shell_do(plink_cmd1)

    listOfFiles = [f'{out_path}.log']
    concat_logs(step, out_path, listOfFiles)
    
    if os.path.isfile(f'{out_path}.mindrem.id'):
        irem = pd.read_csv(f'{out_path}.mindrem.id', sep='\s+', header=None, names=['FID','IID'])
        irem.to_csv(outliers_out, sep='\t', header=True, index=False)

        outlier_count = sum(1 for line in open(f'{outliers_out}'))
        
    else:
        outlier_count = 0
        
    process_complete = True
    
    outfiles_dict = {
        'pruned_samples': f'{outliers_out}',
        'plink_out': f'{out_path}',
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
    sex_fails = f"{out_path}.outliers"

    # check sex 2 methods
    plink_cmd1 = f"{plink_exec} --bfile {geno_path} --check-sex 0.25 0.75 --maf 0.05 --out {sex_tmp1}"
    plink_cmd2 = f"{plink_exec} --bfile {geno_path} --chr 23 --from-bp 2699520 --to-bp 154931043 --maf 0.05 --geno 0.05 --hwe 1E-5 --check-sex  0.25 0.75 --out {sex_tmp2}"

    cmds = [plink_cmd1, plink_cmd2]
    for cmd in cmds:
        shell_do(cmd)

    listOfFiles = [f'{sex_tmp1}.log', f'{sex_tmp2}.log']
    concat_logs(step, out_path, listOfFiles)

    if os.path.isfile(f'{sex_tmp1}.sexcheck') and os.path.isfile(f'{sex_tmp2}.sexcheck'):

        # grab fails from .sexcheck files
        sex1 = pd.read_csv(f'{sex_tmp1}.sexcheck', sep='\s+')
        sex_fail1 = sex1[sex1.STATUS=='PROBLEM']

        sex2 = pd.read_csv(f'{sex_tmp2}.sexcheck', sep='\s+')
        sex_fail2 = sex2[sex2.STATUS=='PROBLEM']

        # combine and output
        sex_fail_df = sex_fail1.append(sex_fail2)
        sex_fail_ids = sex_fail_df.loc[:,['FID','IID']].drop_duplicates(subset=['FID','IID'])
        sex_fail_count = sex_fail_ids.shape[0]
        sex_fail_ids.to_csv(sex_fails, sep='\t', header=True, index=False)

        # remove sex fail samples from geno
        plink_cmd3 = f"{plink2_exec} --bfile {geno_path} --remove {sex_fails} --make-bed --out {out_path}"
        
        shell_do(plink_cmd3)

        listOfFiles = [f'{out_path}.log']
        concat_logs(step, out_path, listOfFiles)
        
        process_complete = True

        # log outputs
        outfiles_dict = {
            'pruned_samples': sex_fails,
            'plink_out': out_path
        }

        metrics_dict = {
            'outlier_count': sex_fail_count
        }
    
    else:
        print('Sex Prune Failed!')
        print(f'Check {sex_tmp1}.log and {sex_tmp2}.log for more information')

        process_complete = False

        outfiles_dict = {
            'pruned_samples': 'Sex Prune Failed!',
            'plink_out': out_path
        }

        metrics_dict = {
            'outlier_count': 0
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
    
    plink_cmd1 = f"{plink2_exec} --bfile {geno_path} --geno 0.01 --maf 0.05 --indep-pairwise 50 5 0.5 --out {het_tmp}"
    plink_cmd2 = f"{plink2_exec} --bfile {geno_path} --extract {het_tmp}.prune.in --make-bed --out {het_tmp2}"
    plink_cmd3 = f"{plink2_exec} --bfile {het_tmp2} --het --out {het_tmp3}"

    cmds1 = [plink_cmd1, plink_cmd2, plink_cmd3]

    for cmd in cmds1:
        shell_do(cmd)

    listOfFiles = [f'{het_tmp}.log', f'{het_tmp2}.log', f'{het_tmp3}.log']
    concat_logs(step, out_path, listOfFiles)
    
    # check if het_tmp3 is created. if not, skip this.
    # NOTE: there may be legitimate reasons for this, for example, if only one sample in genotype file (sometimes happens in ancestry split method)
    hetpath = f'{het_tmp3}.het'
    if os.path.isfile(hetpath):
        het = pd.read_csv(hetpath, sep='\s+')
        het_outliers = het[((het.F <= -0.25) | (het.F >= 0.25))]
        outlier_count = het_outliers.shape[0]
        het_outliers.to_csv(f'{outliers_out}', sep='\t', header=True, index=False)
    
        plink_cmd4 = f"{plink2_exec} --bfile {geno_path} --remove {outliers_out} --make-bed --out {out_path}"

        shell_do(plink_cmd4)

        listOfFiles = [f'{out_path}.log']
        concat_logs(step, out_path, listOfFiles)

        if os.path.isfile(f'{out_path}.bed'):
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
            print(f'Check {out_path}.log for more information')

            outfiles_dict = {
                'pruned_samples': 'Heterozygosity Pruning Failed!',
                'plink_out': [het_tmp, het_tmp2, het_tmp3, out_path]
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


def related_prune(geno_path, out_path, related_cutoff=0.0884, duplicated_cutoff=0.354, prune_related=True, prune_duplicated=True):

    # what step are we running?
    step = "related_prune"
    print()
    print(f"RUNNING: {step}")
    print()

    # make filenames    
    grm1 = f"{out_path}_total_grm"
    grm2 = f"{out_path}_related_grm"
    grm3 = f"{out_path}_duplicated_grm"

    related_out = f"{out_path}_pairs.related"
    related_pruned_out = f"{out_path}.pruned"

    # create pfiles
    king_cmd1 = f'{plink2_exec} --bfile {geno_path} --hwe 0.0001 --mac 2 --make-pgen --out {grm1}'
    # create table of related pairs
    king_cmd2 = f'{plink2_exec} --pfile {grm1} --make-king-table --king-table-filter {related_cutoff} --out {out_path}_pairs'
    # see if any samples are related (includes duplicates)
    king_cmd3 = f'{plink2_exec} --pfile {grm1} --king-cutoff {related_cutoff} --out {grm2}' 
    # see if any samples are duplicated (grm cutoff >= 0.354)
    king_cmd4 = f'{plink2_exec} --pfile {grm1} --king-cutoff {duplicated_cutoff} --out {grm3}' 

    cmds = [king_cmd1, king_cmd2, king_cmd3, king_cmd4]
    for cmd in cmds:
        shell_do(cmd)

    listOfFiles = [f'{grm1}.log', f'{out_path}_pairs.log', f'{grm2}.log', f'{grm3}.log']
    concat_logs(step, out_path, listOfFiles)

    if os.path.isfile(f'{out_path}_pairs.kin0') and os.path.isfile(f'{grm2}.king.cutoff.out.id') and os.path.isfile(f'{grm3}.king.cutoff.out.id'):

        # create .related related pair sample files
        shutil.copy(f'{out_path}_pairs.kin0',f'{out_path}_pairs.related')

        # create .related and .duplicated single sample files
        shutil.copy(f'{grm2}.king.cutoff.out.id',f'{grm2}.related')
        related_count = sum(1 for line in open(f'{grm2}.related'))

        shutil.copy(f'{grm3}.king.cutoff.out.id',f'{grm3}.duplicated')
        duplicated_count = sum(1 for line in open(f'{grm3}.duplicated'))

        related_count = related_count - duplicated_count
        duplicated = pd.read_csv(f'{grm3}.duplicated', '\s+')

        # append duplicated sample ids to related sample ids, drop_duplicates(keep='last) because all duplicated would also be considered related
        if prune_related and prune_duplicated:
            plink_cmd1 = f'{plink2_exec} --pfile {grm1} --remove {grm2}.king.cutoff.out.id --make-bed --out {out_path}'
            shell_do(plink_cmd1) 

            related = pd.read_csv(f'{grm2}.related', '\s+')
            grm_pruned = related.append(duplicated)

            if '#FID' in grm_pruned:
                grm_pruned.rename(columns={"#FID": "FID", "#IID": "IID"}, inplace = True)
            else:
                grm_pruned['FID'] = 0
                grm_pruned.rename(columns={"#IID": "IID"}, inplace = True)
            grm_pruned.drop_duplicates(subset=['FID','IID'], keep='last', inplace=True)
            grm_pruned.to_csv(related_pruned_out, sep='\t', header=True, index=False)
            process_complete = True
        
        if prune_duplicated and not prune_related:
            plink_cmd1 = f'{plink2_exec} --pfile {grm1} --remove {grm3}.king.cutoff.out.id --make-bed --out {out_path}'
            shell_do(plink_cmd1) 

            grm_pruned = duplicated
            
            if '#FID' in grm_pruned:
                grm_pruned.rename(columns={"#FID": "FID", "#IID": "IID"}, inplace = True)
            else:
                grm_pruned['FID'] = 0
                grm_pruned.rename(columns={"#IID": "IID"}, inplace = True)
            grm_pruned.drop_duplicates(subset=['FID','IID'], keep='last', inplace=True)
            grm_pruned.to_csv(related_pruned_out, sep='\t', header=True, index=False)
            process_complete = True
            
        if not prune_related and not prune_duplicated:
            plink_cmd1 = f'echo prune_related and prune_duplicated set to False. Pruning passed'
            shell_do(plink_cmd1)

            process_complete = True
            related_pruned_out = None
            process_complete = True
            
        if not prune_duplicated and prune_related:
            print('This option is invalid. Cannot prune related without also pruning duplicated')
            process_complete = False


        if os.path.isfile('{out_path}.log'):
            listOfFiles = [f'{out_path}.log']
            concat_logs(step, out_path, listOfFiles)

        # remove intermediate files
        os.remove(f'{grm1}.pgen')
        os.remove(f'{grm1}.pvar')
        os.remove(f'{grm1}.psam')
        os.remove(f'{out_path}_pairs.kin0')

        outfiles_dict = {
            'pruned_samples': related_pruned_out,
            'related_samples': related_out,
            'plink_out': out_path
        }

        metrics_dict = {
            'related_count': related_count,
            'duplicated_count': duplicated_count
        }

    else:
        print(f'Relatedness pruning failed!')
        # print(f'Check {grm1}.log, {grm2}.log, {grm3}.log, or {out_path}_pairs.log for more information')
        print(f'Check all_plink_logs.gtlog for more information')
        process_complete = False

        outfiles_dict = {
            'pruned_samples': 'Related Pruning Failed',
            'related_samples': 0, 
            'plink_out': [grm1, grm2, grm3]
        }

        metrics_dict = {
            'related_count': 0,
            'duplicated_count': 0
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
    plink_cmd1 = f"{plink2_exec} --bfile {geno_path} --geno 0.05 --make-bed --out {geno_tmp1}"
    shell_do(plink_cmd1)

    listOfFiles = [f'{geno_tmp1}.log']
    concat_logs(step, out_path, listOfFiles)
    
    # geno pruned count
    geno_snp_count = count_file_lines(f'{geno_tmp1}.bim')
    geno_rm_count = initial_snp_count - geno_snp_count
    
    
    fam = pd.read_csv(f'{geno_path}.fam', sep='\s+', header=None, usecols=[5], names=['case'])
    # check if this contains both cases and controls
    if  all(x in fam['case'].unique() for x in [1, 2]):
        #missingness by case control (--test-missing), using P > 1E-4
        plink_cmd2 = f"{plink_exec} --bfile {geno_tmp1} --test-missing --out {mis_tmp1}"
        shell_do(plink_cmd2)

        listOfFiles = [f'{mis_tmp1}.log']
        concat_logs(step, out_path, listOfFiles)

        mis1 = f'{mis_tmp1}.missing'
        if os.path.isfile(mis1):

            mis = pd.read_csv(mis1, sep='\s+')
            exclude = mis[mis.P <= 0.0001].loc[:,'SNP']
            exclude.to_csv(f'{mis_tmp1}.exclude', sep='\t', header=False, index=False)

            plink_cmd3 = f"{plink2_exec} --bfile {geno_tmp1} --exclude {mis_tmp1}.exclude --make-bed --out {mis_tmp2}"
            shell_do(plink_cmd3)

            listOfFiles = [f'{mis_tmp2}.log']
            concat_logs(step, out_path, listOfFiles)

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
        plink_cmd4 = f"{plink_exec} --bfile {mis_tmp2} --test-mishap --out {hap_tmp1}"
        shell_do(plink_cmd4)

        # read .missing.hap file and grab flanking snps for P <= 0.0001. write flanking snps to file to exclude w bash
        mis_hap = pd.read_csv(f'{hap_tmp1}.missing.hap', sep='\s+')
        mis_hap_snps = list(mis_hap[mis_hap.P <= 0.0001].loc[:,'FLANKING'].str.split('|'))
        snp_ls_df = pd.DataFrame({'snp':[rsid for ls in mis_hap_snps for rsid in ls]})
        snp_ls_df['snp'].to_csv(f'{hap_tmp1}.exclude',sep='\t', header=False, index=False)

        plink_cmd5 = f"{plink2_exec} --bfile {mis_tmp2} --exclude {hap_tmp1}.exclude --make-bed --out {hap_tmp2}"

        # HWE from controls only using P > 1E-4
        plink_cmd6 = f"{plink_exec} --bfile {hap_tmp2} --filter-controls --hwe 1E-4 --write-snplist --out {hwe_tmp1}"
        plink_cmd7 = f"{plink2_exec} --bfile {hap_tmp2} --extract {hwe_tmp1}.snplist --make-bed --out {out_path}"

        cmds3 = [plink_cmd5, plink_cmd6, plink_cmd7]
        for cmd in cmds3:
            shell_do(cmd)

        listOfFiles = [f'{hap_tmp1}.log', f'{hap_tmp2}.log', f'{hwe_tmp1}.log', f'{out_path}.log']
        concat_logs(step, out_path, listOfFiles)

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


def miss_rates(geno_path, out_path, max_threshold=0.05):
    plink_miss_cmd = f'{plink2_exec} --bfile {geno_path} --missing --out {out_path}'

    shell_do(plink_miss_cmd)

    listOfFiles = [f'{out_path}.log']
    concat_logs(step, out_path, listOfFiles)

    # get average call rate
    lmiss = pd.read_csv(f'{out_path}.lmiss', sep='\s+')
    imiss = pd.read_csv(f'{out_path}.imiss', sep='\s+')
    avg_lmiss = lmiss.F_MISS.mean()
    avg_imiss = imiss.F_MISS.mean()
    # print(f'Average Missing Call Rate (lmiss): {avg_lmiss}')
    # print(f'Average Missing Genotyping Rate (imiss): {avg_imiss}')

    i_total = imiss.shape[0]
    thresh_list = np.arange(0.0, max_threshold+0.01, 0.01)
    
    # suggest most-stringent threshold which retains >= 90% of samples
    accept_list = []
    
    for thresh in thresh_list:
        
        i_pass = imiss.loc[imiss.F_MISS<=thresh]
        pass_prop = i_pass.shape[0]/i_total

        if pass_prop < 0.9:
            pass
        else:
            accept_list.append(thresh)
    
    if len(accept_list) > 0:
        suggested_threshold = min(accept_list)
    else:
        print('No acceptable threshold found! Try a less-stringent max_threshold')
        suggested_threshold = None
        
    metrics = {
        'avg_lmiss': avg_lmiss,
        'avg_imiss': avg_imiss,
        'suggested_threshold': suggested_threshold
    }

    return metrics


def plink_pca(geno_path, out_path, build='hg38'):

    step = 'plink_pca'
    print()
    print(f"RUNNING: {step}")
    print()

    hg19_ex_regions = """
    5 44000000 51500000 r1
    6 25000000 33500000 r2
    8 8000000 12000000 r3
    11 45000000 57000000 r4
    """

    hg38_ex_regions = """
    1   47534328    51534328    r1
    2   133742429   137242430   r2
    2   182135273   189135274   r3
    3   47458510    49962567    r4
    3   83450849    86950850    r5
    5   98664296    101164296   r6
    5   129664307   132664308   r7
    5   136164311   139164311   r8
    6   24999772    35032223    r9
    6   139678863   142178863   r10
    8   7142478 13142491    r11
    8   110987771   113987771   r12
    11  87789108    90766832    r13
    12  109062195   111562196   r14
    20  33412194    35912078    r15
    """

    # if build == 'hg19' write hg19_ex_regions exclusion regions to file named by out_path
    if build == 'hg19':
        exclusion_file = f'{out_path}_hg19.txt'
        with open(exclusion_file, 'w') as f:
            f.write(hg19_ex_regions)
    if build == 'hg38':
        exclusion_file = f'{out_path}_hg38.txt'
        with open(exclusion_file, 'w') as f:
            f.write(hg38_ex_regions)
        

    # Filter data
    filter_cmd = f"{plink2_exec} --bfile {geno_path} --maf 0.01 --geno 0.01 --hwe 5e-6 --autosome --exclude {exclusion_file} --make-bed --out {out_path}_tmp"
    shell_do(filter_cmd)
    
    # Prune SNPs
    prune_cmd = f"{plink2_exec} --bfile {out_path}_tmp --indep-pairwise 1000 10 0.02 --autosome --out {out_path}_pruned"
    shell_do(prune_cmd)

    listOfFiles = [f'{out_path}_tmp.log', f'{out_path}_pruned.log']
    concat_logs(step, out_path, listOfFiles)

    # Check if prune.in file exists
    if os.path.isfile(f'{out_path}_pruned.prune.in'):
        # Extract pruned SNPs
        extract_cmd = f"{plink2_exec} --bfile {out_path}_tmp --extract {out_path}_pruned.prune.in --make-bed --out {out_path}"
        shell_do(extract_cmd)
        
        # Calculate/generate PCs
        pca_cmd = f"{plink2_exec} --bfile {out_path} --pca --out {out_path}"
        shell_do(pca_cmd)

        # Remove intermediate files
        # os.remove(f"{out_path}_pruned.log")
        os.remove(f"{out_path}_pruned.prune.in")
        os.remove(f"{out_path}_pruned.prune.out")

        listOfFiles = [f'{out_path}.log']
        concat_logs(step, out_path, listOfFiles)
    
    # Otherwise throw an error (less than 50 samples = bad LD)
    else:
        print()
        print('PCA calculation failed!')
        print(f'Check {out_path}_pruned.log for more information.')
        print('Likely there are <50 samples for this ancestry leading to bad LD calculations.')
        print()

    # Remove intermediate files
    os.remove(f"{out_path}_tmp.bed")
    os.remove(f"{out_path}_tmp.bim")
    os.remove(f"{out_path}_tmp.fam")
    
    # os.remove(f"{out_path}_tmp.log")
    # remove exclusion file
    os.remove(exclusion_file)
