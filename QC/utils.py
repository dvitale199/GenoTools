import subprocess
import sys
import matplotlib.pyplot as plt
from matplotlib import cm 
import numpy as np
import os
import shutil
import pandas as pd

from utils.dependencies import check_plink, check_plink2

plink_exec = check_plink()
plink2_exec = check_plink2()

def shell_do(command, log=False, return_log=False):
    print(f'Executing: {(" ").join(command.split())}', file=sys.stderr)

    res=subprocess.run(command.split(), stdout=subprocess.PIPE)

    if log:
        print(res.stdout.decode('utf-8'))
    if return_log:
        return(res.stdout.decode('utf-8'))
    

def merge_genos(geno_path1, geno_path2, out_name):
    # attempt 1 at merging genos
    bash1 = f"{plink_exec} --bfile {geno_path1} --allow-no-sex --bmerge {geno_path2} --out {out_name} --make-bed"
    shell_do(bash1)

    # if {outname}-merge.missnp file created, snps need to be flipped and merge tried again
    if os.path.isfile(f'{out_name}-merge.missnp'):
        bash2 = f"{plink_exec} --bfile {geno_path1} --allow-no-sex --flip {out_name}-merge.missnp --make-bed --out {geno_path1}_flip"
        bash3 = f"{plink_exec} --bfile {geno_path1}_flip --allow-no-sex --bmerge {geno_path2} --out {out_name}_flip --make-bed"

        cmds1 = [bash2, bash3]

        for cmd in cmds1:
            shell_do(cmd)

        #if another -merge.missnp file is created, these are likely triallelic positions and must be excluded and then try merge again
        if os.path.isfile(f'{out_name}_flip-merge.missnp'):
            bash4 = f"{plink_exec} --bfile {geno_path1}_flip --allow-no-sex --exclude {out_name}_flip-merge.missnp --out {geno_path1}_flip_pruned --make-bed"
            bash5 = f"{plink_exec} --bfile {geno_path1}_flip_pruned --allow-no-sex --bmerge {geno_path2} --out  {out_name} --make-bed"

            cmds2 = [bash4, bash5]

            for cmd in cmds2:
                shell_do(cmd)

        # if second attempt at merge is successful, there are no triallelic snps and we can go ahead and move the _flip files to our _merged_ref_panel filenames 
        # for further processing
        else:
            suffix_list = ['bed','bim','fam']
            for suffix in suffix_list:
                shutil.copy(f'{out_name}_flip.{suffix}',f'{out_name}.{suffix}')

    # if first merge works, there are no alleles in need of flipping and no triallelic positions and we can proceed!
    else:
        pass

        
def ld_prune(geno_path, out_name, window_size=1000, step_size=50, rsq_thresh=0.05):
    # now prune for LD
    ld_prune1 = f'{plink_exec} --bfile {geno_path} --allow-no-sex --indep-pairwise {window_size} {step_size} {rsq_thresh} --autosome --out {geno_path}_pruned_data'
    ld_prune2 = f'{plink_exec} --bfile {geno_path} --allow-no-sex --extract {geno_path}_pruned_data.prune.in --make-bed --out {out_name}'

    ld_prune_cmds = [ld_prune1, ld_prune2]

    for cmd in ld_prune_cmds:
        shell_do(cmd)

        
def random_sample_snps(geno_path, out_name, n=10000):
    rand_samp_snplist = f'{geno_path}_rand_samp.snplist'
    bim = pd.read_csv(f'{geno_path}.bim', sep='\t', header=None)
    ref_panel_random_sample = bim.sample(n, random_state=123)
    ref_panel_random_sample.to_csv(rand_samp_snplist, header=False, index=False, sep='\t')

    rand_sample_cmd = f'{plink_exec} --bfile {geno_path} --allow-no-sex --extract {rand_samp_snplist} --autosome --make-bed --out {out_name}'
    
    shell_do(rand_sample_cmd)
    

def get_common_snps(geno_path1, geno_path2, out_name):  
    
    """
    Gets common snps between 2 genotype files and extracts from geno_path1. outputs plink bed/bim/fam file 
    for geno_path1 with only matching snps from geno_path2
    """
   
    print('Getting Common SNPs')	
 
    bim1 = pd.read_csv(f'{geno_path1}.bim', sep='\t', header=None)
    bim1.columns = ['chr', 'rsid', 'kb', 'pos', 'a1', 'a2']
    bim2 = pd.read_csv(f'{geno_path2}.bim', sep='\t', header=None)
    bim2.columns = ['chr', 'rsid', 'kb', 'pos', 'a1', 'a2']
    
    bim1['rsid'].to_csv(f'{geno_path1}.snplist', sep='\t', header=None, index=None)

    bim1['merge_id'] = bim1['chr'].astype(str) + ':' + bim1['pos'].astype(str) + ':' + bim1['a2'] + ':' + bim1['a1']
    bim2['merge_id1'] = bim2['chr'].astype(str) + ':' + bim2['pos'].astype(str) + ':' + bim2['a2'] + ':' + bim2['a1']
    bim2['merge_id2'] = bim2['chr'].astype(str) + ':' + bim2['pos'].astype(str) + ':' + bim2['a1'] + ':' + bim2['a2']

    common_snps1 = bim2[['rsid','merge_id1','a1','a2']].merge(bim1, how='inner', left_on=['merge_id1'], right_on=['merge_id'])
    common_snps2 = bim2[['rsid','merge_id2','a1','a2']].merge(bim1, how='inner', left_on=['merge_id2'], right_on=['merge_id'])	
    common_snps = pd.concat([common_snps1, common_snps2], axis=0)
    
    flip_cmd = f'{plink_exec} --bfile {geno_path1} --flip {geno_path1}.snplist --make-bed --out {geno_path1}_flip'
    shell_do(flip_cmd)

    bim1_flip = pd.read_csv(f'{geno_path1}_flip.bim', sep='\t', header=None)
    bim1_flip.columns = ['chr', 'rsid', 'kb', 'pos', 'a1', 'a2']

    bim1_flip['merge_id'] = bim1_flip['chr'].astype(str) + ':' + bim1_flip['pos'].astype(str) + ':' + bim1_flip['a2'] + ':' + bim1_flip['a1']
    common_snps1 = bim2[['rsid','merge_id1','a1','a2']].merge(bim1_flip, how='inner', left_on=['merge_id1'], right_on=['merge_id'])
    common_snps2 = bim2[['rsid','merge_id2','a1','a2']].merge(bim1_flip, how='inner', left_on=['merge_id2'], right_on=['merge_id'])	

    common_snps = pd.concat([common_snps, common_snps1, common_snps2], axis=0)
    common_snps = common_snps.drop_duplicates(subset=['rsid_y'], ignore_index=True)

    common_snps_file = f'{out_name}.common_snps'
    common_snps['rsid_y'].to_csv(f'{common_snps_file}', sep='\t', header=False, index=False)
    
    common_snps['rsid_x'].to_csv(f'{out_name}_geno.common_snps', sep='\t', header=False, index=False)
    
    ext_snps_cmd = f'{plink2_exec} --bfile {geno_path1} --extract {common_snps_file} --make-bed --out {out_name}'
    
    shell_do(ext_snps_cmd)

    outfiles = {
        'common_snps': common_snps_file,
        'bed': out_name
    }
    return outfiles


def rm_tmps(tmps, suffixes=None):
    
    if suffixes:
        suffixes=suffixes
    else:
        suffixes=[
            'hh','log','nosex','bed','bim','fam',
            'prune.in','prune.out','sexcheck','het',
            'grm.bim','grm.id','grm.N.bim',
            'missing','missing.hap','exclude','snplist'
        ]

    print()
    print("REMOVING TEMPORARY FILES")
    for tmp in tmps:
        for suf in suffixes:
            tmpfile = f'{tmp}.{suf}'
            try:
                os.remove(tmpfile)
            except OSError:
                pass
            # old method below... remove eventually
            # if os.path.isfile(tmpfile):
            #     os.remove(tmpfile)
            #     print(f"REMOVED: {tmpfile}")
            # else:
            #     pass
    print()


def count_file_lines(file_path):
    
    return sum(1 for line in open(file_path))

