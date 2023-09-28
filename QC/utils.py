import subprocess
import sys
import os
import shutil
import pandas as pd
import warnings 
import numpy as np

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


def upfront_check(geno_path):
    fam = pd.read_csv(f'{geno_path}.fam', header=None, sep = '\s+', 
                      names = ['FID', 'IID', 'Paternal_ID', 'Maternal_ID', 'Sex', 'Phenotype'])
    bim = pd.read_csv(f'{geno_path}.bim', header = None, sep = '\s+', low_memory = False,
                      names = ['CHR', 'SNP_ID', 'Genetic_Distance', 'Position', 'A1', 'A2'])

    sex_counts = fam['Sex'].value_counts().to_dict()
    pheno_counts = fam['Phenotype'].value_counts().to_dict()
    chr_counts = bim['CHR'].value_counts().to_dict()

    # print breakdown of data 
    print("Your data has the following breakdown:")
    print("- Genetic Sex:")
    for sex in sex_counts.keys():
        if sex == 1:
            print(f'{sex_counts[sex]} Males \n')
        if sex == 2:
            print(f'{sex_counts[sex]} Females \n')
        if sex == 0 or sex == -9:
            print(f'{sex_counts[sex]} Unknown \n')

    print("- Phenotypes:")
    for pheno in pheno_counts.keys():
        if pheno == 2:
            print(f'{pheno_counts[pheno]} Cases \n')
        if pheno == 1:
            print(f'{pheno_counts[pheno]} Controls \n')
        if pheno == 0 or pheno == -9:
            print(f'{pheno_counts[pheno]} Missing \n')
            # if pheno_counts[pheno] > 50:  # may want to createe a threshold for how many missing phenos allowed
            #     warnings.warn("You are missing too many phenotypes (over 50).")

    # check for items necessary to the pipeline
    if len(fam) < 50:
        warnings.warn("You do not have enough total participants in your data (less than 50).")
    elif 1 not in pheno_counts.keys() and 2 not in pheno_counts.keys():
        warnings.warn("You are missing at least 1 phenotype class in your data.")
    elif "23" not in chr_counts and "X" not in chr_counts:
        warnings.warn("No sign of sex chromosome in your bim file.")
    elif 1 in sex_counts.keys() or 2 in sex_counts.keys():
        print("You are all set to continue in the pipeline!")


def replace_all(text, dict):
    # replaces any applicable dictionary elements in text
    for i, j in dict.items():
        text = text.replace(i, j)
    return text


def process_log(out_dir, concat_log):
    # exclude lines containing this information from log file
    exclude = ['Hostname', 'Working directory', 'Intel', 'Start time', 'Random number seed', 'RAM detected', 'threads', 'thread', 
    'written to', 'done.', 'End time:', 'Writing', '.bed', '.bim', '.fam', '.id', '.hh', '.sexcheck', '.psam', '-bit',
    '.pvar', '.pgen', '.in', '.out', '.het', '.missing', '.snplist', '.kin0', '.eigenvec', '.eigenval', '(--maf/', 'Step:']

    # save all indices in log file where these instances occur
    step_indices = [i for i, s in enumerate(concat_log) if 'Step:' in s]
    process_indices = [i for i, s in enumerate(concat_log) if 'Process:' in s]
    bfile_indices = [i for i, s in enumerate(concat_log) if '--bfile' in s]
    pfile_indices = [i for i, s in enumerate(concat_log) if '--pfile' in s]

    # combine bfile and pfile indices
    bfile_indices.extend(pfile_indices)
    bfile_indices.sort()

    # add final index of log to traverse entire log
    step_indices.append(len(concat_log))
    start = 0
    stop = 1
    
    # list ancestry only for the following steps
    ancestry_steps = ['related', 'het', 'variant', 'pca']
    # exclude/replace from text
    fillers = ['and', '.']
    replace = {'loaded from': 'loaded', '(see': '', ');': ';'}

    # write final processed log
    with open(f"{out_dir}/cleaned_genotools.log", "w") as f:
        while start < len(step_indices)-1:
            # list step and process names
            step_line = concat_log[step_indices[start]]
            process_line = concat_log[process_indices[start]]
            process_name = process_line.split('_')[0].replace('Process: ', '')
            bfile_line = concat_log[bfile_indices[start]]
            bfile_name = bfile_line.split()[1].replace(out_dir, "")

            # split by second part of plink_pca process name
            if process_name == 'plink':
                process_name = 'pca'

            # split full step name by process name to get better label components & add back in
            step_name = step_line.split(f'{process_name}')[-1]
            step_name = (process_name + step_name).replace('.log', '')

            # find ancestry acronym in bfile input line
            if bfile_name.split("_")[-1].isupper():
                ancestry_check = bfile_name.split("_")[-1]

            # write final labels for concise step name & ancestry of focus
            f.write(f'Step: {step_name}')
            if process_name in ancestry_steps:
                f.write(f'Ancestry: {ancestry_check}\n')
            
            # rewrite concatenated log section by section with exclusion criteria
            for i in range(step_indices[start], step_indices[stop]):
                if "error" in concat_log[i].strip('\n'):
                    f.write(concat_log[i])
                elif concat_log[i].strip('\n') in fillers:
                    pass
                elif len(concat_log[i]) == 1:
                    pass
                elif not any([x in concat_log[i].strip('\n') for x in exclude]):
                    content = replace_all(concat_log[i], replace)
                    f.write(content)

            f.write('\n')
            f.write('\n')
            start += 1
            stop += 1


def concat_logs(step, out_path, listOfFiles):
    # stores concat log in processing directory
    out_dir = os.path.dirname(os.path.abspath(out_path))

    # combine log files into 1 file
    # when transition to Classes: clear log on every new run
    with open(f'{out_dir}/all_plink_logs.log', "a+") as new_file:
        for name in listOfFiles:
            with open(name) as file:
                new_file.write(f'Step: {name}\n')
                new_file.write(f'Process: {step}\n')
                for line in file:
                    new_file.write(line)
        
                new_file.write("\n")

    # remove intermediate log files 
    for files in listOfFiles:
        os.remove(files)

    with open(f'{out_dir}/all_plink_logs.log', 'r') as file:
        process_log(out_dir, file.readlines())


def label_bim_with_genes(bim_file, gene_reference=None, locus_size=1000000):
    """Label SNPs with the gene they are in."""
    # Check if the gene reference file exists
    if gene_reference is None:
        # Get the directory of the current function
        function_dir = os.path.dirname(os.path.abspath(__file__))
        # Get the directory of the file
        ref_dir = os.path.join(function_dir, '..', 'ref')
        # Create the full file path
        my_file_path = os.path.join(ref_dir, 'glist-hg38')
    if not os.path.exists(gene_reference):
        raise FileNotFoundError(f"{gene_reference} not found")

    # Load SNP data from bim file
    snps = pd.read_table(bim_file, sep='\s+', header=None, names=['chr', 'snp_id', 'cm_pos', 'pos', 'a1', 'a2'], dtype={'chr': str})
    # Load gene information
    glist = pd.read_table(gene_reference, sep='\s+', header=None, names=['chr', 'start', 'end', 'name'], dtype={'chr': str})
    
    # convert glist chr X to 23, Y to 24, XY to 25
    glist.loc[glist.chr=='X', 'chr'] = '23'
    glist.loc[glist.chr=='Y', 'chr'] = '24'
    glist.loc[glist.chr=='XY', 'chr'] = '25'

    # Add a new column to hold the gene labels
    snps['gene'] = 'NA'

    # Loop through each row in the gene list
    for _, row in glist.iterrows():

        start = row['start'] - locus_size
        stop = row['end'] + locus_size
        # Find the positions that fall within the current start and stop values
        include_snps = snps.loc[(snps['chr'] == row['chr']) & 
                                (snps['pos'] >= start) & 
                                (snps['pos'] <= stop)].copy()

        # Assign gene name to included SNPs
        include_snps.loc[:, 'gene'] = row['name']

        # Update the label for the included SNPs
        snps.update(include_snps)

    return snps


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
    

def get_common_snps(geno_path1, geno_path2, out_name):  
    
    """
    Gets common snps between 2 genotype files and extracts from geno_path1. outputs plink bed/bim/fam file 
    for geno_path1 with only matching snps from geno_path2
    """
   
    print('Getting Common SNPs')	
    
    # read both bim files
    bim1 = pd.read_csv(f'{geno_path1}.bim', sep='\t', header=None)
    bim1.columns = ['chr', 'rsid', 'kb', 'pos', 'a1', 'a2']
    bim2 = pd.read_csv(f'{geno_path2}.bim', sep='\t', header=None)
    bim2.columns = ['chr', 'rsid', 'kb', 'pos', 'a1', 'a2']
    
    # write bim 1 ids to snplist
    bim1['rsid'].to_csv(f'{geno_path1}.snplist', sep='\t', header=None, index=None)

    # creating merge ids
    bim1['merge_id'] = bim1['chr'].astype(str) + ':' + bim1['pos'].astype(str) + ':' + bim1['a2'] + ':' + bim1['a1']
    bim2['merge_id1'] = bim2['chr'].astype(str) + ':' + bim2['pos'].astype(str) + ':' + bim2['a2'] + ':' + bim2['a1']
    bim2['merge_id2'] = bim2['chr'].astype(str) + ':' + bim2['pos'].astype(str) + ':' + bim2['a1'] + ':' + bim2['a2']

    # two merges and concatenation
    common_snps1 = bim2[['rsid','merge_id1','a1','a2']].merge(bim1, how='inner', left_on=['merge_id1'], right_on=['merge_id'])
    common_snps2 = bim2[['rsid','merge_id2','a1','a2']].merge(bim1, how='inner', left_on=['merge_id2'], right_on=['merge_id'])	
    common_snps = pd.concat([common_snps1, common_snps2], axis=0)
    
    # flip and merge again
    flip_cmd = f'{plink_exec} --bfile {geno_path1} --flip {geno_path1}.snplist --make-bed --out {geno_path1}_flip'
    shell_do(flip_cmd)

    bim1_flip = pd.read_csv(f'{geno_path1}_flip.bim', sep='\t', header=None)
    bim1_flip.columns = ['chr', 'rsid', 'kb', 'pos', 'a1', 'a2']

    bim1_flip['merge_id'] = bim1_flip['chr'].astype(str) + ':' + bim1_flip['pos'].astype(str) + ':' + bim1_flip['a2'] + ':' + bim1_flip['a1']
    common_snps1 = bim2[['rsid','merge_id1','a1','a2']].merge(bim1_flip, how='inner', left_on=['merge_id1'], right_on=['merge_id'])
    common_snps2 = bim2[['rsid','merge_id2','a1','a2']].merge(bim1_flip, how='inner', left_on=['merge_id2'], right_on=['merge_id'])	

    # concat merges and drop duplicates
    common_snps = pd.concat([common_snps, common_snps1, common_snps2], axis=0)
    common_snps = common_snps.drop_duplicates(subset=['chr','pos'], ignore_index=True)

    # write snps to txt and extract
    common_snps_file = f'{out_name}.common_snps'
    common_snps['rsid_y'].to_csv(f'{common_snps_file}', sep='\t', header=False, index=False)
    
    ext_snps_cmd = f'{plink2_exec} --bfile {geno_path1} --extract {common_snps_file} --make-bed --out {out_name}'
    shell_do(ext_snps_cmd)

    # return outfiles
    outfiles = {
        'common_snps': common_snps_file,
        'bed': out_name
    }
    return outfiles


def rm_tmps(tmps, suffixes=None):
    # OVERHAUL
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


def plink_pca(geno_path, out_path, build='hg38'):

    step = 'plink_pca'

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
    try:
        os.remove(f"{out_path}_tmp.bed")
        os.remove(f"{out_path}_tmp.bim")
        os.remove(f"{out_path}_tmp.fam")

        # Remove exclusion file
        os.remove(exclusion_file)

    except OSError:
        pass


def miss_rates(geno_path, out_path, max_threshold=0.05):

    plink_miss_cmd = f'{plink2_exec} --bfile {geno_path} --missing --out {out_path}'

    shell_do(plink_miss_cmd)

    # listOfFiles = [f'{out_path}.log']
    # concat_logs(step, out_path, listOfFiles)

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
