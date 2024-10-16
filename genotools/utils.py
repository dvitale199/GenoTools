# Copyright 2023 The GenoTools Authors. All Rights Reserved.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
# ==============================================================================


import subprocess
import sys
import os
import shutil
import pandas as pd
import warnings
import numpy as np
from scipy.stats import norm
from genotools.dependencies import check_plink, check_plink2

plink_exec = check_plink()
plink2_exec = check_plink2()

def gt_header():

    header = """
     ██████╗ ███████╗███╗  ██╗ █████╗ ████████╗ █████╗  █████╗ ██╗      ██████╗
    ██╔════╝ ██╔════╝████╗ ██║██╔══██╗╚══██╔══╝██╔══██╗██╔══██╗██║     ██╔════╝
    ██║  ██╗ █████╗  ██╔██╗██║██║  ██║   ██║   ██║  ██║██║  ██║██║     ╚█████╗
    ██║  ╚██╗██╔══╝  ██║╚████║██║  ██║   ██║   ██║  ██║██║  ██║██║      ╚═══██╗
    ╚██████╔╝███████╗██║ ╚███║╚█████╔╝   ██║   ╚█████╔╝╚█████╔╝███████╗██████╔╝
    ╚═════╝ ╚══════╝╚═╝  ╚══╝ ╚════╝    ╚═╝    ╚════╝  ╚════╝ ╚══════╝╚═════╝
    """
    return header


def shell_do(command, print_cmd=False, log=False, return_log=False, err=False):
    if print_cmd:
        print(f'Executing: {(" ").join(command.split())}', file=sys.stderr)

    res=subprocess.run(command.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    output = res.stdout.decode('utf-8') + res.stderr.decode('utf-8')

    if log:
        print(output)
    if return_log:
        return output
    if err:
        return res.stderr.decode('utf-8')


def bfiles_to_pfiles(bfile_path=None, pfile_path=None):
    # check if both are none
    if not (bfile_path or pfile_path):
        print()
        print('ERROR: Need either PLINK1.9 or PLINK2 binaries!')
        print()

    elif bfile_path and pfile_path:
        print()
        print('ERROR: Cannot accept both PLINK1.9 and PLINK2 binaries simulaneously!')
        print()

    elif bfile_path and (not pfile_path):
        if not os.path.isfile(f'{bfile_path}.bed'):
            raise FileNotFoundError(f'{bfile_path} does not exist.')

        convert_cmd = f'{plink2_exec} --bfile {bfile_path} --make-pgen psam-cols=fid,parents,sex,pheno1,phenos --out {bfile_path}'
        shell_do(convert_cmd)

        if os.path.isfile(f'{bfile_path}.log'):
            os.remove(f'{bfile_path}.log')

    else:
        if not os.path.isfile(f'{pfile_path}.pgen'):
            raise FileNotFoundError(f'{pfile_path} does not exist.')

        convert_cmd = f'{plink2_exec} --pfile {pfile_path} --make-bed --out {pfile_path}'
        shell_do(convert_cmd)

        if os.path.isfile(f'{pfile_path}.log'):
            os.remove(f'{pfile_path}.log')


def vcf_to_pfiles(vcf_path):
    if not os.path.isfile(vcf_path):
        raise FileNotFoundError(f'{vcf_path} does not exist.')

    prefix = vcf_path.split('.vcf')[0]

    convert_cmd1 = f'{plink2_exec} --vcf {vcf_path} --make-bed --out {prefix}'
    shell_do(convert_cmd1)

    if not os.path.isfile(f'{prefix}.bed'):
        raise FileNotFoundError(f'{prefix} bed/bim/fam files do not exist. Conversion from VCF failed')

    bfiles_to_pfiles(bfile_path=prefix)

    if os.path.isfile(f'{prefix}.pgen'):
        os.remove(f'{prefix}.bed')
        os.remove(f'{prefix}.bim')
        os.remove(f'{prefix}.fam')
    else:
        raise FileNotFoundError(f'{prefix} pgen/pvar/psam files do not exist. Conversion from bed/bim/fam failed.')



def upfront_check(geno_path, args):
    if os.path.isfile(f'{args["out"]}_all_logs.log') and not args['skip_fails']:
        raise ValueError(f'{args["out"]}_all_logs.log, which means the pipeline has previously been run on this output file!\n \
                         Please rerun with "--skip_fails" True flag to ignore this, or write output to a new file name.')

    if not os.path.isfile(f'{geno_path}.pgen'):
        raise FileNotFoundError(f"{geno_path} does not exist.")

    # if no pgen present, but bed is present, and skip fails is True, convert to pgen
    if not os.path.isfile(f'{geno_path}.pgen') and os.path.isfile(f'{geno_path}.bed') and (not args['skip_fails']):
        warnings.warn(f'{geno_path} exists but it is in PLINK1.9 binary format. Converting to PLINK2 binaries...', stacklevel=2)
        bfiles_to_pfiles(bfile_path=geno_path)

    sam = pd.read_csv(f'{geno_path}.psam', sep = '\s+')
    var = pd.read_csv(f'{geno_path}.pvar', delimiter='\t', comment='#', header=None, names=['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'FILTER', 'INFO'], low_memory=False)

    if 'SEX' not in sam.columns:
        raise KeyError(f'{geno_path}.psam is missing SEX column. Even if no SEX information is present, GenoTools requires a SEX column.')

    if 'PHENO1' not in sam.columns:
        raise KeyError(f'{geno_path}.psam is missing PHENO1 column. Even if no PHENO1 information is present, GenoTools requires a PHENO1 column.')

    sex_counts = sam['SEX'].value_counts().to_dict()
    pheno_counts = sam['PHENO1'].value_counts().to_dict()
    chr_counts = var['#CHROM'].value_counts().to_dict()

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

    if args['skip_fails']:
        return args

    else:
        # skip sex check when no sex in fam or no X chromosome
        if args['sex'] is not None:
            if (1 not in sex_counts.keys()) and (2 not in sex_counts.keys()):
                warnings.warn('You tried calling sex prune but no sample sex data is available. Skipping...', stacklevel=2)
                args['sex'] = None
            elif ('23' not in chr_counts.keys()) and ('X' not in chr_counts.keys()):
                warnings.warn('You tried calling sex prune but no X chromosome data is available. Skipping...', stacklevel=2)
                args['sex'] = None

        # change hwe prune to be run without controls filtered when no controls present
        if (args['hwe'] is not None) and (args['filter_controls'] == True) and (1 not in pheno_counts.keys()):
            warnings.warn('You tried calling hwe prune with controls filtered but no controls are available. Skipping...', stacklevel=2)
            args['filter_controls'] = False

        # skip case control when called without cases or controls present
        if (args['case_control'] is not None) and ((1 not in pheno_counts.keys()) or (2 not in pheno_counts.keys())):
            warnings.warn('You tried calling case-control prune but only cases or controls are available, not both. Skipping...', stacklevel=2)
            args['case_control'] = None

        # skip het prune if less than 50 samples are present
        if (args['het'] is not None) and (var.shape[0] < 50):
            warnings.warn('You tried calling het prune with less than 50 samples. Skipping...', stacklevel=2)
            args['het'] = None

        return args


def replace_all(text, dict):
    # replaces any applicable dictionary elements in text
    for i, j in dict.items():
        text = text.replace(i, j)
    return text


def process_log(out_path, concat_log):
    out_dir = os.path.dirname(os.path.abspath(out_path))

    # ancestry labels
    ancestries = ['AFR', 'SAS', 'EAS', 'EUR', 'AMR', 'AJ', 'CAS', 'MDE', 'FIN', 'AAC']

    # exclude lines containing this information from log file
    exclude = ['Hostname', 'Working directory', 'Intel', 'Start time', 'Random number seed', 'RAM detected', 'threads', 'thread',
    'written to', 'done.', 'End time:', 'Writing', '.bed', '.bim', '.fam', '.id', '.hh', '.sexcheck', '.psam', '-bit', 'from',
    '.pvar', '.pgen', '.in', '.out', '.het', '.missing', '.snplist', '.kin0', '.eigenvec', '.eigenval', '(--maf/', '--make-bed to','+']

    # save all indices in log file where these instances occur
    step_indices = [i for i, s in enumerate(concat_log) if 'Log:' in s]
    out_indices = [i for i, s in enumerate(concat_log) if any(x in s for x in ('--prefix', '--out'))]
    ancestry_ran = [True if 'split_cohort_ancestry' in line else False for line in concat_log]

    # add final index of log to traverse entire log
    step_indices.append(len(concat_log))
    start = 0
    stop = 1

    # exclude/replace from text
    fillers = ['and', '.']
    replace = {'loaded from': 'loaded', '(see': '', ');': ';'}

    # write final processed log
    with open(f"{out_path}_cleaned_logs.log", "w") as f:
        header = gt_header()
        f.write(header)
        f.write("\n")

        while start < len(step_indices)-1:
            # list step and process names
            out_line =  concat_log[out_indices[start]]
            out_name = out_line.split()[1].replace(out_dir, "")

            # write final labels for concise step name & ancestry of focus
            if any(ancestry_ran):
                # find ancestry acronym in output line
                ancestry_tokens = [s for s in out_name.split("_") if s in ancestries]
                if len(ancestry_tokens) >= 1:
                    f.write(f'Ancestry: {ancestry_tokens[-1]}\n')

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


def concat_logs(step, out_path, listOfFiles, default_dir=None):
    if not os.path.isabs(out_path):
        out_path = os.path.abspath(out_path)

    if '_tmp/' in out_path:
        out_dir = os.path.dirname(os.path.dirname(out_path))
    else:
        out_dir = os.path.dirname(out_path)
    
    if not out_dir:
        if default_dir:
            out_dir = default_dir
        else:
            raise ValueError("Output directory is empty and no default directory provided.")

    if not os.path.exists(out_dir):
        raise FileNotFoundError(f"Output directory does not exist: {out_dir}")

    log_paths = []
    for file in os.listdir(out_dir):
        if file.endswith("_all_logs.log"):
            log_paths.append(file)

    if len(log_paths) == 0:
        log_path = os.path.join(out_dir, os.path.split(out_path)[1] + "_all_logs.log")

    elif len(log_paths) == 1:
        log_path = os.path.join(out_dir, log_paths[0])

    else:
        mtimes = {path: os.path.getmtime(os.path.join(out_dir, path)) for path in log_paths}
        most_recent_log = max(mtimes, key=mtimes.get)
        log_path = os.path.join(out_dir, most_recent_log)

    with open(log_path, "a+") as new_file:
        for name in listOfFiles:
            with open(name) as file:
                new_file.write(f'Log: {name}\n')
                new_file.write(f'Process: {step}\n')
                for line in file:
                    new_file.write(line)
                new_file.write("\n")
                
    for files in listOfFiles:
        os.remove(files)

    with open(log_path, 'r') as file:
        out_path = log_path.replace('_all_logs.log', '')
        process_log(out_path, file.readlines())


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


def miss_rates(geno_path, out_path, max_threshold=0.05):

    plink_miss_cmd = f'{plink2_exec} --pfile {geno_path} --missing --out {out_path}'

    shell_do(plink_miss_cmd)

    # listOfFiles = [f'{out_path}.log']
    # concat_logs(step, out_path, listOfFiles)

    # get average call rate
    vmiss = pd.read_csv(f'{out_path}.vmiss', sep='\s+')
    smiss = pd.read_csv(f'{out_path}.smiss', sep='\s+')
    avg_vmiss = vmiss.F_MISS.mean()
    avg_smiss = smiss.F_MISS.mean()
    # print(f'Average Missing Call Rate (lmiss): {avg_lmiss}')
    # print(f'Average Missing Genotyping Rate (imiss): {avg_imiss}')

    s_total = smiss.shape[0]
    thresh_list = np.arange(0.0, max_threshold+0.01, 0.01)

    # suggest most-stringent threshold which retains >= 90% of samples
    accept_list = []

    for thresh in thresh_list:

        s_pass = smiss.loc[smiss.F_MISS<=thresh]
        pass_prop = s_pass.shape[0]/s_total

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
        'avg_lmiss': avg_vmiss,
        'avg_imiss': avg_smiss,
        'suggested_threshold': suggested_threshold
    }

    return metrics


def zscore_pval_conversion(zscores=None, pvals=None, stats=None):

    # neither zscore or pvals provided provided
    if zscores is None and pvals is None:
        print('Conversion Failed!')
        print('Either p-values or z-scores must be provided')

    # both zscores and pvals provided
    elif zscores is not None and pvals is not None:
        print('Conversion Failed!')
        print('Provide only p-values or z-scores, not both')

    # pvals provided but stats not provided to determine sign of zscore
    elif pvals is not None and stats is None:
        print('Conversion Failed!')
        print('Stats must be provided when going from p-values to z-scores')

    else:
        # convert pvals to zscores using stats to get proper sign
        if zscores is None:
            z = np.where(stats > 0, norm.isf(pvals/2), -norm.isf(pvals/2))
            return z
        # convert zscores to pvals
        if pvals is None:
            p = 2*norm.sf(abs(zscores))
            return p