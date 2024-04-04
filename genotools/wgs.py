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


import pandas as pd
import numpy as np
import os
import shutil
import platform
from genotools.utils import shell_do, count_file_lines, concat_logs, bfiles_to_pfiles
from genotools.dependencies import check_plink, check_plink2, check_king

plink_exec = check_plink()
plink2_exec = check_plink2()
king_exec = check_king()

class WholeGenomeSeqQC:

    def __init__(self, geno_path=None, out_path=None):
        self.geno_path = geno_path
        self.out_path = out_path


    def run_depth_prune(self, depth_min=10):

        """
        Prunes SNPs based on depth.

        Parameters:
        - depth_min (int, optional): Depth threshold where all genotype calls below are excluded.

        Returns:
        - dict: A structured dictionary containing:
            * 'pass': Boolean indicating the successful completion of the process.
            * 'step': The label for this procedure ('wgs_depth_prune').
            * 'metrics': Metrics associated with the pruning, such as 'depth_rm_count'.
            * 'output': Dictionary containing paths to the generated output files.
        """

        geno_path = self.geno_path
        out_path = self.out_path

        step = "wgs_depth_prune"

        # get initial snps count
        initial_snp_count = count_file_lines(f'{geno_path}.pvar') - 1

        # remove variants where DP field is less than desired depth
        plink_cmd = f"{plink2_exec} --pfile {geno_path} --extract-if-info DP >= {depth_min} --make-pgen --out {out_path}"
        # plink_cmd = f"{plink2_exec} --pfile {geno_path} --vcf-min-dp {depth_min} --out {out_path}"
        shell_do(plink_cmd)

        listOfFiles = [f'{out_path}.log']
        concat_logs(step, out_path, listOfFiles)

        # filter depth count
        depth_snp_count = count_file_lines(f'{out_path}.pvar') - 1
        depth_rm_count = initial_snp_count - depth_snp_count

        process_complete = True

        outfiles_dict = {
            'plink_out': out_path
        }

        metrics_dict = {
            'depth_rm_count': depth_rm_count
        }

        out_dict = {
            'pass': process_complete,
            'step': step,
            'metrics': metrics_dict,
            'output': outfiles_dict
        }

        return out_dict


    def run_filter_prune(self):

        """
        Prunes SNPs based on whether or not they PASS all input filters.

        Returns:
        - dict: A structured dictionary containing:
            * 'pass': Boolean indicating the successful completion of the process.
            * 'step': The label for this procedure ('wgs_filter_prune').
            * 'metrics': Metrics associated with the pruning, such as 'filter_rm_count'.
            * 'output': Dictionary containing paths to the generated output files.
        """

        geno_path = self.geno_path
        out_path = self.out_path

        step = "wgs_filter_prune"

        # get initial snp count
        initial_snp_count = count_file_lines(f'{geno_path}.pvar') - 1

        # variants that failed one or more filters tracked by FILTER field
        plink_cmd1 = f"{plink2_exec} --pfile {geno_path} --var-filter --make-pgen --out {out_path}"
        shell_do(plink_cmd1)

        listOfFiles = [f'{out_path}.log']
        concat_logs(step, out_path, listOfFiles)

        # filter pruned count
        filter_snp_count = count_file_lines(f'{out_path}.pvar') - 1
        filter_rm_count = initial_snp_count - filter_snp_count

        process_complete = True

        outfiles_dict = {
            'plink_out': out_path
        }

        metrics_dict = {
            'filter_rm_count': filter_rm_count
        }

        out_dict = {
            'pass': process_complete,
            'step': step,
            'metrics': metrics_dict,
            'output': outfiles_dict
        }

        return out_dict


    def run_freemix_check(self, preBqsr_path):
        # preBqsr_selfSM.tsv file?
        geno_path = self.geno_path
        out_path = self.out_path

        step = "wgs_filter_prune"

        # Does ZH manually check something by creating a 'CHECK' col?
        # preBqsr['CHECK'] = preBqsr.FREELK1 - preBqsr.FREELK0

        # create filenames
        freemix_fails = f'{out_path}.freemix_fails'

        preBqsr = pd.read_csv(f'{preBqsr_path}', sep='\s+')
        freemix_fail_df = preBqsr[preBqsr.FREEMIX >= 0.03]
        freemix_fail_ids = freemix_fail_df['sample_id']
        freemix_fail_count = freemix_fail_ids.shape[0]
        freemix_fail_ids.to_csv(freemix_fails, sep='\t', header=False, index=False)

        plink_cmd = f'{plink2_exec} --pfile {geno_path} --remove {freemix_fails} --make-pgen psam-cols=fid,parents,sex,pheno1,phenos --out {out_path}'
        shell_do(plink_cmd)

        listOfFiles = [f'{out_path}.log']
        concat_logs(step, out_path, listOfFiles)

        process_complete = True

        outfiles_dict = {
            'pruned_samples': f'{freemix_fails}',
            'plink_out': f'{out_path}',
        }

        metrics_dict = {
            'outlier_count': freemix_fail_count
        }

        out_dict = {
            'pass': process_complete,
            'step': step,
            'metrics': metrics_dict,
            'output': outfiles_dict
        }

        return out_dict


    def run_coverage_check(self, wgs_metrics_path, min_mean_coverage=25):
        # wgs_metrics.tsv file?
        geno_path = self.geno_path
        out_path = self.out_path

        step = "coverage_check"

        # create filenames
        coverage_fails = f'{out_path}.coverage_fails'

        wgs_metrics = pd.read_csv(f'{wgs_metrics_path}', sep='\s+')
        coverage_fail_df = wgs_metrics[wgs_metrics.MEAN_COVERAGE < min_mean_coverage]
        coverage_fail_ids = coverage_fail_df['sample_id']
        coverage_fail_count = coverage_fail_ids.shape[0]
        coverage_fail_ids.to_csv(coverage_fails, sep='\t', header=False, index=False)

        plink_cmd = f'{plink2_exec} --pfile {geno_path} --remove {coverage_fails} --make-pgen psam-cols=fid,parents,sex,pheno1,phenos --out {out_path}'
        shell_do(plink_cmd)

        listOfFiles = [f'{out_path}.log']
        concat_logs(step, out_path, listOfFiles)

        process_complete = True

        outfiles_dict = {
            'pruned_samples': f'{coverage_fails}',
            'plink_out': f'{out_path}',
        }

        metrics_dict = {
            'outlier_count': coverage_fail_count
        }

        out_dict = {
            'pass': process_complete,
            'step': step,
            'metrics': metrics_dict,
            'output': outfiles_dict
        }

        return out_dict


    def run_duplication_check(self):
        # ZH doesn't prune anything here
        # is she manually checking for something?
        pass

    def run_titv_check(self, variant_calling_summary_metrics_path):
        # variant_calling_summary_metrics.tsv file?
        geno_path = self.geno_path
        out_path = self.out_path

        step = "titv_check"

        # create filenames
        titv_ratio_fails = f'{out_path}.titv_ratio_fails'

        vc_summary_metrics = pd.read_csv(variant_calling_summary_metrics_path, sep='\t')
        titv_ratio_fail_df = vc_summary_metrics[vc_summary_metrics.DBSNP_TITV < 2]
        titv_ratio_fail_ids = titv_ratio_fail_df['sample_id']
        titv_ratio_fail_count = titv_ratio_fail_ids.shape[0]
        titv_ratio_fail_ids.to_csv(titv_ratio_fails, sep='\t', header=False, index=False)

        plink_cmd = f'{plink2_exec} --pfile {geno_path} --remove {titv_ratio_fails} --make-pgen psam-cols=fid,parents,sex,pheno1,phenos --out {out_path}'
        shell_do(plink_cmd)

        listOfFiles = [f'{out_path}.log']
        concat_logs(step, out_path, listOfFiles)

        process_complete = True

        outfiles_dict = {
            'pruned_samples': f'{titv_ratio_fails}',
            'plink_out': f'{out_path}',
        }

        metrics_dict = {
            'outlier_count': titv_ratio_fail_count
        }

        out_dict = {
            'pass': process_complete,
            'step': step,
            'metrics': metrics_dict,
            'output': outfiles_dict
        }

        return out_dict

    def run_het_check(self):
        # can we use the method in qc.py file?
        # or do we want to do per chr method here?
        # TAKES IN GP2_1KG_BFILES_PATH and GP2_BFILES_PATH? whats the difference between
        geno_path = self.geno_path
        out_path = self.out_path

        step = "wgs_het_check"

        # create filenames
        het_var_list = f'{out_path}.het_variants_list'
        het_vars = f'{out_path}.het_variants'
        het_test = f'{out_path}.het_test'
        het_outliers = f'{out_path}.het_outliers'

        # create heterozygous variant list
        # ZH specifies number of threads here? --threads 2
        plink_cmd1 = f'{plink_exec} --bfile {geno_path} --geno 0.01 --maf 0.05 --hwe 1E-4 --indep-pairwise 50 5 0.5 --out {het_var_list}'

        # extract heterozygous variants from input geno files
        # ZH specifies memory and threads here? --memory 20000 --threads 4
        plink_cmd2 = f'{plink_exec} --bfile {geno_path} --extract {het_var_list}.prune.in --make-bed --out {het_vars}'

        # generate heterozygosity test results
        # ZH specifies threads here? --threads 2
        plink_cmd3 = f'{plink_exec} --bfile {geno_path} --het --out {het_test}'

        cmds = [plink_cmd1, plink_cmd2, plink_cmd3]
        for cmd in cmds:
            shell_do(cmd)

        listOfFiles = [f'{het_var_list}.log', f'{het_vars}.log', f'{het_test}.log']
        concat_logs(step, out_path, listOfFiles)

        # ZH does some manual checking/file creation here with awk
        # outputs het outliers outside of (-0.15, 0.15) to 1kg_HETEROZYGOSITY_OUTLIERS1.txt
        # outputs het outliers outside of (-0.25, 0.25) to 1kg_HETEROZYGOSITY_OUTLIERS2.txt

        het_df = pd.read_csv(f'{het_test}.het', sep='\s+')
        check1 = het_df['F'].mean() + (3*het_df['F'].std())
        check2 = het_df['F'].mean() - (3*het_df['F'].std())
        het_outlier_df = het_df.loc[(het_df['F'] >= check1) | (het_df['F'] <= check2)]
        het_outlier_ids = het_outlier_df['IID']
        het_outlier_count = het_outlier_ids.shape[0]
        het_outlier_ids.to_csv(het_outliers, sep='\t', header=False, index=False)

        # do we want to remove het outliers?
        plink_cmd = f'{plink2_exec} --pfile {geno_path} --remove {het_outliers} --make-pgen psam-cols=fid,parents,sex,pheno1,phenos --out {out_path}'
        shell_do(plink_cmd)

        listOfFiles = [f'{out_path}.log']
        concat_logs(step, out_path, listOfFiles)

        process_complete = True

        outfiles_dict = {
            'pruned_samples': f'{het_outliers}',
            'plink_out': f'{out_path}',
        }

        metrics_dict = {
            'outlier_count': het_outlier_count
        }

        out_dict = {
            'pass': process_complete,
            'step': step,
            'metrics': metrics_dict,
            'output': outfiles_dict
        }

        return out_dict


    def run_check_callrate(self):
        # maybe can use the callrate prune method in qc.py file? ZH seems to not split by chr?

        # split by chr here (or do we expect the input geno to already be split by chromosome?)
        # run call rate prune
        # average callrate across each chr for each sample
        # then remerge the chrs for each sample?
        geno_path = self.geno_path
        out_path = self.out_path

        step = "callrate_check"

        # create filenames
        call_rates = f'{geno_path}_call_rates'
        call_rate_fail = f'{geno_path}.callrate_fail'

        # split by chr (if n autosomes, n+1 is the X chromosome, n+2 is Y, n+3 is XY, and n+4 is MT per plink2)
        per_chr_callrate = dict()
        for chr in range(1,27):
            plink_cmd1 = f'{plink2_exec} --pfile {geno_path} --chr {chr} --make-pgen psam-cols=fid,parents,sex,pheno1,phenos --out {geno_path}_{chr}'

            # ZH specifies threads here? --threads 4
            plink_cmd2 = f'{plink2_exec} --pfile {geno_path}_{chr} --missing --out {call_rates}_{chr}'

            cmds = [plink_cmd1, plink_cmd2]
            for cmd in cmds:
                shell_do(cmd)

            listOfFiles = [f'{geno_path}_{chr}.log', f'{call_rates}_{chr}.log']
            concat_logs(step, out_path, listOfFiles)

            # store per sample callrates by chromosome
            # why use .smiss file from --missing instead of --mind (like in qc.py)?
            callrate_df = pd.read_csv(f'{call_rates}_{chr}.smiss', sep='\s+')
            callrate_df[f'callrate_{chr}'] = 1 - callrate_df.F_MISS
            per_chr_callrate[chr] = callrate_df[['IID', f'callrate_{chr}']].set_index('IID')

        # merge all callrates/chr together
        # get average across all chromosomes per sample
        callrate_dfs = list(per_chr_callrate.values())
        one_callrate_df = callrate_dfs.pop(0)
        all_callrate_df = one_callrate_df.join(callrate_dfs, how='outer')
        all_callrate_df['callrate'] = all_callrate_df.mean()

        # count sample as fail if average callrate across all chr is less than 0.95
        call_rate_fail_ids = all_callrate_df[all_callrate_df.callrate < 0.95]['IID']
        call_rate_fail_count = call_rate_fail_ids.shape[0]
        call_rate_fail_ids.to_csv(call_rate_fail, sep='\t', header=False, index=False)

        process_complete = True

        outfiles_dict = {
            'pruned_samples': f'{call_rate_fail}',
            'plink_out': f'{out_path}',
        }

        metrics_dict = {
            'outlier_count': call_rate_fail_count
        }

        out_dict = {
            'pass': process_complete,
            'step': step,
            'metrics': metrics_dict,
            'output': outfiles_dict
        }

        return out_dict


    def run_sex_check(self):
        # can we use the method in qc.py file?
        # but we want to extract the chr23 to check just that for sex?
        # ZH does both automated and manual sex check?
        # -> automated sex check just calls maf=0.05 and check-sex [0.25, 0.75]
        # -> manual sex check singles out chr23:bp=[2781479, 155701383] with maf=0.05, geno=0.05, hwe=1E-5 and then check-sex [0.25, 0.75]
        geno_path = self.geno_path
        out_path = self.out_path

        step = "wgs_sex_check"
        pass

    def run_genetic_dups_check(self, preBqsr_path):
        # utilizes KING software
        # runs --related (like in qc.py)
        # then drops the sample with lower depth (using data from preBqsr file)
        pass



