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
from collections import defaultdict
from genotools.utils import shell_do, count_file_lines, concat_logs, bfiles_to_pfiles
from genotools.dependencies import check_plink, check_plink2, check_king

plink_exec = check_plink()
plink2_exec = check_plink2()
king_exec = check_king()

class WholeGenomeSeqQC:

    def __init__(self, geno_path=None, out_path=None, keep_all=True):
        self.geno_path = geno_path
        self.out_path = out_path
        self.keep_all = keep_all


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
        preBqsr = pd.read_csv(preBqsr_path, sep='\s+')
        geno_samples = pd.read_csv(f'{geno_path}.psam', sep='\s+')
        preBqsr['CHECK'] = preBqsr.FREELK1 - preBqsr.FREELK0
        preBqsr = pd.merge(preBqsr, geno_samples, how='inner', left_on='sample_id', right_on='IID')
        if preBqsr.shape[0] != geno_samples.shape[0]:
            print('Warning: metrics file does not contain info for all samples in geno files!')
        preBqsr = preBqsr[['sample_id', 'AVG_DP', 'FREEMIX', 'FREELK1', 'FREELK0', 'CHECK']]

        # create filenames
        freemix_flagged = f'{out_path}.freemix_flagged'
        freemix_fails = f'{out_path}.freemix_fails'

        preBqsr = pd.read_csv(f'{preBqsr_path}', sep='\s+')
        freemix_flagged_df = preBqsr[preBqsr.FREEMIX >= 0.03]
        freemix_flagged_ids = freemix_flagged_df['sample_id']
        freemix_flagged_count = freemix_flagged_ids.shape[0]
        freemix_flagged_df.to_csv(freemix_flagged, sep='\t', header=False, index=False)

        # PER ZH: 'Problem samples MAY be failures and should be removed if CHECK value is large (need to define limits)'
        # TODO: define CHECK limit (arbitrarily set as 1e6 for now)
        freemix_fail_df = freemix_flagged_df[freemix_flagged_df.CHECK >= 1e6]
        freemix_fail_ids = freemix_fail_df['sample_id']
        freemix_fail_count = freemix_fail_ids.shape[0]
        freemix_fail_ids.to_csv(freemix_fails, sep='\t', header=False, index=False)
        if not self.keep_all:
            plink_cmd = f'{plink2_exec} --pfile {geno_path} --remove {freemix_fails} --make-pgen psam-cols=fid,parents,sex,pheno1,phenos --out {out_path}'
            shell_do(plink_cmd)

            listOfFiles = [f'{out_path}.log']
            concat_logs(step, out_path, listOfFiles)

        process_complete = True

        outfiles_dict = {
            'flagged_count': f'{freemix_flagged}',
            'pruned_count': f'{freemix_fails}',
            'plink_out': f'{out_path}',
        }

        metrics_dict = {
            'flagged_count': freemix_flagged_count,
            'removed_count': freemix_fail_count
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
        geno_samples = pd.read_csv(f'{geno_path}.psam', sep='\s+')
        wgs_metrics = pd.merge(wgs_metrics, geno_samples, how='inner', left_on='sample_id', right_on='IID')
        if wgs_metrics.shape[0] != geno_samples.shape[0]:
            print('Warning: metrics file does not contain info for all samples in geno files!')
        coverage_fail_df = wgs_metrics[wgs_metrics.MEAN_COVERAGE < min_mean_coverage]
        coverage_fail_ids = coverage_fail_df['sample_id']
        coverage_fail_count = coverage_fail_ids.shape[0]
        coverage_fail_ids.to_csv(coverage_fails, sep='\t', header=False, index=False)

        out_dict_key = 'flagged_samples'

        # PER ZH: 'Failed samples should be noted'
        if not self.keep_all:
            plink_cmd = f'{plink2_exec} --pfile {geno_path} --remove {coverage_fails} --make-pgen psam-cols=fid,parents,sex,pheno1,phenos --out {out_path}'
            shell_do(plink_cmd)

            listOfFiles = [f'{out_path}.log']
            concat_logs(step, out_path, listOfFiles)

            out_dict_key = 'pruned_samples'

        process_complete = True

        outfiles_dict = {
            out_dict_key: f'{coverage_fails}',
            'plink_out': f'{out_path}'
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


    def run_titv_check(self, variant_calling_summary_metrics_path):
        geno_path = self.geno_path
        out_path = self.out_path

        step = "titv_check"

        # create filenames
        titv_ratio_fails = f'{out_path}.titv_ratio_fails'

        vc_summary_metrics = pd.read_csv(variant_calling_summary_metrics_path, sep='\s+')
        geno_samples = pd.read_csv(f'{geno_path}.psam', sep='\s+')
        vc_summary_metrics = pd.merge(vc_summary_metrics, geno_samples, how='inner', left_on='sample_id', right_on='IID')
        if vc_summary_metrics.shape[0] != geno_samples.shape[0]:
            print('Warning: metrics file does not contain info for all samples in geno files!')
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
        geno_path = self.geno_path
        out_path = self.out_path

        step = "het_check"

        # create filenames
        out_filename = os.path.basename(out_path)
        het_dir = f'{out_path}/het_dir'
        if not os.path.exists(f'{het_dir}'):
            os.makedirs(f'{het_dir}')
        het_tmp = f"{het_dir}/{out_filename}_tmp"
        het_tmp2 = f"{het_dir}/{out_filename}_tmp2"
        het_tmp3 = f"{het_dir}/{out_filename}_tmp3"

        plink_cmd1 = f"{plink2_exec} --pfile {geno_path} --geno 0.01 --maf 0.05 --indep-pairwise 50 5 0.5 --out {het_tmp}"
        plink_cmd2 = f"{plink2_exec} --pfile {geno_path} --extract {het_tmp}.prune.in --make-pgen psam-cols=fid,parents,sex,pheno1,phenos --out {het_tmp2}"
        plink_cmd3 = f"{plink2_exec} --pfile {het_tmp2} --het --out {het_tmp3}"

        cmds1 = [plink_cmd1, plink_cmd2, plink_cmd3]

        for cmd in cmds1:
            shell_do(cmd)

        listOfFiles = [f'{het_tmp}.log', f'{het_tmp2}.log', f'{het_tmp3}.log']
        concat_logs(step, out_path, listOfFiles)

        hetpath = f'{het_tmp3}.het'
        if os.path.isfile(hetpath):
            os.remove(f'{het_tmp}.prune.in')
            os.remove(f'{het_tmp}.prune.out')
            os.remove(f'{het_tmp2}.pgen')
            os.remove(f'{het_tmp2}.pvar')
            os.remove(f'{het_tmp2}.psam')
        else:
            print(f'Heterozygosity pruning failed!')
            print(f'Check {het_tmp}.log, {het_tmp2}.log, or {het_tmp3}.log for more information')


    def run_check_callrate(self):
        geno_path = self.geno_path
        out_path = self.out_path

        step = "callrate_check"

        # create filenames
        out_filename = os.path.basename(out_path)
        smiss_dir = f'{out_path}/smiss_dir'
        if not os.path.exists(f'{smiss_dir}'):
            os.makedirs(f'{smiss_dir}')
        callrates = f'{smiss_dir}/{out_filename}_callrates'
        callrate_fail = f'{out_path}.callrate_fails'

        plink_cmd = f'{plink2_exec} --pfile {geno_path} --missing --out {callrates}'
        shell_do(plink_cmd)

        listOfFiles = [f'{call_rates}.log']
        concat_logs(step, out_path, listOfFiles)

        # callrate_df = pd.read_csv(f'{callrates}.smiss', sep='\s+')
        # callrate_df['CALLRATE'] = 1 - callrate_df.F_MISS

        # # count sample as fail if average callrate across all chr is less than 0.95
        # callrate_fail_ids = callrate_df[callrate_df.callrate < 0.95]['IID']
        # callrate_fail_count = callrate_fail_ids.shape[0]
        # callrate_fail_ids.to_csv(callrate_fail, sep='\t', header=False, index=False)

        # if not self.keep_all:
        #     plink_cmd = f'{plink2_exec} --pfile {geno_path} --remove {callrate_fail} --make-pgen psam-cols=fid,parents,sex,pheno1,phenos --out {out_path}'
        #     shell_do(plink_cmd)

        #     listOfFiles = [f'{out_path}.log']
        #     concat_logs(step, out_path, listOfFiles)

        # process_complete = True

        # outfiles_dict = {
        #     'pruned_samples': f'{callrate_fail}',
        #     'plink_out': f'{out_path}',
        # }

        # metrics_dict = {
        #     'outlier_count': callrate_fail_count
        # }

        # out_dict = {
        #     'pass': process_complete,
        #     'step': step,
        #     'metrics': metrics_dict,
        #     'output': outfiles_dict
        # }

        # return out_dict


    def run_sex_check(self, check_sex=[0.25,0.75]):
        geno_path = self.geno_path
        out_path = self.out_path

        step = "wgs_sex_check"

        # create filenames
        sex_tmp = f'{out_path}_tmp'
        sex_fails = f'{out_path}.sex_fails'

        # convert to bfiles
        bfiles_to_pfiles(pfile_path=geno_path)

        plink_cmd = f'{plink_exec} --bfile {geno_path} --chr 23 --from-bp2781479 --to-bp 155701383 --maf 0.05 --geno 0.05 --hwe 1E-5 --check-sex  0.25 0.75 --out {sex_tmp}'
        shell_do(plink_cmd)

        listOfFiles = [f'{sex_tmp}.log']
        concat_logs(step, out_path, listOfFiles)

        if os.path.isfile(f'{sex_tmp}.sexcheck'):
            sex = pd.read_csv(f'{sex_tmp}.sexcheck', sep='\s+')
            outliers = sex[sex.STATUS=='PROBLEM']

            sex_fail_ids = outliers.loc[:,['FID','IID']]
            sex_fail_count = sex_fail_ids.shape[0]
            sex_fail_ids.to_csv(sex_fails, sep='\t', header=True, index=False)

            process_complete = True

            for file in [f'{sex_tmp}.hh',f'{sex_tmp}.sexcheck']:
                if os.path.isfile(file):
                    os.remove(file)

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
            print(f'Check {sex_tmp}.log for more information')

            process_complete = False

            outfiles_dict = {
                'pruned_samples': 'Sex Prune Failed!',
                'plink_out': out_path
            }

            metrics_dict = {
                'outlier_count': 0
            }

        os.remove(f'{geno_path}.bed')
        os.remove(f'{geno_path}.bim')
        os.remove(f'{geno_path}.fam')

        out_dict = {
            'pass': process_complete,
            'step': step,
            'metrics': metrics_dict,
            'output': outfiles_dict
        }

        return out_dict

    # TODO: MAYBE HAVE CHECKING FOR DUPLICATES/RELATED WITHIN the wgsWrapper CLASS INSTEAD?
    # def run_genetic_dups_check(self, preBqsr_path):
    #     # utilizes KING software
    #     # runs --related (like in qc.py)
    #     # then drops the sample with lower depth (using data from preBqsr file)
    #     geno_path = self.geno_path
    #     out_path = self.out_path

    #     step = "wgs_related_check"

    #     # create filenames
    #     pre_king = f'{geno_path}_pre_king'
    #     kin_est = f'{geno_path}_kin_est'
    #     dup_samples = f'{geno_path}.dup_to_exclude'

    #     plink_cmd = f'{plink2_exec} --pfile {geno_path} --geno 0.01 --maf 0.05 --make-bed --out {pre_king}'

    #     king_cmd = f'{king_exec} --prefix {kin_est} -b {pre_king}.bed --related --degree 4 --duplicate'

    #     cmds = [plink_cmd, king_cmd]
    #     for cmd in cmds:
    #         shell_do(cmd)

    #     listOfFiles = [f'{pre_king}.log']
    #     # in qc.py sam manually wrote the king log file; should sam do that here as well?
    #     concat_logs(step, out_path, listOfFiles)

    #     # ZH just prints out the within family and across family relationships found
    #     # -> doesn't look like she removes relatives
    #     # -> only removes duplicates
    #     #    - duplicates within same cohort removed
    #     #    - duplicates in difference cohorts: core cohort/participants with associated RNAseq prioritized; otherwise higher coverage kept

    #     # .con file created from KING relatedness run?
    #     # can't find documentation regarding it online
    #     if os.path.isfile(f'{kin_est}.con'):
    #         process_complete = True
    #         dup_table = pd.read_csv(f'{kin_est}.con', sep='\s+')
    #         if dup_table.shape[0] != 0:
    #             preBqsr_df = pd.read_csv(preBqsr_path, sep='\s+')
    #             temp_depth_df = pd.merge(dup_table, preBqsr_df.rename(columns={'AVG_DP':'AVD_DP1'}), left_on='IID1', right_on='participant_id', how='left')
    #             depth_df = pd.merge(temp_depth_df, preBqsr_df.rename(columns={'AVD_DP':'AVD_DP2'}), left_on='IID2', right_on='participant_id', how='left')
    #             dup = pd.DataFrame(depth_df, columns=['IID1', 'AVG_DP1', 'IID2', 'AVG_DP2'])
    #             # get ids of those we want to exclude
    #             dup['lower_depth'] = dup[['AVG_DP1', 'AVG_DP2']].idxmin(axis=1).str.replace('AVG_DP', '')
    #             dup.loc[dup['lower_depth']=='1', 'dup_samples_to_exclude'] = dup.loc[dup['lower_depth']=='1', 'IID1']
    #             dup.loc[dup['lower_depth']=='2', 'dup_samples_to_exclude'] = dup.loc[dup['lower_depth']=='2', 'IID2']
    #             dup_sample_ids = dup['dup_samples_to_exclude']
    #             dup_sample_ids.to_csv(dup_samples, header=False, index=False)
    #             dup_sample_count = dup_sample_ids.shape[0]

    #             plink_cmd = f'{plink2_exec} --pfile {geno_path} --remove {dup_samples} --make-pgen psam-cols=fid,parents,sex,pheno1,phenos --out {out_path}'
    #             shell_do(plink_cmd)

    #             listOfFiles = [f'{out_path}.log']
    #             concat_logs(step, out_path, listOfFiles)

    #             outfiles_dict = {
    #                 'pruned_samples': f'{dup_samples}',
    #                 'plink_out': f'{out_path}',
    #             }

    #             metrics_dict = {
    #                 'outlier_count': dup_sample_count
    #             }

    #             out_dict = {
    #                 'pass': process_complete,
    #                 'step': step,
    #                 'metrics': metrics_dict,
    #                 'output': outfiles_dict
    #             }

    #         else:

    #             outfiles_dict = {
    #                 'pruned_samples': None,
    #                 'plink_out': f'{out_path}',
    #             }

    #             metrics_dict = {
    #                 'outlier_count': 0
    #             }

    #             out_dict = {
    #                 'pass': process_complete,
    #                 'step': step,
    #                 'metrics': metrics_dict,
    #                 'output': outfiles_dict
    #             }
    #     else:
    #         process_complete = False
    #         outfiles_dict = {
    #             'pruned_samples': None,
    #             'plink_out': f'{out_path}',
    #         }

    #         metrics_dict = {
    #             'outlier_count': 0
    #         }

    #         out_dict = {
    #             'pass': process_complete,
    #             'step': step,
    #             'metrics': metrics_dict,
    #             'output': outfiles_dict
    #         }

    #     return out_dict




class WholeGenomeSeqWrapper:

    def __init__(self, shards_dir=None, out_path=None, keep_all=True, swarm=False, shard_key=None, preBqsr_path=None, wgs_metrics_path=None, variant_calling_summary_metrics_path=None):
        self.shards_dir = shards_dir
        self.out_path = out_path
        self.keep_all = keep_all

        # wgs qc metrics files
        self.shard_key = shard_key
        self.preBqsr = preBqsr_path
        self.wgs_metrics = wgs_metrics_path
        self.var_calling_summary_metrics_path = variant_calling_summary_metrics_path

        # create swarm file dir if using swarm
        if swarm:
            os.makedirs(f'{self.out_path}/swarm_scripts')


    '''
    create our own pipeline method within pipeline.py?

    execute_wgs_pipeline(steps, steps_dict, shards_dir, out_path, wgs_qc, ancestry, assoc, args, tmp_dir)
        if step in wgs_steps:
            wgs_wrapper = WholeGenomeSeqWrapper(shards_dir, out, keep, wgs_qc_metrics_files)
            for each shard in shards_dir:
                wgs_qc = WholeGenomeSeqQC(geno, out, keep) # creates one instance of WholeGenomeSeqQC per shard
                if freemix:
                    wgs_qc.run_freemix_check()
                if coverage:
                    wgs_qc.run_coverage_check()
                if titv_ratio:
                    wgs_qc.run_titv_ratio_check()
                if het:
                    wgs_qc.run_het()
                    (generate output to out_path/het_dir)
                if callrate:
                    wgs_qc.run_callrate()
                    (generate output to out_path/callrate_dir)

            if freemix:
                wgs_wrapper.freemix_out_metrics(): combine outputs for each shard into one output? to generate out_dict
            if coverage:
                wgs_wrapper.coverage_out_metrics(): combine outputs for each shard into one output? to generate out_dict
            if titv_ratio:
                wgs_wrapper.titv_ratio_metrics(): combine outputs for each shard into one output? to generate out_dict
            if sex:
                x_geno = wgs_wrapper.find_x_shard()
                wgs_qc = WholeGenomeSeqQC(x_geno, out, keep)
                wgs_qc.run_sex_check()
            if het:
                wgs_wrapper.merge_sample_het()
            if callrate:
                wgs_wrapper.merge_sample_callrate()

            if related/duplicate:
                wgs_wrapper.extract_ref_snps()
                wgs_wrapper.run_wgs_related_duplicated()
                wgs_wrapper.drop_rel_dups_less_coverage()

    Within each method:
        if swarm:
            create swarm file call through shell_do()
    '''

    def freemix_out_metrics(self):
        '''
        loop through shards freemix out
        compile all information into one output file
        cleanup freemix out per shard
        '''pass

    def coverage_out_metrics(self):
        '''
        loop through shards coverage out
        compile all information into one output file
        cleanup coverage out per shard
        '''
        pass

    def titv_ratio_metrics(self):
        '''
        loop through shards titv ratio out
        compile all information into one output file
        cleanup titv ratio out per shard
        '''
        pass

    def find_x_shard(self):
        '''
        use shard key to find shard with interval on x chromosome for sex check
        '''
        shard_key = self.shard_key


    def merge_sample_het(self, shards_dir, het_filter=[-0.25,0.25]):
        '''
        -> point to dir where plink --het was outputted for each shard
        '''
        out_path = self.out_path

        step = "heterozygosity check"

        # create filenames
        het_fails = f'{out_path}.het_fails'

        heterosygosities = defaultdict(list)
        for filename in os.listdir(shards_dir):
            if filename.endswith('.het'):
                het = pd.read_csv(f'{filename}.het', setp='\s+')
                het['F'] = het['F'].map(lambda x:[x])
                shard_hets = pd.Series(het.F.values, index=het.IID).to_dict()
                for k,v in shard_hets.items():
                    heterosygosities[k].extend(v)

        for k,v in heterozygosities.items():
            heterosygosities[k] = sum(v) / len(v)

        heterozygosities = pd.DataFrame.from_dict(heterozygosities, orient='index').reset_index()
        heterozygosities = heterozygosities.rename({columns={'index':'IID', 0:'F'}})

        outliers = heterosygosities[((heterosygosities.F <= het_filter[0]) | (heterosygosities.F >= het_filter[1]))]
        het_fail_ids = outliers['IID']
        het_fail_count = het_fail_ids.shape[0]
        het_fail_ids.to_csv(f'{het_fails}', sep='\t', header=False, index=False)

        process_complete = True

        outfiles_dict = {
            'pruned_samples': f'{het_fails}',
            'plink_out': f'{out_path}' # TODO: if removing include plink_out
        }

        metrics_dict = {
            'outlier_count': het_fail_count
        }

        out_dict = {
            'pass': process_complete,
            'step': step,
            'metrics': metrics_dict,
            'output': outfiles_dict
        }

        return out_dict


    def merge_sample_callrate(self, shards_dir):
        '''
        -> point to directory where plink --missing outputted for each shard
        '''
        out_path = self.out_path

        step = "callrate_check"

        # create filenames
        callrate_fails = f'{out_path}.callrate_fails'

        callrates = defaultdict(list)
        for filename in os.listdir(shards_dir):
            if filename.endswith('.smiss'):
                smiss = pd.read_csv(f'{filename}.smiss', sep='\s+')
                smiss['CALLRATE'] = smiss['F_MISS'].map(lambda x: [1-x])
                shard_callrates = pd.Series(smiss.CALLRATE.values, index=smiss.IID).to_dict()
                for k,v in shard_callrates.items():
                    callrates[k].extend(v)

        for k,v in callrates.items():
            callrates[k] = sum(v) / len(v)

        callrates = pd.DataFrame.from_dict(callrates, orient='index').reset_index()
        callrates = callrates.rename({columns={'index':'IID', 0:'CALLRATE'}})

        outliers = callrates[callrates.CALLRATE < 0.95]
        callrate_fail_ids = outliers['IID']
        callrate_fail_count = callrate_fail_ids.shape[0]
        callrate_fail_ids.to_csv(f'{callrate_fails}', sep='\t', header=False, index=False)

        process_complete = True

        outfiles_dict = {
            'pruned_samples': f'{callrate_fails}',
            'plink_out': f'{out_path}' # TODO: if removing include plink_out
        }

        metrics_dict = {
            'outlier_count': callrate_fail_count
        }

        out_dict = {
            'pass': process_complete,
            'step': step,
            'metrics': metrics_dict,
            'output': outfiles_dict
        }

        return out_dict

    def extract_ref_snps(self, ref_panel_path):
        '''
        get chr:pos of all snps in ref panel
        use snp map to get chr:pos range of each shard
        for each snp in ref panel
            get shard that contains that snp, add to set of shards for ref panel extraction

        for each shard containing a ref panel snp
            extract ref panel snp

        merge all extracted files together (output path to merged)
        '''
        shards_dir = self.shards_dir
        shards_key = self.shards_key

        # format shards key
        shards_key['start'] = shards_key['interval'].str.split(':', expand=True)[1].str.split('-', expand=True)[0].astype(int)
        shards_key['end'] = shards_key['interval'].str.split(':', expand=True)[1].str.split('-', expand=True)[1].astype(int)

        # load in ref panel snps
        # TODO: assume ref panel snps in repo?
        ref_panel = pd.read_csv(f'{ref_panel_path}', sep='\s+', header=None, names=['snpID'])
        ref_panel['CHR'] = ref_panel['snpID'].str.split(':', expand=True)[0]
        ref_panel['POS'] = ref_panel['snpID'].str.split(':', expand=True)[1].astype(int)

        # find shard that contains each ref snp
        ref_shard_merge = ref_panel.merge(shards_key, left_on='CHR', right_on='chr', how='outer')
        ref_shard_interval = ref_shard_merge[(ref_shard_merge['POS']>=ref_shard_merge['start']) & (ref_shard_merge['POS']<=ref_shard_merge['end'])]
        ref_panel_with_shard = ref_shard_interval[['snpID', 'CHR', 'POS', 'shard']]

        # make dir for extracted ref panel
        if not os.path.exists(f'{shards_dir}/ref_overlap'):
            os.makedirs(f'{shards_dir}/ref_overlap')

        ref_overlap_dir = f'{shards_dir}/ref_overlap'

        # extract ref snps from each corresponding shard
        shards_to_extract_from = set(zip(ref_panel_with_shard.CHR, ref_panel_with_shard.shard))
        for chrom, shard in shards_to_extract_from:
            shard_path = f'{shards_dir}/{chrom}_{shard}'
            plink_extract = f'{plink2_exec} --pfile {shard_path} --extract {ref_panel_path} --make-pgen --out {ref_overlap_dir}/{chrom}_{shard}_ref_overlap'
            shell_do(plink_extract)

        # merge all the extracted ref_panel snps
        with open(f'{ref_overlap_dir}/to_merge.txt', 'w') as f:
            for chrom, shard in shards_to_extract_from:
                path = f'{ref_overlap_dir}/{chrom}_{shard}_ref_overlap'
                f.write(f'{path}\n')
        f.close()

        plink_merge = f'{plink2_exec} --pmerge-list {ref_overlap_dir}/to_merge.txt --make-pgen --out {ref_overlap_dir}/merged_ref_overlap'
        shell_do(plink_merge)

        # return path to merged ref panel overalp pfiles
        return f'{ref_overlap_dir}/merged_ref_overlap'


    def run_wgs_related_duplicates(self):
        '''
        call KING related on path to merged ref panel snps
        '''
        pass

    def drop_rel_dups_less_coverage(self):
        '''
        with output from KING, check coverage for each pair
        remove the sample with less coverage
        '''
        pass


    # ZH's final wgs qc output that combines the qc metrics into one file to get pass/fail samples
    # def wgs_qc_output(self, wgs_metrics_path, sex_check_path, preBqsr_path, dup_exclude_path):
    #     # compiling all the wgs qc metrics into one file
    #     # ZH does this in final step of smk with wgs_qc.py
    #     out_path = self.out_path

    #     step = 'wgs_qc_metrics_file'

    #     wgs_qc_out = f'{out_path}.wgs_qc_out'

    #     cov = pd.read_csv(f'{wgs_metrics_path}', sep='\s+')
    #     cov = cov[['sample_id', 'MEAN_COVERAGE']]
    #     sex = pd.read_csv(f'{sex_check_path}', sep='\s+')
    #     con = pd.read_csv(f'{preBqsr_path}', sep='\s+')
    #     wgs = pd.merge(pd.merge(cov, con, on='sample_id', how='outer'), sex, left_on='sample_id', right_on='IID', how='outer')
    #     wgs['sample_id'] = wgs['sample_id'].fillna(wgs['IID'])
    #     wgs.loc[~wgs['IID'].isnull(), 'called'] = '1'

    #     out = wgs.drop(columns=['FID', 'IID'])
    #     out['Contamination_QC'] = np.where((out['FREEMIX'] >= 0.03), 'FAIL', 'PASS')
    #     out['Coverage_QC'] = np.where((out['MEAN_COVERAGE'] < 25), 'FAIL', 'PASS')
    #     out['Sex_QC'] = np.where((out['STATUS'] == 'PROBLEM'), 'FAIL', 'PASS')
    #     out['Contamination_QC'].mask(out['STATUS'].isna(), np.nan, inplace=True)
    #     out['Coverage_QC'].mask(out['Status'].isna(), np.nan, inplace=True)
    #     out['Sex_QC'].mask(out['STATUS'].isna(), np.nan, inplace=True)

    #     cols = ['Contamination_QC', 'Coverage_QC', 'Sex_QC']
    #     out.loc[(out[cols]=='PASS').all(axis=1), 'Overall_SampleQC'] = 'PASS'
    #     out.loc[(out['Overall_SampleQC']!='PASS') & (~out['STATUS'].isnull()), 'Overall_SampleQC'] = 'FAIL'

    #     cols = ['Contamination_QC', 'Sex_QC']
    #     out.loc[(out[cols]=='PASS').all(axis=1), 'Included in analysis'] = '1'
    #     out.loc[(out['Included in analysis']!='1') & (~out['STATUS'].isnull()), 'Included in analysis'] = '0'
    #     out = out.drop('STATUS', axis=1)
    #     out.rename(columns = {'F': 'F(chrX_inbreeding_coefficient)',
    #                           'PEDSEX':'Reported_SEX',
    #                           'SNPSEX':'Inferred_Sex'
    #                           }, inplace=True)

    #     dup_samples_to_exclude = pd.read_csv(f'{dup_exclude_path}', sep='\s+', header=None, names=['ID'])
    #     dup_samples_to_exclude = list(dup_samples_to_exclude['IID'])

    #     out['cohort_included_analysis'] = out['Included in analysis']
    #     out.loc[out['sample_id'].isin(dup_samples_to_exclude)]
    #     out.to_csv(f'{wgs_qc_out}', sep='\t', index=False)

    #     outfiles_dict = {
    #         'wgs_qc_out': wgs_qc_out
    #     }

    #     out_dict = {
    #         'process_complete': True,
    #         'step': step,
    #         'output': outfiles_dict
    #     }

    #     return out_dict