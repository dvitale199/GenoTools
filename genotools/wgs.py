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
        preBqsr['CHECK'] = preBqsr.FREELK1 - preBqsr.FREELK0
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


    def run_duplication_check(self, dup_metrics_path):
        # ZH doesn't prune anything here
        # is she manually checking for something?
        # duplicate_metrics.tsv?

        # duplication_df = pd.read_csv(f'{dup_metrics_path}', sep='\s+')
        # dup_flagged = duplication_df.loc[duplication_df.sample_id.isin()]
        pass

    def run_titv_check(self, variant_calling_summary_metrics_path):
        # variant_calling_summary_metrics.tsv file?
        geno_path = self.geno_path
        out_path = self.out_path

        step = "titv_check"

        # create filenames
        titv_ratio_fails = f'{out_path}.titv_ratio_fails'

        vc_summary_metrics = pd.read_csv(variant_calling_summary_metrics_path, sep='\s+')
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

    '''
    ZH plans on using the sampleQC het method instead of this one
    => call sampleQC's het within wgs pipeline
    '''
    # def run_het_check(self):
    #     # can we use the method in qc.py file?
    #     # or do we want to do per chr method here?
    #     # TAKES IN GP2_1KG_BFILES_PATH and GP2_BFILES_PATH? whats the difference between
    #     # ZH has three different criteria for checking het
    #     # 1) F cutoff > 0.15 and < -0.15
    #     # 2) F cutoff > 0.25 and < -0.25
    #     # 3) mean(F) + 3SD(F)
    #     geno_path = self.geno_path
    #     out_path = self.out_path

    #     step = "wgs_het_check"

    #     # create filenames
    #     # are some of these just temp files? should we use temp file names and then remove?
    #     het_var_list = f'{out_path}_het_variant_list'
    #     het_vars = f'{out_path}_het_variants'
    #     het_test = f'{out_path}_het_test'
    #     het_outliers1 = f'{out_path}.het_outliers1'
    #     het_outliers2 = f'{out_path}.het_outliers2'
    #     het_out = f'{out_path}.het_out'

    #     # create heterozygous variant list
    #     # ZH specifies number of threads here? --threads 2
    #     plink_cmd1 = f'{plink_exec} --bfile {geno_path} --geno 0.01 --maf 0.05 --hwe 1E-4 --indep-pairwise 50 5 0.5 --out {het_var_list}'

    #     # extract heterozygous variants from input geno files
    #     # ZH specifies memory and threads here? --memory 20000 --threads 4
    #     plink_cmd2 = f'{plink_exec} --bfile {geno_path} --extract {het_var_list}.prune.in --make-bed --out {het_vars}'

    #     # generate heterozygosity test results
    #     # ZH specifies threads here? --threads 2
    #     plink_cmd3 = f'{plink_exec} --bfile {geno_path} --het --out {het_test}'

    #     cmds = [plink_cmd1, plink_cmd2, plink_cmd3]
    #     for cmd in cmds:
    #         shell_do(cmd)

    #     listOfFiles = [f'{het_var_list}.log', f'{het_vars}.log', f'{het_test}.log']
    #     concat_logs(step, out_path, listOfFiles)

    #     # ZH does some manual checking/file creation here with awk
    #     # outputs het outliers outside of (-0.15, 0.15) to 1kg_HETEROZYGOSITY_OUTLIERS1.txt
    #     # outputs het outliers outside of (-0.25, 0.25) to 1kg_HETEROZYGOSITY_OUTLIERS2.txt
    #     het_df = pd.read_csv(f'{het_test}.het', sep='\s+')
    #     het_outliers1_df = het_df[(het_df.F <= -0.15) & (het_df.F >= 0.15)]
    #     het_outliers1_ids = het_outliers1_df[['FID', 'IID']]
    #     het_outliers1_count = het_outliers1_df.shape[0]
    #     het_outliers1_ids.to_csv(f'{het_outliers1}', sep='\t', header=False, index=False)
    #     het_outliers2_df = het_df[(het_df.F <= -0.25) & (het_df.F >= 0.25)]
    #     het_outliers2_ids = het_outliers2_df[['FID', 'IID']]
    #     het_outliers2_count = het_outliers2_df.shape[0]
    #     het_outliers2_ids.to_csv(f'{het_outliers2}', sep='\t', header=False, index=False)

    #     check1 = het_df['F'].mean() + (3*het_df['F'].std())
    #     check2 = het_df['F'].mean() - (3*het_df['F'].std())
    #     het_outlier_df = het_df.loc[(het_df['F'] >= check1) | (het_df['F'] <= check2)]
    #     het_outlier_ids = het_outlier_df['IID']
    #     het_outlier_count = het_outlier_ids.shape[0]
    #     het_outlier_ids.to_csv(het_out, sep='\t', header=False, index=False)

    #     out_dict_key = 'flagged_samples'

    #     # PER ZH: 'heterozygosity outliers MAY be failures and should be examined for ravial bias'
    #     # do we want to remove het outliers?
    #     if not self.keep_all:
    #         plink_cmd = f'{plink2_exec} --pfile {geno_path} --remove {het_out} --make-pgen psam-cols=fid,parents,sex,pheno1,phenos --out {out_path}'
    #         shell_do(plink_cmd)

    #         listOfFiles = [f'{out_path}.log']
    #         concat_logs(step, out_path, listOfFiles)

    #         out_dict_key = 'pruned_samples'

    #     process_complete = True

    #     outfiles_dict = {
    #         out_dict_key: f'{het_out}',
    #         'outliers1' : f'{het_outliers1}',
    #         'outliers2' : f'{het_outliers2}',
    #         'plink_out': f'{out_path}',
    #     }

    #     metrics_dict = {
    #         'outlier_count': het_outlier_count,
    #         'outliers1_count': het_outliers1_count,
    #         'outliers2_count': het_outliers2_count
    #     }

    #     out_dict = {
    #         'pass': process_complete,
    #         'step': step,
    #         'metrics': metrics_dict,
    #         'output': outfiles_dict
    #     }

    #     return out_dict


    def run_check_callrate(self):
        # maybe can use the callrate prune method in qc.py file? ZH seems to not split by chr?

        # INPUT IS ALREADY SPLIT INTO CHUNKS DO NOT SPLIT
        # run call rate prune
        # need an output with average callrate per sample for each input (chunk)
        # then need a method that takes in these outputs for each sample, and sums all the averages
        # thus getting per sample callrates
        # TODO: update cr method to include an output with per sample per chunk info
        # TODO: create a meta method that takes in this output and creates per sample callrate
        geno_path = self.geno_path
        out_path = self.out_path

        step = "callrate_check"

        # create filenames
        call_rates = f'{geno_path}_call_rates'
        call_rate_fail = f'{geno_path}.callrate_fail'

        plink_cmd = f'{plink2_exec} --pfile {geno_path} --missing --out {call_rates}'
        shell_do(plink_cmd)

        listOfFiles = [f'{call_rates}.log']
        concat_logs(step, out_path, listOfFiles)

        callrate_df = pd.read_csv(f'{call_rates}.smiss', sep='\s+')
        callrate_df['callrate'] = 1 - callrate_df.F_MISS

        # count sample as fail if average callrate across all chr is less than 0.95
        call_rate_fail_ids = callrate_df[callrate_df.callrate < 0.95]['IID']
        call_rate_fail_count = call_rate_fail_ids.shape[0]
        call_rate_fail_ids.to_csv(call_rate_fail, sep='\t', header=False, index=False)

        if not self.keep_all:
            plink_cmd = f'{plink2_exec} --pfile {geno_path} --remove {call_rate_fail} --make-pgen psam-cols=fid,parents,sex,pheno1,phenos --out {out_path}'
            shell_do(plink_cmd)

            listOfFiles = [f'{out_path}.log']
            concat_logs(step, out_path, listOfFiles)

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

    def run_merge_sample_callrates(self, list_of_shards):
        pass

    def run_sex_check(self):
        # can we use the method in qc.py file?
        # ZH does both automated and manual sex check (like method in qc.py):
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
        geno_path = self.geno_path
        out_path = self.out_path

        step = "wgs_related_check"

        # create filenames
        pre_king = f'{geno_path}_pre_king'
        kin_est = f'{geno_path}_kin_est'
        dup_samples = f'{geno_path}.dup_to_exclude'

        plink_cmd = f'{plink2_exec} --pfile {geno_path} --geno 0.01 --maf 0.05 --make-bed --out {pre_king}'

        king_cmd = f'{king_exec} --prefix {kin_est} -b {pre_king}.bed --related --degree 4 --duplicate'

        cmds = [plink_cmd, king_cmd]
        for cmd in cmds:
            shell_do(cmd)

        listOfFiles = [f'{pre_king}.log']
        # in qc.py sam manually wrote the king log file; should sam do that here as well?
        concat_logs(step, out_path, listOfFiles)

        # ZH just prints out the within family and across family relationships found
        # -> doesn't look like she removes relatives
        # -> only removes duplicates
        #    - duplicates within same cohort removed
        #    - duplicates in difference cohorts: core cohort/participants with associated RNAseq prioritized; otherwise higher coverage kept

        # .con file created from KING relatedness run?
        # can't find documentation regarding it online
        if os.path.isfile(f'{kin_est}.con'):
            process_complete = True
            dup_table = pd.read_csv(f'{kin_est}.con', sep='\s+')
            if dup_table.shape[0] != 0:
                preBqsr_df = pd.read_csv(preBqsr_path, sep='\s+')
                temp_depth_df = pd.merge(dup_table, preBqsr_df.rename(columns={'AVG_DP':'AVD_DP1'}), left_on='IID1', right_on='participant_id', how='left')
                depth_df = pd.merge(temp_depth_df, preBqsr_df.rename(columns={'AVD_DP':'AVD_DP2'}), left_on='IID2', right_on='participant_id', how='left')
                dup = pd.DataFrame(depth_df, columns=['IID1', 'AVG_DP1', 'IID2', 'AVG_DP2'])
                # get ids of those we want to exclude
                dup['lower_depth'] = dup[['AVG_DP1', 'AVG_DP2']].idxmin(axis=1).str.replace('AVG_DP', '')
                dup.loc[dup['lower_depth']=='1', 'dup_samples_to_exclude'] = dup.loc[dup['lower_depth']=='1', 'IID1']
                dup.loc[dup['lower_depth']=='2', 'dup_samples_to_exclude'] = dup.loc[dup['lower_depth']=='2', 'IID2']
                dup_sample_ids = dup['dup_samples_to_exclude']
                dup_sample_ids.to_csv(dup_samples, header=False, index=False)
                dup_sample_count = dup_sample_ids.shape[0]

                plink_cmd = f'{plink2_exec} --pfile {geno_path} --remove {dup_samples} --make-pgen psam-cols=fid,parents,sex,pheno1,phenos --out {out_path}'
                shell_do(plink_cmd)

                listOfFiles = [f'{out_path}.log']
                concat_logs(step, out_path, listOfFiles)

                outfiles_dict = {
                    'pruned_samples': f'{dup_samples}',
                    'plink_out': f'{out_path}',
                }

                metrics_dict = {
                    'outlier_count': dup_sample_count
                }

                out_dict = {
                    'pass': process_complete,
                    'step': step,
                    'metrics': metrics_dict,
                    'output': outfiles_dict
                }

            else:

                outfiles_dict = {
                    'pruned_samples': None,
                    'plink_out': f'{out_path}',
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
        else:
            process_complete = False
            outfiles_dict = {
                'pruned_samples': None,
                'plink_out': f'{out_path}',
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


    def wgs_qc_output(self, wgs_metrics_path, sex_check_path, preBqsr_path, dup_exclude_path):
        # compiling all the wgs qc metrics into one file
        # ZH does this in final step of smk with wgs_qc.py
        out_path = self.out_path

        step = 'wgs_qc_metrics_file'

        wgs_qc_out = f'{out_path}.wgs_qc_out'

        cov = pd.read_csv(f'{wgs_metrics_path}', sep='\s+')
        cov = cov[['sample_id', 'MEAN_COVERAGE']]
        sex = pd.read_csv(f'{sex_check_path}', sep='\s+')
        con = pd.read_csv(f'{preBqsr_path}', sep='\s+')
        wgs = pd.merge(pd.merge(cov, con, on='sample_id', how='outer'), sex, left_on='sample_id', right_on='IID', how='outer')
        wgs['sample_id'] = wgs['sample_id'].fillna(wgs['IID'])
        wgs.loc[~wgs['IID'].isnull(), 'called'] = '1'

        out = wgs.drop(columns=['FID', 'IID'])
        out['Contamination_QC'] = np.where((out['FREEMIX'] >= 0.03), 'FAIL', 'PASS')
        out['Coverage_QC'] = np.where((out['MEAN_COVERAGE'] < 25), 'FAIL', 'PASS')
        out['Sex_QC'] = np.where((out['STATUS'] == 'PROBLEM'), 'FAIL', 'PASS')
        out['Contamination_QC'].mask(out['STATUS'].isna(), np.nan, inplace=True)
        out['Coverage_QC'].mask(out['Status'].isna(), np.nan, inplace=True)
        out['Sex_QC'].mask(out['STATUS'].isna(), np.nan, inplace=True)

        cols = ['Contamination_QC', 'Coverage_QC', 'Sex_QC']
        out.loc[(out[cols]=='PASS').all(axis=1), 'Overall_SampleQC'] = 'PASS'
        out.loc[(out['Overall_SampleQC']!='PASS') & (~out['STATUS'].isnull()), 'Overall_SampleQC'] = 'FAIL'

        cols = ['Contamination_QC', 'Sex_QC']
        out.loc[(out[cols]=='PASS').all(axis=1), 'Included in analysis'] = '1'
        out.loc[(out['Included in analysis']!='1') & (~out['STATUS'].isnull()), 'Included in analysis'] = '0'
        out = out.drop('STATUS', axis=1)
        out.rename(columns = {'F': 'F(chrX_inbreeding_coefficient)',
                              'PEDSEX':'Reported_SEX',
                              'SNPSEX':'Inferred_Sex'
                              }, inplace=True)

        dup_samples_to_exclude = pd.read_csv(f'{dup_exclude_path}', sep='\s+', header=None, names=['ID'])
        dup_samples_to_exclude = list(dup_samples_to_exclude['IID'])

        out['cohort_included_analysis'] = out['Included in analysis']
        out.loc[out['sample_id'].isin(dup_samples_to_exclude)]
        out.to_csv(f'{wgs_qc_out}', sep='\t', index=False)

        outfiles_dict = {
            'wgs_qc_out': wgs_qc_out
        }

        out_dict = {
            'process_complete': True,
            'step': step,
            'output': outfiles_dict
        }

        return out_dict