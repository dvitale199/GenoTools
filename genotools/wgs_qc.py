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
import glob
import time
from collections import defaultdict
from genotools.utils import shell_do, concat_logs, bfiles_to_pfiles, count_file_lines
from genotools.dependencies import check_plink, check_plink2, check_king
from genotools.qc import SampleQC

plink_exec = check_plink()
plink2_exec = check_plink2()
king_exec = check_king()

class WholeGenomeSeqQC:

    '''
    execute_wgs_pipeline(steps, steps_dict, shards_dir, out_path, wgs_qc, ancestry, assoc, args, tmp_dir):
        # all functionality within WGS_qc class
        if step in wgs_steps:
            wgs_qc = WholeGenomeSeqQC(shards_dir, out_path, slurm, keep, wgs_qc_metrics_files)
            # take the first shard (assuming all shards have all samples)
            # and run freemix, coverage, titv_ratio on one shard
            if freemix:
                wgs_qc.run_freemix_check(shard_geno_path, check_limit)
            if coverage:
                wgs_qc.run_coverage_check(shard_geno_path, min_mean_coverage)
            if titv_ratio:
                wgs_qc.run_titv_ratio_check(shard_geno_path)
            if het:
                wgs_qc.run_het()
                # run het on each shard
                # merge outputs from each shard
            if callrate:
                wgs_qc.run_callrate()
                # run callrate on each shard
                # merge outputs from each shard
            if sex:
                wgs_qc.run_sex_check()
                # within sex_check() method:
                # find shards with chrX range, merge, run --sex on specific chrX range
            if related/duplicate:
                wgs_qc.run_related_dups()
                # within related() method:
                # extract ref snps, merge, run relatedness
    '''

    def __init__(self, shards_dir=None, out_path=None, keep_all=True, slurm=False, slurm_user=None, \
                 shard_key_path=None, \
                 preBqsr_path=None, \
                 wgs_metrics_path=None, \
                 variant_calling_summary_metrics_path=None, \
                 ref_panel_path=None):

        self.shards_dir = shards_dir
        self.out_path = out_path
        self.keep_all = keep_all

        # wgs qc metrics files
        self.shard_key = shard_key_path
        self.preBqsr = preBqsr_path
        self.wgs_metrics = wgs_metrics_path
        self.var_calling_summary_metrics_path = variant_calling_summary_metrics_path

        self.ref_panel_path = ref_panel_path

        # create list of shard filenames
        self.shard_filenames = list(set([f"{self.shards_dir}/{f.split('.')[0]}" for f in os.listdir(self.shards_dir)]))

        self.slurm = slurm
        self.slurm_user = slurm_user
        self.slurm_scripts = f'{self.out_path}_slurm_scripts'
        if self.slurm:
            os.makedirs(f'{self.slurm_scripts}', exist_ok=True)
            os.makedirs(f'{self.slurm_scripts}/logs', exist_ok=True)

    def run_freemix_check(self, geno_path=None, check_limit=1e6):
        out_path = self.out_path
        keep_all = self.keep_all
        preBqsr_path = self.preBqsr
        shard_filenames = self.shard_filenames

        step = "wgs_freemix_filter_prune"

        if geno_path is None:
            geno_path = shard_filenames[0]

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

        freemix_flagged_df = preBqsr[preBqsr.FREEMIX >= 0.03]
        freemix_flagged_ids = freemix_flagged_df['sample_id']
        freemix_flagged_count = freemix_flagged_ids.shape[0]
        freemix_flagged_df.to_csv(freemix_flagged, sep='\t', header=False, index=False)

        # PER ZH: 'Problem samples MAY be failures and should be removed if CHECK value is large (need to define limits)'
        # TODO: define CHECK limit (arbitrarily set default as 1e6 for now)
        freemix_fail_df = freemix_flagged_df[freemix_flagged_df.CHECK >= check_limit]
        freemix_fail_ids = freemix_fail_df['sample_id']
        freemix_fail_count = freemix_fail_ids.shape[0]
        freemix_fail_ids.to_csv(freemix_fails, sep='\t', header=False, index=False)
        if not keep_all:
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


    def run_coverage_check(self, geno_path=None, min_mean_coverage=25):
        out_path = self.out_path
        keep_all = self.keep_all
        wgs_metrics_path = self.wgs_metrics
        shard_filenames = self.shard_filenames

        step = "coverage_check"

        if geno_path is None:
            geno_path = shard_filenames[0]

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

        # PER ZH: 'Failed samples should be noted (not removed)'
        if not keep_all:
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


    def run_titv_check(self, geno_path=None):
        out_path = self.out_path
        variant_calling_summary_metrics_path = self.var_calling_summary_metrics_path
        # titv ratio fails should be removed regardless
        # keep_all = self.keep_all
        shard_filenames = self.shard_filenames

        step = "titv_check"

        if geno_path is None:
            geno_path = shard_filenames[0]

        # create filenames
        titv_ratio_fails = f'{out_path}.titv_ratio_fails'

        vc_summary_metrics = pd.read_csv(variant_calling_summary_metrics_path, sep='\s+')
        geno_samples = pd.read_csv(f'{geno_path}.psam', sep='\s+')
        vc_summary_metrics = pd.merge(vc_summary_metrics, geno_samples, how='inner', left_on='sample_id', right_on='IID')
        if vc_summary_metrics.shape[0] != geno_samples.shape[0]:
            print('Warning: metrics file does not contain info for all samples in geno files!')

        titv_ratio_fail_df = vc_summary_metrics[vc_summary_metrics.DBSNP_TITV < 2] # variable? or hard-code?
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


    def run_het_check(self, geno_path, het_out):
        '''
        called within merge_sample_het() on each shard
        '''
        out_path = self.out_path
        slurm = self.slurm

        step = "het_check"

        if slurm:
            # create filenames
            shard_geno_path = f'{os.path.dirname(geno_path)}/${{shard_name}}'
            het_tmp = f'{het_out}_${{shard_name}}_tmp'
            het_tmp2 = f'{het_out}_${{shard_name}}_tmp2'
            het_tmp3 = f'{het_out}_${{shard_name}}_tmp3'

            plink_cmd1 = f'{plink2_exec} --pfile {shard_geno_path} --geno 0.01 --maf 0.05 --indep-pairwise 50 5 0.5 --out {het_tmp}'
            plink_cmd2 = f'{plink2_exec} --pfile {shard_geno_path} --extract {het_tmp}.prune.in --make-pgen psam-cols=fid,parents,sex,pheno1,phenos --out {het_tmp2}'
            plink_cmd3 = f'{plink2_exec} --pfile {het_tmp2} --het --out {het_tmp3}'

            plink_cmds = [plink_cmd1, plink_cmd2, plink_cmd3]

        else:
            # create filenames
            het_tmp = f"{het_out}_tmp"
            het_tmp2 = f"{het_out}_tmp2"
            het_tmp3 = f"{het_out}_tmp3"

            plink_cmd1 = f"{plink2_exec} --pfile {geno_path} --geno 0.01 --maf 0.05 --indep-pairwise 50 5 0.5 --out {het_tmp}"
            plink_cmd2 = f"{plink2_exec} --pfile {geno_path} --extract {het_tmp}.prune.in --make-pgen psam-cols=fid,parents,sex,pheno1,phenos --out {het_tmp2}"
            plink_cmd3 = f"{plink2_exec} --pfile {het_tmp2} --het --out {het_tmp3}"
            plink_cmds = [plink_cmd1, plink_cmd2, plink_cmd3]

            for cmd in plink_cmds:
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

        return plink_cmds


    def merge_sample_het(self, het_filter=[-0.25,0.25]):
        '''
        running het on all shards
        merging all of the shards' .het output for sample het qc
        '''
        out_path = self.out_path
        shard_filenames = self.shard_filenames
        slurm = self.slurm
        slurm_user = self.slurm_user
        slurm_scripts = self.slurm_scripts
        keep_all = self.keep_all

        step = "het_check"

        # create het_out dir
        out_filename = os.path.basename(out_path)
        het_dir = f'{os.path.dirname(out_path)}/{out_filename}_het'
        if not os.path.exists(f'{het_dir}'):
            os.makedirs(f'{het_dir}')

        het_out = f'{het_dir}/{out_filename}'
        if slurm:
            # create config file for submitting job array
            with open(f'{slurm_scripts}/sample_het.config', 'w') as f:
                f.write(f'ArrayTaskID\tshard_name\n')
                for i in range(len(shard_filenames)):
                    if i<len(shard_filenames)-1:
                        f.write(f'{i+1}\t{os.path.basename(shard_filenames[i])}\n')
                    else:
                        f.write(f'{i+1}\t{os.path.basename(shard_filenames[i])}')
            f.close()

            # create sbatch file
            with open(f'{slurm_scripts}/sample_het.sh', 'w') as f:
                f.write(f'#!/usr/bin/env bash\n\n')
                f.write(f'config={slurm_scripts}/sample_het.config\n')
                f.write(f"shard_name=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {{print $2}}' $config)\n\n")
                cmds = self.run_het_check(geno_path=shard_filenames[0], het_out=f'{het_out}')
                for cmd in cmds:
                    f.write(f'{cmd}\n\n')
            f.close()

            # run slurm
            slurm_cmd = f'sbatch --cpus-per-task=4 --mem=32g \
                            --array=1-{len(shard_filenames)} \
                            --output={slurm_scripts}/logs/sample_het_%A_%a.out \
                            --time=0:30:0 {slurm_scripts}/sample_het.sh'

            job_id = shell_do(slurm_cmd, return_log=True).strip()

            # ensure all jobs are 'COMPLETED' before concatting logs
            sacct_out = list()
            job_complete = dict()

            # only checks of 'COMPLETED' or 'FAILED' job status
            while (not all(status in ['COMPLETED', 'FAILED'] for status in job_complete.values())) or len(job_complete)!=len(shard_filenames):
                sacct_cmd = f'sacct --user {slurm_user} -X -b -j {job_id}'
                sacct_out = shell_do(sacct_cmd, return_log=True)
                sacct_list = [val.split() for val in sacct_out.split('\n') if job_id in val]
                ids = list()
                statuses = list()
                for val in sacct_list:
                    ids.append(val[0])
                    statuses.append(val[1])
                job_complete = dict(zip(ids, statuses))

            # concat logs
            listOfFiles = list()
            for shard in shard_filenames:
                shard_name = shard.split('/')[-1]
                het_log = f'{het_out}_{shard_name}_tmp.log'
                het_log2 = f'{het_out}_{shard_name}_tmp2.log'
                het_log3 = f'{het_out}_{shard_name}_tmp3.log'
                listOfFiles.append(het_log)
                listOfFiles.append(het_log2)
                listOfFiles.append(het_log3)
            concat_logs(step, out_path, listOfFiles)

        else:
            # loop through shards and run het
            for shard in shard_filenames:
                shard_name = shard.split('/')[-1]
                het_outpath = f'{het_out}_{shard_name}'
                self.run_het_check(shard, het_outpath)

        # create filenames
        het_fails = f'{out_path}.het_fails'

        if len(glob.glob1(het_dir,'*.het')) == len(shard_filenames):
            heterozygosities = defaultdict(list)
            for filename in os.listdir(het_dir):
                if filename.endswith('.het'):
                    het = pd.read_csv(f'{callrate_dir}/filename', sep='\s+')
                    het['F'] = het['F'].map(lambda x:[x])
                    shard_hets = pd.Series(het.F.values, index=het.IID).to_dict()
                    for k,v in shard_hets.items():
                        heterozygosities[k].extend(v)

            for k,v in heterozygosities.items():
                heterozygosities[k] = sum(v) / len(v)

            heterozygosities = pd.DataFrame.from_dict(heterozygosities, orient='index').reset_index()
            heterozygosities = heterozygosities.rename(columns={'index':'IID', 0:'F'})

            outliers = heterozygosities[((heterozygosities.F <= het_filter[0]) | (heterozygosities.F >= het_filter[1]))]
            het_fail_ids = outliers['IID']
            het_fail_count = het_fail_ids.shape[0]
            het_fail_ids.to_csv(het_fails, sep='\t', header=False, index=False)

            process_complete = True

            outfiles_dict = {
                'pruned_samples': f'{het_fails}',
                'plink_out': f'{out_path}' # TODO: if removing include plink_out
            }

            metrics_dict = {
                'outlier_count': het_fail_count
            }
        else:
            process_complete = False
            outfiles_dict = {
                'pruned_samples': 'Heterozygosity Pruning Failed!',
                'plink_out': f'{out_path}' # TODO: if removing include plink_out
            }

            metrics_dict = {
                'outlier_count': 0
            }

            print(f'At least one file failed WGS heterozygosity pruning!')
            print(f'Check the {out_path}.log for more information')

        out_dict = {
            'pass': process_complete,
            'step': step,
            'metrics': metrics_dict,
            'output': outfiles_dict
        }

        return out_dict


    def run_check_callrate(self, geno_path, callrate_out):
        '''
        called within merge_callrate_check() on each shard
        '''
        out_path = self.out_path
        slurm = self.slurm

        step = "callrate_check"

        if slurm:
            plink_cmd = f'{plink2_exec} --pfile {os.path.dirname(geno_path)}/${{shard_name}} --missing --out {callrate_out}_${{shard_name}}_callrates'

        else:
            # create filenames
            callrates = f'{callrate_out}_callrates'

            plink_cmd = f'{plink2_exec} --pfile {geno_path} --missing --out {callrates}'
            shell_do(plink_cmd)

            listOfFiles = [f'{callrates}.log']
            concat_logs(step, out_path, listOfFiles)

        return plink_cmd


    def merge_sample_callrate(self, callrate_threshold=0.95):
        '''
        running callrate on all shards
        merging all the shards' .callrate output for sample callrate qc
        '''
        out_path = self.out_path
        shard_filenames = self.shard_filenames
        slurm = self.slurm
        slurm_user = self.slurm_user
        slurm_scripts = self.slurm_scripts
        keep_all = self.keep_all

        step = "callrate_check"

        # create callrate_out dir
        out_filename = os.path.basename(out_path)
        callrate_dir = f'{os.path.dirname(out_path)}/{out_filename}_callrate'
        if not os.path.exists(f'{callrate_dir}'):
            os.makedirs(f'{callrate_dir}')

        callrate_out = f'{callrate_dir}/{out_filename}'
        if slurm:
            # create config file for submitting job array
            with open(f'{slurm_scripts}/sample_callrate.config', 'w') as f:
                f.write(f'ArrayTaskID\tshard_name\n')
                for i in range(len(shard_filenames)):
                    if i<len(shard_filenames)-1:
                        f.write(f'{i+1}\t{os.path.basename(shard_filenames[i])}\n')
                    else:
                        f.write(f'{i+1}\t{os.path.basename(shard_filenames[i])}')
            f.close()

            # create sbatch file
            with open(f'{slurm_scripts}/sample_callrate.sh', 'w') as f:
                f.write(f'#!/usr/bin/env bash\n\n')
                f.write(f'config={slurm_scripts}/sample_callrate.config\n')
                f.write(f"shard_name=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {{print $2}}' $config)\n\n")
                cmd = self.run_check_callrate(geno_path=shard_filenames[0], callrate_out=f'{callrate_out}')
                f.write(f'{cmd}\n')
            f.close()

            # run slurm
            slurm_cmd = f'sbatch --cpus-per-task=4 --mem=32g \
                            --array=1-{len(shard_filenames)} \
                            --output={slurm_scripts}/logs/sample_callrate_%A_%a.out \
                            --time=0:30:0 {slurm_scripts}/sample_callrate.sh'

            job_id = shell_do(slurm_cmd, return_log=True).strip()

            # ensure all jobs are 'COMPLETED' before concatting logs
            sacct_out = list()
            job_complete = dict()
            while (not all(status in ['COMPLETED', 'FAILED'] for status in job_complete.values())) or len(job_complete)!=len(shard_filenames):
                sacct_cmd = f'sacct --user {slurm_user} -X -b -j {job_id}'
                sacct_out = shell_do(sacct_cmd, return_log=True)
                sacct_list = [val.split() for val in sacct_out.split('\n') if job_id in val]
                ids = list()
                statuses = list()
                for val in sacct_list:
                    ids.append(val[0])
                    statuses.append(val[1])
                job_complete = dict(zip(ids, statuses))

            # concat logs
            listOfFiles = list()
            for shard in shard_filenames:
                shard_name = shard.split('/')[-1]
                callrate_log = f'{callrate_out}_{shard_name}_callrates.log'
                listOfFiles.append(callrate_log)
            concat_logs(step, out_path, listOfFiles)
        else:
            # loop through each shard in shards_dir and run callrate check
            for shard in shard_filenames:
                shard_name = shard.split('/')[-1]
                callrate_outpath = f'{callrate_out}_{shard_name}'
                self.run_check_callrate(shard, callrate_outpath)

        # create filenames
        callrate_fails = f'{out_path}.callrate_fails'

        if len(glob.glob1(callrate_dir,'*.smiss')) == len(shard_filenames):
            callrates = defaultdict(list)
            for filename in os.listdir(callrate_dir):
                if filename.endswith('.smiss'):
                    smiss = pd.read_csv(f'{callrate_dir}/{filename}', sep='\s+')
                    smiss['CALLRATE'] = smiss['F_MISS'].map(lambda x: [1-x])
                    shard_callrates = pd.Series(smiss.CALLRATE.values, index=smiss.IID).to_dict()
                    for k,v in shard_callrates.items():
                        callrates[k].extend(v)

            for k,v in callrates.items():
                callrates[k] = sum(v) / len(v)

            callrates = pd.DataFrame.from_dict(callrates, orient='index').reset_index()
            callrates = callrates.rename(columns={'index':'IID', 0:'CALLRATE'})

            outliers = callrates[callrates.CALLRATE < callrate_threshold]
            callrate_fail_ids = outliers['IID']
            callrate_fail_count = callrate_fail_ids.shape[0]
            callrate_fail_ids.to_csv(f'{callrate_fails}', sep='\t', header=False, index=False)

            process_complete = True

            outfiles_dict = {
                'pruned_samples': f'{callrate_fails}',
                'plink_out': f'{out_path}' # TODO: if removing flagged samples, diff plink_out
            }

            metrics_dict = {
                'outlier_count': callrate_fail_count
            }

        else:
            print('At least one file failed WGS callrate pruning!')
            print(f'Check {out_path}.log for more information')

            process_complete = False

            outfiles_dict = {
                'pruned_samples': 'WGS Callrate Pruning Failed!',
                'plink_out': f'{out_path}'
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


    def find_x_shard(self):
        '''
        use shard key to find shard with interval on x chromosome for sex check
        called within sex check
        '''
        shards_dir = self.shards_dir
        shard_key = self.shard_key
        out_path = self.out_path
        shard_key = pd.read_csv(shard_key, dtype={'shard':str})
        shard_key['start'] = shard_key['interval'].str.split(':', expand=True)[1].str.split('-', expand=True)[0].astype(int)
        shard_key['end'] = shard_key['interval'].str.split(':', expand=True)[1].str.split('-', expand=True)[1].astype(int)
        shard_key_x = shard_key[shard_key['chr']=='chrX']

        x_start = 2781479
        x_end = 155701383

        shard_key_x_range = shard_key_x[((shard_key_x['start']>=x_start) & (shard_key_x['start']<=x_end)) \
                                | ((shard_key_x['end']>=x_start) & (shard_key_x['end']<=x_end)) \
                                | ((shard_key_x['start']<=x_start) & (shard_key_x['end']>=x_end))]

        # make dir to keep merged X chr
        out_filename = os.path.basename(out_path)
        x_dir = f'{os.path.dirname(out_path)}/{out_filename}_chrX'
        if not os.path.exists(f'{x_dir}'):
            os.makedirs(f'{x_dir}')

        # convert all X chr shards into bfiles
        x_range_shards = sorted(set(zip(shard_key_x_range.chr, shard_key_x_range.shard)))
        for chrom, shard in x_range_shards:
            # requires files to be named 'chr_shard'
            bfiles_to_pfiles(pfile_path=f'{shards_dir}/{chrom}_{shard}')

        # merge X chr together to run sex check
        with open(f'{x_dir}/to_merge.txt', 'w') as f:
            for chrom, shard in x_range_shards:
                path = f'{shards_dir}/{chrom}_{shard}'
                f.write(f'{path}\n')
        f.close()

        plink_merge = f'{plink_exec} --merge-list {x_dir}/to_merge.txt --make-bed --out {x_dir}/x_range'
        shell_do(plink_merge)

        for chrom, shard in x_range_shards:
            os.remove(f'{shards_dir}/{chrom}_{shard}.bed')
            os.remove(f'{shards_dir}/{chrom}_{shard}.bim')
            os.remove(f'{shards_dir}/{chrom}_{shard}.fam')

        # return out path to X chr plink files containing sex check range
        return f'{x_dir}/x_range'


    def run_sex_check(self, check_sex=[0.25,0.75]):
        out_path = self.out_path

        step = "wgs_sex_check"

        # generate chrX range geno file
        geno_path = self.find_x_shard()

        # create filenames
        sex_tmp = f'{out_path}_tmp'
        sex_fails = f'{out_path}.sex_fails'

        plink_cmd = f'{plink_exec} --bfile {geno_path} --chr 23 --from-bp 2781479 --to-bp 155701383 \
                    --maf 0.05 --geno 0.05 --hwe 1E-5 --check-sex {check_sex[0]} {check_sex[1]} --out {sex_tmp}'
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
        shard_key = self.shard_key
        out_path = self.out_path

        # format shards key
        shard_key = pd.read_csv(shard_key, dtype={'shard':str})
        shard_key['start'] = shard_key['interval'].str.split(':', expand=True)[1].str.split('-', expand=True)[0].astype(int)
        shard_key['end'] = shard_key['interval'].str.split(':', expand=True)[1].str.split('-', expand=True)[1].astype(int)

        # load in ref panel snps
        # TODO: assume ref panel snps in repo?
        ref_panel = pd.read_csv(f'{ref_panel_path}', sep='\s+', header=None, names=['snpID'])
        ref_panel['CHR'] = ref_panel['snpID'].str.split(':', expand=True)[0]
        ref_panel['POS'] = ref_panel['snpID'].str.split(':', expand=True)[1].astype(int)

        # find shard that contains each ref snp
        ref_shard_merge = ref_panel.merge(shard_key, left_on='CHR', right_on='chr', how='outer')
        ref_shard_interval = ref_shard_merge[(ref_shard_merge['POS']>=ref_shard_merge['start']) & (ref_shard_merge['POS']<=ref_shard_merge['end'])]
        ref_panel_with_shard = ref_shard_interval[['snpID', 'CHR', 'POS', 'shard']]

        # make dir for extracted ref panel
        out_filename = os.path.basename(out_path)
        ref_overlap_dir = f'{os.path.dirname(out_path)}/{out_filename}_ref_overlap'
        if not os.path.exists(f'{ref_overlap_dir}'):
            os.makedirs(f'{ref_overlap_dir}')

        # extract ref snps from each corresponding shard
        shards_to_extract_from = set(zip(ref_panel_with_shard.CHR, ref_panel_with_shard.shard))
        for chrom, shard in shards_to_extract_from:
            # requires files to be named 'chr_shard'
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


    def run_relatedness(self):
        '''
        check for both duplicates and relatedness
        '''
        out_path = self.out_path
        ref_panel_path = self.ref_panel_path

        step = "relatedness_check"

        subset_geno = self.extract_ref_snps(ref_panel_path)
        sampleQC = SampleQC(subset_geno, out_path)
        related_out_dict = sampleQC.run_related_prune(prune_related=False, prune_duplicated=False)
        related_out_files = related_out_dict['output']
        related_samples_path = related_out_files['related_samples']

        if os.path.isfile(f'{related_samples_path}'):
            rel_samples_to_exclude = self.drop_related(related_samples_path)
            rel_samples_to_exclude_count = count_file_lines(rel_samples_to_exclude)
            # TODO: should we remove these samples from each shard?
            process_complete = True

            outfiles_dict = {
                'samples_to_exclude': rel_samples_to_exclude,
                'plink_out': out_path
            }

            metrics_dict = {
                'outlier_count': rel_samples_to_exclude_count
            }

        else:
            print('WGS Relatedness Failed!')
            print(f'Check all_plink_logs.gtlog for more information')

            process_complete = False

            outfiles_dict = {
                'samples_to_exclude': None,
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


    def drop_related(self, related_samples_path):
        '''
        drop duplicates (and related pairs) with less coverage
        '''
        preBqsr_path = self.preBqsr
        out_path = self.out_path

        step = "relatedness_check"

        preBqsr = pd.read_csv(preBqsr_path, sep='\s+')
        related_samples = pd.read_csv(related_samples_path)

        if related_samples.shape[0] != 0:
            depth_df = pd.merge(related_samples, preBqsr.rename(columns={"AVG_DP":"AVG_DP1"}), left_on='IID1', right_on='sample_id', how='left')
            depth_df = pd.merge(depth_df, preBqsr.rename(columns={"AVG_DP":"AVG_DP2"}), left_on='IID2', right_on='sample_id', how='left')
            related = pd.DataFrame(depth_df, columns=['IID1', 'AVG_DP1', 'IID2', 'AVG_DP2'])
            # related.to_csv(f'{related_samples_path}.depths', sep='\t', index=False)
            # TODO: what to do when one of the samples doesn't have depth info?
            related['lower_depth'] = related[['AVG_DP1', 'AVG_DP2']].idxmin(axis=1).str.replace("AVG_DP","")
            related.loc[related['lower_depth']=='1','rel_samples_to_exclude'] = related.loc[related['lower_depth']=='1','IID1']
            related.loc[related['lower_depth']=='2','rel_samples_to_exclude'] = related.loc[related['lower_depth']=='2','IID2']
            rel_samples_to_exclude = set(related['rel_samples_to_exclude'])

        with open(f'{out_path}.rel_samples_to_exclude', 'w') as f:
            for sample in rel_samples_to_exclude:
                f.write(f'{sample}\n')

        return f'{out_path}.rel_samples_to_exclude'