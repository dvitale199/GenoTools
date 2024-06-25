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
from collections import defaultdict
from genotools.utils import shell_do, concat_logs, bfiles_to_pfiles
from genotools.dependencies import check_plink, check_plink2, check_king

plink_exec = check_plink()
plink2_exec = check_plink2()
king_exec = check_king()

class WholeGenomeSeqQC:

    '''
    execute_wgs_pipeline(steps, steps_dict, shards_dir, out_path, wgs_qc, ancestry, assoc, args, tmp_dir):
        # all functionality within WGS_qc class
        if step in wgs_steps:
            wgs_qc = WholeGenomeSeqQC(shards_dir, out_path, swarm, keep, wgs_qc_metrics_files)
            # take the first shard (assuming all shards have all samples)
            # and run freemix, coverage, titv_ratio on one shard
            if freemix:
                wgs_qc.run_freemix_check(shard_geno_path, qc_file, check_limit)
            if coverage:
                wgs_qc.run_coverage_check(shard_geno_path, qc_file, min_mean_coverage)
            if titv_ratio:
                wgs_qc.run_titv_ratio_check(shard_geno_path, qc_file)
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

    def __init__(self, shards_dir=None, out_path=None, keep_all=True, swarm=False, shard_key_path=None, preBqsr_path=None, wgs_metrics_path=None, variant_calling_summary_metrics_path=None):
        self.shards_dir = shards_dir
        self.out_path = out_path
        self.keep_all = keep_all

        # wgs qc metrics files
        self.shard_key = shard_key_path
        self.preBqsr = preBqsr_path
        self.wgs_metrics = wgs_metrics_path
        self.var_calling_summary_metrics_path = variant_calling_summary_metrics_path

        # create list of shard filenames
        self.shard_filenames = list(set([f.split('.')[0] for f in os.listdir(self.shards_dir)]))

        self.swarm = swarm
        if self.swarm:
            os.makedirs(f'{self.out_path}/swarm_scripts')

    def run_freemix_check(self, geno_path, preBqsr_path, check_limit=1e6):
        out_path = self.out_path
        keep_all = self.keep_all

        step = "wgs_freemix_filter_prune"

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


    def run_coverage_check(self, geno_path, wgs_metrics_path, min_mean_coverage=25):
        out_path = self.out_path
        keep_all = self.keep_all

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


    def run_titv_check(self, geno_path, variant_calling_summary_metrics_path):
        out_path = self.out_path
        # titv ratio fails should be removed regardless
        # keep_all = self.keep_all

        step = "titv_check"

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
        step = "het_check"

        # create filenames
        het_tmp = f"{het_out}_tmp"
        het_tmp2 = f"{het_out}_tmp2"
        het_tmp3 = f"{het_out}_tmp3"

        plink_cmd1 = f"{plink2_exec} --pfile {geno_path} --geno 0.01 --maf 0.05 --indep-pairwise 50 5 0.5 --out {het_tmp}"
        plink_cmd2 = f"{plink2_exec} --pfile {geno_path} --extract {het_tmp}.prune.in --make-pgen psam-cols=fid,parents,sex,pheno1,phenos --out {het_tmp2}"
        plink_cmd3 = f"{plink2_exec} --pfile {het_tmp2} --het --out {het_tmp3}"

        cmds = [plink_cmd1, plink_cmd2, plink_cmd3]

        for cmd in cmds:
            shell_do(cmd)

        listOfFiles = [f'{het_tmp}.log', f'{het_tmp2}.log', f'{het_tmp3}.log']
        concat_logs(step, het_out, listOfFiles)

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


    def merge_sample_het(self, het_filter=[-0.25,0.25]):
        '''
        merging all of the shards' .het output for sample het qc
        '''
        out_path = self.out_path
        shard_filenames = self.shard_filenames
        swarm = self.swarm
        keep_all = self.keep_all

        step = "het_check"

        # create het_out dir
        out_filename = os.path.basename(out_path)
        het_dir = f'{os.path.dirname(out_path)}/het_out'
        if not os.path.exists(f'{het_dir}'):
            os.makedirs(f'{het_dir}')

        het_out = f'{het_dir}/{out_filename}'
        if swarm:
            # TODO: create swarm file and run het
            pass
        else:
            # loop through shards and run het
            for shard in shard_filenames:
                self.run_het_check(shard, het_out)

        # create filenames
        het_fails = f'{out_path}.het_fails'

        heterosygosities = defaultdict(list)
        for filename in os.listdir(het_dir):
            if filename.endswith('.het'):
                het = pd.read_csv(filename, sep='\s+')
                het['F'] = het['F'].map(lambda x:[x])
                shard_hets = pd.Series(het.F.values, index=het.IID).to_dict()
                for k,v in shard_hets.items():
                    heterosygosities[k].extend(v)

        for k,v in heterozygosities.items():
            heterosygosities[k] = sum(v) / len(v)

        heterozygosities = pd.DataFrame.from_dict(heterozygosities, orient='index').reset_index()
        heterozygosities = heterozygosities.rename(columns={'index':'IID', 0:'F'})

        outliers = heterosygosities[((heterosygosities.F <= het_filter[0]) | (heterosygosities.F >= het_filter[1]))]
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

        step = "callrate_check"

        # create filenames
        callrates = f'{callrate_out}_callrates'

        plink_cmd = f'{plink2_exec} --pfile {geno_path} --missing --out {callrates}'
        shell_do(plink_cmd)

        listOfFiles = [f'{callrates}.log']
        concat_logs(step, out_path, listOfFiles)


    def merge_sample_callrate(self):
        '''
        merging all the shards' .callrate output for sample callrate qc
        '''
        out_path = self.out_path
        shard_filenames = self.shard_filenames
        swarm = self.swarm
        keep_all = self.keep_all

        step = "callrate_check"

        # create callrate_out dir
        out_filename = os.path.basename(out_path)
        callrate_dir = f'{os.path.dirname(out_path)}/callrate_dir'
        if not os.path.exists(f'{callrate_dir}'):
            os.makedirs(f'{callrate_dir}')

        callrate_out = f'{callrate_dir}/{out_filename}'
        if swarm:
            # TODO: create swarm file and run callrate check on each shard
            pass
        else:
            # loop through each shard in shards_dir and run callrate check
            for shard in shard_filenames:
                self.run_check_callrate(shard, callrate_out)

        # create filenames
        callrate_fails = f'{out_path}.callrate_fails'

        callrates = defaultdict(list)
        for filename in os.listdir(callrate_dir):
            if filename.endswith('.smiss'):
                smiss = pd.read_csv(f'{filename}.smiss', sep='\s+')
                smiss['CALLRATE'] = smiss['F_MISS'].map(lambda x: [1-x])
                shard_callrates = pd.Series(smiss.CALLRATE.values, index=smiss.IID).to_dict()
                for k,v in shard_callrates.items():
                    callrates[k].extend(v)

        for k,v in callrates.items():
            callrates[k] = sum(v) / len(v)

        callrates = pd.DataFrame.from_dict(callrates, orient='index').reset_index()
        callrates = callrates.rename(columns={'index':'IID', 0:'CALLRATE'})

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


    def find_x_shard(self):
        '''
        use shard key to find shard with interval on x chromosome for sex check
        called within sex check
        '''
        shards_dir = self.shards_dir
        shard_key = self.shard_key
        shard_key = pd.read_csv(shard_key, sep='\s+', dtype={'shard':str})
        shard_key['start'] = shard_key['interval'].str.split(':', expand=True)[1].str.split('-', expand=True)[0].astype(int)
        shard_key['end'] = shard_key['interval'].str.split(':', expand=True)[1].str.split('-', expand=True)[1].astype(int)
        shard_key_x = shard_key[shard_key['chr']=='chrX']

        x_start = 2781479
        x_end = 155701383

        shard_key_x_range = shard_key_x[((shard_key_x['start']>=x_start) & (shard_key_x['start']<=x_end)) \
                                | ((shard_key_x['end']>=x_start) & (shard_key_x['end']<=x_end)) \
                                | ((shard_key_x['start']<=x_start) & (shard_key_x['end']>=x_end))]

        # make dir to keep merged X chr
        x_dir = f'{shards_dir}/x_range'
        if not os.path.exists(x_dir):
            os.makedirs(x_dir)

        # merge X chr together to run sex check
        x_range_shards = sorted(set(zip(shard_key_x_range.chr, shard_key_x_range.shard)))
        with open(f'{x_dir}/to_merge.txt', 'w') as f:
            for chrom, shard in x_range_shards:
                path = f'{shards_dir}/{chrom}_{shard}'
                f.write(f'{path}\n')
        f.close()

        plink_merge = f'plink2 --pmerge-list {shards_dir}/x_range/to_merge.txt --make-pgen --out {x_dir}/x_range'
        shell_do(plink_merge)

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

        # convert to bfiles
        bfiles_to_pfiles(pfile_path=geno_path)

        plink_cmd = f'{plink_exec} --bfile {geno_path} --chr 23 --from-bp2781479 --to-bp 155701383 --maf 0.05 --geno 0.05 --hwe 1E-5 --check-sex  {check_sex[0]} {check_sex[1]} --out {sex_tmp}'
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

        # format shards key
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
        ref_overlap_dir = f'{shards_dir}/ref_overlap'
        if not os.path.exists(f'{ref_overlap_dir}'):
            os.makedirs(f'{ref_overlap_dir}')

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


    def run_relatedness(self):
        '''
        check for both duplicates and relatedness
        '''
        # TODO
        pass

    def drop_duplicates(self):
        '''
        drop duplicates (and related pairs) with less coverage
        '''
        # TODO
        pass
