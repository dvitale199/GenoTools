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


    def run_freemix_check():
        # oreVqsr_selfSM.tsv file?
        pass

    def run_coverage_check():
        # wgs_metrics.tsv file?
        pass

    def run_duplication_check():
        # ZH doesn't prune anything here
        pass

    def run_titv_check():
        # variant_valling_summary_metrics.tsv file?
        pass

    def run_het_check():
        # can we use the method in qc.py file?
        # or do we want to do per chr method here?
        pass

    def run_check_callrate():
        # again can we use the method in qc.py file?
        # or do we want to do per chr method here?
        # some call_rates.imiss file?
        pass

    def run_sex_check():
        # again can we use the method in qc.py file?
        # but we want to extract the chr23 to check just that for sex?
        pass

    def run_genetic_dups_check():
        # utilizes KING software
        # runs --related (like in qc.py)
        # then drops the sample with lower depth (using data from preBqsr file)
        pass

    

