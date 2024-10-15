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

class SampleQC:

    def __init__(self, geno_path=None, out_path=None):
        self.geno_path = geno_path
        self.out_path = out_path

    def run_callrate_prune(self, mind=0.02):

        """
        Execute call rate pruning on genotype data using PLINK.

        This function performs call rate pruning based on a user-specified threshold for
        the missing genotype rate. If an individual's missing genotype rate surpasses the
        given threshold, the individual will be excluded from the dataset.

        Parameters:
        - mind (float, optional): Threshold for permissible missing genotype rate. Default is set to 0.02, which corresponds to 2%.

        Returns:
        - dict: A structured dictionary containing:
            * 'pass': Boolean indicating the successful completion of the process.
            * 'step': The label for this procedure ('callrate_prune').
            * 'metrics': Metrics associated with the pruning, such as the 'outlier_count'.
            * 'output': Dictionary containing paths to the generated output files like 'pruned_samples' and 'plink_out'.
        """

        geno_path = self.geno_path
        out_path = self.out_path

        # Check that paths are set
        if geno_path is None or out_path is None:
            raise ValueError("Both geno_path and out_path must be set before calling this method.")

        # Check path validity
        if not os.path.exists(f'{geno_path}.pgen'):
            raise FileNotFoundError(f"{geno_path} does not exist.")

        # Check type of mind
        if not isinstance(mind, (int, float)):
            raise TypeError("mind should be of type int or float.")

        # Check valid range for mind
        if mind < 0 or mind > 1:
            raise ValueError("mind should be between 0 and 1.")

        step = "callrate_prune"

        outliers_out = f'{out_path}.outliers'

        plink_cmd1 = f"{plink2_exec} --pfile {geno_path} --mind {mind} --make-pgen psam-cols=fid,parents,sex,pheno1,phenos --out {out_path}"

        shell_do(plink_cmd1)

        listOfFiles = [f'{out_path}.log']
        concat_logs(step, out_path, listOfFiles)

        if os.path.isfile(f'{out_path}.mindrem.id'):
            irem = pd.read_csv(f'{out_path}.mindrem.id', sep='\s+')
            irem.to_csv(outliers_out, sep='\t', header=True, index=False)

            outlier_count = sum(1 for line in open(f'{outliers_out}')) - 1

            # remove mindrem.id file
            os.remove(f'{out_path}.mindrem.id')

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


    def run_sex_prune(self, check_sex=[0.25,0.75]):

        """
        Execute sex-based pruning on genotype data using PLINK.

        Parameters:
        - check_sex (list of float, optional): Two values indicating the lower and upper bounds for sex discrepancy checks. Default values are [0.25, 0.75].

        Returns:
        - dict: A structured dictionary containing:
            * 'pass': Boolean indicating the successful completion of the process.
            * 'step': The label for this procedure ('sex_prune').
            * 'metrics': Metrics associated with the pruning, such as 'outlier_count'.
            * 'output': Dictionary containing paths to the generated output files.
        """

        geno_path = self.geno_path
        out_path = self.out_path

        # Check that paths are set
        if geno_path is None or out_path is None:
            raise ValueError("Both geno_path and out_path must be set before calling this method.")

        # Check path validity
        if not os.path.exists(f'{geno_path}.pgen'):
            raise FileNotFoundError(f"{geno_path} does not exist.")

        # Check check_sex is two values
        if len(check_sex) != 2:
            raise ValueError("check_sex should contain two values.")

        # Check type of check_sex
        if (not isinstance(check_sex[0], float)) or (not isinstance(check_sex[1], float)) or (not isinstance(check_sex, list)):
            raise TypeError("check_sex values should be a list of two floats.")

        # Check check_sex is within bounds
        if (check_sex[0] < 0 or check_sex[0] > 1) or (check_sex[1] < 0 or check_sex[1] > 1):
            raise ValueError("check_sex values should be between 0 and 1.")

        step = "sex_prune"

        # create filenames
        sex_tmp1 = f"{out_path}_tmp1"
        sex_tmp2 = f"{out_path}_tmp2"
        sex_fails = f"{out_path}.outliers"

        # convert to bfiles
        bfiles_to_pfiles(pfile_path=geno_path)

        # check sex 2 methods
        plink_cmd1 = f"{plink_exec} --bfile {geno_path} --check-sex 0.25 0.75 --maf 0.05 --out {sex_tmp1}"
        plink_cmd2 = f"{plink_exec} --bfile {geno_path} --chr 23 --from-bp 2781479 --to-bp 155701383 --maf 0.05 --geno 0.05 --hwe 1E-5 --check-sex  0.25 0.75 --out {sex_tmp2}"

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
            sex_fail_df = pd.concat([sex_fail1, sex_fail2], ignore_index=True)
            sex_fail_ids = sex_fail_df.loc[:,['FID','IID']].drop_duplicates(subset=['FID','IID'])
            sex_fail_ids = sex_fail_ids.rename({'FID':'#FID'}, axis=1)
            sex_fail_count = sex_fail_ids.shape[0]
            sex_fail_ids.to_csv(sex_fails, sep='\t', header=True, index=False)

            # remove sex fail samples from geno
            plink_cmd3 = f"{plink2_exec} --pfile {geno_path} --remove {sex_fails} --make-pgen psam-cols=fid,parents,sex,pheno1,phenos --out {out_path}"

            shell_do(plink_cmd3)

            listOfFiles = [f'{out_path}.log']
            concat_logs(step, out_path, listOfFiles)

            process_complete = True

            for file in [f'{sex_tmp1}.hh',f'{sex_tmp1}.sexcheck',f'{sex_tmp2}.hh',f'{sex_tmp2}.sexcheck']:
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
            print(f'Check {sex_tmp1}.log and {sex_tmp2}.log for more information')

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


    def run_het_prune(self, het_filter=[-0.15,0.15]):

        """
        Execute heterozygosity-based pruning on genotype data using PLINK.

        Parameters:
        - het_filter (list of float, optional): Two values indicating the bounds for heterozygosity outlier checks. Default values are [-0.15, 0.15].

        Returns:
        - dict: A structured dictionary containing:
            * 'pass': Boolean indicating the successful completion of the process.
            * 'step': The label for this procedure ('het_prune').
            * 'metrics': Metrics associated with the pruning, such as 'outlier_count'.
            * 'output': Dictionary containing paths to the generated output files.
        """

        geno_path = self.geno_path
        out_path = self.out_path

        # Check that paths are set
        if geno_path is None or out_path is None:
            raise ValueError("Both geno_path and out_path must be set before calling this method.")

        # Check path validity
        if not os.path.exists(f'{geno_path}.pgen'):
            raise FileNotFoundError(f"{geno_path} does not exist.")

        # Check het_filter is two values
        if len(het_filter) != 2:
            raise ValueError("check_sex should contain two values.")

        # Check type of het_filter
        if (not isinstance(het_filter[0], float)) or (not isinstance(het_filter[1], float)) or (not isinstance(het_filter, list)):
            raise TypeError("check_sex values should be of type int or float.")

        step = "het_prune"

        het_tmp = f"{out_path}_tmp"
        het_tmp2 = f"{out_path}_tmp2"
        het_tmp3 = f"{out_path}_tmp3"
        outliers_out = f"{out_path}.outliers"

        # variant(maf=0.05, geno=0.01, indep_pairwise=[50,5,0.5])

        plink_cmd1 = f"{plink2_exec} --pfile {geno_path} --geno 0.01 --maf 0.05 --indep-pairwise 50 5 0.5 --out {het_tmp}"
        plink_cmd2 = f"{plink2_exec} --pfile {geno_path} --extract {het_tmp}.prune.in --make-pgen psam-cols=fid,parents,sex,pheno1,phenos --out {het_tmp2}"
        plink_cmd3 = f"{plink2_exec} --pfile {het_tmp2} --het --out {het_tmp3}"

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

            if (het_filter[0] == -1) and (het_filter[1] == -1):
                het['HET'] = (het['OBS_CT']-het['O(HOM)'])/het['OBS_CT']
                het_mean = het['HET'].mean()
                het_std = het['HET'].std()
                het_low = het_mean - (3*het_std)
                het_high = het_mean + (3*het_std)
                het_outliers = het[(het['HET'] < het_low) | (het['HET'] > het_high)]
            
            else:
                het_outliers = het[((het.F <= het_filter[0]) | (het.F >= het_filter[1]))]

            outlier_count = het_outliers.shape[0]
            het_outliers.to_csv(f'{outliers_out}', sep='\t', header=True, index=False)

            plink_cmd4 = f"{plink2_exec} --pfile {geno_path} --remove {outliers_out} --make-pgen psam-cols=fid,parents,sex,pheno1,phenos --out {out_path}"

            shell_do(plink_cmd4)

            listOfFiles = [f'{out_path}.log']
            concat_logs(step, out_path, listOfFiles)

            outfiles_dict = {
                'pruned_samples': outliers_out,
                'plink_out': out_path
            }

            metrics_dict = {
                'outlier_count': outlier_count
            }

            process_complete = True

            # remove intermediate files
            os.remove(f'{het_tmp}.prune.in')
            os.remove(f'{het_tmp}.prune.out')
            os.remove(f'{het_tmp2}.pgen')
            os.remove(f'{het_tmp2}.pvar')
            os.remove(f'{het_tmp2}.psam')
            os.remove(f'{het_tmp3}.het')

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


    def run_related_prune(self, related_cutoff=0.0884, duplicated_cutoff=0.354, prune_related=True, prune_duplicated=True):

        """
        Execute pruning based on relatedness and duplication checks on genotype data using PLINK and KING.

        Parameters:
        - related_cutoff (float, optional): Threshold for relatedness check. Default value is 0.0884.
        - duplicated_cutoff (float, optional): Threshold for duplication check. Default value is 0.354.
        - prune_related (bool, optional): Whether to prune related samples. Default is True.
        - prune_duplicated (bool, optional): Whether to prune duplicated samples. Default is True.

        Returns:
        - dict: A structured dictionary containing:
            * 'pass': Boolean indicating the successful completion of the process.
            * 'step': The label for this procedure ('related_prune').
            * 'metrics': Metrics associated with the pruning.
            * 'output': Dictionary containing paths to the generated output files.
        """

        geno_path = self.geno_path
        out_path = self.out_path

        step = "related_prune"

        # Check that paths are set
        if geno_path is None or out_path is None:
            raise ValueError("Both geno_path and out_path must be set before calling this method.")

        # Check path validity
        if not os.path.exists(f'{geno_path}.pgen'):
            raise FileNotFoundError(f"{geno_path} does not exist.")

        # Check type of related and duplicated cutoff
        if not (isinstance(related_cutoff, float) and isinstance(duplicated_cutoff, float)):
            raise TypeError("related_cutoff and duplicated_cutoff should be of type int or float.")

        # Check valid range for geno_threshold
        if related_cutoff < 0 or related_cutoff > 1 or duplicated_cutoff < 0 or duplicated_cutoff > 1:
            raise ValueError("related_cutoff and duplicated_cutoff should be between 0 and 1.")

        king1 = f"{out_path}_related_king"
        king2 = f"{out_path}_duplicated_king"

        related_pairs = f"{out_path}_pairs"
        related_out = f"{related_pairs}.related"
        related_pruned_out = f"{out_path}.pruned"

        # create pfiles
        # create table of related pairs
        king_cmd1 = f'{plink2_exec} --pfile {geno_path} --make-king-table --make-king triangle bin --king-table-filter {related_cutoff} --out {related_pairs}'
        # see if any samples are related (includes duplicates)
        king_cmd2 = f'{plink2_exec} --pfile {geno_path} --king-cutoff {related_pairs} {related_cutoff} --out {king1}'
        # see if any samples are duplicated (grm cutoff >= 0.354)
        king_cmd3 = f'{plink2_exec} --pfile {geno_path} --king-cutoff {related_pairs} {duplicated_cutoff} --out {king2}'

        cmds = [king_cmd1, king_cmd2, king_cmd3]
        for cmd in cmds:
            shell_do(cmd)

        listOfFiles = [f'{related_pairs}.log', f'{king1}.log', f'{king2}.log']
        concat_logs(step, out_path, listOfFiles)

        if os.path.isfile(f'{related_pairs}.kin0') and os.path.isfile(f'{king1}.king.cutoff.out.id') and os.path.isfile(f'{king2}.king.cutoff.out.id'):

            # create .related related pair sample files
            kinship = pd.read_csv(f'{related_pairs}.kin0', sep='\s+')
            kinship['REL'] = pd.cut(x=kinship['KINSHIP'], bins=[-np.inf, 0.0884, 0.177, 0.354, np.inf], labels=['unrel', 'second_deg', 'first_deg', 'duplicate'])
            kinship.to_csv(f'{related_pairs}.related', index=False)

            # create .related and .duplicated single sample files
            shutil.copy(f'{king1}.king.cutoff.out.id',f'{king1}.related')
            related_count = sum(1 for line in open(f'{king1}.related'))-1

            shutil.copy(f'{king2}.king.cutoff.out.id',f'{king2}.duplicated')
            duplicated_count = sum(1 for line in open(f'{king2}.duplicated'))-1

            related_count = related_count - duplicated_count
            duplicated = pd.read_csv(f'{king2}.duplicated', sep = '\s+')

            # concat duplicated sample ids to related sample ids, drop_duplicates(keep='last) because all duplicated would also be considered related
            if prune_related and prune_duplicated:
                plink_cmd1 = f'{plink2_exec} --pfile {geno_path} --remove {king1}.king.cutoff.out.id --make-pgen psam-cols=fid,parents,sex,pheno1,phenos --out {out_path}'
                shell_do(plink_cmd1)

                related = pd.read_csv(f'{king1}.related', sep = '\s+')
                grm_pruned = pd.concat([related, duplicated], ignore_index=True)

                if '#FID' in grm_pruned:
                    grm_pruned.rename(columns={"#IID": "IID"}, inplace = True)
                else:
                    grm_pruned['#FID'] = 0
                    grm_pruned.rename(columns={"#IID": "IID"}, inplace = True)
                grm_pruned.drop_duplicates(subset=['#FID','IID'], keep='last', inplace=True)
                grm_pruned.to_csv(related_pruned_out, sep='\t', header=True, index=False)
                process_complete = True

            if prune_duplicated and not prune_related:
                plink_cmd1 = f'{plink2_exec} --pfile {geno_path} --remove {king2}.king.cutoff.out.id --make-pgen psam-cols=fid,parents,sex,pheno1,phenos --out {out_path}'
                shell_do(plink_cmd1)

                grm_pruned = duplicated

                if '#FID' in grm_pruned:
                    grm_pruned.rename(columns={"#IID": "IID"}, inplace = True)
                else:
                    grm_pruned['#FID'] = 0
                    grm_pruned.rename(columns={"#IID": "IID"}, inplace = True)
                grm_pruned.drop_duplicates(subset=['#FID','IID'], keep='last', inplace=True)
                grm_pruned.to_csv(related_pruned_out, sep='\t', header=True, index=False)
                process_complete = True
                related_count = 0

            if not prune_related and not prune_duplicated:
                plink_cmd1 = f'echo prune_related and prune_duplicated set to False. Pruning passed'
                shell_do(plink_cmd1)

                process_complete = True
                related_pruned_out = None

                related_count = 0
                duplicated_count = 0

            if not prune_duplicated and prune_related:
                print('This option is invalid. Cannot prune related without also pruning duplicated')
                process_complete = False
                related_pruned_out = None

                related_count = 0
                duplicated_count = 0


            if os.path.isfile(f'{out_path}.log'):
                listOfFiles = [f'{out_path}.log']
                concat_logs(step, out_path, listOfFiles)

            # remove intermediate files
            os.remove(f'{king1}.king.cutoff.in.id')
            os.remove(f'{king1}.king.cutoff.out.id')
            os.remove(f'{king1}.related')
            os.remove(f'{king2}.duplicated')
            os.remove(f'{king2}.king.cutoff.in.id')
            os.remove(f'{king2}.king.cutoff.out.id')
            os.remove(f'{related_pairs}.king.bin')
            os.remove(f'{related_pairs}.king.id')
            os.remove(f'{related_pairs}.kin0')

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
            print(f'Check all_plink_logs.gtlog for more information')
            process_complete = False

            outfiles_dict = {
                'pruned_samples': 'Related Pruning Failed',
                'related_samples': None,
                'plink_out': [king1, king2]
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


    def run_confirming_kinship(self):

        """
        Find samples with discordant FIDs using PLINK and KING.

        Parameters:
        - None

        Returns:
        - dict: A structured dictionary containing:
            * 'pass': Boolean indicating the successful completion of the process.
            * 'step': The label for this procedure ('confirming_kinship').
            * 'metrics': Metrics associated with the counts.
            * 'output': Dictionary containing paths to the generated output files.
        """
        geno_path = self.geno_path
        out_path = self.out_path

        step = 'confirming_kinship'

        duplicated_cutoff = 0.354 # or greater
        first_deg_cutoff = 0.177 # to 0.354
        second_deg_cutoff = 0.0884 # to 0.177
        third_deg_cutoff = 0.0442 # to 0.0884

        # create output files
        temp = f'{out_path}_temp'
        same_fid_unrelated = f'{out_path}_same_fid.unrelated'
        diff_fid_related = f'{out_path}_diff_fid.related'
        num_same_fid_unrelated = 0
        num_diff_fid_related = 0

        # check OS is NOT macOS
        # Will only run this method if os is linux or windows
        if platform.system() != 'Linux':
        # if (platform.system() != 'Linux') or (platform.system() != 'Windows'):
            print('Relatedness Assessment can only run on a linux OS!')
            process_complete = False

            outfiles_dict = {
                'same_fid_unrelated': None,
                'diff_fid_related': None
            }

            metrics_dict = {
                'same_fid_unrelated_count': 0,
                'diff_fid_related_count': 0
            }

            out_dict = {
            'pass': process_complete,
            'step': step,
            'metrics': metrics_dict,
            'output': outfiles_dict
            }

            return out_dict


        # if data is bfiles, proceed
        if not os.path.isfile(f'{geno_path}.bed'):
            # # if data is vcf, convert to pfiles
            # if os.path.isfile(f'{geno_path}.vcf'):
            #     vfiles_to_pfiles(vcf_path=geno_path)
            # if data is in pfiles, convert to bfiles for KING
            if os.path.isfile(f'{geno_path}.pgen'):
                bfiles_to_pfiles(pfile_path=geno_path)

        print('If data does NOT contain pedigree (PAT/MAT) info, the Error column will be 1.0 for any found relationships')

        # run KING
        king_cmd = f'{king_exec} -b {geno_path}.bed --related --build --degree 3 --prefix {temp}'
        shell_do(king_cmd)

        if (os.path.isfile(f'{temp}.kin0')) or (os.path.isfile(f'{temp}.kin')):
            # pairs with different FID but found to be related
            if os.path.isfile(f'{temp}.kin0'):
                os.rename(f'{temp}.kin0', diff_fid_related)
                num_diff_fid_related = sum(1 for _ in open(diff_fid_related)) - 1

            # pairs with same FID but found to be unrelated
            if os.path.isfile(f'{temp}.kin'):
                kin = pd.read_csv(f'{temp}.kin', sep='\s+')
                unrelated = kin[kin['Kinship']<=third_deg_cutoff]
                unrelated.to_csv(same_fid_unrelated, sep='\t', header=True, index=False)
                num_same_fid_unrelated = sum(1 for _ in open(same_fid_unrelated)) - 1

            # create log file
            with open(f'{out_path}.log', 'w') as log:
                log.write('KING relatedness options in effect:\n')
                log.write('  --related\n')
                log.write('  --build\n')
                log.write('  --degree 3\n')
                log.write(f'  --prefix {temp}\n')
                log.write(f'  -b {geno_path}.bed\n')
                log.write('Note: if no pedigree info given, Error column will be 1.0 for any found relationship\n')
                log.write(f'{num_same_fid_unrelated + num_diff_fid_related} total noncongruent relationships found\n')
                # log.write('Potential pedigree is as follows:\n')
                # with open(f'{temp}build.log', 'r') as build:
                #     for line in build:
                #         if line != '\n':
                #         log.write(line)

            shutil.copyfile(f'{out_path}.log', f'{out_path}_test.log')

            if os.path.isfile(f'{out_path}.log'):
                listOfFiles = [f'{out_path}.log']
                concat_logs(step, out_path, listOfFiles)

            # remove intermediate files
            os.remove(f'{temp}.kin')
            os.remove(f'{temp}.seg')
            os.remove(f'{temp}X.kin')
            os.remove(f'{temp}X.kin0')
            os.remove(f'{temp}.segments.gz')
            os.remove(f'{temp}allsegs.txt')
            os.remove(f'{temp}build.log')
            os.remove(f'{temp}updateids.txt')
            os.remove(f'{temp}updateparents.txt')

            process_complete = True

            outfiles_dict = {
                'same_fid_unrelated': same_fid_unrelated,
                'diff_fid_related': diff_fid_related
            }

            metrics_dict = {
                'same_fid_unrelated_count': num_same_fid_unrelated,
                'diff_fid_related_count': num_diff_fid_related
            }

        else:
            print('Relatedness Assessment has failed!')
            process_complete = False

            outfiles_dict = {
                'same_fid_unrelated': None,
                'diff_fid_related': None
            }

            metrics_dict = {
                'same_fid_unrelated_count': 0,
                'diff_fid_related_count': 0
            }

        out_dict = {
            'pass': process_complete,
            'step': step,
            'metrics': metrics_dict,
            'output': outfiles_dict
        }

        return out_dict


class VariantQC:

    def __init__(self, geno_path=None, out_path=None):
        self.geno_path = geno_path
        self.out_path = out_path

    def run_geno_prune(self, geno_threshold=0.05):

        """
        Prunes SNPs based on missing genotype data.

        Parameters:
        - geno_threshold (float, optional): Maximum allowable proportion of missing genotyping for a SNP. Default is 0.05 (5%).

        Returns:
        - dict: A structured dictionary containing:
            * 'pass': Boolean indicating the successful completion of the process.
            * 'step': The label for this procedure ('geno_prune').
            * 'metrics': Metrics associated with the pruning, such as 'geno_removed_count'.
            * 'output': Dictionary containing paths to the generated output files.
        """

        geno_path = self.geno_path
        out_path = self.out_path

        # Check that paths are set
        if geno_path is None or out_path is None:
            raise ValueError("Both geno_path and out_path must be set before calling this method.")

        # Check path validity
        if not os.path.exists(f'{geno_path}.pgen'):
            raise FileNotFoundError(f"{geno_path} does not exist.")

        # Check type of geno_treshold
        if not isinstance(geno_threshold, float):
            raise TypeError("geno_threshold should be of type int or float.")

        # Check valid range for geno_threshold
        if geno_threshold < 0 or geno_threshold > 1:
            raise ValueError("geno_threshold should be between 0 and 1.")

        step = "geno_prune"

        # get initial snp count
        initial_snp_count = count_file_lines(f'{geno_path}.pvar') - 1

        # variant missingness
        plink_cmd1 = f"{plink2_exec} --pfile {geno_path} --geno 0.05 --make-pgen psam-cols=fid,parents,sex,pheno1,phenos --out {out_path}"
        shell_do(plink_cmd1)

        listOfFiles = [f'{out_path}.log']
        concat_logs(step, out_path, listOfFiles)

        # geno pruned count
        geno_snp_count = count_file_lines(f'{out_path}.pvar') - 1
        geno_rm_count = initial_snp_count - geno_snp_count

        process_complete = True

        outfiles_dict = {
            'plink_out': out_path
        }

        metrics_dict = {
            'geno_removed_count': geno_rm_count
        }

        out_dict = {
            'pass': process_complete,
            'step': step,
            'metrics': metrics_dict,
            'output': outfiles_dict
        }

        return out_dict


    def run_case_control_prune(self, p_threshold=1e-4):

        """
        Prunes SNPs based on missing genotype data, stratified by case-control status.

        Parameters:
        - p_threshold (float, optional): P-value threshold for testing the difference in missing rates between cases and controls. Default is 1e-4.

        Returns:
        - dict: A structured dictionary containing:
            * 'pass': Boolean indicating the successful completion of the process.
            * 'step': The label for this procedure ('case_control_missingness_prune').
            * 'metrics': Metrics associated with the pruning, such as 'mis_removed_count'.
            * 'output': Dictionary containing paths to the generated output files.
        """

        geno_path = self.geno_path
        out_path = self.out_path

        # Check that paths are set
        if geno_path is None or out_path is None:
            raise ValueError("Both geno_path and out_path must be set before calling this method.")

        # Check path validity
        if not os.path.exists(f'{geno_path}.pgen'):
            raise FileNotFoundError(f"{geno_path} does not exist.")

        # Check type of p_threshold
        if not isinstance(p_threshold, float):
            raise TypeError("p_threshold should be of type int or float.")

        # Check valid range for p_threshold
        if p_threshold < 0 or p_threshold > 1:
            raise ValueError("p_threshold should be between 0 and 1.")

        step = "case_control_missingness_prune"
        mis_tmp = f'{out_path}_mis_tmp'

        # convert to bfiles
        bfiles_to_pfiles(pfile_path=geno_path)

        # get initial snp count
        initial_snp_count = count_file_lines(f'{geno_path}.bim')

        fam = pd.read_csv(f'{geno_path}.fam', sep='\s+', header=None, usecols=[5], names=['case'])
        # check if this contains both cases and controls
        if  all(x in fam['case'].unique() for x in [1, 2]):
            # missingness by case control (--test-missing), using P > 1E-4
            plink_cmd1 = f"{plink_exec} --bfile {geno_path} --test-missing --out {mis_tmp}"
            shell_do(plink_cmd1)

            listOfFiles = [f'{mis_tmp}.log']
            concat_logs(step, out_path, listOfFiles)

            mis1 = f'{mis_tmp}.missing'
            if os.path.isfile(mis1):

                mis = pd.read_csv(mis1, sep='\s+')
                exclude = mis[mis.P <= p_threshold].loc[:,'SNP']
                exclude.to_csv(f'{mis_tmp}.exclude', sep='\t', header=False, index=False)

                plink_cmd2 = f"{plink2_exec} --bfile {geno_path} --exclude {mis_tmp}.exclude --make-pgen psam-cols=fid,parents,sex,pheno1,phenos --out {out_path}"
                shell_do(plink_cmd2)

                listOfFiles = [f'{out_path}.log']
                concat_logs(step, out_path, listOfFiles)

                # mis pruned count
                mis_snp_count = count_file_lines(f'{out_path}.pvar') - 1
                mis_rm_count = initial_snp_count - mis_snp_count

                process_complete = True

                for file in [f'{mis_tmp}.exclude', f'{mis_tmp}.hh', f'{mis_tmp}.missing']:
                    if os.path.isfile(file):
                        os.remove(file)

            else:
                print(f'Case/Control Missingness pruning failed!')
                print(f'Check {mis_tmp}.log for more information')

                process_complete = False
                mis_snp_count = count_file_lines(f'{geno_path}.bim')
                mis_rm_count = 0
        else:
            print(f'Case/control pruning failed! May be missing Controls!')
            print(f'Check number of Cases/Controls')

            process_complete = False

            mis_rm_count = 0

        os.remove(f'{geno_path}.bed')
        os.remove(f'{geno_path}.bim')
        os.remove(f'{geno_path}.fam')

        outfiles_dict = {
            'plink_out': out_path
        }

        metrics_dict = {
            'mis_removed_count': mis_rm_count
        }

        out_dict = {
            'pass': process_complete,
            'step': step,
            'metrics': metrics_dict,
            'output': outfiles_dict
        }

        return out_dict


    def run_haplotype_prune(self, p_threshold=1e-4):

        """
        Prunes SNPs based on missing haplotype data.

        Parameters:
        - p_threshold (float, optional): P-value threshold for missingness by haplotype. Default is 1e-4.

        Returns:
        - dict: A structured dictionary containing:
            * 'pass': Boolean indicating the successful completion of the process.
            * 'step': The label for this procedure ('haplotype_prune').
            * 'metrics': Metrics associated with the pruning, such as 'haplotype_removed_count'.
            * 'output': Dictionary containing paths to the generated output files.
        """

        geno_path = self.geno_path
        out_path = self.out_path

        # Check that paths are set
        if geno_path is None or out_path is None:
            raise ValueError("Both geno_path and out_path must be set before calling this method.")

        # Check path validity
        if not os.path.exists(f'{geno_path}.pgen'):
            raise FileNotFoundError(f"{geno_path} does not exist.")

        # Check type of p_threshold
        if not isinstance(p_threshold, float):
            raise TypeError("p_threshold should be of type int or float.")

        # Check valid range for p_threshold
        if p_threshold < 0 or p_threshold > 1:
            raise ValueError("p_threshold should be between 0 and 1.")

        step = "haplotype_prune"

        # make tmp names
        # missingness by haplotype
        hap_tmp = f'{out_path}_hap_tmp'

        # convert to bfiles
        bfiles_to_pfiles(pfile_path=geno_path)

        # get initial snp count
        initial_snp_count = count_file_lines(f'{geno_path}.bim')


        # if sample size is over 10k, correct and make P more stringent
        sample_size = count_file_lines(f'{geno_path}.fam')
        if sample_size > 10000:
            p_threshold = 0.05/sample_size


        # missingness by haplotype (--test-mishap), using P > 1E-4
        plink_cmd1 = f"{plink_exec} --bfile {geno_path} --maf 0.05 --test-mishap --out {hap_tmp}"
        shell_do(plink_cmd1)

        listOfFiles = [f'{hap_tmp}.log']
        concat_logs(step, out_path, listOfFiles)

        if os.path.isfile(f'{hap_tmp}.missing.hap'):

            # read .missing.hap file and grab flanking snps for P <= 0.0001. write flanking snps to file to exclude w bash
            mis_hap = pd.read_csv(f'{hap_tmp}.missing.hap', sep='\s+')
            mis_hap_snps = list(mis_hap[mis_hap.P <= p_threshold].loc[:,'FLANKING'].str.split('|'))
            snp_ls_df = pd.DataFrame({'snp':[rsid for ls in mis_hap_snps for rsid in ls]})
            snp_ls_df['snp'].to_csv(f'{hap_tmp}.exclude',sep='\t', header=False, index=False)

            plink_cmd2 = f"{plink2_exec} --bfile {geno_path} --exclude {hap_tmp}.exclude --make-pgen psam-cols=fid,parents,sex,pheno1,phenos --out {out_path}"
            shell_do(plink_cmd2)

            listOfFiles = [f'{out_path}.log']
            concat_logs(step, out_path, listOfFiles)

            # hap pruned count
            hap_snp_count = count_file_lines(f'{out_path}.pvar') - 1
            hap_rm_count = initial_snp_count - hap_snp_count

        else:
            print(f'Haplotype pruning failed!')
            print(f'Check {hap_tmp}.log for more information')
            process_complete = False
            hap_snp_count = count_file_lines(f'{geno_path}.bim')
            hap_rm_count = 0

        outfiles_dict = {
            'plink_out': out_path
        }

        metrics_dict = {
            'haplotype_removed_count': hap_rm_count
        }

        process_complete = True

        for file in [f'{hap_tmp}.hh']:
            if os.path.isfile(file):
                os.remove(file)

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


    def run_hwe_prune(self, hwe_threshold=1e-4, filter_controls=False):

        """
        Prunes SNPs based on Hardy-Weinberg equilibrium (HWE) p-values.

        Parameters:
        - hwe_threshold (float, optional): P-value threshold for SNP filtering based on HWE. Default is 1e-4.
        - filter_controls (bool, optional): If set to True, perform HWE test only on control samples. Default is False.

        Returns:
        - dict: A structured dictionary containing:
            * 'pass': Boolean indicating the successful completion of the process.
            * 'step': The label for this procedure ('hwe_prune').
            * 'metrics': Metrics associated with the pruning, such as 'hwe_removed_count'.
            * 'output': Dictionary containing paths to the generated output files.
        """

        geno_path = self.geno_path
        out_path = self.out_path

        # Check that paths are set
        if geno_path is None or out_path is None:
            raise ValueError("Both geno_path and out_path must be set before calling this method.")

        # Check path validity
        if not os.path.exists(f'{geno_path}.pgen'):
            raise FileNotFoundError(f"{geno_path} does not exist.")

        # Check type of hwe_threshold
        if not isinstance(hwe_threshold, float):
            raise TypeError("p_threshold should be of type int or float.")

        # Check type of filter_controls
        if not isinstance(filter_controls, bool):
            raise TypeError("filter_controls should be of type boolean.")

        step = "hwe_prune"

        hwe_tmp = f'{out_path}_hwe_tmp'

        # convert to bfiles
        bfiles_to_pfiles(pfile_path=geno_path)

        # get initial snp count
        initial_snp_count = count_file_lines(f'{geno_path}.bim')

        if filter_controls:
            # HWE from controls only using P > 1E-4
            plink_cmd1 = f"{plink_exec} --bfile {geno_path} --filter-controls --hwe {hwe_threshold} --write-snplist --out {hwe_tmp}"
        else:
            # HWE using P > 1E-4
            plink_cmd1 = f"{plink_exec} --bfile {geno_path} --hwe {hwe_threshold} --write-snplist --out {hwe_tmp}"

        plink_cmd2 = f"{plink2_exec} --bfile {geno_path} --extract {hwe_tmp}.snplist --make-pgen psam-cols=fid,parents,sex,pheno1,phenos --out {out_path}"

        cmds = [plink_cmd1, plink_cmd2]
        for cmd in cmds:
            shell_do(cmd)

        listOfFiles = [f'{hwe_tmp}.log', f'{out_path}.log']
        concat_logs(step, out_path, listOfFiles)

        if os.path.isfile(f'{out_path}.pvar'):
            process_complete = True
            # hwe pruned count
            final_snp_count = count_file_lines(f'{out_path}.pvar') - 1
            hwe_rm_count = initial_snp_count - final_snp_count
        else:
            print(f'HWE pruning failed!')
            print(f'Check {hwe_tmp}.log or {out_path}.log for more information')
            process_complete = False
            hwe_rm_count = 0

        outfiles_dict = {
            'plink_out': out_path
        }

        metrics_dict = {
            'hwe_removed_count': hwe_rm_count
        }

        for file in [f'{hwe_tmp}.hh',f'{hwe_tmp}.snplist']:
            if os.path.isfile(file):
                os.remove(file)

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


    def run_ld_prune(self, window_size=50, step_size=5, r2_threshold=0.5):

        """
        Prunes SNPs based on Linkage Disequilibrium

        Parameters:
        - window_size (int, optional): Size of sliding window over the genome in variant count.
        - step_size (int, optional): Step of sliding window over the genome.
        - r2_threshold (float, optional): r^2 threshold to include a variant.

        Returns:
        - dict: A structured dictionary containing:
            * 'pass': Boolean indicating the successful completion of the process.
            * 'step': The label for this procedure ('ld_prune').
            * 'metrics': Metrics associated with the pruning, such as 'ld_removed_count'.
            * 'output': Dictionary containing paths to the generated output files.
        """

        geno_path = self.geno_path
        out_path = self.out_path

        # Check that paths are set
        if geno_path is None or out_path is None:
            raise ValueError("Both geno_path and out_path must be set before calling this method.")

        # Check path validity
        if not os.path.exists(f'{geno_path}.pgen'):
            raise FileNotFoundError(f"{geno_path} does not exist.")

        # Check type of window_size
        if not isinstance(window_size, int):
            raise TypeError("window_size should be of type int.")

        # Check type of step_size
        if not isinstance(step_size, int):
            raise TypeError("step_size should be of type int.")

        # Check type of r2_threshold
        if not isinstance(r2_threshold, float):
            raise ValueError("r2_threshold should be of type float.")

        step = "ld_prune"

        # temp file
        ld_temp = f'{out_path}_ld_temp'

        # get initial snp count
        initial_snp_count = count_file_lines(f'{geno_path}.pvar') - 1

        # get list of SNPs to be extracted
        plink_cmd1 = f"{plink2_exec} --pfile {geno_path} --indep-pairwise {window_size} {step_size} {r2_threshold} --out {ld_temp}"
        # and extract
        plink_cmd2 = f"{plink2_exec} --pfile {geno_path} --extract {ld_temp}.prune.in --make-pgen psam-cols=fid,parents,sex,pheno1,phenos --out {out_path}"

        cmds = [plink_cmd1, plink_cmd2]
        for cmd in cmds:
            shell_do(cmd)

        listOfFiles = [f'{ld_temp}.log', f'{out_path}.log']
        concat_logs(step, out_path, listOfFiles)

        if os.path.isfile(f'{out_path}.pvar'):
            # ld pruned count
            ld_snp_count = count_file_lines(f'{out_path}.pvar') - 1
            ld_rm_count = initial_snp_count - ld_snp_count
            os.remove(f'{ld_temp}.prune.in')
            os.remove(f'{ld_temp}.prune.out')
            process_complete = True

        else:
            print(f'LD pruning failed!')
            print(f'Check {ld_temp}.log or {out_path}.log for more information')
            process_complete = False
            ld_rm_count = 0

        # process_complete = True

        # os.remove(f'{ld_temp}.prune.in')
        # os.remove(f'{ld_temp}.prune.out')

        outfiles_dict = {
            'plink_out': out_path
        }

        metrics_dict = {
            'ld_removed_count': ld_rm_count
        }

        out_dict = {
            'pass': process_complete,
            'step': step,
            'metrics': metrics_dict,
            'output': outfiles_dict
        }

        return out_dict
