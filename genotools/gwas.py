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


import os
import pandas as pd
import numpy as np
import warnings
from scipy.stats import ncx2
from genotools.utils import shell_do, concat_logs
from genotools.dependencies import check_plink, check_plink2


plink_exec = check_plink()
plink2_exec = check_plink2()

class Assoc:
    def __init__(self, geno_path=None, out_path=None, pca=10, build='hg38', gwas=True, pheno_name='PHENO1', covar_path=None, covar_names=None, maf_lambdas=False):
        self.geno_path = geno_path
        self.out_path = out_path
        self.pca = pca
        self.build = build
        self.gwas = gwas
        self.pheno_name = pheno_name
        self.covar_path = covar_path
        self.covar_names = covar_names
        self.maf_lambdas = maf_lambdas


    def write_exclusion_file(self):
        if self.build == 'hg19':
            exclusion_regions = """
        5 44000000 51500000 r1
        6 25000000 33500000 r2
        8 8000000 12000000 r3
        11 45000000 57000000 r4
        """
        elif self.build == 'hg38':
            exclusion_regions = """
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
        else:
            raise ValueError("Invalid build specified.")

        exclusion_file = f'{self.out_path}_{self.build}.exclusion'
        with open(exclusion_file, 'w') as f:
            f.write(exclusion_regions)

        return exclusion_file


    def run_pca_pruning(self, maf=0.01, geno=0.01, hwe=5e-6, indep_pairwise=[1000, 10, 0.02]):

        step = 'plink_pruning'

        exclusion_file = self.write_exclusion_file()

        # Filter data
        filter_cmd = f"{plink2_exec} --pfile {self.geno_path} --maf {maf} --geno {geno} --hwe {hwe} --autosome --exclude {exclusion_file} --make-pgen psam-cols=fid,parents,sex,pheno1,phenos --out {self.out_path}_tmp"
        shell_do(filter_cmd)

        # Prune SNPs
        prune_cmd = f"{plink2_exec} --pfile {self.out_path}_tmp --indep-pairwise {indep_pairwise[0]} {indep_pairwise[1]} {indep_pairwise[2]} --autosome --out {self.out_path}_pruned"
        shell_do(prune_cmd)

        # Check if prune.in file exists
        if os.path.isfile(f'{self.out_path}_pruned.prune.in'):
            # Extract pruned SNPs
            extract_cmd = f"{plink2_exec} --pfile {self.out_path}_tmp --extract {self.out_path}_pruned.prune.in --make-pgen psam-cols=fid,parents,sex,pheno1,phenos --out {self.out_path}"
            shell_do(extract_cmd)

            listOfFiles = [f'{self.out_path}_tmp.log', f'{self.out_path}_pruned.log']
            concat_logs(step, self.out_path, listOfFiles)

            process_complete = True
        else:
            print(f'Check {self.out_path}_pruned.log for more information.')
            print('Likely there are <50 samples for this ancestry leading to bad LD calculations.')
            print()

            process_complete = False

        # Remove intermediate files
        try:
            os.remove(f"{self.out_path}_tmp.pgen")
            os.remove(f"{self.out_path}_tmp.pvar")
            os.remove(f"{self.out_path}_tmp.psam")
            os.remove(f"{self.out_path}_pruned.prune.in")
            os.remove(f"{self.out_path}_pruned.prune.out")
            os.remove(exclusion_file)
        except OSError:
            pass


        out_dict={
            'pass':process_complete,
            'step': step,
            'output': f'{self.out_path}'
            }

        return out_dict


    def run_plink_pca(self):

        step = 'plink_pca'

        # Calculate/generate PCs
        pca_pruned = self.run_pca_pruning()
        pca_cmd = f"{plink2_exec} --pfile {pca_pruned['output']} --pca {self.pca} --out {self.out_path}"
        shell_do(pca_cmd)

        listOfFiles = [f'{self.out_path}.log']
        concat_logs(step, self.out_path, listOfFiles)

        if os.path.isfile(f'{self.out_path}.eigenvec'):
            process_complete = True

            outfiles_dict = {
                'eigenvec': f'{self.out_path}.eigenvec',
                'eigenval': f'{self.out_path}.eigenval'
            }
        else:
            print(f'PCA failed!')
            print(f'Check {self.out_path}.log for more information')

            process_complete = False

            outfiles_dict = {
                'eigenvec': None,
                'eigenval': None
            }

        out_dict = {
            'pass':process_complete,
            'step': step,
            'output': outfiles_dict
        }

        return out_dict

    @staticmethod
    def calculate_inflation(pval_array, normalize=False, ncases=None, ncontrols=None):

            step = 'lambda_calculation'

            # need cases and controls if normalizing
            if normalize and (not ncases or not ncontrols):
                print(f'Inflation Calculation failed!')
                print(f'If normalizing, please add ncases and ncontrols')

                process_complete = False
                inflation = 0
            else:
                # calculate inflation
                num = ncx2.ppf(1-pval_array, 1, nc=0)
                denom = ncx2.ppf(0.5, 1, nc=0)
                inflation = np.nanmedian(num) / denom

                # normalize when necessary
                if normalize:
                    inflation1000 = 1 + (inflation - 1) * (1/ncases + 1/ncontrols) / (1/1000 + 1/1000)
                    inflation = inflation1000

                process_complete = True

            metrics_dict = {
                'inflation': inflation
            }

            out_dict = {
                'pass': process_complete,
                'step': step,
                'metrics': metrics_dict
            }

            return out_dict


    def run_gwas(self, covars=True):

        step = 'GWAS'

        # make pheno file from .psam
        psam = pd.read_csv(f'{self.geno_path}.psam', sep='\s+')
        print('psam produced')

        # check if psam has #FID and IID
        if '#FID' in psam.columns:
            psam[['#FID', 'IID', 'PHENO1']].to_csv(f'{self.geno_path}.pheno', index=False)
        else:
            psam[['#IID', 'PHENO1']].to_csv(f'{self.geno_path}.pheno', index=False)

        glm_options = (
            'hide-covar '
            'firth-fallback '
            'no-x-sex '
            'cols=+a1freq,'
            '+a1freqcc,'
            '+a1count,'
            '+totallele,'
            '+a1countcc,'
            '+totallelecc,'
            '+err'
            )

        if covars:
            # covar names are column names of covariate file minus #FID and IID unless specified
            if self.covar_path and not self.covar_names:
                covar_cols = pd.read_csv(f'{self.covar_path}', sep='\s+', header=0, nrows=0).columns.tolist()
                if '#FID' in covar_cols:
                    # account for 1st two cols being #FID and IID
                    covar_names = ','.join(covar_cols[2:])
                else:
                    covar_names = ','.join(covar_cols[1:])
            else:
                # if no covar path then default is PCA
                covar_names = self.covar_names

            gwas_cmd = (
                f"{plink2_exec} --pfile {self.geno_path} "
                f"--glm {glm_options} "
                f"--pheno-name {self.pheno_name} "
                f"--pheno {self.geno_path}.pheno "
                f"--covar {self.covar_path} "
                f"--covar-variance-standardize "
                f"--covar-name {covar_names} "
                f"--out {self.out_path}"
                )

        else:
            gwas_cmd = (
                f"{plink2_exec} --pfile {self.geno_path} "
                f"--glm {glm_options} allow-no-covars "
                f"--pheno-name {self.pheno_name} "
                f"--pheno {self.geno_path}.pheno "
                f"--out {self.out_path}"
                )

        shell_do(gwas_cmd)

        if os.path.isfile(f'{self.out_path}.PHENO1.glm.logistic.hybrid'):
            print('logistic')

            # calculate inflation
            gwas_df = pd.read_csv(f'{self.out_path}.PHENO1.glm.logistic.hybrid', sep='\s+', dtype={'#CHROM': str})
            psam = pd.read_csv(f'{self.geno_path}.psam', sep='\s+', dtype={'PHENO1':str})
            ncontrols = psam.PHENO1.value_counts()['1']
            ncases = psam.PHENO1.value_counts()['2']

            # add pruning step here (pre lambdas)
            gwas_df_add = gwas_df.loc[gwas_df.TEST=='ADD']

            if self.maf_lambdas:
                gwas_df_maf = gwas_df_add[gwas_df_add['A1_FREQ']>=0.01] 
                gwas_df_maf = gwas_df_maf[gwas_df_maf['A1_FREQ']<=0.99]
                gwas_df_maf = gwas_df_maf[~gwas_df_maf['P'].isna()]

                lambda_maf_dict = Assoc.calculate_inflation(gwas_df_maf.P, normalize=False)
                lambda1000_maf_dict = Assoc.calculate_inflation(gwas_df_maf.P, normalize=True, ncases=ncases, ncontrols=ncontrols)

            # calculate inflation
            lambda_dict = Assoc.calculate_inflation(gwas_df_add.P, normalize=False)
            lambda1000_dict = Assoc.calculate_inflation(gwas_df_add.P, normalize=True, ncases=ncases, ncontrols=ncontrols)

            if self.maf_lambdas:
                metrics_dict = {
                    'lambda': lambda_dict['metrics']['inflation'],
                    'lambda1000': lambda1000_dict['metrics']['inflation'],
                    'lambda_maf': lambda_maf_dict['metrics']['inflation'],
                    'lambda1000_maf': lambda1000_maf_dict['metrics']['inflation'],
                    'cases': ncases,
                    'controls': ncontrols
                }
            else:
                metrics_dict = {
                    'lambda': lambda_dict['metrics']['inflation'],
                    'lambda1000': lambda1000_dict['metrics']['inflation'],
                    'cases': ncases,
                    'controls': ncontrols
                    }

            process_complete = True

            outfiles_dict = {
                'gwas_output': f'{self.out_path}.PHENO1.glm.logistic.hybrid'
            }

        elif os.path.isfile(f'{self.out_path}.PHENO1.glm.linear'):
            print('linear')

            # calculate inflation
            gwas_df = pd.read_csv(f'{self.out_path}.PHENO1.glm.linear', sep='\s+', dtype={'#CHROM': str})

            # add pruning step here (pre lambdas)
            gwas_df_add = gwas_df.loc[gwas_df.TEST=='ADD']

            # calculate inflation
            lambda_dict = self.calculate_inflation(gwas_df_add.P, normalize=False)

            metrics_dict = {
                'lambda': lambda_dict['metrics']['inflation'],
                'lambda1000': np.nan,
                'cases': np.nan,
                'controls': np.nan
                }
            process_complete = True

            outfiles_dict = {
                'gwas_output': f'{self.out_path}.PHENO1.glm.linear'
            }

        else:
            print(f'Association failed!')
            print(f'Check {self.out_path}.log for more information')

            process_complete = False

            metrics_dict = {
                'lambda': None,
                'lambda1000': None,
                'cases': None,
                'controls': None
                }

            outfiles_dict = {
                'gwas_output': None
            }

        out_dict = {
            'pass': process_complete,
            'step': step,
            'metrics': metrics_dict,
            'output': outfiles_dict
        }

        return out_dict


    def run_prs(self):
        print('COMING SOON')


    def run_association(self):
        # Check that paths are set
        if self.geno_path is None or self.out_path is None:
            raise ValueError("Both geno_path and out_path must be set before calling this method.")

        # Check path validity
        if not os.path.exists(f'{self.geno_path}.pgen'):
            raise FileNotFoundError(f"{self.geno_path} does not exist.")

        assoc = {}

        # calculate pcs
        if self.pca:
            assoc['pca'] = self.run_plink_pca()

        # run gwas
        if self.gwas:
            # if pca called and covars passed, warn and default to passed covars
            if os.path.isfile(f'{self.out_path}.eigenvec') and (self.covar_path is not None):
                warnings.warn('PCA ran and Covars passed! Defaulting to passed covars! Recommend merging the PCAs with the inputted Covars and rerunning.')
                assoc['gwas'] = self.run_gwas(covars=True)

            # if pca called and no covars pass, defult to pca
            elif os.path.isfile(f'{self.out_path}.eigenvec') and (self.covar_path is None):
                self.covar_path = f'{self.out_path}.eigenvec'
                assoc['gwas'] = self.run_gwas(covars=True)

            # otherwise run GWAS with no covars
            else:
                assoc['gwas'] = self.run_gwas(covars=False)
        
        # create an overall assoc pass flag
        if self.pca and self.gwas:
            assoc['pass'] = True if ((assoc['pca']['pass']) and (assoc['gwas']['pass'])) else False
        elif self.pca:
            assoc['pass'] = True if (assoc['pca']['pass']) else False
        else:
            assoc['pass'] = True if (assoc['gwas']['pass']) else False

        return assoc