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
    def __init__(self, geno_path=None, out_path=None, pca=10, build='hg38', gwas=True, covar_path=None, covar_names=None):
        self.geno_path = geno_path
        self.out_path = out_path
        self.pca = pca
        self.build = build
        self.gwas = gwas
        self.covar_path = covar_path
        self.covar_names = covar_names
    

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
        
        exclusion_file = f'{self.out_path}_{self.build}.txt'
        with open(exclusion_file, 'w') as f:
            f.write(exclusion_regions)
        
        return exclusion_file
    
    
    def run_pca_pruning(self, maf=0.01, geno=0.01, hwe=5e-6, indep_pairwise=[1000, 10, 0.02]):
        
        step = 'plink_pruning'

        exclusion_file = self.write_exclusion_file()

        # Filter data
        filter_cmd = f"{plink2_exec} --pfile {self.geno_path} --maf {maf} --geno {geno} --hwe {hwe} --autosome --exclude {exclusion_file} --make-pgen psam-cols=fid,parents,sex,phenos --out {self.out_path}_tmp"
        # shell_do(filter_cmd)
        
        # Prune SNPs
        prune_cmd = f"{plink2_exec} --pfile {self.out_path}_tmp --indep-pairwise {indep_pairwise[0]} {indep_pairwise[1]} {indep_pairwise[2]} --autosome --out {self.out_path}_pruned"
        # shell_do(prune_cmd)

        # Check if prune.in file exists
        if os.path.isfile(f'{self.out_path}_pruned.prune.in'):
            # Extract pruned SNPs
            extract_cmd = f"{plink2_exec} --pfile {self.out_path}_tmp --extract {self.out_path}_pruned.prune.in --make-pgen psam-cols=fid,parents,sex,phenos --out {self.out_path}"
            shell_do(extract_cmd)

            listOfFiles = [f'{self.out_path}_tmp.log', f'{self.out_path}_pruned.log']
            print(self.out_path)
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


    def calculate_inflation(self, pval_array, normalize=False, ncases=None, ncontrols=None):

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
        psam[['#FID', 'IID', 'PHENO1']].to_csv(f'{self.geno_path}_pheno.txt', index=False)

        if covars:
            # covar names are column names of covariate file minus #FID and IID unless specified
            if self.covar_path and not self.covar_names:
                covar_cols = pd.read_csv(f'{self.covar_path}', sep='\s+', header=0, nrows=0).columns.tolist()
                # account for 1st two cols being #FID and IID
                covar_names = ','.join(covar_cols[2:])
            else:
                # if no covar path then default is PCA
                covar_names = self.covar_names

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
            
            gwas_cmd = (
                f"{plink2_exec} --pfile {self.geno_path} "
                f"--glm {glm_options} "
                f"--pheno-name PHENO1 "
                f"--pheno {self.geno_path}_pheno.pheno "
                f"--covar {self.covar_path} "
                f"--covar-variance-standardize "
                f"--covar-name {covar_names} "
                f"--out {self.out_path}"
                )
        
        else:
            gwas_cmd = (
                f"{plink2_exec} --pfile {self.geno_path} "
                f"--glm {glm_options} "
                f"--pheno-name PHENO1 "
                f"--pheno {self.geno_path}_pheno.pheno "
                f"--out {self.out_path}"
                )

            gwas_cmd = f'{plink2_exec} --pfile {self.geno_path} \
                        --glm allow-no-covars firth-fallback no-x-sex cols=+a1freq,+a1freqcc,+a1count,+totallele,+a1countcc,+totallelecc,+err \
                        --pheno-name PHENO1 \
                        --pheno {self.geno_path}_pheno.txt \
                        --out {self.out_path}'
            
        if os.path.isfile(f'{self.out_path}.PHENO1.glm.logistic.hybrid'):
            process_complete = True

            outfiles_dict = {
                'gwas_output': f'{self.out_path}.PHENO1.glm.logistic.hybrid'
            }
        
        else:
            print(f'Heterozygosity pruning failed!')
            print(f'Check {self.out_path}.log for more information')

            process_complete = False

            outfiles_dict = {
                'gwas_output': None
            }

        # will need to calculate inflation using outputs of GWAS
        shell_do(gwas_cmd)

        process_complete = True

        out_dict = {
            'pass': process_complete,
            'step': step,
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


        # calculate pcs
        if self.pca:
            pca = self.run_plink_pca()

        # run gwas
        if self.gwas:
            # if pca called and covars passed, warn and default to passed covars
            if os.path.isfile(f'{self.out_path}.eigenvec') and (self.covar_path is not None):
                warnings.warn('PCA ran and Covar passed! Defaulting to passed covars!')
                gwas = self.run_gwas(covars=True)
            
            # if pca called and no covars pass, defult to pca
            elif os.path.isfile(f'{self.out_path}.eigenvec') and (self.covar_path is None):
                self.covar_path = f'{self.out_path}.eigenvec'
                gwas = self.run_gwas(covars=True)
            
            # otherwise run GWAS with no covars
            else:
                gwas = self.run_gwas(covars=False)
