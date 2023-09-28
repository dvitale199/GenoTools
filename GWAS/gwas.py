import pandas as pd
import numpy as np
from scipy.stats import ncx2
from QC.utils import shell_do
from utils.dependencies import check_plink, check_plink2


plink_exec = check_plink()
plink2_exec = check_plink2()

class GWAS:
    def __init__(self, geno_path, out_path, covar_path = None, covar_names = None):
        self.geno_path = geno_path
        self.out_path = out_path
        self.covar_path = covar_path
        self.covar_names = covar_names

    def calculate_inflation(self, pval_array, normalize=False, ncases=None, ncontrols=None):
        # what step are we running?
        step = 'lambda_calculation'
        print()
        print(f'RUNNING: {step}')
        print()

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

    def gwas(self):
        # what step are we running?
        step = 'GWAS'
        print()
        print(f'RUNNING: {step}')
        print()

        # make pheno file from .psam
        psam = pd.read_csv(f'{self.geno_path}.psam', sep='\s+')
        psam[['#FID', 'IID', 'PHENO1']].to_csv(f'{self.geno_path}_pheno.txt', index=False)

        # covar names are column names of covariate file minus #FID and IID unless specified
        if self.covar_path and not self.covar_names:
            covar_cols = pd.read_csv(f'{self.covar_path}', sep='\s+', header=0, nrows=0).columns.tolist()
            # account for 1st two cols being #FID and IID
            covar_names = ','.join(covar_cols[2:])
        else:
            # if no covar path then default is PCA ?
            covar_names = self.covar_names

        # adjust input to be PLINK2 pgen
        gwas_cmd = f'{plink2_exec} --pfile {self.geno_path} \
                    --glm hide-covar firth-fallback no-x-sex cols=+a1freq,+a1freqcc,+a1count,+totallele,+a1countcc,+totallelecc,+err \
                    --pheno-name PHENO1 \
                    --pheno {self.geno_path}_pheno.txt \
                    --covar {self.covar_path} \
                    --covar-variance-standardize \
                    --covar-name {covar_names} \
                    --out {self.out_path}'

        # will need to calculate inflation using outputs of GWAS
        shell_do(gwas_cmd)