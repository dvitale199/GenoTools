# GenoTools Reference Panel Preparation

## Overview
This documentation provides a detailed description of how to prepare a reference panel for use in the `GenoTools` ancestry module. Preparing a custom reference panel will allow for the prediction of samples within ancestry groups outside of the 11 provided by the default reference panel.

---

### Format
In order to run the proper pruning steps, as well as use the custom reference panel in the ancestry module, please make sure the reference genotypes are in PLINK1.9 binary format (.bed/.bim/.fam): https://www.cog-genomics.org/plink/1.9/formats

---

### Pruning Script
Please run the following script in order to prune the reference panel for use in the `GenoTools` ancestry model:

```python
import subprocess
import sys
import os
import shutil
import pandas as pd
import numpy as np

from genotools.dependencies import check_plink, check_plink2
from genotools.utils import shell_do

def prune_geno(ref_path, ld_path):
    prune_path = f'{ref_path}_maf_geno_hwe'
    final_path = f'{ref_path}_pruned'
    
    # read bim and get palindromes
    bim = pd.read_csv(f'{ref_path}.bim', sep='\s+', header=None)
    bim.columns = ['chr','id','kb','pos','alt','ref']
    palindromes = bim.loc[((bim.alt == 'A') & (bim.ref == 'T')) | ((bim.alt == 'T') & (bim.ref == 'A')) | ((bim.alt == 'C') & (bim.ref == 'G')) | ((bim.alt == 'G') & (bim.ref == 'C'))]
    palindromes['id'].to_csv(f'{ref_path}_palindromes.snplist', header=False, index=False, sep='\t')
    
        
    # maf, geno, hwe command
    plink_cmd1 = f'{plink_exec} --bfile {ref_path}\
        --maf 0.05\
        --geno 0.01\
        --hwe 0.0001\
        --indep-pairwise 1000 10 0.02\
        --exclude {geno_path}_palindromes.snplist\
        --autosome\
        --allow-extra-chr\
        --make-bed\
        --out {prune_path}'

    shell_do(plink_cmd1)

    # ld command
    plink_cmd2 = f'{plink2_exec} --bfile {prune_path}\
        --exclude range {ld_path}\
        --autosome\
        --allow-extra-chr\
        --make-bed\
        --out {final_path}'
    
    shell_do(plink2_cmd2)

if __name__ == '__main__':
    ref_path = '/path/to/reference/reference/genotypes'
    ld_path = '/path/to/ld/exclusion/regions'

    prune_geno(ref_path, ld_path)
```

---

### Pruning Steps
1. Find palindrome SNPs from .bim file, and write to ```python
{ref_path}_palindromes.snplist```