# GenoTools Reference Panel Preparation

## Overview
This documentation provides a detailed description of how to prepare a reference panel for use in the `GenoTools` ancestry module. Preparing a custom reference panel will allow for the prediction of samples within ancestry groups outside of the 11 included in the provided reference panel.

---

### Format
In order to run the proper pruning steps, as well as use the custom reference panel in the ancestry module, please make sure the reference genotypes are in PLINK1.9 binary format (.bed/.bim/.fam): https://www.cog-genomics.org/plink/1.9/formats. In addition, you will need to create a tab-delimited file containing the sample IDs as well as the reference ancestry groups (see Preparing the `--ref_labels` File section).

---

### Pruning Script
Please run the following script in order to prune the custom reference panel for use in the `GenoTools` ancestry model:

```python
import subprocess
import sys
import os
import shutil
import pandas as pd
import numpy as np

from genotools.dependencies import check_plink, check_plink2
from genotools.utils import shell_do

def prune_ref_panel(ref_path, ld_path):
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
    ld_path = '/path/to/ld/exclusion_regions.txt'

    prune_ref_panel(ref_path, ld_path)
```

### Parameters
- **ref_path**: Path to PLINK 1.9 format reference panel genotype file (before the *.bed/bim/fam).
- **ld_path**: Path to exclusion regions file in tab-delimited .txt format (see LD Exclusion Region File Format section).

---

### Description of Pruning Steps
1. Find palindrome SNPs from .bim file, and write to ```{ref_path}_palindromes.snplist```.
2. Prune reference genotypes for ```--maf```, ```--geno```, ```--hwe```, and linkage disequilibrium, as well as removing the previously identified palindrome SNPs.
3. Exclude high-LD regions from reference genotypes.

---

### LD Exclusion Region File Format (hg38)
Copy this to a .txt file and pass to the ```prune_ref_panel()``` function in the provided Python script:
```
6	24999772	33532223 r1
6	26726947	26726981 r2
8	8142478	12142491 r3
20	28534739	28534795 r4
11	50809938	50809993 r5
2	94531009	94531062 r6
11	54525074	55026708 r7
12	34504916	34504974 r8
11	55104000	55104074 r9
```

---

### Preparing the `--ref_labels` File
In order to render predictions using the new reference panel, the `GenoTools` ancestry module needs information on which ancestry group each sample belongs to. The file passed to the `--ref_labels` flag when running the ancestry module must be tab-delimited, and contain three columns corresponding to `FID`, `IID`, and `label`. Please note this file should not have a header row and the IDs must match those in the `.fam` file for the process to run correctly. For more information on how to use the `--ref_labels` flag, please see the `train_new_model.md` documentation. Below is an example of how the `--ref_labels` file should look:
```
FID1    IID1    ANCESTRY1
FID2    IID2    ANCESTRY2
FID3    IID3    ANCESTRY2
FID4    IID4    ANCESTRY3
```
