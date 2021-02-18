import subprocess
import argparse
import pandas as pd
import os
import glob
import shutil
import sys

# local imports
from QC.utils import shell_do

def call_rate_pruning(geno, out):

    step = "PRUNING FOR CALL RATE"
    print(step)
    bash1 = f"awk '{{print $1,$2,$6}}' {geno}.fam > {geno}.phenos"
    bash2 = f"plink --bfile {geno} --mind 0.05 --make-bed --out {out}_call_rate"
    bash3 = f"mv {out}_call_rate.irem {out}_CALL_RATE_OUTLIERS.txt"

    cmds = [bash1, bash2, bash3]

    for cmd in cmds:
        shell_do(cmd)
    
    outliers_count = sum(1 for line in open(f'{out}_CALL_RATE_OUTLIERS.txt'))

    out_dict = {
        'step': step,
        'outliers_count': outliers_count,
        'outliers_path': f'{out}_CALL_RATE_OUTLIERS.txt',
        'geno_outpath': [f'{out}_call_rate.{suffix}' for suffix in {'bed','bim','fam'}],
        'phenos_path': f'{geno}.phenos'
    }

    return out_dict