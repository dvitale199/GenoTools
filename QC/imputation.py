import subprocess
import argparse
import os
import requests
import json
import time
import glob

# local imports
from QC.utils import shell_do, rm_tmps, count_file_lines

def impute_prep_data(geno_path, out_path):
    
    
    workdir = os.getcwd()
    os.chdir(out_path)

    plink1 = f'plink --bfile {geno_path} --freq --out {geno_path}'
    check_bim_cmd = f'perl {ref_path}/HRC-1000G-check-bim.pl -b {geno_path}.bim -f {geno_path}.frq -r {ref_path}/HRC.r1-1.GRCh37.wgs.mac5.sites.tab -h'
    bash1 = 'sh Run-plink.sh'

    cmds = [plink1, check_bim_cmd, bash1]

    for cmd in cmds:
        shell_do(cmd)

    os.chdir(workdir)
