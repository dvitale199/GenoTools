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

    geno_path2 = geno_path.replace(out_path, '')
    GWAS_workdir = os.getcwd()
    step1 = "PREP PLINK FILES FOR IMPUTATION"

    # download file to check
    bash1 = "wget http://www.well.ox.ac.uk/~wrayner/tools/HRC-1000G-check-bim.v4.2.5.zip -P " + out_path
    bash2 = "unzip " + out_path + "HRC-1000G-check-bim.v4.2.5.zip -d " + out_path
    bash3 = "wget ftp://ngs.sanger.ac.uk/production/hrc/HRC.r1-1/HRC.r1-1.GRCh37.wgs.mac5.sites.tab.gz -P " + out_path
    bash4 = "gunzip " + out_path + "HRC.r1-1.GRCh37.wgs.mac5.sites.tab.gz"
    # make .frq file
    bash5 = "plink --bfile " + geno_path + " --freq --out " + geno_path
    
    cmds1 = [bash1, bash2, bash3, bash4, bash5]
    #run first set of commands
    for cmd in cmds1:
        shell_do(cmd)
    
    bash6 = "perl HRC-1000G-check-bim.pl -b " + geno_path2 + ".bim -f " + geno_path2 + ".frq -r HRC.r1-1.GRCh37.wgs.mac5.sites.tab -h"
    # then run to fix data
    bash7 = "sh Run-plink.sh"
    
    driver2 = Driver(geno_path2)
    #change directory for second set of commands
    os.chdir(out_path)
    
    cmds2 = [bash6, bash7]
    step2 = "Run HRC-1000G-check-bim.pl and Run-plink.sh to fix data for HRC"
    driver2.run_cmds(cmds2, step2)
    
    # now change back to working dir
    os.chdir(GWAS_workdir)