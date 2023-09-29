import pandas as pd
import argparse
import shutil

# local imports
from QC.imputation import *
# from QC.qc import callrate_prune, het_prune, sex_prune, related_prune, variant_prune, avg_miss_rates
# from Ancestry.ancestry import run_ancestry, split_cohort_ancestry
# from QC.utils import shell_do


parser = argparse.ArgumentParser(description='Arguments for Genotyping Imputation submission (data in Plink .bim/.bam/.fam format)')

parser.add_argument('--geno', type=str, default='nope', help='Genotype: (string file path). Path to PLINK format genotype file, everything before the *.bed/bim/fam [default: nope].')
# parser.add_argument('--temp', type=str, default='nope', help='Prefix for output (including path)')
parser.add_argument('--out', type=str, default='nope', help='Full path to directory for imputation server output')
parser.add_argument('--token', type=str, default='nope', help='API token for TOPMed Imputation server')
# parser.add_argument('--ref_panel', type=str, default='nope', help='Imputation Server Reference Panel (more information: https://www.well.ox.ac.uk/~wrayner/tools/)')
# parser.add_argument('--check_bim_pl', type=str, default='nope', help='Imputation Server HRC-1000G-check-bim.pl script (more information: https://www.well.ox.ac.uk/~wrayner/tools/)')

# *****FUTURE FIX******
# having issues running imputation_data_prep() in swarm job. for now, take geno_path as output prefix from imputation_data_prep() and create vcf_list
# geno_path is everything before "_pre_impute_chr{str(i)}.vcf.gz"

args = parser.parse_args()

geno_path = args.geno
# temp_path = args.temp
token = args.token
# ref_panel = args.ref_panel
# check_bim_pl = args.check_bim_pl
out_path = args.out

vcf_list = [f'{geno_path}_pre_impute_chr{str(i)}.vcf.gz' for i in range(1,24)]

# impute = run_auto_imputation(geno_path=geno_path, temp_path=temp_path, out_path=out_path, token=token, ref_panel=ref_panel, check_bim_pl=check_bim_pl, password='imputer')
impute = run_auto_imputation(vcf_list=vcf_list, out_path=out_path, token=token, password='imputer')

print(impute)