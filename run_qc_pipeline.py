import pandas as pd
import argparse

# local imports
from QC.qc import callrate_prune, het_prune, sex_prune, related_prune, variant_prune, avg_miss_rates
from Ancestry.ancestry import run_ancestry, split_cohort_ancestry
from QC.utils import shell_do


parser = argparse.ArgumentParser(description='Arguments for Genotyping QC (data in Plink .bim/.bam/.fam format)')
parser.add_argument('--geno', type=str, default='nope', help='Genotype: (string file path). Path to PLINK format genotype file, everything before the *.bed/bim/fam [default: nope].')
parser.add_argument('--ref', type=str, default='nope', help='Genotype: (string file path). Path to PLINK format reference genotype file, everything before the *.bed/bim/fam.')
parser.add_argument('--ref_labels', type=str, default='nope', help='tab-separated plink-style IDs with ancestry label (FID  IID label) with no header')
parser.add_argument('--out', type=str, default='out', help='Prefix for output (including path)')
# parser.add_argument('--rare', default=False, action="store_true", help='Pruning toggle for rare variants. If --rare is used, final MAF pruning (0.01) will not be conducted, otherwise, rare variants will be pruned')

args = parser.parse_args()

geno_path = args.geno
ref_panel = args.ref
ref_labels = args.ref_labels
out = args.out

# geno_path = '/data/CARD/PD/GP2/raw_genotypes/shulman_ny/plink/shulman'
# hard code in reference for now--- definitely change later
# ref_dir_path = '/data/LNG/vitaled2/1kgenomes'
# ref_panel = f'{ref_dir_path}/1kg_ashkj_ref_panel_gp2_pruned'
# ref_labels = f'{ref_dir_path}/ref_panel_ancestry.txt'

steps = []

avg_miss = avg_miss_rates(geno_path, f'{geno_path}_missing')

# prune callrate
callrate_out = f'{geno_path}_callrate'
callrate = callrate_prune(geno_path, callrate_out)
steps.append(callrate)

# prune sex
sex_out = f'{callrate_out}_sex'
sex = sex_prune(callrate_out, sex_out)
steps.append(sex)

# ancestry estimation
ancestry_out = f'{sex_out}_ancestry'
ancestry = run_ancestry(geno_path=sex_out, out_path=ancestry_out, ref_panel=ref_panel, ref_labels=ref_labels)

# split cohort by ancestry group
pred_labels_path = ancestry['output']['predicted_labels']['labels_outpath']
cohort_split = split_cohort_ancestry(geno_path=sex_out, labels_path=pred_labels_path, out_path=ancestry_out)

# run 
het_dict = dict()
related_dict = dict()
variant_dict = dict()
avg_miss_dict = dict()

for geno, label in zip(cohort_split['paths'], cohort_split['labels']):

    # het
    het_out = f'{geno}_het'
    het = het_prune(geno, het_out)
    het_dict[label] = het
    # steps.append(het)

    # related
    related_out = f'{het_out}_related'
    related = related_prune(het_out, related_out)
    related_dict[label] = related
    # steps.append(related)
    # variant
    variant_out = f'{related_out}_variant'
    variant = variant_prune(related_out, variant_out)
    variant_dict[label] = variant
    # steps.append(variant)
    avg_miss_out = f'{variant_out}_missing'
    avg_miss = avg_miss_rates(variant_out, avg_miss_out)
    avg_miss_dict[label] = avg_miss

    steps2 = [het_dict, related_dict, variant_dict]

with open(f'{out}_QC_REPORT.txt', 'w') as f:
    f.write("QC REPORT\n")
    f.write("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n")
    f.write("$$$$$$$$$$$$$$ Whole cohort steps $$$$$$$$$$$$$$$\n")
    f.write("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n")
    f.write("\n")
    f.write(f"Missingness per Ancestry BEFORE QC: {avg_miss}\n")
    f.write("\n")

    for step in steps:
        f.write("\n")
        f.write(f"STEP: {step['step']}\n")
        f.write(f"Metrics: {step['metrics']}\n")
        f.write("\n")
    
    f.write("\n")
    f.write("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n")
    f.write("$$$$$$$$$$ Individual Ancestry steps $$$$$$$$$$$$\n")
    f.write("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n")
    f.write("\n")
    
    for item in steps2:
        for key, step in item.items():
            f.write("\n")
            f.write(f"Ancestry Group: {key}\n")
            f.write(f"STEP: {step['step']}\n")
            f.write(f"{step['metrics']}\n")
            f.write("\n")

    for key, value in avg_miss_dict.items():
        f.write("\n")
        f.write(f"Ancestry Group: {key}\n")
        f.write(f"Missingness per Ancestry AFTER QC: {value}\n")
        f.write("\n")
    
    f.write("\n")
    f.write("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n")
    f.write("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n")
    f.write("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n")
    f.write("\n")