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
out_path = args.out

# sample-level pruning and metrics
avg_miss = avg_miss_rates(geno_path, f'{geno_path}_missing')
avg_miss

callrate_out = f'{geno_path}_callrate'
callrate = callrate_prune(geno_path, callrate_out)

sex_out = f'{callrate_out}_sex'
sex = sex_prune(callrate_out, sex_out)


# run ancestry methods
ancestry_out = f'{sex_out}_ancestry'
ancestry = run_ancestry(geno_path=sex_out, out_path=ancestry_out, ref_panel=ref_panel, ref_labels=ref_labels)

# get ancestry counts to add to output .h5 later
ancestry_counts_df = pd.DataFrame(ancestry['metrics']['predicted_counts']).reset_index()
ancestry_counts_df.columns = ['label', 'count']


# split cohort into individual ancestry groups
pred_labels_path = ancestry['output']['predicted_labels']['labels_outpath']
cohort_split = split_cohort_ancestry(geno_path=sex_out, labels_path=pred_labels_path, out_path=ancestry_out)


# ancestry-specific pruning steps
het_dict = dict()
related_dict = dict()
variant_dict = dict()

for geno, label in zip(cohort_split['paths'], cohort_split['labels']):

    # related
    related_out = f'{geno}_related'
    related = related_prune(geno, related_out)
    related_dict[label] = related
    
    # het
    het_out = f'{related_out}_het'
    het = het_prune(geno, het_out)
    het_dict[label] = het
    
    # variant
    variant_out = f'{het_out}_variant'
    if het['pass']:
        variant = variant_prune(het_out, variant_out)
        variant_dict[label] = variant
    else:
        variant = variant_prune(related_out, variant_out)
        variant_dict[label] = variant


# build report- eventually make this an individual method
steps = [callrate, sex]
steps2 = [het_dict, related_dict, variant_dict]
metrics_df = pd.DataFrame()

for item in steps:
    
    step = item['step']
    pf = item['pass']
    level = 'sample'
    ancestry_label = 'all'
    
    for metric, value in item['metrics'].items():
        tmp_metrics_df = pd.DataFrame({'step':[step], 'pruned_count':[value], 'metric':[metric], 'ancestry':[ancestry_label], 'level':[level], 'pass': [pf]})
        metrics_df = metrics_df.append(tmp_metrics_df)


for item in steps2:
    for ancestry_label, metrics in item.items():
        
        step = metrics['step']
        pf = metrics['pass']
        
        if step in ['het_prune','related_prune']:
            level = 'sample'
        else:
            level = 'variant'

        for metric, value in metrics['metrics'].items():
            tmp_metrics_df = pd.DataFrame({'step':[step], 'pruned_count':[value], 'metric':[metric], 'ancestry':[ancestry_label], 'level':[level], 'pass': [pf]})
            metrics_df = metrics_df.append(tmp_metrics_df)

metrics_df.reset_index(drop=True, inplace=True)


# build output hdf
metrics_outfile = f'{out_path}.QC.metrics.h5'

le = ancestry['data']['label_encoder']
confusion_matrix = ancestry['data']['confusion_matrix']
conf_mat_df = pd.DataFrame(confusion_matrix)
conf_mat_df.columns = le.inverse_transform([i for i in range(8)])
conf_mat_df.index = le.inverse_transform([i for i in range(8)])

ref_pcs = ancestry['data']['ref_pcs']
projected_pcs = ancestry['data']['projected_pcs']
total_umap = ancestry['data']['total_umap']
ref_umap = ancestry['data']['ref_umap']
new_samples_umap = ancestry['data']['new_samples_umap']
pred_ancestry_labels = ancestry['data']['predict_data']['ids']

metrics_df.to_hdf(metrics_outfile, key='QC', mode='w')
ancestry_counts_df.to_hdf(metrics_outfile, key='ancestry_counts')
pred_ancestry_labels.to_hdf(metrics_outfile, key='ancestry_labels')
conf_mat_df.to_hdf(metrics_outfile, key='confusion_matrix', index=True)
ref_pcs.to_hdf(metrics_outfile, key='ref_pcs')
projected_pcs.to_hdf(metrics_outfile, key='projected_pcs')
total_umap.to_hdf(metrics_outfile, key='total_umap')
ref_umap.to_hdf(metrics_outfile, key='ref_umap')
new_samples_umap.to_hdf(metrics_outfile, key='new_samples_umap')







# geno_path = '/data/CARD/PD/GP2/raw_genotypes/shulman_ny/plink/shulman'
# hard code in reference for now--- definitely change later
# ref_dir_path = '/data/LNG/vitaled2/1kgenomes'
# ref_panel = f'{ref_dir_path}/1kg_ashkj_ref_panel_gp2_pruned'
# ref_labels = f'{ref_dir_path}/ref_panel_ancestry.txt'

# steps = []

# avg_miss = avg_miss_rates(geno_path, f'{geno_path}_missing')

# # prune callrate
# callrate_out = f'{geno_path}_callrate'
# callrate = callrate_prune(geno_path, callrate_out)
# steps.append(callrate)

# # prune sex
# sex_out = f'{callrate_out}_sex'
# sex = sex_prune(callrate_out, sex_out)
# steps.append(sex)

# # ancestry estimation
# ancestry_out = f'{sex_out}_ancestry'
# ancestry = run_ancestry(geno_path=sex_out, out_path=ancestry_out, ref_panel=ref_panel, ref_labels=ref_labels)

# # split cohort by ancestry group
# pred_labels_path = ancestry['output']['predicted_labels']['labels_outpath']
# cohort_split = split_cohort_ancestry(geno_path=sex_out, labels_path=pred_labels_path, out_path=ancestry_out)

# # run 
# het_dict = dict()
# related_dict = dict()
# variant_dict = dict()
# avg_miss_dict = dict()

# for geno, label in zip(cohort_split['paths'], cohort_split['labels']):

#     # het
#     het_out = f'{geno}_het'
#     het = het_prune(geno, het_out)
#     het_dict[label] = het
#     # steps.append(het)

#     # related
#     related_out = f'{het_out}_related'
#     related = related_prune(het_out, related_out)
#     related_dict[label] = related
#     # steps.append(related)
#     # variant
#     variant_out = f'{related_out}_variant'
#     variant = variant_prune(related_out, variant_out)
#     variant_dict[label] = variant
#     # steps.append(variant)
#     avg_miss_out = f'{variant_out}_missing'
#     avg_miss = avg_miss_rates(variant_out, avg_miss_out)
#     avg_miss_dict[label] = avg_miss

#     steps2 = [het_dict, related_dict, variant_dict]

# with open(f'{out}_QC_REPORT.txt', 'w') as f:
#     f.write("QC REPORT\n")
#     f.write("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n")
#     f.write("$$$$$$$$$$$$$$ Whole cohort steps $$$$$$$$$$$$$$$\n")
#     f.write("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n")
#     f.write("\n")
#     f.write(f"Missingness per Ancestry BEFORE QC: {avg_miss}\n")
#     f.write("\n")

#     for step in steps:
#         f.write("\n")
#         f.write(f"STEP: {step['step']}\n")
#         f.write(f"Metrics: {step['metrics']}\n")
#         f.write("\n")
    
#     f.write("\n")
#     f.write("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n")
#     f.write("$$$$$$$$$$ Individual Ancestry steps $$$$$$$$$$$$\n")
#     f.write("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n")
#     f.write("\n")
    
#     for item in steps2:
#         for key, step in item.items():
#             f.write("\n")
#             f.write(f"Ancestry Group: {key}\n")
#             f.write(f"STEP: {step['step']}\n")
#             f.write(f"{step['metrics']}\n")
#             f.write("\n")

#     for key, value in avg_miss_dict.items():
#         f.write("\n")
#         f.write(f"Ancestry Group: {key}\n")
#         f.write(f"Missingness per Ancestry AFTER QC: {value}\n")
#         f.write("\n")
    
#     f.write("\n")
#     f.write("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n")
#     f.write("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n")
#     f.write("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n")
#     f.write("\n")