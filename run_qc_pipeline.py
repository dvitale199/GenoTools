import pandas as pd
import argparse
import shutil
import os

# local imports
from QC.qc import callrate_prune, het_prune, sex_prune, related_prune, variant_prune
from Ancestry.ancestry import run_ancestry, split_cohort_ancestry
from QC.utils import shell_do


parser = argparse.ArgumentParser(description='Arguments for Genotyping QC (data in Plink .bim/.bam/.fam format)')
parser.add_argument('--geno', type=str, default='nope', help='Genotype: (string file path). Path to PLINK format genotype file, everything before the *.bed/bim/fam [default: nope].')
parser.add_argument('--ref', type=str, default='nope', help='Genotype: (string file path). Path to PLINK format reference genotype file, everything before the *.bed/bim/fam.')
parser.add_argument('--ref_labels', type=str, default='nope', help='tab-separated plink-style IDs with ancestry label (FID  IID label) with no header')
parser.add_argument('--model', type=str, default=None, help='Path to pickle file with trained ancestry model for passed reference panel')
parser.add_argument('--callrate', type=float, default=0.02, help='Minimum Callrate threshold for QC')
parser.add_argument('--out', type=str, default='nope', help='Prefix for output (including path)')

args = parser.parse_args()

geno_path = args.geno
ref_panel = args.ref
ref_labels = args.ref_labels
model_path = args.model
callrate = args.callrate
out_path = args.out

# sample-level pruning and metrics
callrate_out = f'{geno_path}_callrate'
callrate = callrate_prune(geno_path, callrate_out, mind=callrate)

sex_out = f'{callrate_out}_sex'
sex = sex_prune(callrate_out, sex_out)


# run ancestry methods
ancestry_out = f'{sex_out}_ancestry'
ancestry = run_ancestry(geno_path=sex_out, out_path=ancestry_out, ref_panel=ref_panel, ref_labels=ref_labels, model_path=model_path)

# get ancestry counts to add to output .h5 later
ancestry_counts_df = pd.DataFrame(ancestry['metrics']['adjusted_counts']).reset_index()
ancestry_counts_df.columns = ['label', 'count']


# split cohort into individual ancestry groups
pred_labels_path = ancestry['output']['adjusted_labels']['labels_outpath']
cohort_split = split_cohort_ancestry(geno_path=sex_out, labels_path=pred_labels_path, out_path=ancestry_out)

# ancestry-specific pruning steps
het_dict = dict()
related_dict = dict()
variant_dict = dict()

for geno, label in zip(cohort_split['paths'], cohort_split['labels']):

    # related
    related_out = f'{geno}_related'
    related = related_prune(geno, related_out, prune_related=False)
    related_dict[label] = related
    
    # het
    het_out = f'{related_out}_het'
    het = het_prune(related_out, het_out)
    het_dict[label] = het
    
    # variant
    variant_out = f'{het_out}_variant'
    if het['pass']:
        variant = variant_prune(het_out, variant_out)
        variant_dict[label] = variant
    else:
        variant = variant_prune(related_out, variant_out)
        variant_dict[label] = variant



# copy output to out_path
for label, data in variant_dict.items():
    if data['pass']:
        for suffix in ['bed','bim','fam','log']:
            plink_file = f"{data['output']['plink_out']}.{suffix}"
            plink_outfile = f'{out_path}_{label}.{suffix}'
            shutil.copyfile(src=plink_file, dst=plink_outfile)

# copy list of related samples to out_path
for label, data in related_dict.items():
    if data['pass']:
        related_file = f"{data['output']['related_samples']}"
        related_outfile = f"{out_path}_{label}.related"
        shutil.copyfile(src=related_file, dst=related_outfile)

# build report- eventually make this an individual method
steps = [callrate, sex]
steps2 = [het_dict, related_dict, variant_dict]
metrics_df = pd.DataFrame()
pruned_samples_df = pd.DataFrame()

for item in steps:
    
    step = item['step']
    pf = item['pass']
    level = 'sample'
    ancestry_label = 'all'
    
    for metric, value in item['metrics'].items():
        tmp_metrics_df = pd.DataFrame({'step':[step], 'pruned_count':[value], 'metric':[metric], 'ancestry':[ancestry_label], 'level':[level], 'pass': [pf]})
        metrics_df = metrics_df.append(tmp_metrics_df)
    
    samplefile = item['output']['pruned_samples']
    if os.path.isfile(samplefile):
        pruned = pd.read_csv(samplefile, sep='\t')
        if pruned.shape[0] > 0:
            pruned.loc[:,'step'] = step
            pruned_samples_df = pruned_samples_df.append(pruned[['FID','IID','step']])
        
for item in steps2:
    for ancestry_label, metrics in item.items():
        
        step = metrics['step']
        pf = metrics['pass']
        
        if step in ['het_prune','related_prune']:
            level = 'sample'

            samplefile = metrics['output']['pruned_samples']
            if os.path.isfile(samplefile):
                pruned = pd.read_csv(samplefile, sep='\t', header=0, usecols=[0,1], names=['FID','IID'])
                if pruned.shape[0] > 0:
                    pruned.loc[:,'step'] = step
                    pruned_samples_df = pruned_samples_df.append(pruned[['FID','IID','step']])
            
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
pred_ancestry_labels = ancestry['data']['admix_data']['ids']

metrics_df.to_hdf(metrics_outfile, key='QC', mode='w')
pruned_samples_df.to_hdf(metrics_outfile, key='pruned_samples')
ancestry_counts_df.to_hdf(metrics_outfile, key='ancestry_counts')
pred_ancestry_labels.to_hdf(metrics_outfile, key='ancestry_labels')
conf_mat_df.to_hdf(metrics_outfile, key='confusion_matrix', index=True)
ref_pcs.to_hdf(metrics_outfile, key='ref_pcs')
projected_pcs.to_hdf(metrics_outfile, key='projected_pcs')
total_umap.to_hdf(metrics_outfile, key='total_umap')
ref_umap.to_hdf(metrics_outfile, key='ref_umap')
new_samples_umap.to_hdf(metrics_outfile, key='new_samples_umap')

