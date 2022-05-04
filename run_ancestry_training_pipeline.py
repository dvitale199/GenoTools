import pandas as pd
import argparse
import shutil
import os

from Ancestry.ancestry import get_raw_files, munge_training_data, calculate_pcs, train_umap_classifier
from QC.utils import shell_do

parser = argparse.ArgumentParser(description='Arguments for Training Ancestry Model')
parser.add_argument('--geno', type=str, default='nope', help='Genotype: (string file path). Path to PLINK format genotype file, everything before the *.bed/bim/fam')
parser.add_argument('--ref', type=str, default='nope', help='Genotype: (string file path). Path to PLINK format reference genotype file, everyting before the *.bed/bim/fam')
parser.add_argument('--ref_labels', type=str, default='nope', help='tab-separated plink-style IDs with ancestry label (FID  IID label) with no header')
parser.add_argument('--out', type=str, default='nope', help='Prefix for model and common SNPs output (including path)')

args = parser.parse_args()

geno_path = args.geno
ref_panel = args.ref
ref_labels = args.ref_labels
out_path = args.out

outdir = os.path.dirname(out_path)
plot_dir = f'{outdir}/plot_ancestry'

raw = get_raw_files(geno_path=geno_path, ref_path=ref_panel, labels_path=ref_labels, out_path=out_path, train=True)

train_split = munge_training_data(labeled_ref_raw=raw['raw_ref'])

calc_pcs = calculate_pcs(
        X_train=train_split['X_train'], 
        X_test=train_split['X_test'],
        y_train=train_split['y_train'],
        y_test=train_split['y_test'],
        train_ids=train_split['train_ids'],
        test_ids=train_split['test_ids'],
        raw_geno=raw['raw_geno'],
        label_encoder=train_split['label_encoder'],
        out=out_path,
        plot_dir=plot_dir
    )

trained_clf = train_umap_classifier(
        X_train=calc_pcs['X_train'],
        X_test=calc_pcs['X_test'],
        y_train=train_split['y_train'],
        y_test=train_split['y_test'],
        label_encoder=train_split['label_encoder'],
        out=out_path,
        plot_dir=plot_dir
)

# output - not really sure what we want here
print('Path to common SNPs between genotype file and reference panel:')
print(f"{raw['out_paths']['geno_common_snps_bed']}.txt")
print()
print('Path to trained model pickle file:')
print(trained_clf['model_path'])