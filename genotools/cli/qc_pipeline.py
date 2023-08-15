import pandas as pd
import argparse
import shutil
import os

# local imports
from utils.dependencies import check_plink, check_plink2
from QC.qc import QC_processor
from Ancestry.ancestry import Ancestry
from QC.utils import QC_util

plink_exec = check_plink()
plink2_exec = check_plink2()


def main(args):

    geno_path = args.geno
    ref_panel = args.ref
    ref_labels = args.ref_labels
    model_path = args.model
    callrate = args.callrate
    out_path = args.out

    # sample-level pruning and metrics
    qc_processor = QC_processor(plink_exec, plink2_exec)
    callrate_out = f'{geno_path}_callrate'
    qc_processor.callrate_prune(geno_path, callrate_out, mind=callrate)

    sex_out = f'{qc_processor.get_recent_file()}_sex'
    qc_processor.sex_prune(callrate_out, sex_out)


    # run ancestry methods
    ancestry_out = f'{qc_processor.get_recent_file()}_ancestry'
    ancestry = Ancestry()
    ancestry.run_ancestry(geno_path, ancestry_out, ref_panel, ref_labels, model_path)

    ### more class usage examples

    # qc_processor.het_prune(geno_path, out_path)
    # qc_processor.related_prune(geno_path, out_path, related_cutoff=0.0884, duplicated_cutoff=0.354, prune_related=True, prune_duplicated=True)
    # qc_processor.variant_prune(geno_path, out_path)
    # qc_processor.miss_rates(geno_path, out_path, max_threshold=0.05)
    # qc_processor.plink_pca(geno_path, out_path, build='hg38')

    # qc_util = QC_util()
    # geno_path1 = "path_to_genotype_file1"
    # geno_path2 = "path_to_genotype_file2"
    # out_name = "output_name"
    # qc_util.merge_genos(geno_path1, geno_path2, out_name)

    ### remaining steps in pipeline
    # get ancestry counts to add to output .h5 later
    # split cohort into individual ancestry groups
    # ancestry-specific pruning steps
    # copy output to out_path
    # copy list of related samples to out_path
    # build output hdf
    

