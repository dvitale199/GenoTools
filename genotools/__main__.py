import pandas as pd
import argparse
import shutil
import os

# local imports
from genotools.utils import dependencies
from genotools.cli import qc_pipeline


def qc_pipeline():
    parser = argparse.ArgumentParser(description='Arguments for Genotyping QC (data in Plink .bim/.bam/.fam format)')
    parser.add_argument('--geno', type=str, default='nope', help='Genotype: (string file path). Path to PLINK format genotype file, everything before the *.bed/bim/fam [default: nope].')
    parser.add_argument('--ref', type=str, default='nope', help='Genotype: (string file path). Path to PLINK format reference genotype file, everything before the *.bed/bim/fam.')
    parser.add_argument('--ref_labels', type=str, default='nope', help='tab-separated plink-style IDs with ancestry label (FID  IID label) with no header')
    parser.add_argument('--model', type=str, default=None, help='Path to pickle file with trained ancestry model for passed reference panel')
    parser.add_argument('--callrate', type=float, default=0.02, help='Minimum Callrate threshold for QC')
    parser.add_argument('--out', type=str, default='nope', help='Prefix for output (including path)')

    args = parser.parse_args()

    dependencies.check_dependencies()
    qc_pipeline.main(args) # can eventually adjust to *args for more general inputs

if __name__ == "__main__":
    # call different flag handlers
    qc_pipeline()