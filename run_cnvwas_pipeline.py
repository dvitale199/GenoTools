from CNV.cnv import CNV_WAS
import argparse

parser = argparse.ArgumentParser(description='Arguments for Running CNV-WAS Pipeline.')    
parser.add_argument('--cnv_dosage_file', type=str, default='Nope.', help='Path to CNV Dosage file (output from run_cnv_dosage_pipeline.py')
parser.add_argument('--out_path', type=str, default='Nope.', help='Path to output summary stats.')
parser.add_argument('--pheno', type=str, default='Nope.', help='path to phenotype file')
parser.add_argument('--covar', type=str, default='Nope.', help='path to covariate file')

args = parser.parse_args()

cnv_dosage_file = args.cnv_dosage_file
pheno = args.pheno
covar = args.covar
out_path = args.out_path

CNV_WAS(cnv_dosage_file, pheno, covar, out_path)