import os
from genotools.gwas_copy import *
from genotools.utils import *


geno_path = '/data/GP2/projects/20240109_AOO_GWAS_/age_subset2/GP2_release6_10k'
out_path = '/data/GP2/projects/20240109_AOO_GWAS_/logs_age/GP2_release6_10k_out'

my_project = Assoc(geno_path=geno_path, out_path=out_path, pca=10, build='hg38', gwas=True, covar_path=None, covar_names=None)

my_project.write_exclusion_file()

print('exclusion file wriiten')

my_project.run_pca_pruning()

print('data pruned')

my_project.run_association()
d