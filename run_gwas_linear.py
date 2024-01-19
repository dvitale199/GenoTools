import os
import pandas as pd
import numpy as np
from genotools.gwas import *
from genotools.utils import *


geno_path = '/data/GP2/projects/20240109_AOO_GWAS_/gp2_release6_subset_1k'
out_path = '/data/GP2/projects/20240109_AOO_GWAS_/gp2_release6_subset_1k_out'

my_project = Assoc(geno_path=geno_path, out_path=out_path, pca=10, build='hg38', gwas=True, covar_path=None, covar_names=None)

my_project.run_plink_pca()