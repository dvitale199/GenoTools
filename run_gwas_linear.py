import os
print(os.getcwd())
print(sys.path)
import pandas as pd
import numpy as np
import warnings
from scipy.stats import ncx2
from genotools.dependencies import check_plink, check_plink2
from genotools.gwas import *
from genotools.utils import *


project_path = '/data/GP2/projects/20240109_AOO_GWAS_/gp2_release6_subset_1k'
pca = run_plink_pca(project_path)
run_gwas(project_path, covars=True, linear=True)