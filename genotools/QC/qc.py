import subprocess
import argparse
import pandas as pd
import numpy as np
import os
import glob
import shutil
import sys

# local imports
from QC.utils import shell_do, rm_tmps, count_file_lines, concat_logs

from utils.dependencies import check_plink, check_plink2

plink_exec = check_plink()
plink2_exec = check_plink2()

class QC_processor:
    def __init__(self, plink_exec, plink2_exec):
        self.plink_exec = plink_exec
        self.plink2_exec = plink2_exec
    
    def callrate_prune(self, geno_path, out_path, mind=0.02):
        # Method implementation
        self.latest_file = out_path
        pass

    def sex_prune(self, geno_path, out_path, check_sex=[0.25, 0.75]):
        # Method implementation
        self.latest_file = out_path
        
        pass

    def het_prune(self, geno_path, out_path):
        # Method implementation
        pass

    def related_prune(self, geno_path, out_path, related_cutoff=0.0884, duplicated_cutoff=0.354, prune_related=True, prune_duplicated=True):
        # Method implementation
        pass

    def variant_prune(self, geno_path, out_path):
        # Method implementation
        pass

    def miss_rates(self, geno_path, out_path, max_threshold=0.05):
        # Method implementation
        pass

    def plink_pca(self, geno_path, out_path, build='hg38'):
        # Method implementation
        pass

    def get_recent_file(self):
        return self.latest_file