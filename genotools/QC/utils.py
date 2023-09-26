import subprocess
import sys
import matplotlib.pyplot as plt
from matplotlib import cm 
import numpy as np
import os
import shutil
import pandas as pd

class QC_util:
    def __init__(self):
        self.plink_exec = self.check_plink()
        self.plink2_exec = self.check_plink2()

    def check_plink(self):
        # Implementation of check_plink function
        pass

    def check_plink2(self):
        # Implementation of check_plink2 function
        pass

    def shell_do(self, command, log=False, return_log=False):
        # Implementation of shell_do function
        pass

    def replace_all(self, text, dict):
        # Implementation of replace_all function
        pass

    def process_log(self, out_dir, concat_log):
        # Implementation of process_log function
        pass

    def concat_logs(self, step, out_path, listOfFiles):
        # Implementation of concat_logs function
        pass

    def label_bim_with_genes(self, bim_file, gene_reference=None, locus_size=1000000):
        # Implementation of label_bim_with_genes function
        pass

    def merge_genos(self, geno_path1, geno_path2, out_name):
        # Implementation of merge_genos function
        pass

    def ld_prune(self, geno_path, out_name, window_size=1000, step_size=50, rsq_thresh=0.05):
        # Implementation of ld_prune function
        pass

    def random_sample_snps(self, geno_path, out_name, n=10000):
        # Implementation of random_sample_snps function
        pass

    def get_common_snps(self, geno_path1, geno_path2, out_name):
        # Implementation of get_common_snps function
        pass

    def rm_tmps(self, tmps, suffixes=None):
        # Implementation of rm_tmps function
        pass

    def count_file_lines(self, file_path):
        # Implementation of count_file_lines function
        pass