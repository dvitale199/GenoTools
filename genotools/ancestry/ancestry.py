import pandas as pd
import subprocess
import sys
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import seaborn as sns
from matplotlib import cm 
import numpy as np
import os
import shutil
from umap import UMAP
from sklearn import preprocessing, metrics
from sklearn.impute import SimpleImputer
from sklearn.decomposition import PCA
from sklearn.model_selection import train_test_split, GridSearchCV, StratifiedKFold
from sklearn.pipeline import Pipeline
from sklearn.svm import LinearSVC
import plotly.express as px
import plotly
import joblib
import pickle as pkl

#local imports
from QC.utils import shell_do, get_common_snps, rm_tmps, merge_genos

from utils.dependencies import check_plink, check_plink2

class Ancestry:
    def __init__(self):
        self.plink_exec = self.check_plink()
        self.plink2_exec = self.check_plink2()

    def check_plink(self):
        # Your implementation for checking plink executable goes here
        pass

    def check_plink2(self):
        # Your implementation for checking plink2 executable goes here
        pass

    def get_raw_files(self, geno_path, ref_path, labels_path, out_path, train):
        # Your implementation for get_raw_files goes here
        pass

    def munge_training_data(self, labeled_ref_raw):
        # Your implementation for munge_training_data goes here
        pass

    def calculate_pcs(self, X_train, X_test, y_train, y_test, train_ids, test_ids, raw_geno, label_encoder, out, plot_dir):
        # Your implementation for calculate_pcs goes here
        pass

    def train_umap_classifier(self, X_train, X_test, y_train, y_test, label_encoder, out, plot_dir, input_param_grid=None):
        # Your implementation for train_umap_classifier goes here
        pass

    def load_umap_classifier(self, pkl_path, X_test, y_test):
        # Your implementation for load_umap_classifier goes here
        pass

    def predict_ancestry_from_pcs(self, projected, pipe_clf, label_encoder, out):
        # Your implementation for predict_ancestry_from_pcs goes here
        pass

    def umap_transform_with_fitted(self, ref_pca, X_new, y_pred, classifier=None):
        # Your implementation for umap_transform_with_fitted goes here
        pass

    def split_cohort_ancestry(self, geno_path, labels_path, out_path):
        # Your implementation for split_cohort_ancestry goes here
        pass

    def run_ancestry(self, geno_path, out_path, ref_panel, ref_labels, model_path, train_param_grid=None):
        # Your implementation for run_ancestry goes here
        pass
