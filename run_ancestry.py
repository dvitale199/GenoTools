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
import umap.plot
from sklearn import preprocessing, metrics
from sklearn.model_selection import train_test_split, GridSearchCV, StratifiedKFold
from sklearn.pipeline import Pipeline
from sklearn.svm import LinearSVC
import joblib
import plotly.express as px
import plotly
# import sklearn.externals.joblib as extjoblib
import joblib
import argparse

#local imports
# from gwas.qc import QC
from gwas.utils import shell_do, merge_genos, ld_prune, random_sample_snps, get_common_snps
from gwas.ancestry import ancestry_prune, flash_pca, pca_projection, plot_3d

def calculate_pcs(geno, ref, ref_labels, out):

    geno_name = geno
    out_path = out
    ref_panel = ref

    ref_panel_common_snps = f'{out_path}/ref_panel_common_snps'

    # create directories for output files
    if os.path.exists(f'{out_path}/plot_ancestry'):
        pass
    else:
        os.mkdir(f'{out_path}/plot_ancestry')

    if os.path.exists(f'{out_path}/models'):
        pass
    else:
        os.mkdir(f'{out_path}/models')

    # get common snps between ref panel and geno
    get_common_snps(ref_panel, geno_name, ref_panel_common_snps)

    # run pca on ref panel
    flash_pca(ref_panel_common_snps, ref_panel_common_snps, dim=50)

    # read ancestry file
    ancestry = pd.read_csv(f'{ref_labels}', sep='\t', header=None, names=['FID','IID','label'])
    ref_fam = pd.read_csv(f'{ref_panel}.fam', sep=' ', header=None)
    ref_labeled = ref_fam.merge(ancestry, how='left', left_on=[0,1], right_on=['FID','IID'])

    pca = pd.read_csv(f'{ref_panel_common_snps}.pcs', sep='\t')

    # combined_labels
    labeled_pca = pca.merge(ref_labeled, how='left', on=['FID','IID'])
    labeled_pca.drop(columns=[0,1,2,3,4,5],inplace=True)
    print(labeled_pca.label.value_counts())

    # plot it!
    plot_3d(labeled_pca, color='label', title='Reference Panel PCA', plot_out=f'{out_path}/plot_ancestry/plot_ref_panel_pcs', x='PC1', y='PC2', z='PC3')


    # get reference alleles from ref_panel_common_snps
    ref_panel_common_snps_ref_alleles = f'{ref_panel_common_snps}.ref_allele'
    ref_panel_common_snps_bim = pd.read_csv(f'{ref_panel_common_snps}.bim', header=None, sep='\t')
    ref_panel_common_snps_bim.columns = ['chr', 'rsid', 'kb', 'pos', 'a1', 'a2']
    ref_panel_common_snps_bim[['rsid','a1']].to_csv(ref_panel_common_snps_ref_alleles, sep='\t', header=False, index=False)

    geno_common_snps = f'{geno_name}_common_snps'
    common_snps = f'{ref_panel_common_snps}.common_snps'

    ext_snps_cmd = f'plink --bfile {geno_name} --extract {common_snps} --reference-allele {ref_panel_common_snps_ref_alleles} --make-bed --out {geno_common_snps}'

    shell_do(ext_snps_cmd)

    # project new samples onto ref pcs
    pca_projection(geno_common_snps, f'{ref_panel_common_snps}.meansd', f'{ref_panel_common_snps}.loadings', f'{geno_common_snps}.projections')

    projected = pd.read_csv(f'{geno_common_snps}.projections', sep='\t')
    projected['label'] = 'new'
    total_pca = labeled_pca.append(projected)

    plot_3d(total_pca, color='label', title='New Samples Projected on Reference Panel', plot_out=f'{out_path}/plot_ancestry/plot_PCA_projected_new_samples', x='PC1', y='PC2', z='PC3')


def train_umap_classifier_pipeline(pca_loadings):


    # grab X and Y
    X = labeled_pca.drop(columns=['label'])
    y = labeled_pca.label

    #encode labels
    le = preprocessing.LabelEncoder()
    y = le.fit_transform(y)

    # train/test split 1kg pca data
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=123)

    IDs_train = X_train[['FID', 'IID']]
    IDs_test = X_test[['FID', 'IID']]
    X_train = X_train.drop(columns=['FID','IID'])
    X_test = X_test.drop(columns=['FID','IID'])

    # Transformation with UMAP followed by classification with svc
    umap = UMAP(random_state=123)
    svc = LinearSVC(dual=False, random_state=123)
    pipeline = Pipeline([("umap", umap), ("svc", svc)])
    param_grid = {
        "umap__n_neighbors": [5, 20],
        "umap__n_components": [15, 25],
        "umap__a":[0.75, 1.0, 1.5],
        "umap__b": [0.25, 0.5, 0.75],
        "svc__C": [10**i for i in range(-3,2)],
    }

    cross_validation = StratifiedKFold(n_splits=5, shuffle=True, random_state=123)
    pipe_grid = GridSearchCV(pipeline, param_grid, cv=cross_validation, scoring='balanced_accuracy')
    pipe_grid.fit(X_train, y_train)

    print(pipe_grid.best_score_)
    print(pipe_grid.best_params_)

    pipe_clf = pipe_grid.best_estimator_
    print(f"Accuracy on the test set with UMAP transformation: {pipe_clf.score(X_test, y_test)}")
    pipe_clf_pred = pipe_clf.predict(X_test)

    pipe_clf_c_matrix = metrics.confusion_matrix(y_test, pipe_clf_pred)
    fig, ax = plt.subplots(figsize=(10,10))
    sns.heatmap(pipe_clf_c_matrix, annot=True, fmt='d',
                xticklabels=le.inverse_transform([i for i in range(8)]), yticklabels=le.inverse_transform([i for i in range(8)]))
    plt.ylabel('Actual')
    plt.xlabel('Predicted')
    plt.show()
    fig.savefig(f'{out_path}/plot_ancestry/plot_umap_linearsvc_ancestry_conf_matrix.png')

    # dump best estimator to pkl
    joblib.dump(pipe_clf, f'{out_path}/models/umap_linearsvc_ancestry_model.pkl')





# set new samples aside for labeling after training the model
X_new = projected.drop(columns=['FID','IID','label'])
X_new_ids = projected[['FID','IID']]
X_new


# predict new samples
ancestry_pred = pipe_clf.predict(X_new)
X_new_ids.loc[:,'label'] = le.inverse_transform(ancestry_pred)
print('predicted:\n', X_new_ids.label.value_counts())
X_new_ids.to_csv(f'{geno_name}_umap_linearsvc_predicted_labels.txt', sep='\t')


# return pipe_grid from train_umap_classifier_pipeline()
def plot_umap():
    a = pipe_grid.best_params_['umap__a']
    b = pipe_grid.best_params_['umap__b']

    n_components = pipe_grid.best_params_['umap__n_components']
    n_neighbors = pipe_grid.best_params_['umap__n_neighbors']
    X_ = X.drop(columns=['FID','IID'])
    umapper = UMAP(random_state=123, n_components=n_components, n_neighbors=n_neighbors, a=a, b=b).fit(X_)

    ref_umap = pd.DataFrame(umapper.transform(X_))
    ref_umap.loc[:,'label'] = le.inverse_transform(y)
    new_samples_umap = pd.DataFrame(umapper.transform(X_new))
    new_samples_umap.loc[:,'label'] = le.inverse_transform(ancestry_pred)

    ref_umap.loc[:,'dataset'] = 'ref'
    new_samples_umap.loc[:, 'dataset'] = 'predicted'
    total_umap = ref_umap.append(new_samples_umap)

    plot_3d(total_umap, color='label', symbol='dataset', plot_out=f'{out_path}/plot_ancestry/plot_total_umap', x=0,y=1,z=2)
    plot_3d(ref_umap, color='label', symbol='dataset', plot_out=f'{out_path}/plot_ancestry/plot_ref_umap', x=0,y=1,z=2)
    plot_3d(new_samples_umap, color='label', symbol='dataset', plot_out=f'{out_path}/plot_ancestry/plot_predicted_samples_umap', x=0,y=1,z=2)

