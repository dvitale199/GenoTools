import pandas as pd
import subprocess
import sys
import pickle as pkl
import numpy as np
from sklearn import metrics
import json
import os


if __name__=='__main__':
    X_test = pd.read_csv(f'/app/X_test.txt', sep='\t')
    y_test = pd.read_csv(f'/app/y_test.txt', sep='\t', header=None)
    projected = pd.read_csv(f'/app/projected.txt', sep='\t')

    pkl_path = f'/app/GP2_merge_APRIL_2023_ready_genotools_callrate_sex_ancestry_umap_linearsvc_ancestry_model.pkl'
    pkl_in = open(pkl_path, 'rb')
    pipe_clf = pkl.load(pkl_in)
    pkl_in.close()

    test_acc = pipe_clf.score(X_test, y_test)
    margin_of_error = 1.96 * np.sqrt((test_acc * (1-test_acc)) / np.shape(y_test)[0])
    low_ci95 = test_acc-margin_of_error
    high_ci95 = test_acc+margin_of_error
    acc_dict = {'test_acc':test_acc, 'margin_of_error':margin_of_error}

    pipe_clf_pred = pipe_clf.predict(X_test)

    X_new = projected.drop(columns=['FID','IID','label'])
    y_pred = pipe_clf.predict(X_new)
    params = pipe_clf.get_params()

    params_dict = {'umap__a':params['umap__a'],'umap__b':params['umap__b'],
                   'umap__n_components':params['umap__n_components'],'umap__n_neighbors':params['umap__n_neighbors']}

    with open(f'/app/accuracy.json', 'w') as f:
        json.dump(acc_dict, f)
        f.close()

    pd.Series(pipe_clf_pred).to_csv(f'/app/pipe_clf_pred.txt', sep='\t', header=None, index=None)

    pd.Series(y_pred).to_csv(f'/app/predicted_labels.txt', sep='\t', header=None, index=None)

    with open(f'/app/params.json', 'w') as f:
        json.dump(params_dict, f)
        f.close()