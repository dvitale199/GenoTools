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
from sklearn.model_selection import train_test_split, GridSearchCV, StratifiedKFold
from sklearn.pipeline import Pipeline
from sklearn.svm import LinearSVC
import plotly.express as px
import plotly
import joblib

#local imports
from QC.utils import shell_do, get_common_snps, rm_tmps


def ancestry_prune(geno_path, out_path=None):
    '''Pruning of --maf 0.05, --geno 0.01, --hwe 0.0001, palindrome snps, and high-LD regions for ancestry methods.
    
    Parameters: 
    geno_path (string): path to plink genotype (everything before .bed/.bim/.fam).
    out_path (string): path to output plink genotype (everything before .bed/.bim/.fam).
    
    Returns:
    None.
    '''

    geno_ancestry_prune_tmp = f'{geno_path}_ancestry_prune_tmp'
    
    if out_path:
        geno_ancestry_prune_out = out_path
    else:
        geno_ancestry_prune_out = f'{geno_path}_ancestry_prune'

    # prune geno_path for geno, maf, hwe, and palindromes
    geno_bim = pd.read_csv(f'{geno_path}.bim', sep='\t', header=None)

    # find and drop palindromes in geno_path .bim file
    geno_bim.columns = ['chr', 'rsid', 'kb', 'pos', 'a1', 'a2']
    palindromes = geno_bim.loc[((geno_bim.a1 == 'A') & (geno_bim.a2 == 'T')) | ((geno_bim.a1 == 'T') & (geno_bim.a2 == 'A')) | ((geno_bim.a1 == 'C') & (geno_bim.a2 == 'G')) | ((geno_bim.a1 == 'G') & (geno_bim.a2 == 'C'))]
    palindromes['rsid'].to_csv(f'{geno_path}_palindromes.snplist', header=False, index=False, sep='\t')

    plink_cmd1 = f'plink --bfile {geno_path}\
     --maf 0.05\
     --geno 0.01\
     --hwe 0.0001\
     --autosome\
     --allow-no-sex\
     --exclude {geno_path}_palindromes.snplist\
     --make-bed\
     --out {geno_ancestry_prune_tmp}' 

    # exclude high-LD regions
    plink_cmd2 = f'plink --bfile {geno_ancestry_prune_tmp}\
     --exclude range ref_data/hg19_exclusion_regions.txt\
     --autosome\
     --allow-no-sex\
     --make-bed\
     --out {geno_ancestry_prune_out}'

    cmds = [plink_cmd1, plink_cmd2]

    for cmd in cmds:
        shell_do(cmd)
    
    rm_tmps([f'{geno_ancestry_prune_tmp}'], ['bed','bim','fam','log'])

        
def flash_pca(geno_path, out_name, dim=20):

    """Runs flashpca on input genotype.

    Parameters:
    geno_path (string): path to plink genotype (everything before .bed/.bim/.fam).
    out_path (string): path to flashpca output (everything before .pcs/.vec/.loadings).

    Returns:
    dict: output paths to flashpca output files.
    """

    flashpca_cmd = f'\
    flashpca --bfile {geno_path}\
     -d {dim}\
     --outpc {out_name}.pcs\
     --outvec {out_name}.vec\
     --outval {out_name}.val\
     --outpve {out_name}.pve\
     --outload {out_name}.loadings\
     --outmeansd {out_name}.meansd'
    
    shell_do(flashpca_cmd)

    out_dict = {
        'outpc':f'{out_name}.pcs',
        'outvec': f'{out_name}.vec',
        'outval': f'{out_name}.val',
        'outpve': f'{out_name}.pve',
        'outload': f'{out_name}.loadings',
        'outmeansd': f'{out_name}.meansd'
        }
        
    return out_dict


    
def pca_projection(geno_path, inmeansd, inload, outproj):
    
    '''Project new samples on to PCs in flashpca format. Must use the EXACT same snps (use get_common_snps() below).

    Parameters:
    geno_path (string): path to plink genotype (everything before .bed/.bim/.fam).
    inmeansd (string): path to flashpca meansd file
    inload (string): path to flashpca loadings file
    outproj (string): path to output projected samples
    
    Returns:
    dict: paths to output files
    '''
    
    project_geno_pcs_cmd = f'\
    flashpca --bfile {geno_path}\
     --project\
     --inmeansd {inmeansd}\
     --inload {inload}\
     --outproj {outproj}\
     -v'

    shell_do(project_geno_pcs_cmd)

    out_dict = {
        'outpc':f'{outproj}.pcs',
        'outvec': f'{outproj}.vec',
        'outval': f'{outproj}.val',
        'outload': f'{outproj}.loadings',
        }
    return out_dict


def plot_3d(labeled_df, color, symbol=None, plot_out=None, x='PC1', y='PC2', z='PC3', title=None, x_range=None, y_range=None, z_range=None):
    '''
    Parameters: 
    labeled_df (Pandas dataframe): labeled ancestry dataframe
    color (string): color of ancestry label. column name containing labels for ancestry in labeled_pcs_df
    symbol (string): symbol of secondary label (for example, predicted vs reference ancestry). default: None
    plot_out (string): filename to output filename for .png and .html plotly images
    x (string): column name of x-dimension
    y (string): column name of y-dimension
    z (string): column name of z-dimension
    title (string, optional): title of output scatterplot
    x_range (list of floats [min, max], optional): range for x-axis
    y_range (list of floats [min, max], optional): range for y-axis
    z_range (list of floats [min, max], optional): range for z-axis

    Returns:
    3-D scatterplot (plotly.express.scatter_3d). If plot_out included, will write .png static image and .html interactive to plot_out filename
        
    '''    
    fig = px.scatter_3d(
        labeled_df,
        x=x,
        y=y,
        z=z,
        color=color,
        symbol=symbol,
        title=title,
        color_discrete_sequence=px.colors.qualitative.Bold,
        range_x=x_range,
        range_y=y_range,
        range_z=z_range
    )

    fig.show()

    if plot_out:
        fig.write_image(f'{plot_out}.png', width=1980, height=1080)
        fig.write_html(f'{plot_out}.html')


def calculate_pcs(geno, ref, labels, out, plot_dir, keep_temp=True):

    '''Calculate Principal Components for reference (.bed/.bim/.fam) and project new samples. Plot PCs.

    Parameters:
    geno (string): path to plink genotype (everything before .bed/.bim/.fam).
    ref (string): path to reference genotype (everything before .bed/.bim/.fam)
    labels (string): path to tab-separated file containing labels (FID,IID,label)
    plotdir (string): path to directory where plots will be written to
    keep_temp (boolean): Default=True. if false, delete intermediate files
    
    Returns:
    out_dict = {
        'labeled_ref_pca': labeled reference pca df (pandas dataframe),
        'new_samples_projected': projected sample pca df (pandas dataframe),
        'outpaths': paths to output files (dictionary of strings),
        'temp_paths': paths to temp files (dictionary of strings)
        }
    '''

    step = "calculate_pcs"
    print()
    print(f"RUNNING: {step}")
    print()

    outdir = os.path.dirname(out)
    out_paths = {}
    temp_paths = {}

    # ref_labeled_pca_out = f'{out}/'

    ref_common_snps = f'{outdir}/ref_common_snps'

    # get common snps between ref panel and geno
    common_snps_files = get_common_snps(ref, geno, ref_common_snps)

    # add common_snps_files output paths to temp_paths
    temp_paths = {**temp_paths, **common_snps_files}

    # run pca on ref panel
    ref_common_snps_pca = f'{ref_common_snps}'
    ref_pca = flash_pca(ref_common_snps, ref_common_snps_pca, dim=50)

    # add ref_pca to out_paths
    out_paths = {**out_paths, **ref_pca}

    # read ancestry file with reference labels 
    ancestry = pd.read_csv(f'{labels}', sep='\t', header=None, names=['FID','IID','label'])
    ref_fam = pd.read_csv(f'{ref}.fam', sep=' ', header=None)
    ref_labeled = ref_fam.merge(ancestry, how='left', left_on=[0,1], right_on=['FID','IID'])

    
    pca = pd.read_csv(ref_pca['outpc'], sep='\t')

    # combined_labels
    labeled_pca = pca.merge(ref_labeled, how='left', on=['FID','IID'])
    labeled_pca.drop(columns=[0,1,2,3,4,5],inplace=True)

    print()
    print()
    print("Labeled Reference Ancestry Counts:")
    print(labeled_pca.label.value_counts())
    print()
    print()

    # plot it!
#     plot_3d(labeled_pca, color='label', title='Reference Panel PCA', plot_out=f'{plot_dir}/plot_ref_pcs', x='PC1', y='PC2', z='PC3')

    # get reference alleles from ref_common_snps
    ref_common_snps_ref_alleles = f'{ref_common_snps}.ref_allele'
    ref_common_snps_bim = pd.read_csv(f'{ref_common_snps}.bim', header=None, sep='\t')
    ref_common_snps_bim.columns = ['chr', 'rsid', 'kb', 'pos', 'a1', 'a2']
    ref_common_snps_bim[['rsid','a1']].to_csv(ref_common_snps_ref_alleles, sep='\t', header=False, index=False)
    temp_paths['ref_alleles'] = ref_common_snps_ref_alleles

    geno_common_snps = f'{out}_common_snps'
    common_snps = f'{ref_common_snps}.common_snps'
    out_paths['common_snps'] = common_snps

    ext_snps_cmd = f'plink --bfile {geno} --extract {common_snps} --reference-allele {ref_common_snps_ref_alleles} --make-bed --out {geno_common_snps}'

    shell_do(ext_snps_cmd)

    print(geno_common_snps, f'{ref_common_snps}.meansd', f'{ref_common_snps}.loadings', f'{geno_common_snps}.projections')

    # project new samples onto ref pcs
    projection = pca_projection(geno_common_snps, f'{ref_common_snps}.meansd', f'{ref_common_snps}.loadings', f'{geno_common_snps}.projections')
    out_paths = {**out_paths, **projection}

    projected = pd.read_csv(f'{geno_common_snps}.projections', sep='\t')
    projected['label'] = 'new'
    total_pca = labeled_pca.append(projected)

#     plot_3d(total_pca, color='label', title='New Samples Projected on Reference Panel', plot_out=f'{plot_dir}/plot_PCA_projected_new_samples', x='PC1', y='PC2', z='PC3')

    # write output files
    labeled_pca.to_csv(f'{out}_labeled_ref_pca.txt', sep='\t', index=None)
    projected.to_csv(f'{out}_projected_new_pca.txt', sep='\t', index=None)
    out_paths['ref_pca'] = f'{out}_labeled_ref_pca.txt'
    out_paths['projected_pca'] = f'{out}_projected_new_pca.txt'

    out_dict = {
        'labeled_ref_pca': labeled_pca,
        'new_samples_projected': projected,
        'outpaths': out_paths,
        'temp_paths': temp_paths
        }

    return out_dict


def munge_training_pca_loadings(labeled_pca):
    """ Train/test split and label encode labeled PCA data for training
    
    
    """
    step = "munge_pca_loadings"
    print()
    print(f"RUNNING: {step}")
    print()

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

    out_dict = {
        'X_train': X_train,
        'X_test': X_test,
        'y_train': y_train,
        'y_test': y_test,
        'train_ids': IDs_train,
        'test_ids': IDs_test,
        'label_encoder': le,
        'X_all':X,
        'y_all':y
        }
    
    return out_dict


def train_umap_classifier(X_train, X_test, y_train, y_test, label_encoder, plot_dir, model_dir, input_param_grid=None):
    
    """Train UMAP classifier pipeline
    
    
    """
    step = "train_umap_classifier"
    print()
    print(f"RUNNING: {step}")
    print()

    if input_param_grid:
        param_grid = input_param_grid
    else:
        param_grid = {
        "umap__n_neighbors": [5, 20],
        "umap__n_components": [15, 25],
        "umap__a":[0.75, 1.0, 1.5],
        "umap__b": [0.25, 0.5, 0.75],
        "svc__C": [10**i for i in range(-3,3)],
    }

    le = label_encoder

    # Transformation with UMAP followed by classification with svc
    umap = UMAP(random_state=123)
    svc = LinearSVC(dual=False, random_state=123)
    pipeline = Pipeline([("umap", umap), ("svc", svc)])
    
    cross_validation = StratifiedKFold(n_splits=5, shuffle=True, random_state=123)
    pipe_grid = GridSearchCV(pipeline, param_grid, cv=cross_validation, scoring='balanced_accuracy')
    pipe_grid.fit(X_train, y_train)

    train_acc = pipe_grid.best_score_
    print(f'Training Balanced Accuracy: {train_acc}')
    print(f'Best Parameters: {pipe_grid.best_params_}')

    pipe_clf = pipe_grid.best_estimator_
    test_acc = pipe_clf.score(X_test, y_test)
    print(f"Balanced Accuracy on Test Set: {test_acc}")
    pipe_clf_pred = pipe_clf.predict(X_test)

    pipe_clf_c_matrix = metrics.confusion_matrix(y_test, pipe_clf_pred)
    
    # eventually make this separate function
    # need to get x and y tick labels from 
#     fig, ax = plt.subplots(figsize=(10,10))
#     sns.heatmap(pipe_clf_c_matrix, annot=True, fmt='d',
#                 xticklabels=le.inverse_transform([i for i in range(8)]), yticklabels=le.inverse_transform([i for i in range(8)]))
#     plt.ylabel('Actual')
#     plt.xlabel('Predicted')
#     plt.show()
#     fig.savefig(f'{plot_dir}/plot_umap_linearsvc_ancestry_conf_matrix.png')

    # dump best estimator to pkl
    joblib.dump(pipe_clf, f'{model_dir}/umap_linearsvc_ancestry_model.pkl')

    out_dict = {
        'classifier': pipe_clf,
        'label_encoder' : le,
        'best_params': pipe_grid.best_params_,
        'confusion_matrix': pipe_clf_c_matrix,
        'fitted_pipe_grid': pipe_grid,
        'train_accuracy': train_acc,
        'test_accuracy': test_acc
    }
    
    return out_dict


def predict_ancestry_from_pcs(projected, pipe_clf, label_encoder, out):
    
    le = label_encoder

    # set new samples aside for labeling after training the model
    X_new = projected.drop(columns=['FID','IID','label'])
    # X_new_ids = projected[['FID','IID']]

    # predict new samples
    y_pred = pipe_clf.predict(X_new)
    # X_new_ids.loc[:,'label'] = le.inverse_transform(ancestry_pred)
    ancestry_pred = le.inverse_transform(y_pred)
    projected.loc[:,'label'] = ancestry_pred

    print()
    print('predicted:\n', projected.label.value_counts())
    print()

    projected[['FID','IID','label']].to_csv(f'{out}_umap_linearsvc_predicted_labels.txt', sep='\t', index=False)

    data_out = {
        'ids': projected.loc[:,['FID','IID','label']],
        'X_new': X_new,
        'y_pred': ancestry_pred,
        'label_encoder': le
    }

    # metrics_dict = {
    #     'predicted_labels_counts': projected.label.value_counts(),
    # }

    outfiles_dict = {
        'labels_outpath': f'{out}_umap_linearsvc_predicted_labels.txt'
    }

    out_dict = {
        'data': data_out,
        'metrics': projected.label.value_counts(),
        'output': outfiles_dict
    }

    return out_dict


def umap_transform_with_fitted(X_new, X_ref, y_pred, y_ref, label_encoder, fitted_pipe_grid=None):
    
    le = label_encoder
    pipe_grid = fitted_pipe_grid
    # pipe_clf = pipe_grid.best_estimator_

    # X_new = X_new_pred.drop(columns=['FID','IID','label'])
    X_ = X_ref.drop(columns=['FID','IID'])

    # if fitted_pipe_grid provided, use those params else, use UMAP defaults
    if fitted_pipe_grid:
        a = pipe_grid.best_params_['umap__a']
        b = pipe_grid.best_params_['umap__b']

        n_components = pipe_grid.best_params_['umap__n_components']
        n_neighbors = pipe_grid.best_params_['umap__n_neighbors']

        umapper = UMAP(random_state=123, n_components=n_components, n_neighbors=n_neighbors, a=a, b=b).fit(X_)
    
    else:
        umapper = UMAP(random_state=123).fit(X_)

    ref_umap = pd.DataFrame(umapper.transform(X_))
    ref_umap.loc[:,'label'] = le.inverse_transform(y_ref)

    new_samples_umap = pd.DataFrame(umapper.transform(X_new))
    # pred = pipe_clf.predict(X_new)
    # new_samples_umap.loc[:,'label'] = le.inverse_transform(pred)
    
    # y_pred_labels = le.inverse_transform(y_pred)
    new_samples_umap.loc[:,'label'] = y_pred

    ref_umap.loc[:,'dataset'] = 'ref'
    new_samples_umap.loc[:, 'dataset'] = 'predicted'
    total_umap = ref_umap.append(new_samples_umap)

    out_dict = {
        'total_umap': total_umap,
        'ref_umap': ref_umap,
        'new_samples_umap': new_samples_umap
    }

    return out_dict


def split_cohort_ancestry(geno_path, labels_path, out_path):
    pred_labels = pd.read_csv(labels_path, sep='\t')
    labels_list = list()
    outfiles = list()
    for label in pred_labels.label.unique():
        labels_list.append(label)
        outname = f'{out_path}_{label}'
        outfiles.append(outname)
        ancestry_group_outpath = f'{outname}.samples'
        pred_labels[pred_labels.label == label][['FID','IID']].to_csv(ancestry_group_outpath, index=False, header=False, sep='\t')

        plink_cmd = f'\
plink --bfile {geno_path} \
--keep {ancestry_group_outpath} \
--make-bed \
--out {outname}'

        shell_do(plink_cmd)
    
    output_dict = {
        'labels': labels_list,
        'paths': outfiles
    }
    
    return output_dict


def run_ancestry(geno_path, out_path, ref_panel, ref_labels, train_param_grid=None):
    step = "predict_ancestry"
    print()
    print(f"RUNNING: {step}")
    print()

    outdir = os.path.dirname(out_path)
    plot_dir = f'{outdir}/plot_ancestry'
    model_dir = f'{outdir}/models'
    pc_dir = f'{outdir}/pcs'

    # create directories if not already in existence
    os.makedirs(plot_dir, exist_ok=True)
    os.makedirs(model_dir, exist_ok=True)
    os.makedirs(pc_dir, exist_ok=True)

    calc_pcs = calculate_pcs(
        geno=geno_path,
        ref=ref_panel,
        labels=ref_labels,
        out=pc_dir,
        plot_dir=plot_dir,
        keep_temp=True
    )

    train_split = munge_training_pca_loadings(calc_pcs['labeled_ref_pca'])

    trained_clf = train_umap_classifier(
        X_train=train_split['X_train'],
        X_test=train_split['X_test'],
        y_train=train_split['y_train'],
        y_test=train_split['y_test'],
        label_encoder=train_split['label_encoder'],
        plot_dir=plot_dir,
        model_dir=model_dir,
        input_param_grid=train_param_grid
    )

    pred = predict_ancestry_from_pcs(
        projected=calc_pcs['new_samples_projected'],
        pipe_clf=trained_clf['classifier'],
        label_encoder=train_split['label_encoder'],
        out=out_path
    )

    umap_transforms = umap_transform_with_fitted(
        X_new=pred['data']['X_new'],
        X_ref=train_split['X_all'],
        y_pred=pred['data']['y_pred'],
        y_ref=train_split['y_all'],
        label_encoder=train_split['label_encoder'],
        fitted_pipe_grid=trained_clf['fitted_pipe_grid']
    )
    
#     x_min, x_max = min(umap_transforms['total_umap'].iloc[:,0]), max(umap_transforms['total_umap'].iloc[:,0])
#     y_min, y_max = min(umap_transforms['total_umap'].iloc[:,1]), max(umap_transforms['total_umap'].iloc[:,1])
#     z_min, z_max = min(umap_transforms['total_umap'].iloc[:,2]), max(umap_transforms['total_umap'].iloc[:,2])

#     x_range = [x_min-5, x_max+5]
#     y_range = [y_min-5, y_max+5]
#     z_range = [z_min-5, z_max+5]

#     plot_3d(
#         umap_transforms['total_umap'],
#         color='label',
#         symbol='dataset',
#         plot_out=f'{plot_dir}/plot_total_umap',
#         title='UMAP of New and Reference Samples',
#         x=0,
#         y=1,
#         z=2,
#         x_range=x_range,
#         y_range=y_range,
#         z_range=z_range
#     )

#     plot_3d(
#         umap_transforms['ref_umap'],
#         color='label',
#         symbol='dataset',
#         plot_out=f'{plot_dir}/plot_ref_umap',
#         title="UMAP of Reference Samples",
#         x=0,
#         y=1,
#         z=2,
#         x_range=x_range,
#         y_range=y_range,
#         z_range=z_range
#     )

#     plot_3d(
#         umap_transforms['new_samples_umap'],
#         color='label',
#         symbol='dataset',
#         plot_out=f'{plot_dir}/plot_predicted_samples_umap',
#         title='UMAP of New Samples',
#         x=0,
#         y=1,
#         z=2,
#         x_range=x_range,
#         y_range=y_range,
#         z_range=z_range
#     )

    # return more stuff as needed but for now, just need predicted labels, predicted labels out path, and predicted counts
    data_dict = {
        'predict_data': pred['data'],
        'confusion_matrix': trained_clf['confusion_matrix'],
        'ref_pcs': calc_pcs['labeled_ref_pca'],
        'projected_pcs': calc_pcs['new_samples_projected'],
        'total_umap': umap_transforms['total_umap'],
        'ref_umap': umap_transforms['ref_umap'],
        'new_samples_umap': umap_transforms['new_samples_umap'],
        'label_encoder': train_split['label_encoder']
        }

    metrics_dict = {
        'predicted_counts': pred['metrics'],
        'train_accuracy': trained_clf['train_accuracy'],
        'test_accuracy': trained_clf['test_accuracy']
        }

    outfiles_dict = {
        'predicted_labels': pred['output']
    }
    
    out_dict = {
        'step': step,
        'data': data_dict,
        'metrics': metrics_dict,
        'output': outfiles_dict
    }

    return out_dict