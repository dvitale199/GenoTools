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

#local imports
from QC.utils import shell_do, get_common_snps, rm_tmps, merge_genos
plink2 = f'{os.path.dirname(os.path.abspath(__file__))}/../exec/plink2'

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


def get_raw_files(geno_path, ref_path, labels_path, out_path):
    step = "get_raw_files"
    print()
    print(f"RUNNING: {step}")
    print()

    outdir = os.path.dirname(out_path)
    out_paths = {}

    # callrate prune ref panel and geno before getting common snps
    ref_prune_path = f'{outdir}/ref_callrate_pruned'
    ref_prune_cmd = f'plink --bfile {ref_path} --geno 0.1 --make-bed --out {ref_prune_path}'
    shell_do(ref_prune_cmd)
    out_paths['ref_pruned_bed'] = ref_prune_path

    geno_prune_path = f'{out_path}_callrate_pruned'
    geno_prune_cmd = f'plink --bfile {geno_path} --geno 0.1 --make-bed --out {geno_prune_path}'
    shell_do(geno_prune_cmd)
    out_paths['geno_pruned_bed'] = geno_prune_path

    # get common snps between ref panel and geno
    ref_common_snps = f'{outdir}/ref_common_snps'
    common_snps_files = get_common_snps(ref_prune_path, geno_prune_path, ref_common_snps)

    # add common_snps_files output paths to out_paths
    out_paths = {**out_paths, **common_snps_files}

    # get raw version of common snps - reference panel
    raw_ref_cmd = f'{plink2} --bfile {ref_common_snps} --recode A --out {ref_common_snps}'
    shell_do(raw_ref_cmd)

    # read in raw common snps
    ref_raw = pd.read_csv(f'{ref_common_snps}.raw', sep='\s+')

    # separate IDs and SNPs
    ref_ids = ref_raw[['FID','IID']]
    ref_snps = ref_raw.drop(columns=['FID', 'IID', 'PAT', 'MAT', 'SEX', 'PHENOTYPE'], axis=1)
    
    col_names = ['FID','IID'] + list(ref_snps.columns)

    # mean imputation for missing SNPs data
    mean_imp = SimpleImputer(missing_values=np.nan, strategy='mean')
    ref_snps = mean_imp.fit_transform(ref_snps)
    ref_snps = pd.DataFrame(ref_snps)

    ref_raw = pd.concat([ref_ids,ref_snps], axis=1)
    ref_raw.columns = col_names

    # read ancestry file with reference labels 
    ancestry = pd.read_csv(f'{labels_path}', sep='\t', header=None, names=['FID','IID','label'])
    ref_fam = pd.read_csv(f'{ref_path}.fam', sep=' ', header=None)
    ref_labeled = ref_fam.merge(ancestry, how='left', left_on=[0,1], right_on=['FID','IID'])

    # combined_labels
    labeled_ref_raw = ref_raw.merge(ref_labeled, how='left', on=['FID','IID'])
    labeled_ref_raw.drop(columns=[0,1,2,3,4,5],inplace=True)

    print()
    print()
    print("Labeled Reference Ancestry Counts:")
    print(labeled_ref_raw.label.value_counts())
    print()
    print()

    # get reference alleles from ref_common_snps
    ref_common_snps_ref_alleles = f'{ref_common_snps}.ref_allele'
    ref_common_snps_bim = pd.read_csv(f'{ref_common_snps}.bim', header=None, sep='\t')
    ref_common_snps_bim.columns = ['chr', 'rsid', 'kb', 'pos', 'a1', 'a2']
    ref_common_snps_bim[['rsid','a1']].to_csv(ref_common_snps_ref_alleles, sep='\t', header=False, index=False)
    out_paths['ref_alleles'] = ref_common_snps_ref_alleles

    geno_common_snps = f'{out_path}_common_snps'
    common_snps = f'{ref_common_snps}.common_snps'
    ref_common_snps_bim[['rsid']].to_csv(f'{geno_common_snps}.txt', sep='\t', header=False, index=False)
    out_paths['geno_common_snps_bed'] = geno_common_snps

    # extracting common snps
    ext_snps_cmd = f'{plink2} --bfile {geno_prune_path} --extract {common_snps} --alt1-allele {ref_common_snps_ref_alleles} --make-bed --out {geno_common_snps}'

    shell_do(ext_snps_cmd)

    # getting raw version of common snps - genotype
    raw_geno_cmd = f'{plink2} --bfile {geno_common_snps} --recode A --out {geno_common_snps}'
    shell_do(raw_geno_cmd)

    # read in raw genotypes
    raw_geno = pd.read_csv(f'{geno_common_snps}.raw', sep='\s+')

    # separate IDs and SNPs
    geno_ids = raw_geno[['FID','IID']]
    geno_snps = raw_geno.drop(columns=['FID', 'IID', 'PAT', 'MAT', 'SEX', 'PHENOTYPE'], axis=1)

    # mean imputation for missing SNPs data
    geno_snps = mean_imp.fit_transform(geno_snps)
    geno_snps = pd.DataFrame(geno_snps)

    raw_geno = pd.concat([geno_ids, geno_snps], axis=1)
    raw_geno.columns = col_names
    raw_geno['label'] = 'new'

    # merge reference panel with genotypes
    merge_common_snps = f'{out_path}_merge_ref'
    merge_genos(ref_common_snps, geno_common_snps, merge_common_snps)
    out_paths['merge_bed'] = merge_common_snps

    out_dict = {
        'raw_ref': labeled_ref_raw,
        'raw_geno': raw_geno,
        'out_paths': out_paths
    }

    return out_dict


def munge_training_data(labeled_ref_raw):
    """ Train/test split and label encode labeled PCA data for training
    
    
    """
    step = "munge_pca_loadings"
    print()
    print(f"RUNNING: {step}")
    print()

    X = labeled_ref_raw.drop(columns=['label'])
    y = labeled_ref_raw.label

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


def calculate_pcs(X_train, X_test, y_train, y_test, train_ids, test_ids, raw_geno, label_encoder, out, plot_dir):
    step = "calculate_pcs"
    print()
    print(f"RUNNING: {step}")
    print()

    out_paths = {}

    train_labels = label_encoder.inverse_transform(y_train)
    test_labels = label_encoder.inverse_transform(y_test)

    # mean and SD for flashPCA style scaling
    # paper: https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0093766
    train_mean = X_train.mean(axis=0)
    train_flash_sd = np.sqrt((train_mean/2)*(1-(train_mean/2)))

    # set up pca
    n_pcs = 50
    sk_pca = PCA(n_components=n_pcs)

    # column names for PCA
    col_names = ['PC'+str(i+1) for i in range(n_pcs)]

    # transform training data
    train_pca = transform(X_train, train_mean, train_flash_sd, sk_pca, col_names, True)
    X_train = train_pca.copy()
    train_pca['label'] = train_labels
    train_ids = train_ids.reset_index(drop=True)
    train_pca = pd.concat([train_ids, train_pca], axis=1)

    # plot_3d(train_pca, color='label', title='Reference Panel PCA - Training', plot_out=f'{plot_dir}/plot_train_skPCA', x='PC1', y='PC2', z='PC3')

    # transform testing data
    test_pca = transform(X_test, train_mean, train_flash_sd, sk_pca, col_names)
    X_test = test_pca.copy()
    test_pca['label'] = test_labels
    test_ids = test_ids.reset_index(drop=True)
    test_pca = pd.concat([test_ids, test_pca], axis=1)

    # get full reference panel pca
    ref_pca = train_pca.append(test_pca)

    # plot_3d(ref_pca, color='label', title='Reference Panel PCA - All', plot_out=f'{plot_dir}/plot_ref_skPCA', x='PC1', y='PC2', z='PC3')

    geno_ids = raw_geno[['FID','IID','label']]
    geno_snps = raw_geno.drop(columns=['FID','IID','label'], axis=1)

    # transform new samples
    projected = transform(geno_snps, train_mean, train_flash_sd, sk_pca, col_names)
    projected['label'] = geno_ids['label']

    # project new samples onto reference panel
    total_pca = ref_pca.append(projected)

    # plot_3d(total_pca, color='label', title='New Samples Projected on Reference Panel', plot_out=f'{plot_dir}/plot_projected_skPCA', x='PC1', y='PC2', z='PC3')

    projected = pd.concat([geno_ids[['FID','IID']], projected], axis=1)
    
    train_pca.to_csv(f'{out}_labeled_train_pca.txt', sep='\t', index=None)
    ref_pca.to_csv(f'{out}_labeled_ref_pca.txt', sep='\t', index=None)
    projected.to_csv(f'{out}_projected_new_pca.txt', sep='\t', index=None)
    out_paths['train_pca'] = f'{out}_labeled_train_pca.txt'
    out_paths['ref_pca'] = f'{out}_labeled_ref_pca.txt'
    out_paths['projected_pca'] = f'{out}_projected_new_pca.txt'

    out_dict = {
        'X_train': X_train,
        'X_test': X_test,
        'labeled_train_pca': train_pca,
        'labeled_ref_pca': ref_pca,
        'new_samples_projected': projected,
        'out_paths': out_paths
    }

    return out_dict


def transform(data, mean, sd, pca, col_names, fit=False):
    # flashPCA-style scaling
    data = (data-mean)/sd

    # fit-transform or transform data with PCA
    if fit:
        data_pca = pca.fit_transform(data)
    else:
        data_pca = pca.transform(data)
    
    # dataframe and named columns
    data_pca = pd.DataFrame(data_pca)
    data_pca.columns = col_names

    return data_pca


def train_umap_classifier(X_train, X_test, y_train, y_test, label_encoder, out, plot_dir, input_param_grid=None):
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
    # fig, ax = plt.subplots(figsize=(10,10))
    # sns.heatmap(pipe_clf_c_matrix, annot=True, fmt='d',
            #   xticklabels=le.inverse_transform([i for i in range(8)]), yticklabels=le.inverse_transform([i for i in range(8)]))
    # plt.ylabel('Actual')
    # plt.xlabel('Predicted')
    # plt.show()
    # fig.savefig(f'{plot_dir}/plot_umap_linearsvc_ancestry_conf_matrix.png')

    # dump best estimator to pkl
    model_path = f'{out}_umap_linearsvc_ancestry_model.pkl'
    joblib.dump(pipe_clf, model_path)

    out_dict = {
        'classifier': pipe_clf,
        'label_encoder' : le,
        'best_params': pipe_grid.best_params_,
        'confusion_matrix': pipe_clf_c_matrix,
        'fitted_pipe_grid': pipe_grid,
        'train_accuracy': train_acc,
        'test_accuracy': test_acc,
        'model_path': model_path
    }
    
    return out_dict


def load_umap_classifier(pkl_path, X_test, y_test):
    step = "load_umap_classifier"
    print()
    print(f"RUNNING: {step}")
    print()

    # load trained umap classifier from pickle file
    pipe_clf = joblib.load(pkl_path)

    # test accuracy
    test_acc = pipe_clf.score(X_test, y_test)
    print(f'Balanced Accuracy on Test Set: {test_acc}')

    # confustion matrix
    pipe_clf_pred = pipe_clf.predict(X_test)
    pipe_clf_c_matrix = metrics.confusion_matrix(y_test, pipe_clf_pred)

    out_dict = {
        'classifier': pipe_clf,
        'confusion_matrix': pipe_clf_c_matrix,
        'test_accuracy': test_acc
    }

    return out_dict


def predict_ancestry_from_pcs(projected, pipe_clf, label_encoder, out):
    step = "predict_ancestry"
    print()
    print(f"RUNNING: {step}")
    print()

    le = label_encoder

    # set new samples aside for labeling after training the model
    X_new = projected.drop(columns=['FID','IID','label'])

    # predict new samples
    y_pred = pipe_clf.predict(X_new)
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

    outfiles_dict = {
        'labels_outpath': f'{out}_umap_linearsvc_predicted_labels.txt'
    }

    out_dict = {
        'data': data_out,
        'metrics': projected.label.value_counts(),
        'output': outfiles_dict
    }

    return out_dict


def run_admixture(merged_geno_path, predicted_labels, train_pca, out_path):
    step = "run_admixture"
    print()
    print(f"RUNNING: {step}")
    print()

    # making sure FID and IID are good for merging
    predicted_labels['FID'] = predicted_labels['FID'].astype(str)
    predicted_labels['IID'] = predicted_labels['IID'].astype(str)

    # change AFR and AAC predictions to '-' for supervised admixture .pop file
    predicted_pop = predicted_labels.copy()
    predicted_pop.loc[predicted_pop['label'] == 'AAC', 'label'] = '-'
    predicted_pop.loc[predicted_pop['label'] == 'AFR', 'label'] = '-'

    # get training labels
    train_labels = train_pca[['FID','IID','label']]

    # change AAC to AFR for supervised admixture with 7 ancestry groups
    train_labels.loc[train_labels['label'] == 'AAC', 'label'] = 'AFR'

    # append train_labels to predictions
    combined_labels = predicted_pop.append(train_labels)
    combined_labels['FID'] = combined_labels['FID'].astype(str)
    combined_labels['IID'] = combined_labels['IID'].astype(str)

    # write to text file
    combined_labels_path = f'{out_path}_train_and_projected_ids.txt'
    combined_labels[['FID','IID']].to_csv(combined_labels_path, sep='\t', index=None, header=None)
    
    # plink command to keep only training set and new ids for running admixture
    keep_out = f'{out_path}_merge_train'
    keep_cmd = f'plink --bfile {merged_geno_path} --keep {combined_labels_path} --make-bed --out {keep_out}'
    shell_do(keep_cmd)

    # read in fam file to match up labels
    fam = pd.read_csv(f'{keep_out}.fam', sep='\s+', header=None)
    fam.columns = ['FID','IID','PAT','MAT','SEX','PHENO']
    fam['FID'] = fam['FID'].astype(str)
    fam['IID'] = fam['IID'].astype(str)

    # merge fam file and labels and write to .pop file for supervised admixture
    merge = fam.merge(combined_labels, how='inner', on=['FID','IID'])
    merge['label'].to_csv(f'{keep_out}.pop', sep='\t', index=None, header=None)

    # run admixture
    out_dir = os.path.split(f'{keep_out}.bed')[0]
    admixture_cmd = f'cd {out_dir} && admixture {keep_out}.bed 7 --supervised'
#     shell_do(admixture_cmd)
    # should grab exit status from here to catch errors. coming soon
    os.system(admixture_cmd)
    
    
    # read admixture results
    # admix_results = f"{keep_out.split('/')[-1]}.7.Q"
    admix_results = f"{keep_out}.7.Q"
    q_df = pd.read_csv(admix_results, sep='\s+', header=None)
    q_df.columns = [f'pop{i}' for i in range(1,8)]

    # get IDs from fam file
    q_df['FID'], q_df['IID'] = fam['FID'], fam['IID']
    q_df['FID'] = q_df['FID'].astype(str)
    q_df['IID'] = q_df['IID'].astype(str)

    # only adjust the new samples that were intially labelled AFR or AAC
    q_pop = q_df.merge(predicted_labels, left_on=['FID','IID'], right_on=['FID','IID'])
    q_pop_aac = q_pop[q_pop['label'] == 'AAC']
    q_pop_afr = q_pop[q_pop['label'] == 'AFR']
    q_pop_afr = q_pop_afr.append(q_pop_aac)
    
    # finding AFR column
    max_val = 0
    max_col = None
    for col in q_pop_afr.columns:
        if col not in ['FID','IID','label']:
            if q_pop_afr[col].mean() > max_val:
                max_col = col
    
    # making admixture adjustment
    q_pop.loc[(q_pop['label'] == 'AAC') & (q_pop[max_col] > 0.9), 'label'] = 'AFR'
    q_pop.loc[(q_pop['label'] == 'AFR') & (q_pop[max_col] < 0.9), 'label'] = 'AAC'

    adjusted_labels = q_pop[['FID','IID','label']]
    adjusted_labels_path = f'{out_path}_adjusted_labels.txt'
    adjusted_labels.to_csv(adjusted_labels_path, sep='\t', index=None)

    print()
    print('adjusted:\n', q_pop.label.value_counts())
    print()

    data_out = {
        'ids': adjusted_labels,
        'admix_proportions': q_pop
    }

    outfiles_dict = {
        'labels_outpath': adjusted_labels_path
    }

    out_dict = {
        'data': data_out,
        'metrics': q_pop.label.value_counts(),
        'output': outfiles_dict
    }

    return out_dict 


def umap_transform_with_fitted(ref_pca, X_new, y_pred, classifier=None):
    # 
    y_ref = ref_pca.loc[:,'label']
    X_ = ref_pca.drop(columns=['FID','IID','label'])

    y_pred = y_pred.drop(columns=['FID','IID'])

    # if classifier provided, use those params else, use UMAP defaults
    if classifier:
        params = classifier.get_params()

        a = params['umap__a']
        b = params['umap__b']

        n_components = params['umap__n_components']
        n_neighbors = params['umap__n_neighbors']

        umapper = UMAP(random_state=123, n_components=n_components, n_neighbors=n_neighbors, a=a, b=b).fit(X_)
    
    else:
        umapper = UMAP(random_state=123).fit(X_)

    # transform and assign labels
    ref_umap = pd.DataFrame(umapper.transform(X_))
    y_ref = y_ref.reset_index(drop=True)
    ref_umap = ref_umap.reset_index(drop=True)
    ref_umap.loc[:,'label'] = y_ref

    new_samples_umap = pd.DataFrame(umapper.transform(X_new))
    new_samples_umap.loc[:,'label'] = y_pred

    # assign dataset and get full UMAP
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


def run_ancestry(geno_path, out_path, ref_panel, ref_labels, model_path, train_param_grid=None):
    step = "predict_ancestry"
    print()
    print(f"RUNNING: {step}")
    print()

    outdir = os.path.dirname(out_path)
    plot_dir = f'{outdir}/plot_ancestry'

    # create directories if not already in existence
    # os.makedirs(plot_dir, exist_ok=True)

    raw = get_raw_files(
        geno_path=geno_path,
        ref_path=ref_panel,
        labels_path=ref_labels,
        out_path=out_path
    )

    train_split = munge_training_data(labeled_ref_raw=raw['raw_ref'])

    calc_pcs = calculate_pcs(
        X_train=train_split['X_train'], 
        X_test=train_split['X_test'],
        y_train=train_split['y_train'],
        y_test=train_split['y_test'],
        train_ids=train_split['train_ids'],
        test_ids=train_split['test_ids'],
        raw_geno=raw['raw_geno'],
        label_encoder=train_split['label_encoder'],
        out=out_path,
        plot_dir=plot_dir
    )

    if model_path:
        trained_clf = load_umap_classifier(
            pkl_path=model_path,
            X_test=calc_pcs['X_test'],
            y_test=train_split['y_test']
        )

    else:
        trained_clf = train_umap_classifier(
            X_train=calc_pcs['X_train'],
            X_test=calc_pcs['X_test'],
            y_train=train_split['y_train'],
            y_test=train_split['y_test'],
            label_encoder=train_split['label_encoder'],
            out=out_path,
            plot_dir=plot_dir
        )

    pred = predict_ancestry_from_pcs(
        projected=calc_pcs['new_samples_projected'],
        pipe_clf=trained_clf['classifier'],
        label_encoder=train_split['label_encoder'],
        out=out_path
    )

    admix = run_admixture(
        merged_geno_path=raw['out_paths']['merge_bed'],
        predicted_labels=pred['data']['ids'],
        train_pca=calc_pcs['labeled_train_pca'],
        out_path=out_path
    )

    umap_transforms = umap_transform_with_fitted(
        ref_pca=calc_pcs['labeled_ref_pca'],
        X_new=pred['data']['X_new'],
        y_pred=admix['data']['ids'],
        classifier=trained_clf['classifier']
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
        'train_pcs': calc_pcs['labeled_train_pca'],
        'ref_pcs': calc_pcs['labeled_ref_pca'],
        'projected_pcs': calc_pcs['new_samples_projected'],
        'admix_data': admix['data'],
        'total_umap': umap_transforms['total_umap'],
        'ref_umap': umap_transforms['ref_umap'],
        'new_samples_umap': umap_transforms['new_samples_umap'],
        'label_encoder': train_split['label_encoder']
    }

    metrics_dict = {
        'predicted_counts': pred['metrics'],
        'adjusted_counts': admix['metrics'],
        'test_accuracy': trained_clf['test_accuracy']
    }

    outfiles_dict = {
        'predicted_labels': pred['output'],
        'adjusted_labels': admix['output']
    }
    
    out_dict = {
        'step': step,
        'data': data_dict,
        'metrics': metrics_dict,
        'output': outfiles_dict
    }

    return out_dict
