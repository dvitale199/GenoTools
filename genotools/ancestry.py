import os
import pandas as pd
import numpy as np
from umap import UMAP
from sklearn import preprocessing, metrics
from sklearn.impute import SimpleImputer
from sklearn.decomposition import PCA
from sklearn.model_selection import train_test_split, GridSearchCV, StratifiedKFold
from sklearn.pipeline import Pipeline
from xgboost import XGBClassifier
import pickle as pkl
import json

from genotools.utils import shell_do, get_common_snps, concat_logs
from genotools.dependencies import check_plink, check_plink2

plink_exec = check_plink()
plink2_exec = check_plink2()

class ancestry:
    def __init__(self, geno_path, ref_panel, ref_labels, out_path, model_path=None, containerized=False, singularity=False, train_param_grid=None):
        # initialize passed variables
        self.geno_path = geno_path
        self.ref_panel = ref_panel
        self.ref_labels = ref_labels
        self.out_path = out_path
        self.model_path = model_path
        self.containerized = containerized
        self.singulatity = singularity
        self.train_param_grid = train_param_grid
        #NOTE: plotting may eventually go somewhere different
        self.plot_dir = f'{os.path.dirname(out_path)}/plot_ancestry'
        os.makedirs(self.plot_dir, exist_ok=True)

        # setting train variable to false if there is a model path or containerized predictions
        if self.model_path or self.containerized:
            self.train = False
        else:
            self.train = True

        # Check that paths are set
        if not all([geno_path, ref_panel, ref_labels, out_path]):
            raise ValueError("Please make sure geno_path, ref_panel, ref_labels, and out_path are all set when initializing this class.")

        # Check path validity
        if not os.path.exists(f'{geno_path}.bed'):
            raise FileNotFoundError(f"{geno_path} does not exist.")
        elif not os.path.exists(f'{ref_panel}.bed'):
            raise FileNotFoundError(f"{ref_panel} does not exist.")
        elif not os.path.exists(f'{ref_labels}'):
            raise FileNotFoundError(f"{ref_labels} does not exist.")
        elif not (model_path or containerized or singularity or train_param_grid):
            raise Exception('Must provide a path to pkl\'d model, request containerized predictions, or train a new model with provided parametrs!')  
        elif model_path and containerized:
            raise Warning('Model path provided and containerized predictions requested! Defaulting to containerized predictions!')
    

    def get_raw_files(self):
        """
        Processes reference and genotype data for prediction, including variant pruning, extracting common SNPs, 
        obtaining raw versions of common SNPs, and preparing data with labels.

        Returns:
        dict: A dictionary containing the following keys and associated values:
            - 'raw_ref': A DataFrame with labeled raw reference data.
            - 'raw_geno': A DataFrame with labeled raw genotype data.
            - 'out_paths': A dictionary with output file paths.
        """

        step = "get_raw_files"

        outdir = os.path.dirname(self.out_path)
        out_paths = {}

        # variant prune geno before getting common snps
        geno_prune_path = f'{self.out_path}_variant_pruned'
        geno_prune_cmd = f'{plink2_exec} --pfile {self.geno_path} --geno 0.1 --make-bed --out {geno_prune_path}'
        shell_do(geno_prune_cmd)
        out_paths['geno_pruned_bed'] = geno_prune_path

        ref_common_snps = f'{outdir}/ref_common_snps'
        common_snps_file = f'{ref_common_snps}.common_snps'

        # during training get common snps between ref panel and geno
        if self.train:
            common_snps_files = get_common_snps(self.ref_panel, geno_prune_path, ref_common_snps)
            # add common_snps_files output paths to out_paths
            out_paths = {**out_paths, **common_snps_files}
        # otherwise extract common snps from training
        else:
            if not os.path.exists(f'{common_snps_file}'):
                raise FileNotFoundError(f"{common_snps_file} does not exist.")
            else:
                extract_cmd = f'{plink2_exec} --bfile {self.ref_panel} --extract {common_snps_file} --make-bed --out {ref_common_snps}'
                shell_do(extract_cmd)

                listOfFiles = [f'{ref_common_snps}.log']
                concat_logs(step, self.out_path, listOfFiles)
    
                # add to out_paths (same as common_snps_files)
                out_paths['common_snps'] = common_snps_file
                out_paths['bed'] = ref_common_snps
                
        if not os.path.exists(f'{ref_common_snps}.bed'):
            raise FileNotFoundError(f"{ref_common_snps} PLINK binaries (bed/bim/fam) do not exist.")

        # get raw version of common snps - reference panel
        raw_ref_cmd = f'{plink2_exec} --bfile {ref_common_snps} --recode A --out {ref_common_snps}'
        shell_do(raw_ref_cmd)

        if not os.path.exists(f'{ref_common_snps}.raw'):
            raise FileNotFoundError(f"{ref_common_snps}.raw does not exist.")

        # read in raw common snps
        ref_raw = pd.read_csv(f'{ref_common_snps}.raw', sep='\s+')

        # separate IDs and snps
        ref_ids = ref_raw[['FID','IID']]
        ref_snps = ref_raw.drop(columns=['FID', 'IID', 'PAT', 'MAT', 'SEX', 'PHENOTYPE'], axis=1)
        
        # change snp column names to avoid sklearn warning/future error
        ref_snps_cols = ref_snps.columns.str.extract('(.*)_')[0]
        ref_snps.columns = ref_snps_cols

        # full col names
        col_names = ['FID','IID'] + list(ref_snps_cols)

        ref_raw = pd.concat([ref_ids,ref_snps], axis=1)
        ref_raw.columns = col_names

        # read ancestry file with reference labels 
        ancestry = pd.read_csv(f'{self.ref_labels}', sep='\t', header=None, names=['FID','IID','label'])
        ref_fam = pd.read_csv(f'{self.ref_panel}.fam', sep='\s+', header=None)
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

        geno_common_snps = f'{self.out_path}_common_snps'

        geno_common_snps_files = get_common_snps(geno_prune_path, ref_common_snps, geno_common_snps)

        # read geno common snps bim file
        geno_common_snps_bim = pd.read_csv(f'{geno_common_snps}.bim', sep='\s+', header=None)
        geno_common_snps_bim.columns = ['chr', 'rsid', 'kb', 'pos', 'a1', 'a2']
        
        # make chr:pos merge ids
        ref_common_snps_bim['merge_id'] = ref_common_snps_bim['chr'].astype(str) + ':' + ref_common_snps_bim['pos'].astype(str)
        geno_common_snps_bim['merge_id'] = geno_common_snps_bim['chr'].astype(str) + ':' + geno_common_snps_bim['pos'].astype(str)

        # merge and write over geno common snps files so snp ids match
        merge_common_snps_bim = geno_common_snps_bim[['merge_id','a1','a2']].merge(ref_common_snps_bim, how='inner', on=['merge_id'])
        merge_common_snps_bim[['chr','rsid','kb','pos','a1_x','a2_x']].to_csv(f'{geno_common_snps}.bim', sep='\t', header=None, index=None)

        # dictionary to switch alleles
        switch = {'A':'T','T':'A','C':'G','G':'C'}

        # finding alleles to be switched
        merge_common_snps_bim['a1_x_switch'] = merge_common_snps_bim['a1_x'].map(switch)
        merge_common_snps_switch = merge_common_snps_bim[(merge_common_snps_bim['a1_y'] != merge_common_snps_bim['a1_x']) & (merge_common_snps_bim['a1_y'] != merge_common_snps_bim['a1_x_switch'])]
        merge_common_snps_switch[['rsid','a2_x']].to_csv(f'{geno_common_snps}_switch.alleles', sep='\t', header=False, index=False)

        # getting raw version of common snps and setting alleles - genotype
        raw_geno_cmd = f'{plink2_exec} --bfile {geno_common_snps} --alt1-allele {geno_common_snps}_switch.alleles --recode A --out {geno_common_snps}'
        shell_do(raw_geno_cmd)

        # read in raw genotypes
        raw_geno = pd.read_csv(f'{geno_common_snps}.raw', sep='\s+')

        # separate IDs and SNPs
        geno_ids = raw_geno[['FID','IID']]
        geno_snps = raw_geno.drop(columns=['FID', 'IID', 'PAT', 'MAT', 'SEX', 'PHENOTYPE'], axis=1)

        # change col names to match ref
        geno_snps.columns = geno_snps.columns.str.extract('(.*)_')[0]

        # adding missing snps when not training
        missing_cols = []
        if not self.train:
            for col in ref_snps.columns:
                if col not in geno_snps.columns:
                    missing_cols += [pd.Series(np.repeat(2, geno_snps.shape[0]), name=col)]
            
            if len(missing_cols) > 0:
                missing_cols = pd.concat(missing_cols, axis=1)
                geno_snps = pd.concat([geno_snps, missing_cols], axis=1)
            # reordering columns to match ref
            geno_snps = geno_snps[ref_snps.columns]

        raw_geno = pd.concat([geno_ids, geno_snps], axis=1)
        raw_geno.columns = col_names
        raw_geno['label'] = 'new'

        # concat logs
        listOfFiles = [f'{geno_prune_path}.log', f'{geno_prune_path}_flip.log', f'{ref_common_snps}.log', f'{geno_common_snps}.log']
        concat_logs(step, self.out_path, listOfFiles)
        
        # remove intermediate files
        extensions = ['bim', 'bed', 'fam', 'hh', 'snplist', 'ref_allele', 'alleles', 'raw', 'common_snps']
        files = [geno_prune_path, ref_common_snps, f'{geno_prune_path}_flip', f'{self.out_path}_common_snps',
                 f'{self.out_path}_common_snps_switch']

        for file in files:
            for ext in extensions:
                file_ext = f'{file}.{ext}'
                if os.path.exists(file_ext):
                    os.remove(file_ext)

        out_dict = {
            'raw_ref': labeled_ref_raw,
            'raw_geno': raw_geno,
            'out_paths': out_paths
        }

        return out_dict


    def munge_training_data(self, labeled_ref_raw):
        """
        Preprocesses labeled raw data for prediction by performing train/test split, label encoding, and data formatting.

        Args:
        labeled_ref_raw (DataFrame): Labeled raw reference data with PCA information.

        Returns:
        dict: A dictionary containing the following keys and associated values:
            - 'X_train': Features of the training set.
            - 'X_test': Features of the test set.
            - 'y_train': Labels of the training set.
            - 'y_test': Labels of the test set.
            - 'train_ids': IDs of samples in the training set.
            - 'test_ids': IDs of samples in the test set.
            - 'label_encoder': A fitted LabelEncoder for label transformation.
            - 'X_all': All features (including both training and test sets).
            - 'y_all': All labels (including both training and test sets).
        """
       
        step = "munge_pca_loadings"

        X = labeled_ref_raw.drop(columns=['label'])
        y = labeled_ref_raw.label

        #encode labels
        le = preprocessing.LabelEncoder()
        y = le.fit_transform(y)

        # train/test split 1kg pca data
        X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=123, stratify=y)

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


    def calculate_pcs(self, X_train, X_test, y_train, y_test, train_ids, test_ids, raw_geno, label_encoder):
        """
        Calculates principal components (PCs) for training and testing datasets, transforms the data,
        and projects new samples onto the reference panel.

        Args:
        X_train (DataFrame): Features of the training set.
        X_test (DataFrame): Features of the test set.
        y_train (Series): Labels of the training set.
        y_test (Series): Labels of the test set.
        train_ids (DataFrame): IDs of samples in the training set.
        test_ids (DataFrame): IDs of samples in the test set.
        raw_geno (DataFrame): Raw genotype data.
        label_encoder (LabelEncoder): A fitted LabelEncoder for label transformation.

        Returns:
        dict: A dictionary containing the following keys and associated values:
            - 'X_train': Transformed features of the training set.
            - 'X_test': Transformed features of the test set.
            - 'labeled_train_pca': Labeled PCA data for the training set.
            - 'labeled_ref_pca': Labeled PCA data for the reference panel.
            - 'new_samples_projected': Projected PCA data for new samples.
            - 'out_paths': Dictionary containing output file paths.
        """

        step = "calculate_pcs"

        out_paths = {}

        train_labels = label_encoder.inverse_transform(y_train)
        test_labels = label_encoder.inverse_transform(y_test)

        mean_imp = SimpleImputer(missing_values=np.nan, strategy='mean')
        X_train = mean_imp.fit_transform(X_train)
        X_test = mean_imp.transform(X_test)

        # mean and SD for flashPCA style scaling
        # paper: https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0093766
        train_mean = X_train.mean(axis=0)
        train_flash_sd = np.sqrt((train_mean/2)*(1-(train_mean/2)))

        # set up pca
        n_pcs = 50
        sk_pca = PCA(n_components=n_pcs, svd_solver='full')

        # column names for PCA
        col_names = ['PC'+str(i+1) for i in range(n_pcs)]

        # transform training data
        train_pca = self.transform(X_train, train_mean, train_flash_sd, sk_pca, col_names, True)
        X_train = train_pca.copy()
        train_pca['label'] = train_labels
        train_ids = train_ids.reset_index(drop=True)
        train_pca = pd.concat([train_ids, train_pca], axis=1)

        # create df for eigenvalues and explained variance ratio
        ev = pd.Series(sk_pca.explained_variance_)
        evr = pd.Series(sk_pca.explained_variance_ratio_)
        ev_df = pd.concat([pd.Series(col_names), ev, evr], axis=1)
        ev_df.columns = ['PC','eigenvalue','explained_variance_ratio']
        ev_df.to_csv(f'{self.out_path}_pca_eigenvalues.txt', sep='\t', index=False)

        # plot_3d(train_pca, color='label', title='Reference Panel PCA - Training', plot_out=f'{plot_dir}/plot_train_skPCA', x='PC1', y='PC2', z='PC3')

        # transform testing data
        test_pca = self.transform(X_test, train_mean, train_flash_sd, sk_pca, col_names)
        X_test = test_pca.copy()
        test_pca['label'] = test_labels
        test_ids = test_ids.reset_index(drop=True)
        test_pca = pd.concat([test_ids, test_pca], axis=1)

        # get full reference panel pca
        ref_pca = pd.concat([train_pca, test_pca], ignore_index=True)

        # plot_3d(ref_pca, color='label', title='Reference Panel PCA - All', plot_out=f'{plot_dir}/plot_ref_skPCA', x='PC1', y='PC2', z='PC3')

        geno_ids = raw_geno[['FID','IID','label']]
        geno_snps = raw_geno.drop(columns=['FID','IID','label'], axis=1)

        geno_snps = mean_imp.transform(geno_snps)

        # transform new samples
        projected = self.transform(geno_snps, train_mean, train_flash_sd, sk_pca, col_names)
        projected['label'] = geno_ids['label']

        # project new samples onto reference panel
        total_pca = pd.concat([ref_pca, projected])

        # plot_3d(total_pca, color='label', title='New Samples Projected on Reference Panel', plot_out=f'{plot_dir}/plot_projected_skPCA', x='PC1', y='PC2', z='PC3')

        projected = pd.concat([geno_ids[['FID','IID']], projected], axis=1)
        
        train_pca.to_csv(f'{self.out_path}_labeled_train_pca.txt', sep='\t', index=None)
        ref_pca.to_csv(f'{self.out_path}_labeled_ref_pca.txt', sep='\t', index=None)
        projected.to_csv(f'{self.out_path}_projected_new_pca.txt', sep='\t', index=None)
        out_paths['train_pca'] = f'{self.out_path}_labeled_train_pca.txt'
        out_paths['ref_pca'] = f'{self.out_path}_labeled_ref_pca.txt'
        out_paths['projected_pca'] = f'{self.out_path}_projected_new_pca.txt'

        out_dict = {
            'X_train': X_train,
            'X_test': X_test,
            'labeled_train_pca': train_pca,
            'labeled_ref_pca': ref_pca,
            'new_samples_projected': projected,
            'out_paths': out_paths
        }

        return out_dict


    def transform(self, data, mean, sd, pca, col_names, fit=False):
        """
        Applies flashPCA-style scaling and PCA transformation to the input data.

        Args:
        data (DataFrame): Input data to be transformed.
        mean (Series): Mean values for flashPCA-style scaling.
        sd (Series): Standard deviation values for flashPCA-style scaling.
        pca (PCA): PCA model for transforming data.
        col_names (list): List of column names for the transformed data.
        fit (bool, optional): If True, fit-transform the data with PCA. If False, only transform. Default is False.

        Returns:
        DataFrame: Transformed data with named columns.
        """

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


    def train_umap_classifier(self, X_train, X_test, y_train, y_test, label_encoder):
        """
        Train a UMAP to linear XGBoost classifier pipeline.

        Args:
        X_train (DataFrame): Training features.
        X_test (DataFrame): Testing features.
        y_train (Series): Training labels.
        y_test (Series): Testing labels.
        label_encoder: LabelEncoder object for encoding and decoding labels.

        Returns:
        dict: Dictionary containing classifier, label encoder, parameters, confusion matrix, fitted grid, train accuracy, test accuracy, and model path.
        """

        step = "train_umap_classifier"

        if self.train_param_grid:
            param_grid = self.train_param_grid
        else:
            param_grid = {
            "umap__n_neighbors": [5, 20],
            "umap__n_components": [15, 25],
            "umap__a":[0.75, 1.0, 1.5],
            "umap__b": [0.25, 0.5, 0.75],
            "xgb__lambda": [10**i for i in range(-3,3)],
        }

        le = label_encoder

        # Transformation with UMAP followed by classification with svc
        umap = UMAP(random_state=123)
        
        xgb = XGBClassifier(booster='gblinear', random_state=123)
        pipeline = Pipeline([("umap", umap), ("xgb", xgb)])
        
        cross_validation = StratifiedKFold(n_splits=5, shuffle=True, random_state=123)
        pipe_grid = GridSearchCV(pipeline, param_grid, cv=cross_validation, scoring='balanced_accuracy')
        pipe_grid.fit(X_train, y_train)

        results_df = pd.DataFrame(pipe_grid.cv_results_)
        top_results = results_df[results_df['rank_test_score'] == 1]

        train_acc = pipe_grid.best_score_
        print(f'Training Balanced Accuracy: {train_acc}')

        interval = 1.96 * float(top_results['std_test_score'].iloc[0])
        print(f'Training Balanced Accuracy; 95% CI: ({train_acc-interval}, {train_acc+interval})')

        print(f'Best Parameters: {pipe_grid.best_params_}')

        pipe_clf = pipe_grid.best_estimator_
        test_acc = pipe_clf.score(X_test, y_test)
        print(f"Balanced Accuracy on Test Set: {test_acc}")

        margin_of_error = 1.96 * np.sqrt((test_acc * (1-test_acc)) / np.shape(y_test)[0])
        print(f"Balanced Accuracy on Test Set, 95% Confidence Interval: ({test_acc-margin_of_error}, {test_acc+margin_of_error})")

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
        model_path = f'{self.out_path}_umap_linearsvc_ancestry_model.pkl'
        model_file = open(model_path, 'wb')
        pkl.dump(pipe_clf, model_file)
        model_file.close()

        out_dict = {
            'classifier': pipe_clf,
            'label_encoder' : le,
            'params': pipe_grid.best_params_,
            'confusion_matrix': pipe_clf_c_matrix,
            'fitted_pipe_grid': pipe_grid,
            'train_accuracy': train_acc,
            'test_accuracy': test_acc,
            'model_path': model_path
        }
        
        return out_dict
    

    def load_umap_classifier(self, X_test, y_test):
        """
        Load a trained UMAP classifier from a pickle file and evaluate its performance on the test set.

        Args:
        X_test (DataFrame): Testing features.
        y_test (Series): Testing labels.

        Returns:
        dict: Dictionary containing classifier, confusion matrix, test accuracy, and model parameters.
        """

        step = "load_umap_classifier"

        # load trained umap classifier from pickle file
        pkl_in = open(self.model_path, 'rb')
        pipe_clf = pkl.load(pkl_in)
        pkl_in.close()

        # test accuracy
        test_acc = pipe_clf.score(X_test, y_test)
        print(f'Balanced Accuracy on Test Set: {test_acc}')

        margin_of_error = 1.96 * np.sqrt((test_acc * (1-test_acc)) / np.shape(y_test)[0])
        print(f"Balanced Accuracy on Test Set, 95% Confidence Interval: ({test_acc-margin_of_error}, {test_acc+margin_of_error})")

        # confustion matrix
        pipe_clf_pred = pipe_clf.predict(X_test)
        pipe_clf_c_matrix = metrics.confusion_matrix(y_test, pipe_clf_pred)

        # parameters
        params = pipe_clf.get_params()

        out_dict = {
            'classifier': pipe_clf,
            'confusion_matrix': pipe_clf_c_matrix,
            'test_accuracy': test_acc,
            'params': params
        }

        return out_dict


    def predict_ancestry_from_pcs(self, projected, pipe_clf, label_encoder):
        """
        Predict ancestry labels for new samples based on their projected principal components.

        Args:
        projected (DataFrame): Dataframe containing projected principal components of new samples.
        pipe_clf: Trained classifier pipeline.
        label_encoder: Label encoder used for encoding ancestry labels.

        Returns:
        dict: Dictionary containing predicted labels, output data, and metrics.
        """

        step = "predict_ancestry"

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

        projected[['FID','IID','label']].to_csv(f'{self.out_path}_umap_linearsvc_predicted_labels.txt', sep='\t', index=False)

        data_out = {
            'ids': projected.loc[:,['FID','IID','label']],
            'X_new': X_new,
            'y_pred': ancestry_pred,
            'label_encoder': le
        }

        outfiles_dict = {
            'labels_outpath': f'{self.out_path}_umap_linearsvc_predicted_labels.txt'
        }

        out_dict = {
            'data': data_out,
            'metrics': projected.label.value_counts(),
            'output': outfiles_dict
        }

        return out_dict


    def get_containerized_predictions(self, X_test, y_test, projected, label_encoder):
        """
        Get predictions using a containerized environment for UMAP and linear SVM classifier.

        Args:
        X_test (DataFrame): Test data.
        y_test (Series): True labels for the test data.
        projected (DataFrame): Projected principal components of new samples.
        label_encoder: Label encoder used for encoding ancestry labels.

        Returns:
        tuple: Two dictionaries containing trained classifier results and prediction results.
        """

        container_dir = f'{os.path.dirname(__file__)}/container'

        # write test data and projections to txt
        X_test.to_csv(f'{container_dir}/X_test.txt', sep='\t', index=False)
        pd.Series(y_test).to_csv(f'{container_dir}/y_test.txt', sep='\t', index=False, header=False)
        projected.to_csv(f'{container_dir}/projected.txt', sep='\t', index=False)

        if self.singulatity:
            shell_do(f'singularity pull {container_dir}/get_predictions.sif docker://mkoretsky1/genotools_ancestry:python3.8')
            shell_do(f'singularity run --bind {container_dir}:/app {container_dir}/get_predictions.sif')
            os.remove(f'{container_dir}/get_predictions.sif')
        else:
            shell_do(f'docker pull mkoretsky1/genotools_ancestry:python3.8')
            shell_do(f'docker run -v {container_dir}:/app --name get_predictions mkoretsky1/genotools_ancestry:python3.8')
            shell_do(f'docker rm get_predictions')

        # test accuracy
        accuracy_dict_path = f'{container_dir}/accuracy.json'
        with open(accuracy_dict_path, 'r') as f:
            accuracy_dict = json.load(f)
            f.close()
        
        test_acc = accuracy_dict['test_acc']
        print(f'Balanced Accuracy on Test Set: {test_acc}')

        margin_of_error = accuracy_dict['margin_of_error']
        print(f"Balanced Accuracy on Test Set, 95% Confidence Interval: ({test_acc-margin_of_error}, {test_acc+margin_of_error})")
        
        # confusion matrix
        y_test = pd.read_csv(f'{container_dir}/y_test.txt', sep='\t', header=None)
        pipe_clf_pred = pd.read_csv(f'{container_dir}/pipe_clf_pred.txt', sep='\t', header=None)
        pipe_clf_c_matrix = metrics.confusion_matrix(y_test, pipe_clf_pred)

        # parameters
        params_path = f'{container_dir}/params.json'
        with open(params_path) as f:
            params = json.load(f)
            f.close()

        trained_clf_out_dict = {
            'confusion_matrix': pipe_clf_c_matrix,
            'test_accuracy': test_acc,
            'params': params
        }

        le = label_encoder

        X_new = projected.drop(columns=['FID','IID','label'])

        # predicted new samples
        predict_path = f'{container_dir}/predicted_labels.txt'
        y_pred = pd.read_csv(predict_path, sep='\t', header=None)
        ancestry_pred = le.inverse_transform(y_pred)
        projected.loc[:,'label'] = ancestry_pred

        print()
        print('predicted:\n', projected.label.value_counts())
        print()

        projected[['FID','IID','label']].to_csv(f'{self.out_path}_umap_linearsvc_predicted_labels.txt', sep='\t', index=False)

        # remove created files
        files = ['X_test.txt','y_test.txt','projected.txt','accuracy.json','pipe_clf_pred.txt','predicted_labels.txt','params.json']

        for file in files:
            os.remove(f'{container_dir}/{file}')

        data_out = {
            'ids': projected.loc[:,['FID','IID','label']],
            'X_new': X_new,
            'y_pred': ancestry_pred,
            'label_encoder': le
        }

        outfiles_dict = {
            'labels_outpath': f'{self.out_path}_umap_linearsvc_predicted_labels.txt'
        }

        pred_out_dict = {
            'data': data_out,
            'metrics': projected.label.value_counts(),
            'output': outfiles_dict
        }

        return trained_clf_out_dict, pred_out_dict


    def umap_transform_with_fitted(self, ref_pca, X_new, y_pred, params=None):
        """
        Transform data using a fitted UMAP model.

        Args:
        ref_pca (DataFrame): Reference PCA data with labels.
        X_new (DataFrame): New samples to be transformed.
        y_pred (DataFrame): Predicted labels for new samples.
        params (dict, optional): UMAP parameters. Defaults to None.

        Returns:
        dict: Dictionary containing UMAP-transformed data.
        """

        step = 'umap_transform'

        y_ref = ref_pca.loc[:,'label']
        X_ = ref_pca.drop(columns=['FID','IID','label'])

        y_pred = y_pred.drop(columns=['FID','IID'])

        # if params provided, use those else, use UMAP defaults
        if params:
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
        total_umap = pd.concat([ref_umap, new_samples_umap], ignore_index=True)

        out_dict = {
            'total_umap': total_umap,
            'ref_umap': ref_umap,
            'new_samples_umap': new_samples_umap
        }

        return out_dict
    

    def split_cohort_ancestry(self, labels_path, subset=False):
        """
        Split a cohort based on predicted ancestries.

        Args:
        labels_path (str): Path to the file containing predicted labels.
        subset (list, optional): List of ancestries to continue analysis for. Defaults to False.

        Returns:
        dict: Dictionary containing labels and paths for each split ancestry.
        """

        step = 'split_cohort_ancestry'

        pred_labels = pd.read_csv(labels_path, sep='\t')
        labels_list = list()
        outfiles = list()

        # subset is a list of ancestries to continue analysis for passed by the user
        if subset:
            split_labels = subset
        else:
            split_labels = pred_labels.label.unique()

        listOfFiles = []
        for label in split_labels:
            labels_list.append(label)
            outname = f'{self.out_path}_{label}'
            outfiles.append(outname)
            ancestry_group_outpath = f'{outname}.samples'
            pred_labels[pred_labels.label == label][['FID','IID']].to_csv(ancestry_group_outpath, index=False, header=False, sep='\t')

            plink_cmd = plink_cmd = f'{plink2_exec} --pfile {self.geno_path} --keep {ancestry_group_outpath} --make-pgen psam-cols=fid,parents,sex,phenos --out {outname}'
            shell_do(plink_cmd)

            listOfFiles.append(f'{outname}.log')
            
        concat_logs(step, self.out_path, listOfFiles)

        output_dict = {
            'labels': labels_list,
            'paths': outfiles
        }


    def run_ancestry(self):
        """
        Run the ancestry prediction pipeline.

        Returns:
        dict: Dictionary containing data, metrics, and output information.
        """

        step = "predict_ancestry"

        raw = self.get_raw_files()

        train_split = self.munge_training_data(labeled_ref_raw=raw['raw_ref'])

        calc_pcs = self.calculate_pcs(
            X_train=train_split['X_train'], 
            X_test=train_split['X_test'],
            y_train=train_split['y_train'],
            y_test=train_split['y_test'],
            train_ids=train_split['train_ids'],
            test_ids=train_split['test_ids'],
            raw_geno=raw['raw_geno'],
            label_encoder=train_split['label_encoder'],
        )

        if self.containerized == False:
            if self.model_path:
                trained_clf = self.load_umap_classifier(
                    X_test=calc_pcs['X_test'],
                    y_test=train_split['y_test']
                )

            else:
                trained_clf = self.train_umap_classifier(
                    X_train=calc_pcs['X_train'],
                    X_test=calc_pcs['X_test'],
                    y_train=train_split['y_train'],
                    y_test=train_split['y_test'],
                    label_encoder=train_split['label_encoder'],
                )

            pred = self.predict_ancestry_from_pcs(
                projected=calc_pcs['new_samples_projected'],
                pipe_clf=trained_clf['classifier'],
                label_encoder=train_split['label_encoder'],
            )
        
        else:
            trained_clf, pred = self.get_containerized_predictions(
                X_test=calc_pcs['X_test'],
                y_test=train_split['y_test'],
                projected=calc_pcs['new_samples_projected'],
                label_encoder=train_split['label_encoder'],
            )

        umap_transforms = self.umap_transform_with_fitted(
            ref_pca=calc_pcs['labeled_ref_pca'],
            X_new=pred['data']['X_new'],
            y_pred=pred['data']['ids'],
            params=trained_clf['params']
        )

        #NOTE: Just copying over here for the sake of having everything, figure out plotting later
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


        data_dict = {
            'predict_data': pred['data'],
            'confusion_matrix': trained_clf['confusion_matrix'],
            'train_pcs': calc_pcs['labeled_train_pca'],
            'ref_pcs': calc_pcs['labeled_ref_pca'],
            'projected_pcs': calc_pcs['new_samples_projected'],
            'total_umap': umap_transforms['total_umap'],
            'ref_umap': umap_transforms['ref_umap'],
            'new_samples_umap': umap_transforms['new_samples_umap'],
            'label_encoder': train_split['label_encoder']
        }

        metrics_dict = {
            'predicted_counts': pred['metrics'],
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