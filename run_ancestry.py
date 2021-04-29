import os
from Ancestry.ancestry import calculate_pcs, munge_training_pca_loadings, train_umap_classifier, predict_ancestry_from_pcs, umap_transform_with_fitted, plot_3d

###### paths for testing ######
# geno_path = '/data/vitaled2/test_data/mcgill/MCGILL_all_call_rate_sex'
# out_path = '/data/vitaled2/test_data/mcgill/MCGILL_all_call_rate_sex_ancestry'

# ref_dir_path = '/data/LNG/vitaled2/1kgenomes'
# ref_panel = f'{ref_dir_path}/1kg_ashkj_ref_panel_gp2_pruned'
# ref_labels = f'{ref_dir_path}/ref_panel_ancestry.txt'
###############################

def run_ancestry(geno_path, out_path, ref_panel, ref_labels):
    outdir = os.path.dirname(out_path)
    plot_dir = f'{outdir}/plot_ancestry'
    model_dir = f'{outdir}/models'
    temp_dir = f'{outdir}/temp'

    # create directories if not already in existence
    os.makedirs(plot_dir, exist_ok=True)
    os.makedirs(model_dir, exist_ok=True)
    os.makedirs(temp_dir, exist_ok=True)

    calc_pcs = calculate_pcs(
        geno=geno_path,
        ref=ref_panel,
        labels=ref_labels,
        out=out_path,
        plot_dir=plot_dir,
        keep_temp=True
    )

    train_split = munge_training_pca_loadings(calc_pcs['labeled_ref_pca'])

    test_param_grid = {
            "umap__n_neighbors": [5],
            "umap__n_components": [15],
            "umap__a":[1.5],
            "umap__b": [0.25],
            "svc__C": [10**-3],
        }

    trained_clf = train_umap_classifier(
        X_train=train_split['X_train'],
        X_test=train_split['X_test'],
        y_train=train_split['y_train'],
        y_test=train_split['y_test'],
        label_encoder=train_split['label_encoder'],
        plot_dir=plot_dir,
        model_dir=model_dir,
        input_param_grid=test_param_grid
    )

    pred = predict_ancestry_from_pcs(
        projected=calc_pcs['new_samples_projected'],
        pipe_clf=trained_clf['classifier'],
        label_encoder=train_split['label_encoder'],
        out=out_path
    )

    umap_transforms = umap_transform_with_fitted(
        X_new=pred['X_new'],
        X_ref=train_split['X_all'],
        y_pred=pred['y_pred'],
        y_ref=train_split['y_all'],
        label_encoder=train_split['label_encoder'],
        fitted_pipe_grid=trained_clf['fitted_pipe_grid']
    )

    plot_3d(
        umap_transforms['total_umap'],
        color='label',
        symbol='dataset',
        plot_out=f'{plot_dir}/plot_total_umap',
        x=0,
        y=1,
        z=2
    )

    plot_3d(
        umap_transforms['ref_umap'],
        color='label',
        symbol='dataset',
        plot_out=f'{plot_dir}/plot_ref_umap',
        x=0,
        y=1,
        z=2
    )
    plot_3d(
        umap_transforms['new_samples_umap'],
        color='label',
        symbol='dataset',
        plot_out=f'{plot_dir}/plot_predicted_samples_umap',
        x=0,
        y=1,
        z=2
    )
