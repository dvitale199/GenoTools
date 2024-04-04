# GenoTools Train New Ancestry Prediction Models

## Overview
This documentation provides a detailed description of how to train a new model with the `GenoTools` ancestry module. There are two main use cases for training a new model:
    1. When none of the available pretrained models are suited for your data but you would like to use the provided reference panel.
    2. When you want to train a model using a new reference panel with different ancestry groups than those available through the provided reference panel.
For the second use case, please see the `prep_reference_panel.md` documentation on how to properly prepare that reference panel for use in the ancestry module.

---

### 1. Train New Model with the Provided Reference Panel
To download the provided reference panel, you can run the following command:
```genotools-download --ref 1kg_30x_hgdp_ashk_ref_panel```
By default, this will be downloaded to ~/.genotools/ref/ref_panel, but can be downloaded to a location of choice with the --destination flag:
```genotools-download --ref 1kg_30x_hgdp_ashk_ref_panel --destination /path/to/desired/download/location```

To train a new model using this downloaded reference panel, use the following command:
```
genotools \
    --pfile /path/to/genotypes/for/ancestry/prediction \
    --out /path/to/ancestry/prediction/output \
    --ancestry \
    --ref_panel /path/to/downloaded/reference/panel \
    --ref_labels /path/to/reference/ancestry/labels \
```
This command will train a model and render predictions for the provided genotypes.
The model will be saved to: `{out}__umap_linearsvc_ancestry_model.pkl` and the SNPs used to train the model will be saved to `{out}__umap_linearsvc_ancestry_model.common_snps`. Please save these files and keep them in the same directory for future use on your genotypes!

---

### 2. Train New Model with a Custom Reference Panel
Once again, please see the `prep_reference_panel.md` documentation for instructions on how to properly prep your reference panel.

To train a new model using the custom reference panel, use the following command:
```
genotools \
    --pfile /path/to/genotypes/for/ancestry/prediction \
    --out /path/to/ancestry/prediction/output \
    --ancestry \
    --ref_panel /path/to/downloaded/reference/panel \
    --ref_labels /path/to/reference/ancestry/labels \
```
This command will train a model and render predictions for the provided genotypes.
The model will be saved to: `{out}__umap_linearsvc_ancestry_model.pkl` and the SNPs used to train the model will be saved to `{out}__umap_linearsvc_ancestry_model.common_snps`. Please save these files and keep them in the same directory for future use on your genotypes!