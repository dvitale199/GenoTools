# GenoTools Train New Ancestry Prediction Models

## Overview
This documentation provides a detailed description of how to train a new model with the `GenoTools` ancestry module. There are two main use cases for training a new model:
    1. When none of the available pretrained models are suited for your data but you would like to use the provided reference panel and
    2. When you want to train a model using a new reference panel with different ancestry groups than those available through the provided reference panel.
For the second use case, please see the `prep_reference_panel.md` documentation on how to properly prepare that reference panel and labels file for use in the ancestry module.

---

### 1. Train New Model with the Provided Reference Panel
To download the provided reference panel, you can run the following command:
```
genotools-download --ref 1kg_30x_hgdp_ashk_ref_panel
```
By default, this will be downloaded to ~/.genotools/ref/ref_panel, but can be downloaded to a location of choice with the --destination flag:
```
genotools-download --ref 1kg_30x_hgdp_ashk_ref_panel --destination /path/to/desired/download/location
```

To train a new model using this downloaded reference panel, use the following command:
```
genotools \
    --pfile /path/to/genotypes/for/ancestry/prediction \
    --out /path/to/ancestry/prediction/output \
    --ancestry \
    --ref-panel /path/to/downloaded/reference/panel \
    --ref-labels /path/to/downloaded/reference/ancestry/labels
```
This command trains a model and predicts ancestry for the provided genotypes. On a full release-scale cohort it takes several hours, and the preprocessing step needs a high-memory machine (see the memory note below), so it is best run as a batch job.

The model is saved to a **directory**, `{out}_ancestry_model/`, holding `pipeline.pkl`, `common_snps.txt`, `metadata.json` and `requirements.txt`. Keep the directory intact and pass it to `--model` for future runs. (1.x wrote a single `{out}_umap_linearsvc_ancestry_model.pkl` plus a sibling `.common_snps`; 2.x cannot load that format.)

---

### 2. Train New Model with a Custom Reference Panel
Once again, please see the `prep_reference_panel.md` documentation for instructions on how to properly prep your reference panel.

To train a new model using the custom reference panel, use the following command:
```
genotools \
    --pfile /path/to/genotypes/for/ancestry/prediction \
    --out /path/to/ancestry/prediction/output \
    --ancestry \
    --ref-panel /path/to/created/reference/panel \
    --ref-labels /path/to/created/reference/ancestry/labels
```
Same output layout and the same caveats as above.

---

### 3. Refitting an existing model under new libraries

Training normally derives its SNP list from the intersection of the panel and
your cohort, which is why it needs the cohort at all. When the SNP list is
already fixed and is not what you want to change — the usual reason being a
library upgrade, since ancestry calls move with umap/sklearn versions — the
cohort is not needed, and the refit takes minutes rather than hours:

```
python tests/scripts/retrain_reference_model.py \
    --ref-panel /path/to/reference/panel \
    --ref-labels /path/to/reference/ancestry/labels \
    --snplist /path/to/existing/model/common_snps.txt \
    --out /path/to/new/model
```

Measured on the GP2 panel (4,008 samples x 43,173 SNPs): **8 minutes** with the
grid search parallelised. The refit produces no predictions and no per-cohort
diagnostics — validate it against a cohort separately — and it is not
bit-identical to the model it replaces, because the hyperparameter grid's top
is a plateau. Expect the substantive parameters to reproduce and the UMAP shape
pair to move.

---

### Memory

Ancestry preprocessing materializes the cohort as a dense 8-byte matrix, and
peak memory is about **4x** that matrix — roughly **167 GiB** at 130,000 samples
x 43,000 SNPs. Full-release *training* or *prediction* therefore needs a
high-memory machine.

2.2.0 cut the peak from 6.78x the matrix to 4.01x (a real 10k run went
26.84 -> 18.03 GiB, predictions bit-identical); 1.x still has the unimproved
path. Bounding the peak independent of cohort size is
[#292](https://github.com/dvitale199/GenoTools/issues/292).
