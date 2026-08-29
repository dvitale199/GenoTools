"""Cohort splitting by predicted ancestry (legacy-free port).

Faithful port of legacy Ancestry.split_cohort_ancestry (legacy.py:907-963),
using core.executors. Genotype output is identical to the legacy version.
"""

import pandas as pd

from genotools.core.executors import run_command, get_plink2, PLINK2_PSAM_COLS


def split_cohort_by_ancestry(
    labels_path: str,
    geno_path: str,
    out_path: str,
    min_samples: int,
    subset: list | None = None,
) -> dict:
    """Split a cohort into per-ancestry pfiles based on predicted labels."""
    plink2 = get_plink2()

    pred_labels = pd.read_csv(labels_path, sep="\t")
    labels_list = []
    outfiles = []

    split_labels = subset if subset else pred_labels.label.unique()
    # Collected and concatenated once at the end. Seeding with an empty
    # DataFrame and concatenating into it would let that frame's dtype-less
    # object columns decide the result's dtypes instead of the data, which
    # differs by pandas version (object vs str under pandas 3).
    pruned_frames: list[pd.DataFrame] = []

    for label in split_labels:
        if pred_labels[pred_labels.label == label].shape[0] >= min_samples:
            labels_list.append(label)
            outname = f"{out_path}_{label}"
            outfiles.append(outname)
            ancestry_group_outpath = f"{outname}.samples"
            pred_labels[pred_labels.label == label][["FID", "IID"]].to_csv(
                ancestry_group_outpath, index=False, header=False, sep="\t"
            )

            run_command(
                f"{plink2} --pfile {geno_path} --keep {ancestry_group_outpath} "
                f"--make-pgen {PLINK2_PSAM_COLS} --out {outname}",
                tool_name="plink2",
            )

        else:
            pruned_samples_label = pred_labels[pred_labels.label == label].copy()
            pruned_samples_label["step"] = "insufficient_ancestry_sample_n"
            pruned_frames.append(
                pruned_samples_label[["FID", "IID", "step", "label"]]
            )

    pruned_samples = (
        pd.concat(pruned_frames, axis=0, ignore_index=True)
        if pruned_frames
        else pd.DataFrame(columns=["FID", "IID", "step", "label"])
    )

    return {
        "labels": labels_list,
        "paths": outfiles,
        "pruned_samples": pruned_samples,
    }
