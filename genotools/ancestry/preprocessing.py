"""Ancestry PLINK preprocessing (legacy-free port of the Ancestry helpers).

Faithful reimplementations of legacy Ancestry.get_raw_files / clean_up and
utils.get_common_snps, using core.executors. Genotype output is identical to
the legacy versions (proven by differential tests).
"""

import os

import numpy as np
import pandas as pd

from genotools.core.executors import run_command, get_plink, get_plink2
from genotools.core.logging import get_logger

logger = get_logger(__name__)


def get_common_snps(geno_path1: str, geno_path2: str, out_name: str) -> dict:
    """Extract SNPs common to two bfiles from geno_path1. Returns output paths."""
    logger.info("Getting Common SNPs")
    plink = get_plink()
    plink2 = get_plink2()

    bim1 = pd.read_csv(f"{geno_path1}.bim", sep="\t", header=None, dtype={0: str}, low_memory=False)
    bim1.columns = ["chr", "rsid", "kb", "pos", "a1", "a2"]
    bim2 = pd.read_csv(f"{geno_path2}.bim", sep="\t", header=None, dtype={0: str}, low_memory=False)
    bim2.columns = ["chr", "rsid", "kb", "pos", "a1", "a2"]

    bim1["rsid"].to_csv(f"{geno_path1}.snplist", sep="\t", header=None, index=None)

    bim1["merge_id"] = bim1["chr"].astype(str) + ":" + bim1["pos"].astype(str) + ":" + bim1["a2"] + ":" + bim1["a1"]
    bim2["merge_id1"] = bim2["chr"].astype(str) + ":" + bim2["pos"].astype(str) + ":" + bim2["a2"] + ":" + bim2["a1"]
    bim2["merge_id2"] = bim2["chr"].astype(str) + ":" + bim2["pos"].astype(str) + ":" + bim2["a1"] + ":" + bim2["a2"]

    common_snps1 = bim2[["rsid", "merge_id1", "a1", "a2"]].merge(bim1, how="inner", left_on=["merge_id1"], right_on=["merge_id"])
    common_snps2 = bim2[["rsid", "merge_id2", "a1", "a2"]].merge(bim1, how="inner", left_on=["merge_id2"], right_on=["merge_id"])
    common_snps = pd.concat([common_snps1, common_snps2], axis=0)

    run_command(
        f"{plink} --bfile {geno_path1} --flip {geno_path1}.snplist --make-bed --out {geno_path1}_flip",
        tool_name="plink",
    )

    bim1_flip = pd.read_csv(f"{geno_path1}_flip.bim", sep="\t", header=None)
    bim1_flip.columns = ["chr", "rsid", "kb", "pos", "a1", "a2"]

    bim1_flip["merge_id"] = bim1_flip["chr"].astype(str) + ":" + bim1_flip["pos"].astype(str) + ":" + bim1_flip["a2"] + ":" + bim1_flip["a1"]
    common_snps1 = bim2[["rsid", "merge_id1", "a1", "a2"]].merge(bim1_flip, how="inner", left_on=["merge_id1"], right_on=["merge_id"])
    common_snps2 = bim2[["rsid", "merge_id2", "a1", "a2"]].merge(bim1_flip, how="inner", left_on=["merge_id2"], right_on=["merge_id"])

    common_snps = pd.concat([common_snps, common_snps1, common_snps2], axis=0)
    common_snps = common_snps.drop_duplicates(subset=["chr", "pos"], ignore_index=True)

    common_snps_file = f"{out_name}.common_snps"
    common_snps["rsid_y"].to_csv(f"{common_snps_file}", sep="\t", header=False, index=False)

    run_command(
        f"{plink2} --bfile {geno_path1} --extract {common_snps_file} --make-bed --out {out_name}",
        tool_name="plink2",
    )

    return {"common_snps": common_snps_file, "bed": out_name}


def clean_up_files(files: list) -> None:
    """Remove intermediate PLINK files by known extensions (port of Ancestry.clean_up)."""
    extensions = ["bim", "bed", "fam", "hh", "snplist", "ref_allele", "alleles", "raw"]
    for file in files:
        for ext in extensions:
            file_ext = f"{file}.{ext}"
            if os.path.exists(file_ext):
                os.remove(file_ext)


def get_raw_files(
    geno_path: str,
    ref_panel: str,
    ref_labels: str,
    out_path: str,
    train: bool,
    common_snps_file: str | None = None,
) -> dict:
    """Process reference + genotype data into labeled raw feature matrices.

    Faithful port of legacy Ancestry.get_raw_files (train + model-inference
    branches; containerized branch intentionally omitted). In inference mode,
    pass the model's common-SNP list path via `common_snps_file`.
    """
    plink2 = get_plink2()
    out_paths = {}

    # variant prune geno before getting common snps
    geno_prune_path = f"{out_path}_variant_pruned"
    run_command(
        f"{plink2} --pfile {geno_path} --geno 0.1 --max-alleles 2 --make-bed --out {geno_prune_path}",
        tool_name="plink2",
    )
    out_paths["geno_pruned_bed"] = geno_prune_path

    if train:
        ref_common_snps = f"{out_path}_umap_linearsvc_ancestry_model"
        common_snps_files = get_common_snps(ref_panel, geno_prune_path, ref_common_snps)
        out_paths = {**out_paths, **common_snps_files}
    else:
        if common_snps_file is None or not os.path.isfile(common_snps_file):
            raise FileNotFoundError(f"{common_snps_file} file does not exist.")
        ref_common_snps = f"{os.path.dirname(out_path)}/" + os.path.basename(common_snps_file).split(".")[0]
        run_command(
            f"{plink2} --bfile {ref_panel} --extract {common_snps_file} --make-bed --out {ref_common_snps}",
            tool_name="plink2",
        )
        out_paths["common_snps"] = common_snps_file
        out_paths["bed"] = ref_common_snps

    if not os.path.exists(f"{ref_common_snps}.bed"):
        raise FileNotFoundError(f"{ref_common_snps} PLINK binaries (bed/bim/fam) do not exist.")

    # get raw version of common snps - reference panel
    run_command(
        f"{plink2} --bfile {ref_common_snps} --recode A --out {ref_common_snps}",
        tool_name="plink2",
    )
    if not os.path.exists(f"{ref_common_snps}.raw"):
        raise FileNotFoundError(f"{ref_common_snps}.raw does not exist.")

    ref_raw = pd.read_csv(f"{ref_common_snps}.raw", sep=r"\s+")
    ref_ids = ref_raw[["FID", "IID"]]
    ref_snps = ref_raw.drop(columns=["FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE"], axis=1)
    ref_snps_cols = ref_snps.columns.str.extract("(.*)_")[0]
    ref_snps.columns = ref_snps_cols
    col_names = ["FID", "IID"] + list(ref_snps_cols)
    ref_raw = pd.concat([ref_ids, ref_snps], axis=1)
    ref_raw.columns = col_names

    # read ancestry file with reference labels
    ancestry = pd.read_csv(f"{ref_labels}", sep="\t", header=None, names=["FID", "IID", "label"])
    ref_fam = pd.read_csv(f"{ref_panel}.fam", sep=r"\s+", header=None)
    ref_labeled = ref_fam.merge(ancestry, how="left", left_on=[0, 1], right_on=["FID", "IID"])

    labeled_ref_raw = ref_raw.merge(ref_labeled, how="left", on=["FID", "IID"])
    labeled_ref_raw.drop(columns=[0, 1, 2, 3, 4, 5], inplace=True)

    logger.info("Labeled Reference Ancestry Counts:")
    logger.info(f"\n{labeled_ref_raw.label.value_counts()}")

    # get reference alleles from ref_common_snps
    ref_common_snps_ref_alleles = f"{ref_common_snps}.ref_allele"
    ref_common_snps_bim = pd.read_csv(f"{ref_common_snps}.bim", header=None, sep="\t")
    ref_common_snps_bim.columns = ["chr", "rsid", "kb", "pos", "a1", "a2"]
    ref_common_snps_bim[["rsid", "a1"]].to_csv(ref_common_snps_ref_alleles, sep="\t", header=False, index=False)
    out_paths["ref_alleles"] = ref_common_snps_ref_alleles

    geno_common_snps = f"{out_path}_common_snps"
    get_common_snps(geno_prune_path, ref_common_snps, geno_common_snps)

    geno_common_snps_bim = pd.read_csv(f"{geno_common_snps}.bim", sep=r"\s+", header=None)
    geno_common_snps_bim.columns = ["chr", "rsid", "kb", "pos", "a1", "a2"]

    ref_common_snps_bim["merge_id"] = ref_common_snps_bim["chr"].astype(str) + ":" + ref_common_snps_bim["pos"].astype(str)
    geno_common_snps_bim["merge_id"] = geno_common_snps_bim["chr"].astype(str) + ":" + geno_common_snps_bim["pos"].astype(str)

    merge_common_snps_bim = geno_common_snps_bim[["merge_id", "a1", "a2"]].merge(ref_common_snps_bim, how="inner", on=["merge_id"])
    merge_common_snps_bim[["chr", "rsid", "kb", "pos", "a1_x", "a2_x"]].to_csv(f"{geno_common_snps}.bim", sep="\t", header=None, index=None)

    switch = {"A": "T", "T": "A", "C": "G", "G": "C"}
    merge_common_snps_bim["a1_x_switch"] = merge_common_snps_bim["a1_x"].map(switch)
    merge_common_snps_switch = merge_common_snps_bim[
        (merge_common_snps_bim["a1_y"] != merge_common_snps_bim["a1_x"])
        & (merge_common_snps_bim["a1_y"] != merge_common_snps_bim["a1_x_switch"])
    ]
    merge_common_snps_switch[["rsid", "a2_x"]].to_csv(f"{geno_common_snps}_switch.alleles", sep="\t", header=False, index=False)

    run_command(
        f"{plink2} --bfile {geno_common_snps} --alt1-allele {geno_common_snps}_switch.alleles --recode A --out {geno_common_snps}",
        tool_name="plink2",
    )

    raw_geno = pd.read_csv(f"{geno_common_snps}.raw", sep=r"\s+")
    geno_ids = raw_geno[["FID", "IID"]]
    geno_snps = raw_geno.drop(columns=["FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE"], axis=1)
    geno_snps.columns = geno_snps.columns.str.extract("(.*)_")[0]

    # adding missing snps when not training
    missing_cols = []
    if not train:
        for col in ref_snps.columns:
            if col not in geno_snps.columns:
                missing_cols += [pd.Series(np.repeat(2, geno_snps.shape[0]), name=col)]
        if len(missing_cols) > 0:
            missing_cols = pd.concat(missing_cols, axis=1)
            geno_snps = pd.concat([geno_snps, missing_cols], axis=1)
        geno_snps = geno_snps[ref_snps.columns]

    raw_geno = pd.concat([geno_ids, geno_snps], axis=1)
    raw_geno.columns = col_names
    raw_geno["label"] = "new"

    # remove intermediate files (concat_logs dropped: structured logging replaces raw .log aggregation)
    files = [geno_prune_path, f"{geno_prune_path}_flip", f"{out_path}_common_snps_switch"]
    clean_up_files(files)

    return {
        "raw_ref": labeled_ref_raw,
        "raw_geno": raw_geno,
        "out_paths": out_paths,
    }
