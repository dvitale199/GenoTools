"""Ancestry PLINK preprocessing (legacy-free port of the Ancestry helpers).

Faithful reimplementations of legacy Ancestry.get_raw_files / clean_up and
utils.get_common_snps, using core.executors. Genotype output is identical to
the legacy versions (proven by differential tests).
"""

import os

import numpy as np
import pandas as pd

from genotools.core.executors import run_command, get_plink, get_plink2


def get_common_snps(geno_path1: str, geno_path2: str, out_name: str) -> dict:
    """Extract SNPs common to two bfiles from geno_path1. Returns output paths."""
    print("Getting Common SNPs")
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
