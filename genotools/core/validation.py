"""Pre-flight validation of pipeline input (legacy-free).

Ports the validation + data-breakdown behavior of legacy utils.upfront_check.
The legacy data-driven step-skip logic is intentionally NOT ported (deferred;
see the design spec) — the current runner already discarded it.
"""

from pathlib import Path
from typing import Union

import pandas as pd

from genotools.core.exceptions import ValidationError


def validate_input(
    geno_path: Union[str, Path],
    out_path: Union[str, Path],
    skip_fails: bool = False,
) -> None:
    """Validate pipeline input and print a data breakdown.

    Raises:
        ValidationError: if the output log already exists (and not skip_fails),
            or the psam is missing SEX/PHENO1.
        FileNotFoundError: if the input pgen is missing.
    """
    geno_path = str(geno_path)
    out_path = str(out_path)

    if Path(f"{out_path}_all_logs.log").is_file() and not skip_fails:
        raise ValidationError(
            f"{out_path}_all_logs.log exists, which means the pipeline has "
            f"previously been run on this output file! Rerun with --skip_fails "
            f"to ignore this, or write output to a new file name."
        )

    if not Path(f"{geno_path}.pgen").is_file():
        raise FileNotFoundError(f"{geno_path} does not exist.")

    sam = pd.read_csv(f"{geno_path}.psam", sep=r"\s+")
    var = pd.read_csv(
        f"{geno_path}.pvar", delimiter="\t", comment="#", header=None,
        usecols=range(5), names=["#CHROM", "POS", "ID", "REF", "ALT"],
        low_memory=False,
    )

    if "SEX" not in sam.columns:
        raise ValidationError(
            f"{geno_path}.psam is missing SEX column. GenoTools requires a SEX column."
        )
    if "PHENO1" not in sam.columns:
        raise ValidationError(
            f"{geno_path}.psam is missing PHENO1 column. GenoTools requires a PHENO1 column."
        )

    sex_counts = sam["SEX"].value_counts().to_dict()
    pheno_counts = sam["PHENO1"].value_counts().to_dict()

    print("Your data has the following breakdown:")
    print("- Genetic Sex:")
    for sex in sex_counts:
        if sex == 1:
            print(f"{sex_counts[sex]} Males \n")
        if sex == 2:
            print(f"{sex_counts[sex]} Females \n")
        if sex in (0, -9):
            print(f"{sex_counts[sex]} Unknown \n")

    print("- Phenotypes:")
    for pheno in pheno_counts:
        if pheno == 2:
            print(f"{pheno_counts[pheno]} Cases \n")
        if pheno == 1:
            print(f"{pheno_counts[pheno]} Controls \n")
        if pheno in (0, -9):
            print(f"{pheno_counts[pheno]} Missing \n")
