"""Pre-flight validation of pipeline input (legacy-free).

Ports the validation + data-breakdown behavior of legacy utils.upfront_check,
including the data-driven step-skip decisions.
"""

from dataclasses import dataclass
from pathlib import Path
from typing import Union

import pandas as pd

from .exceptions import ValidationError
from .logging import get_logger

logger = get_logger(__name__)


@dataclass(frozen=True)
class ValidationDecisions:
    """Data-driven step-skip decisions derived from the input breakdown."""
    skip_sex: bool = False
    skip_case_control: bool = False
    skip_het: bool = False
    disable_filter_controls: bool = False


def guard_output_not_exists(
    out_path: Union[str, Path], skip_fails: bool = False
) -> None:
    """Guard against re-running over a prior run's output.

    Raises ``ValidationError`` if ``{out_path}_all_logs.log`` already exists and
    ``skip_fails`` is False. Decoupled from ``validate_input`` so the runner can
    call it *before* logging is set up (logging setup creates that log file).
    """
    if Path(f"{out_path}_all_logs.log").is_file() and not skip_fails:
        raise ValidationError(
            f"{out_path}_all_logs.log exists, which means the pipeline has "
            f"previously been run on this output file! Rerun with --skip-fails "
            f"to ignore this (the existing log is preserved as "
            f"{Path(out_path).name}_all_logs.log.N), or write output to a new "
            f"file name."
        )


def validate_input(
    geno_path: Union[str, Path],
    out_path: Union[str, Path],
    skip_fails: bool = False,
    *,
    sex_requested: bool = False,
    het_requested: bool = False,
    hwe_requested: bool = False,
    filter_controls: bool = False,
    case_control_requested: bool = False,
) -> ValidationDecisions:
    """Validate pipeline input, log a data breakdown, and derive step-skip decisions.

    The re-run guard (output already exists) lives in ``guard_output_not_exists``
    so it can run before logging is configured; this function is pure validation.

    Raises:
        ValidationError: if the psam is missing SEX/PHENO1.
        FileNotFoundError: if the input pgen is missing.
    """
    geno_path = str(geno_path)
    out_path = str(out_path)

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

    logger.info("Your data has the following breakdown:")
    logger.info("- Genetic Sex:")
    for sex in sex_counts:
        if sex == 1:
            logger.info(f"  {sex_counts[sex]} Males")
        if sex == 2:
            logger.info(f"  {sex_counts[sex]} Females")
        if sex in (0, -9):
            logger.info(f"  {sex_counts[sex]} Unknown")

    logger.info("- Phenotypes:")
    for pheno in pheno_counts:
        if pheno == 2:
            logger.info(f"  {pheno_counts[pheno]} Cases")
        if pheno == 1:
            logger.info(f"  {pheno_counts[pheno]} Controls")
        if pheno in (0, -9):
            logger.info(f"  {pheno_counts[pheno]} Missing")

    chr_counts = var["#CHROM"].value_counts().to_dict()

    if skip_fails:
        return ValidationDecisions()

    skip_sex = False
    if sex_requested:
        if (1 not in sex_counts) and (2 not in sex_counts):
            logger.warning(
                "You tried calling sex prune but no sample sex data is available. "
                "Skipping...")
            skip_sex = True
        elif ("23" not in chr_counts) and ("X" not in chr_counts):
            logger.warning(
                "You tried calling sex prune but no X chromosome data is "
                "available. Skipping...")
            skip_sex = True

    disable_filter_controls = False
    if hwe_requested and filter_controls and (1 not in pheno_counts):
        logger.warning(
            "You tried calling hwe prune with controls filtered but no controls "
            "are available. Skipping...")
        disable_filter_controls = True

    skip_case_control = False
    if case_control_requested and ((1 not in pheno_counts) or (2 not in pheno_counts)):
        logger.warning(
            "You tried calling case-control prune but only cases or controls are "
            "available, not both. Skipping...")
        skip_case_control = True

    skip_het = False
    if het_requested and (var.shape[0] < 50):
        logger.warning(
            "You tried calling het prune with less than 50 samples. Skipping...")
        skip_het = True

    return ValidationDecisions(
        skip_sex=skip_sex,
        skip_case_control=skip_case_control,
        skip_het=skip_het,
        disable_filter_controls=disable_filter_controls,
    )
