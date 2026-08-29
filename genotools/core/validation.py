"""Pre-flight validation of pipeline input (legacy-free).

Ports the validation + data-breakdown behavior of legacy utils.upfront_check,
including the data-driven step-skip decisions.
"""

from dataclasses import dataclass, field
from pathlib import Path
from typing import Mapping, Optional, Union

import pandas as pd

from .exceptions import ValidationError
from .logging import get_logger

logger = get_logger(__name__)

# PLINK2's --indep-pairwise refuses to estimate LD from fewer than 50 samples
# (it exits 13), and the het step LD-prunes before it can call --het. Below this
# floor het cannot run at all, so it is skipped rather than left to fail.
MIN_HET_SAMPLES = 50

# One reason string per decision, shared by the cohort-level checks in
# ``validate_input`` and the per-dataset checks below so the two cannot drift
# into different explanations of the same finding.
SEX_SKIP_REASON = "no sample sex data is available"
NO_X_SKIP_REASON = "no X chromosome data is available"
CASE_CONTROL_SKIP_REASON = "only cases or controls are available, not both"


def het_floor_reason(n_samples: int) -> str:
    """Phrase the het skip for a dataset that cannot support LD estimation."""
    return (
        f"{n_samples} samples is fewer than the {MIN_HET_SAMPLES} "
        f"PLINK requires to estimate LD"
    )


def count_samples(geno_path: Union[str, Path]) -> int:
    """Count samples in a pfile's .psam, excluding the header line."""
    with open(f"{geno_path}.psam") as f:
        return sum(1 for _ in f) - 1


# --- per-dataset step decisions ----------------------------------------------
#
# ``validate_input`` decides against the whole cohort. Under ``--ancestry``
# every group is its own dataset, and the split can change any answer derived
# from samples: the cohort can clear PLINK's LD floor, hold sample sex, and
# hold both phenotypes while an individual group does none of those. The
# functions below re-decide those three against one dataset, and the runner
# applies them per group.
#
# Only sample-derived checks belong here. The X chromosome check cannot change,
# because the split keeps samples and every group inherits the cohort's pvar.
#
# All of them read the psam as it stands *before* the QC chain runs, so a
# precondition broken by an earlier prune is left to the step to raise.


def _psam_column(geno_path: Union[str, Path], column: str) -> Optional[pd.Series]:
    """Read one column of a pfile's .psam, or None if the column is absent.

    A missing *file* raises: the caller resolved this path, and deciding
    nothing would hide a wrong prefix (which is how every ancestry run once
    died). A missing *column* returns None, leaving the decision to the step,
    which raises its own specific error - ``validate_input`` already rejects a
    cohort psam lacking SEX or PHENO1.
    """
    with open(f"{geno_path}.psam") as f:
        header = f.readline().split()
    if column not in header:
        return None
    frame = pd.read_csv(f"{geno_path}.psam", sep=r"\s+", usecols=[column])
    return frame[column]


def het_skip_reason(geno_path: Union[str, Path]) -> Optional[str]:
    """Decide whether het pruning can run on a specific dataset.

    The step LD-prunes before it can call ``--het``, and ``--indep-pairwise``
    refuses fewer than ``MIN_HET_SAMPLES`` samples.
    """
    n_samples = count_samples(geno_path)
    if n_samples < MIN_HET_SAMPLES:
        return het_floor_reason(n_samples)
    return None


def sex_skip_reason(geno_path: Union[str, Path]) -> Optional[str]:
    """Decide whether sex pruning can run on a specific dataset.

    ``--check-sex`` needs recorded sample sex to compare its calls against. A
    dataset with some sex recorded still runs, matching the cohort-level rule.
    """
    sex = _psam_column(geno_path, "SEX")
    if sex is None:
        return None

    sex_counts = sex.value_counts().to_dict()
    if (1 not in sex_counts) and (2 not in sex_counts):
        return SEX_SKIP_REASON
    return None


def case_control_skip_reason(geno_path: Union[str, Path]) -> Optional[str]:
    """Decide whether case-control missingness can run on a specific dataset.

    ``--test-missing`` needs both cases (2) and controls (1); with only one it
    raises inside the step.
    """
    pheno = _psam_column(geno_path, "PHENO1")
    if pheno is None:
        return None

    pheno_counts = pheno.value_counts().to_dict()
    if (1 not in pheno_counts) or (2 not in pheno_counts):
        return CASE_CONTROL_SKIP_REASON
    return None


@dataclass(frozen=True)
class ValidationDecisions:
    """Data-driven step-skip decisions derived from the input breakdown.

    ``skip_reasons`` maps a step key to why it cannot run. Steps named here are
    reported as skipped rather than dropped silently, so a consumer can tell
    "was not requested" from "was requested but could not run".
    """

    skip_reasons: Mapping[str, str] = field(default_factory=dict)
    disable_filter_controls: bool = False

    @property
    def skip_sex(self) -> bool:
        return "sex" in self.skip_reasons

    @property
    def skip_case_control(self) -> bool:
        return "case_control" in self.skip_reasons

    @property
    def skip_het(self) -> bool:
        return "het" in self.skip_reasons


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

    skip_reasons: dict[str, str] = {}

    if sex_requested:
        if (1 not in sex_counts) and (2 not in sex_counts):
            skip_reasons["sex"] = SEX_SKIP_REASON
        elif ("23" not in chr_counts) and ("X" not in chr_counts):
            skip_reasons["sex"] = NO_X_SKIP_REASON

    disable_filter_controls = False
    if hwe_requested and filter_controls and (1 not in pheno_counts):
        logger.warning(
            "You tried calling hwe prune with controls filtered but no controls "
            "are available. Running hwe without filtering controls...")
        disable_filter_controls = True

    if case_control_requested and ((1 not in pheno_counts) or (2 not in pheno_counts)):
        skip_reasons["case_control"] = CASE_CONTROL_SKIP_REASON

    if het_requested and (sam.shape[0] < MIN_HET_SAMPLES):
        # Pre-2.0 this read var.shape[0] (the variant count, never < 50 on real
        # data), so the guard never fired and het was left to fail inside PLINK.
        skip_reasons["het"] = het_floor_reason(sam.shape[0])

    # Not logged here: the runner announces each skip at the point it applies,
    # which under --ancestry is once per group with that group's own numbers.

    return ValidationDecisions(
        skip_reasons=skip_reasons,
        disable_filter_controls=disable_filter_controls,
    )
