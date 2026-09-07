"""Ancestry PLINK preprocessing (legacy-free port of the Ancestry helpers).

Ports of legacy Ancestry.get_raw_files / clean_up and utils.get_common_snps,
using core.executors. `get_raw_files` deliberately diverges from the legacy
version in one place -- how it fills SNPs the cohort does not carry, see
`MISSING_FILL_STRATEGIES` -- and is otherwise identical to it (proven by
differential tests). `get_common_snps` diverges too: it excludes palindromic
sites and breaks position ties by missingness, both of which the legacy version
got wrong. See the helpers below.

This is also where the ancestry run's diagnostics start: the two questions this
module can answer that no later stage can -- how much of the model's SNP list
the cohort actually carried, and whether the matched sites agree on allele
frequency -- are recorded into an `AncestryDiagnostics` here and carried
through to the report.
"""

import os
from typing import Optional

import numpy as np
import pandas as pd

from genotools.ancestry.config import LEGACY_FILL_VALUE, MISSING_FILL_STRATEGIES
from genotools.ancestry.diagnostics import (
    AncestryDiagnostics,
    DEFAULT_MAX_FILL_FRACTION,
    SnpOverlapReport,
    allele_frequency_concordance,
    allele_frequency_warnings,
    bim_compatibility,
    is_low_overlap,
    snp_overlap_warnings,
)
from genotools.core.exceptions import AncestryError
from genotools.core.executors import run_command, get_plink, get_plink2
from genotools.core.logging import get_logger

logger = get_logger(__name__)


_PALINDROMIC = frozenset({"AT", "TA", "CG", "GC"})


def _drop_palindromic(bim: pd.DataFrame) -> pd.DataFrame:
    """Drop strand-ambiguous (palindromic) A/T and C/G sites.

    Complementing such a pair returns the same pair, so the four-route match in
    `get_common_snps` cannot tell which strand the site came from, and the
    allele-switch test in `get_raw_files` ("differs from the reference allele
    and from its complement") can never fire for one either. A palindromic site
    is therefore accepted at whatever orientation it happened to arrive in,
    which silently inverts its dosages when the strands disagree. Reference
    panels built per `docs/prep_reference_panel.md` already exclude these, so
    this is a guard for panels that did not.
    """
    return bim.loc[~(bim["a1"] + bim["a2"]).isin(_PALINDROMIC)]


def _variant_missingness(geno_path: str) -> pd.Series:
    """Per-variant missing-call rate for `geno_path`, keyed by variant ID.

    Reports to a `_missing` prefix rather than `geno_path` itself: plink2 writes
    a `.log` beside whatever `--out` names, and `geno_path` already has one
    belonging to the step that produced it (the variant prune, or the user's own
    reference panel directory). Overwriting that would destroy the audit trail
    for a different command.
    """
    out_prefix = f"{geno_path}_missing"
    run_command(
        f"{get_plink2()} --bfile {geno_path} --missing variant-only --out {out_prefix}",
        tool_name="plink2",
    )
    vmiss = pd.read_csv(f"{out_prefix}.vmiss", sep="\t", usecols=["ID", "F_MISS"])
    return vmiss.set_index("ID")["F_MISS"]


def _select_one_per_position(common_snps: pd.DataFrame, geno_path1: str) -> pd.DataFrame:
    """Keep one candidate row per chr:pos, preferring the lowest missingness.

    A cohort routinely carries the same site under several probe IDs
    (`rs301801`, `IlmnSeq_rs301801`, `seq_rs301801`), and each match route can
    contribute its own row, so a position can arrive here several times. Those
    candidates are *not* redundant copies: on a 10k GP2 cohort only 4 of 1,078
    contested positions had identical calls, missingness differed by up to 6.9%
    and genotypes by up to 49.9%. Keeping whichever row the merge emitted first
    made that choice arbitrary and sensitive to input variant order, so pick the
    best-called probe instead. Exact ties keep the earlier row, which preserves
    the previous behaviour whenever missingness cannot separate the candidates.
    """
    common_snps = common_snps.reset_index(drop=True)
    duplicated = common_snps.duplicated(subset=["chr", "pos"], keep=False)
    # Only positions whose candidates name *different* variants can be decided
    # by row order; where they all name the same one, the choice cannot matter
    # and the missingness pass would be wasted work.
    contested = (
        common_snps.loc[duplicated].groupby(["chr", "pos"])["rsid_y"].nunique().gt(1)
        if duplicated.any()
        else pd.Series(dtype=bool)
    )
    if not contested.any():
        return common_snps.drop_duplicates(subset=["chr", "pos"], ignore_index=True)

    logger.info(
        f"Resolving {int(contested.sum())} position(s) with competing variants by lowest missingness"
    )
    f_miss = _variant_missingness(geno_path1)
    common_snps["_row"] = range(len(common_snps))
    # Unscored variants sort last rather than winning by default.
    common_snps["_f_miss"] = common_snps["rsid_y"].map(f_miss).fillna(1.0)
    return (
        common_snps.sort_values(["_f_miss", "_row"], kind="stable")
        .drop_duplicates(subset=["chr", "pos"])
        .sort_values("_row")
        .drop(columns=["_row", "_f_miss"])
        .reset_index(drop=True)
    )


def _check_overlap(
    common_snps: pd.DataFrame,
    bim1: pd.DataFrame,
    bim2: pd.DataFrame,
    geno_path1: str,
    geno_path2: str,
) -> None:
    """Fail loudly when two files barely match, and say why they might not.

    A near-empty match is almost never two cohorts genuinely differing. It is
    `chr1` against `1`, hg19 against hg38, or a path pointing at the wrong
    file -- and none of those are visible downstream, where the symptom is
    instead a feature matrix made mostly of fill and a cohort that predicts one
    label. `bim_compatibility` names the mismatch; this decides whether to warn
    or to stop.

    Zero overlap raises: `plink2 --extract` on an empty list fails anyway, with
    a message about a file rather than about the genome.
    """
    n_common = len(common_snps)
    n_variants = (len(bim1), len(bim2))
    if not is_low_overlap(n_common, n_variants):
        return

    compat = bim_compatibility(bim1, bim2)
    if n_common == 0:
        raise AncestryError(
            f"No common variants between {geno_path1} and {geno_path2}.\n"
            f"{compat.format_report()}"
        )
    logger.warning(
        f"Only {n_common} variant(s) are common to {geno_path1} and "
        f"{geno_path2}:\n{compat.format_report()}"
    )


def get_common_snps(geno_path1: str, geno_path2: str, out_name: str) -> dict:
    """Extract SNPs common to two bfiles from geno_path1. Returns output paths."""
    logger.info("Getting Common SNPs")
    plink = get_plink()
    plink2 = get_plink2()

    bim1 = pd.read_csv(f"{geno_path1}.bim", sep="\t", header=None, dtype={0: str}, low_memory=False)
    bim1.columns = ["chr", "rsid", "kb", "pos", "a1", "a2"]
    bim2 = pd.read_csv(f"{geno_path2}.bim", sep="\t", header=None, dtype={0: str}, low_memory=False)
    bim2.columns = ["chr", "rsid", "kb", "pos", "a1", "a2"]

    n_before = len(bim1) + len(bim2)
    bim1 = _drop_palindromic(bim1)
    bim2 = _drop_palindromic(bim2)
    n_palindromic = n_before - len(bim1) - len(bim2)
    if n_palindromic:
        logger.info(f"Excluded {n_palindromic} palindromic (strand-ambiguous) variant(s)")

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
    bim1_flip = _drop_palindromic(bim1_flip)

    bim1_flip["merge_id"] = bim1_flip["chr"].astype(str) + ":" + bim1_flip["pos"].astype(str) + ":" + bim1_flip["a2"] + ":" + bim1_flip["a1"]
    common_snps1 = bim2[["rsid", "merge_id1", "a1", "a2"]].merge(bim1_flip, how="inner", left_on=["merge_id1"], right_on=["merge_id"])
    common_snps2 = bim2[["rsid", "merge_id2", "a1", "a2"]].merge(bim1_flip, how="inner", left_on=["merge_id2"], right_on=["merge_id"])

    common_snps = pd.concat([common_snps, common_snps1, common_snps2], axis=0)
    common_snps = _select_one_per_position(common_snps, geno_path1)
    _check_overlap(common_snps, bim1, bim2, geno_path1, geno_path2)

    common_snps_file = f"{out_name}.common_snps"
    common_snps["rsid_y"].to_csv(f"{common_snps_file}", sep="\t", header=False, index=False)

    run_command(
        f"{plink2} --bfile {geno_path1} --extract {common_snps_file} --make-bed --out {out_name}",
        tool_name="plink2",
    )

    return {"common_snps": common_snps_file, "bed": out_name}


def clean_up_files(files: list) -> None:
    """Remove intermediate PLINK files by known extensions (port of Ancestry.clean_up)."""
    extensions = ["bim", "bed", "fam", "hh", "snplist", "ref_allele", "alleles", "raw", "vmiss"]
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
    fill_strategy: str = "ref-mean",
    max_missing_fraction: float = DEFAULT_MAX_FILL_FRACTION,
    diagnostics: Optional[AncestryDiagnostics] = None,
    diagnostics_prefix: Optional[str] = None,
) -> dict:
    """Process reference + genotype data into labeled raw feature matrices.

    Port of legacy Ancestry.get_raw_files (train + model-inference branches;
    containerized branch intentionally omitted). In inference mode, pass the
    model's common-SNP list path via `common_snps_file`.

    Args:
        geno_path: Cohort pfile prefix.
        ref_panel: Reference panel bfile prefix.
        ref_labels: TSV of `FID IID label` for the reference panel.
        out_path: Output prefix for intermediate files.
        train: Whether to derive the common-SNP list (True) or take it from a
            saved model (False).
        common_snps_file: The model's common-SNP list. Required when
            `train` is False.
        fill_strategy: How to fill a model SNP the cohort does not carry. One
            of `MISSING_FILL_STRATEGIES`.
        max_missing_fraction: Refuse to build the matrix when more than this
            share of the model's SNPs had to be filled. A label predicted from
            mostly-filled features describes the fill.
        diagnostics: Collector to record the SNP overlap and allele-frequency
            findings into. One is created if not supplied, and returned either
            way.
        diagnostics_prefix: Prefix for diagnostic artifacts the user keeps.
            Defaults to `out_path`, which for a partial-output run is a
            temporary directory -- the runner passes the final output prefix so
            the evidence outlives the run that produced it.

    Returns:
        Dict with `raw_ref`, `raw_geno`, `out_paths` and `diagnostics`.

    Raises:
        AncestryError: If `fill_strategy` is unknown, or if the cohort matched
            too little of the model's SNP list to predict from.
        FileNotFoundError: If an expected PLINK output is missing.
    """
    if fill_strategy not in MISSING_FILL_STRATEGIES:
        raise AncestryError(
            f"Unknown fill strategy {fill_strategy!r}; expected one of "
            f"{', '.join(MISSING_FILL_STRATEGIES)}"
        )

    plink2 = get_plink2()
    out_paths = {}
    diagnostics = diagnostics if diagnostics is not None else AncestryDiagnostics()
    diagnostics_prefix = diagnostics_prefix or out_path

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
    ref_snps = ref_raw.drop(columns=["FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE"])
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
    geno_snps = raw_geno.drop(columns=["FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE"])
    geno_snps.columns = geno_snps.columns.str.extract("(.*)_")[0]

    # Do the panel and the cohort agree about the sites they matched? Asked
    # before any fill, so the answer is about real genotypes only. Both
    # matrices count the same allele by this point, so disagreement here is a
    # matching or orientation fault, not population difference.
    concordance = allele_frequency_concordance(ref_snps, geno_snps)
    diagnostics.allele_frequency = concordance
    correlation = concordance["correlation"]
    logger.info(
        f"Panel/cohort allele-frequency concordance over "
        f"{concordance['n_snps']} matched SNP(s): "
        f"r={'n/a' if correlation is None else f'{correlation:.4f}'}, "
        f"{concordance['n_discordant']} discordant, "
        f"{concordance['n_flip_signature']} with a swap signature"
    )
    diagnostics.add_warnings(allele_frequency_warnings(concordance))

    # Adding missing snps when not training. The model's feature matrix has a
    # fixed width and order, so a SNP the cohort does not carry has to be
    # invented -- and how much of the matrix was invented is the number that
    # separates "the model is wrong about this cohort" from "the model was
    # never shown this cohort". Recorded, warned about, and above a threshold
    # refused, because none of it is visible in a predicted label.
    if not train:
        absent = [col for col in ref_snps.columns if col not in set(geno_snps.columns)]
        overlap = SnpOverlapReport(
            n_model_snps=len(ref_snps.columns),
            n_matched=len(ref_snps.columns) - len(absent),
            fill_strategy=fill_strategy,
            filled_snps=tuple(absent),
        )
        diagnostics.snp_overlap = overlap
        logger.info(f"Model SNP overlap: {overlap.format_summary()}")

        if absent:
            filled_path = f"{diagnostics_prefix}_filled_snps.txt"
            with open(filled_path, "w") as handle:
                for snp in absent:
                    handle.write(f"{snp}\n")
            diagnostics.filled_snps_path = filled_path
            logger.info(f"Filled SNP IDs written to: {filled_path}")

        diagnostics.add_warnings(snp_overlap_warnings(overlap))
        if overlap.fraction_filled > max_missing_fraction:
            raise AncestryError(
                f"The cohort carries only {overlap.n_matched} of the model's "
                f"{overlap.n_model_snps} SNPs, so "
                f"{100 * overlap.fraction_filled:.2f}% of the feature matrix "
                f"would be fill -- above the "
                f"{100 * max_missing_fraction:.0f}% limit. A label predicted "
                f"from that describes the fill, not the sample. Check that the "
                f"cohort and the reference panel share a genome build and "
                f"chromosome naming; the matched and filled SNPs are listed "
                f"beside this run. Raise --ancestry-max-missing-snps to "
                f"predict anyway."
            )

        if absent:
            # ``constant`` reproduces 1.x exactly, int dosages included;
            # ``ref-mean`` writes NaN for PCAReducer's reference-fitted imputer
            # to fill. See MISSING_FILL_STRATEGIES.
            if fill_strategy == "constant":
                fill = np.repeat(LEGACY_FILL_VALUE, geno_snps.shape[0])
            else:
                fill = np.full(geno_snps.shape[0], np.nan)
            missing_cols = pd.DataFrame(
                {col: fill for col in absent}, index=geno_snps.index
            )
            geno_snps = pd.concat([geno_snps, missing_cols], axis=1)
        geno_snps = geno_snps[ref_snps.columns]

    raw_geno = pd.concat([geno_ids, geno_snps], axis=1)
    raw_geno.columns = col_names
    raw_geno["label"] = "new"

    # remove intermediate files (concat_logs dropped: structured logging replaces raw .log aggregation)
    files = [geno_prune_path, f"{geno_prune_path}_flip", f"{geno_prune_path}_missing",
             f"{out_path}_common_snps_switch"]
    clean_up_files(files)

    return {
        "raw_ref": labeled_ref_raw,
        "raw_geno": raw_geno,
        "out_paths": out_paths,
        "diagnostics": diagnostics,
    }
