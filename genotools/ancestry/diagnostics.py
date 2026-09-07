# Copyright 2023 The GenoTools Authors. All Rights Reserved.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
# ==============================================================================

"""Instruments for the ancestry *prediction* path.

Training accuracy is measured on reference samples that went through the
training branch of `get_raw_files`; a cohort goes through the inference branch,
which aligns itself to a saved model's SNP list and can fail there in ways
training never sees. The failure that motivated this module is a model with
0.98 test accuracy predicting a single label -- `CAH` -- for every sample in a
cohort, because the cohort matched almost none of the model's SNPs and the
absent ones were filled with a constant. Every sample then carried the same
large offset, projected to one point off the reference manifold, and landed
nearer the global centroid than any ancestry centroid, which is exactly the
`CAH` rule.

Each function here answers one question about that path, and the answers
compose into the pipeline's report:

- `SnpOverlapReport` -- how much of the model's SNP list the cohort actually
  carried, and how much was filled in
- `bim_compatibility` -- when a match comes back near-empty, whether the two
  files even describe the same coordinate space
- `allele_frequency_concordance` -- whether matched sites agree on allele
  frequency, or look strand/allele-swapped
- `pc_drift` -- where the cohort landed in PC space relative to the reference
  it is being compared against
- `admixture_summary` -- what the classifier said before the `CAH` override
  replaced it, and by what margin

They are deliberately pure functions over frames so they can be tested against
hand-built input, per the `select_het_outliers` pattern.
"""

from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from genotools.core.logging import get_logger

logger = get_logger(__name__)


#: Fraction of a model's SNPs that may be absent from a cohort before the run
#: warns. Small overlap losses are normal (a cohort legitimately drops variants
#: at `--geno 0.1`); this is the level at which the projection starts to be
#: driven by fill rather than by data.
WARN_FILL_FRACTION = 0.05

#: Fraction above which prediction is refused outright. At half the model's
#: SNPs filled, a predicted label describes the fill, not the sample. Override
#: with `--ancestry-max-missing-snps`.
DEFAULT_MAX_FILL_FRACTION = 0.5

#: A match this small means the two files probably do not share a coordinate
#: space (chromosome naming, genome build) rather than merely disagreeing on
#: content. Expressed against the smaller of the two variant sets.
LOW_OVERLAP_FRACTION = 0.05

#: Cohort centroid displacement, in reference PC standard deviations, that is
#: too large to be ancestry. Roughly: no real population sits 5 SD off the
#: reference panel's own spread on an early PC.
WARN_PC_SHIFT_SD = 5.0

#: Cohort-to-reference PC spread ratios outside this range. A ratio near zero
#: means the cohort collapsed to a point (the constant-fill signature); a large
#: one means the projection is dominated by something other than ancestry.
WARN_PC_SD_RATIO = (0.1, 10.0)

#: Per-SNP allele-frequency difference (on the 0-1 scale) counted as
#: discordant between panel and cohort.
DISCORDANT_FREQ_DELTA = 0.2

#: How close to `1 - ref_freq` a discordant SNP must sit to be called an
#: allele-swap signature rather than ordinary population difference.
FLIP_SIGNATURE_TOLERANCE = 0.05

#: Fraction of samples labeled CAH that is high enough to suspect the override
#: rather than the cohort. A real cohort has admixed samples; it does not
#: consist of them.
WARN_CAH_FRACTION = 0.5


def _f(value: Any) -> Optional[float]:
    """Plain float for JSON, or None for anything not finite."""
    if value is None:
        return None
    try:
        out = float(value)
    except (TypeError, ValueError):
        return None
    return out if np.isfinite(out) else None


# ---------------------------------------------------------------------------
# 1. How much of the model's SNP list did the cohort actually carry?
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class SnpOverlapReport:
    """What the cohort contributed to a model's feature matrix, and what was
    filled in for it.

    The inference branch of `get_raw_files` must hand the model a matrix with
    exactly the model's SNPs, in order. Any SNP the cohort does not carry has
    to be invented, and this records how much of the matrix was invented -- the
    number that separates "the model is wrong about this cohort" from "the
    model was never shown this cohort".

    Attributes:
        n_model_snps: SNPs the model expects (the width of its feature matrix).
        n_matched: Of those, how many came from the cohort.
        fill_strategy: How the remainder were filled. See
            `preprocessing.MISSING_FILL_STRATEGIES`.
        filled_snps: IDs of the filled SNPs, for writing alongside the run.
    """

    n_model_snps: int
    n_matched: int
    fill_strategy: str
    filled_snps: Tuple[str, ...] = ()

    @property
    def n_filled(self) -> int:
        """SNPs the cohort did not carry."""
        return self.n_model_snps - self.n_matched

    @property
    def fraction_filled(self) -> float:
        """Share of the feature matrix that is fill rather than data."""
        if self.n_model_snps <= 0:
            return 0.0
        return self.n_filled / self.n_model_snps

    def to_dict(self) -> Dict[str, Any]:
        """Counts for the JSON report. The ID list goes to its own file."""
        return {
            "n_model_snps": int(self.n_model_snps),
            "n_matched": int(self.n_matched),
            "n_filled": int(self.n_filled),
            "fraction_filled": _f(self.fraction_filled),
            "fill_strategy": self.fill_strategy,
        }

    def format_summary(self) -> str:
        """One line for the log."""
        return (
            f"cohort matched {self.n_matched}/{self.n_model_snps} model SNPs "
            f"({100 * (1 - self.fraction_filled):.2f}%); "
            f"{self.n_filled} filled with '{self.fill_strategy}' "
            f"({100 * self.fraction_filled:.2f}% of the feature matrix)"
        )


def snp_overlap_warnings(
    report: SnpOverlapReport,
    warn_fraction: float = WARN_FILL_FRACTION,
) -> List[str]:
    """Human-readable concerns about a `SnpOverlapReport`, worst first."""
    warnings: List[str] = []
    if report.n_model_snps <= 0:
        return ["the model reports no SNPs, so nothing could be matched"]
    if report.n_matched == 0:
        warnings.append(
            "the cohort matched NONE of the model's SNPs: every feature is "
            "fill, so predictions describe the fill and not the samples. Check "
            "chromosome naming and genome build against the reference panel."
        )
    elif report.fraction_filled >= warn_fraction:
        warnings.append(
            f"{100 * report.fraction_filled:.2f}% of the model's SNPs were "
            f"absent from the cohort and filled with "
            f"'{report.fill_strategy}'. Predictions are correspondingly "
            f"driven by fill rather than by genotypes."
        )
    return warnings


# ---------------------------------------------------------------------------
# 2. Do the two files even describe the same coordinate space?
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class BimCompatibility:
    """Whether two `.bim` files are talking about the same genome.

    Produced when a variant match comes back small or empty, where the useful
    answer is almost never "these cohorts differ" and almost always one of a
    handful of clerical mismatches: `chr1` against `1`, hg19 against hg38, or a
    file that is not the one intended.

    Attributes:
        n_variants: Variant counts, `(file1, file2)`.
        shared_chroms: Chromosome names present in both.
        only_in_1: Chromosome names unique to the first file (truncated).
        only_in_2: Chromosome names unique to the second file (truncated).
        position_ranges: Per shared chromosome, `(min1, max1, min2, max2)`.
        hints: Ordered, plain-language explanations of what looks wrong.
    """

    n_variants: Tuple[int, int]
    shared_chroms: Tuple[str, ...]
    only_in_1: Tuple[str, ...]
    only_in_2: Tuple[str, ...]
    position_ranges: Dict[str, Tuple[int, int, int, int]] = field(
        default_factory=dict
    )
    hints: Tuple[str, ...] = ()

    def to_dict(self) -> Dict[str, Any]:
        """Serializable form for the report."""
        return {
            "n_variants": [int(n) for n in self.n_variants],
            "shared_chroms": list(self.shared_chroms),
            "only_in_1": list(self.only_in_1),
            "only_in_2": list(self.only_in_2),
            "hints": list(self.hints),
        }

    def format_report(self) -> str:
        """Multi-line text for a log or an exception message."""
        lines = [
            f"  variants:          {self.n_variants[0]} vs {self.n_variants[1]}",
            f"  shared chroms:     {', '.join(self.shared_chroms) or '(none)'}",
        ]
        if self.only_in_1:
            lines.append(f"  only in first:     {', '.join(self.only_in_1)}")
        if self.only_in_2:
            lines.append(f"  only in second:    {', '.join(self.only_in_2)}")
        for chrom, (lo1, hi1, lo2, hi2) in sorted(self.position_ranges.items()):
            lines.append(
                f"  chr {chrom} positions: {lo1}-{hi1} vs {lo2}-{hi2}"
            )
        for hint in self.hints:
            lines.append(f"  hint:              {hint}")
        return "\n".join(lines)


def _strip_chr(values: Sequence[str]) -> set:
    return {str(v)[3:] if str(v).lower().startswith("chr") else str(v) for v in values}


def bim_compatibility(
    bim1: pd.DataFrame,
    bim2: pd.DataFrame,
    max_chroms_listed: int = 8,
    max_ranges: int = 3,
) -> BimCompatibility:
    """Compare the coordinate spaces of two `.bim`-shaped frames.

    Args:
        bim1: Frame with `chr` and `pos` columns.
        bim2: Frame with `chr` and `pos` columns.
        max_chroms_listed: Cap on chromosome names reported per side.
        max_ranges: Cap on per-chromosome position ranges reported.

    Returns:
        A `BimCompatibility` whose `hints` name the likely mismatch.
    """
    chroms1 = {str(c) for c in bim1["chr"].unique()}
    chroms2 = {str(c) for c in bim2["chr"].unique()}
    shared = sorted(chroms1 & chroms2)

    hints: List[str] = []
    if not shared:
        if _strip_chr(chroms1) & _strip_chr(chroms2):
            hints.append(
                "chromosome naming differs ('chr1' vs '1'): the two files "
                "share no chromosome name, but do once a 'chr' prefix is "
                "stripped. Rename one side's contigs before matching."
            )
        else:
            hints.append(
                "no chromosome name is present in both files -- check that "
                "both paths point at the intended data."
            )

    ranges: Dict[str, Tuple[int, int, int, int]] = {}
    for chrom in shared[:max_ranges]:
        p1 = bim1.loc[bim1["chr"].astype(str) == chrom, "pos"]
        p2 = bim2.loc[bim2["chr"].astype(str) == chrom, "pos"]
        if len(p1) and len(p2):
            ranges[chrom] = (int(p1.min()), int(p1.max()), int(p2.min()), int(p2.max()))

    if shared:
        shared_positions = 0
        for chrom in shared:
            pos1 = set(bim1.loc[bim1["chr"].astype(str) == chrom, "pos"])
            pos2 = set(bim2.loc[bim2["chr"].astype(str) == chrom, "pos"])
            shared_positions += len(pos1 & pos2)
        if shared_positions == 0:
            hints.append(
                "the files share chromosomes but not a single position, which "
                "is what a genome-build mismatch looks like (hg19 vs hg38). "
                "Lift one side over before matching."
            )

    return BimCompatibility(
        n_variants=(len(bim1), len(bim2)),
        shared_chroms=tuple(shared[:max_chroms_listed]),
        only_in_1=tuple(sorted(chroms1 - chroms2)[:max_chroms_listed]),
        only_in_2=tuple(sorted(chroms2 - chroms1)[:max_chroms_listed]),
        position_ranges=ranges,
        hints=tuple(hints),
    )


def is_low_overlap(
    n_common: int,
    n_variants: Tuple[int, int],
    floor_fraction: float = LOW_OVERLAP_FRACTION,
) -> bool:
    """Whether a match is small enough to be a coordinate-space problem.

    Measured against the *smaller* input: matching a 43k panel into a 1.9M
    cohort legitimately keeps 43k, so a fraction of the larger side would call
    every healthy run suspicious.
    """
    smaller = min(n for n in n_variants if n > 0) if any(n_variants) else 0
    if smaller <= 0:
        return True
    return (n_common / smaller) < floor_fraction


# ---------------------------------------------------------------------------
# 3. Do matched sites agree on allele frequency?
# ---------------------------------------------------------------------------


def allele_frequency_concordance(
    ref_dosages: pd.DataFrame,
    geno_dosages: pd.DataFrame,
    discordant_delta: float = DISCORDANT_FREQ_DELTA,
    flip_tolerance: float = FLIP_SIGNATURE_TOLERANCE,
) -> Dict[str, Any]:
    """Compare per-SNP allele frequency between reference panel and cohort.

    Both matrices count the same allele at the same sites by the time they
    reach the model, so their frequencies should agree to within population
    difference. They will not if the allele-switch step mis-oriented a site, and
    a mis-oriented site has a signature worth naming: its cohort frequency sits
    near `1 - ref_freq` rather than near `ref_freq`. That distinguishes a
    strand/allele bug from a genuinely different population, which the
    correlation alone does not.

    Args:
        ref_dosages: Reference matrix, samples x SNPs, values 0/1/2 (NaN ok).
        geno_dosages: Cohort matrix over the same SNP columns.
        discordant_delta: Frequency difference counted as discordant.
        flip_tolerance: Closeness to `1 - ref_freq` counted as a swap.

    Returns:
        Dict with `n_snps`, `correlation`, `mean_abs_delta`, `n_discordant`,
        `n_flip_signature` and their fractions. Empty overlap yields zeros and
        a None correlation rather than an exception -- this is a diagnostic and
        must never be the thing that fails a run.
    """
    # A cohort can carry the same site under several probe IDs; if two
    # survived into the matrix, the second is not a second opinion worth
    # averaging, and a duplicated label would break the alignment below.
    geno_dosages = geno_dosages.loc[:, ~geno_dosages.columns.duplicated()]
    ref_dosages = ref_dosages.loc[:, ~ref_dosages.columns.duplicated()]

    shared = [c for c in ref_dosages.columns if c in set(geno_dosages.columns)]
    if not shared:
        return {
            "n_snps": 0,
            "correlation": None,
            "mean_abs_delta": None,
            "n_discordant": 0,
            "fraction_discordant": None,
            "n_flip_signature": 0,
            "fraction_flip_signature": None,
        }

    ref_freq = ref_dosages[shared].mean(axis=0, skipna=True) / 2.0
    geno_freq = geno_dosages[shared].mean(axis=0, skipna=True) / 2.0
    usable = ref_freq.notna() & geno_freq.notna()
    ref_freq, geno_freq = ref_freq[usable], geno_freq[usable]

    n_snps = int(len(ref_freq))
    if n_snps == 0:
        correlation = None
    elif ref_freq.std() == 0 or geno_freq.std() == 0:
        # A constant column has no correlation to report; saying so beats
        # emitting the NaN that pandas would.
        correlation = None
    else:
        correlation = _f(ref_freq.corr(geno_freq))

    delta = (geno_freq - ref_freq).abs()
    discordant = delta > discordant_delta
    flipped = discordant & ((geno_freq - (1.0 - ref_freq)).abs() <= flip_tolerance)

    return {
        "n_snps": n_snps,
        "correlation": correlation,
        "mean_abs_delta": _f(delta.mean()) if n_snps else None,
        "n_discordant": int(discordant.sum()),
        "fraction_discordant": _f(discordant.mean()) if n_snps else None,
        "n_flip_signature": int(flipped.sum()),
        "fraction_flip_signature": _f(flipped.mean()) if n_snps else None,
    }


def allele_frequency_warnings(
    concordance: Dict[str, Any],
    min_correlation: float = 0.7,
    max_flip_fraction: float = 0.01,
) -> List[str]:
    """Concerns about an `allele_frequency_concordance` result."""
    warnings: List[str] = []
    if not concordance.get("n_snps"):
        return warnings
    corr = concordance.get("correlation")
    if corr is not None and corr < min_correlation:
        warnings.append(
            f"panel/cohort allele frequencies correlate at only r={corr:.3f} "
            f"across {concordance['n_snps']} matched SNPs; the match itself is "
            f"suspect (wrong variant identity, or wrong file)."
        )
    flip_fraction = concordance.get("fraction_flip_signature")
    if flip_fraction is not None and flip_fraction > max_flip_fraction:
        warnings.append(
            f"{100 * flip_fraction:.2f}% of matched SNPs have a cohort "
            f"frequency near 1 - panel frequency, the signature of a "
            f"reference-allele or strand swap rather than of population "
            f"difference."
        )
    return warnings


# ---------------------------------------------------------------------------
# 4. Where did the cohort land in PC space?
# ---------------------------------------------------------------------------


def pc_columns(frame: pd.DataFrame, n_pcs: Optional[int] = None) -> List[str]:
    """PC column names in numeric order, optionally capped at `n_pcs`.

    Selected by name rather than by dropping known ID columns: a frame carrying
    an extra diagnostic column must not have it silently treated as a PC.
    """
    cols = [c for c in frame.columns if str(c).startswith("PC")]

    def _index(name: str) -> int:
        try:
            return int(str(name)[2:])
        except ValueError:
            return 1 << 30

    cols.sort(key=_index)
    return cols[:n_pcs] if n_pcs else cols


def pc_drift(
    train_pca: pd.DataFrame,
    projected: pd.DataFrame,
    n_pcs: int = 10,
) -> pd.DataFrame:
    """Per-PC comparison of where the reference sits and where the cohort landed.

    The reference PCA defines the space; the cohort is projected into it. Two
    numbers say whether that projection is believable. `shift_sd` is how far the
    cohort's centre sits from the reference's, in reference standard deviations
    -- ancestry moves a cohort a fraction of an SD to a few SD, not tens.
    `sd_ratio` is how much spread the cohort has relative to the reference; near
    zero means the cohort collapsed to a point, which is what a constant fill
    across most of the feature matrix produces.

    Args:
        train_pca: Reference PCA frame with PC columns.
        projected: Cohort PCA frame with the same PC columns.
        n_pcs: How many leading PCs to report.

    Returns:
        Frame with one row per PC: `pc`, `train_mean`, `train_sd`,
        `cohort_mean`, `cohort_sd`, `shift_sd`, `sd_ratio`.
    """
    cols = [c for c in pc_columns(train_pca, n_pcs) if c in projected.columns]
    rows = []
    for col in cols:
        t = pd.to_numeric(train_pca[col], errors="coerce")
        c = pd.to_numeric(projected[col], errors="coerce")
        t_sd = float(t.std())
        rows.append(
            {
                "pc": col,
                "train_mean": _f(t.mean()),
                "train_sd": _f(t_sd),
                "cohort_mean": _f(c.mean()),
                "cohort_sd": _f(c.std()),
                "shift_sd": _f((c.mean() - t.mean()) / t_sd) if t_sd else None,
                "sd_ratio": _f(c.std() / t_sd) if t_sd else None,
            }
        )
    return pd.DataFrame(rows)


def format_pc_drift(drift: pd.DataFrame) -> str:
    """Fixed-width text table of a `pc_drift` frame, for the log."""
    if drift.empty:
        return "  (no PCs in common between reference and cohort)"
    header = (
        f"  {'PC':<5} {'ref mean':>12} {'ref sd':>10} "
        f"{'cohort mean':>12} {'cohort sd':>10} {'shift (sd)':>11} {'sd ratio':>9}"
    )
    lines = [header, "  " + "-" * (len(header) - 2)]
    for row in drift.itertuples(index=False):
        def _fmt(value: Any, width: int, places: int = 3) -> str:
            return f"{'n/a':>{width}}" if value is None else f"{value:>{width}.{places}f}"

        lines.append(
            f"  {row.pc:<5} {_fmt(row.train_mean, 12)} {_fmt(row.train_sd, 10)} "
            f"{_fmt(row.cohort_mean, 12)} {_fmt(row.cohort_sd, 10)} "
            f"{_fmt(row.shift_sd, 11, 2)} {_fmt(row.sd_ratio, 9, 2)}"
        )
    return "\n".join(lines)


def pc_drift_warnings(
    drift: pd.DataFrame,
    max_shift_sd: float = WARN_PC_SHIFT_SD,
    sd_ratio_range: Tuple[float, float] = WARN_PC_SD_RATIO,
) -> List[str]:
    """Concerns about a `pc_drift` frame, one line per offending PC."""
    warnings: List[str] = []
    low, high = sd_ratio_range
    for row in drift.itertuples(index=False):
        if row.shift_sd is not None and abs(row.shift_sd) > max_shift_sd:
            warnings.append(
                f"{row.pc}: the cohort centre sits {row.shift_sd:.1f} reference "
                f"SD from the reference centre -- too far to be ancestry. The "
                f"feature matrix handed to the projection is the thing to check."
            )
        if row.sd_ratio is None:
            continue
        if row.sd_ratio < low:
            warnings.append(
                f"{row.pc}: the cohort has {row.sd_ratio:.3f}x the reference's "
                f"spread -- it has collapsed to a point, which is what a mostly "
                f"filled feature matrix produces."
            )
        elif row.sd_ratio > high:
            warnings.append(
                f"{row.pc}: the cohort has {row.sd_ratio:.1f}x the reference's "
                f"spread, so something other than ancestry dominates the "
                f"projection."
            )
    return warnings


# ---------------------------------------------------------------------------
# 5. What did the classifier say before the CAH override?
# ---------------------------------------------------------------------------


def admixture_summary(decisions: pd.DataFrame) -> Dict[str, Any]:
    """Summarize the CAH override from a per-sample decision frame.

    The override is the last thing to touch a label and it discards what the
    classifier said, so an all-CAH run looks identical to a broken classifier
    from the outside. Keeping both labels separates them: a sensible spread of
    `label_pre_admixture` under an all-CAH `label` localizes the fault to the
    override (or to the projection feeding it), not to the model.

    Args:
        decisions: Frame with `label_pre_admixture`, `label`, and the
            `margin_to_all` written by `AncestryModel._predict_admixed`.

    Returns:
        Dict with the CAH count and fraction, the pre-override label counts,
        the pre-to-post crosstab, and margin quantiles.
    """
    n_total = int(len(decisions))
    if n_total == 0:
        return {"n_samples": 0, "n_cah": 0, "fraction_cah": None}

    final = decisions["label"].astype(str)
    n_cah = int((final == "CAH").sum())

    summary: Dict[str, Any] = {
        "n_samples": n_total,
        "n_cah": n_cah,
        "fraction_cah": _f(n_cah / n_total),
        "final_counts": {str(k): int(v) for k, v in final.value_counts().items()},
    }

    if "label_pre_admixture" in decisions.columns:
        pre = decisions["label_pre_admixture"].astype(str)
        summary["pre_admixture_counts"] = {
            str(k): int(v) for k, v in pre.value_counts().items()
        }
        summary["overridden_from"] = {
            str(k): int(v) for k, v in pre[final == "CAH"].value_counts().items()
        }
        summary["n_distinct_pre_admixture"] = int(pre.nunique())

    if "margin_to_all" in decisions.columns:
        margin = pd.to_numeric(decisions["margin_to_all"], errors="coerce").dropna()
        if len(margin):
            summary["margin_to_all_quantiles"] = {
                q: _f(margin.quantile(v))
                for q, v in (("p05", 0.05), ("p50", 0.5), ("p95", 0.95))
            }
    return summary


def admixture_warnings(
    summary: Dict[str, Any],
    max_cah_fraction: float = WARN_CAH_FRACTION,
) -> List[str]:
    """Concerns about an `admixture_summary` result."""
    warnings: List[str] = []
    fraction = summary.get("fraction_cah")
    if fraction is None or fraction <= max_cah_fraction:
        return warnings

    detail = ""
    n_pre = summary.get("n_distinct_pre_admixture")
    if n_pre and n_pre > 1:
        detail = (
            f" The classifier itself assigned {n_pre} distinct labels before "
            f"the override, so the classifier is not the thing that failed -- "
            f"check the projection and the SNP overlap above."
        )
    warnings.append(
        f"{100 * fraction:.1f}% of samples were labeled CAH by admixture "
        f"detection.{detail} Re-run with --no-admixture-detection to see the "
        f"classifier's own labels."
    )
    return warnings


# ---------------------------------------------------------------------------
# Aggregate
# ---------------------------------------------------------------------------


@dataclass
class AncestryDiagnostics:
    """Everything the instruments found, as one thing to hand to the report.

    Assembled across two layers -- preprocessing owns the SNP overlap and the
    allele frequencies, the model owns the PC drift and the override -- so it is
    mutable and filled in as a run proceeds.
    """

    snp_overlap: Optional[SnpOverlapReport] = None
    filled_snps_path: Optional[str] = None
    allele_frequency: Optional[Dict[str, Any]] = None
    pc_drift: Optional[pd.DataFrame] = None
    admixture: Optional[Dict[str, Any]] = None
    self_test: Optional[Dict[str, Any]] = None
    plots: List[str] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)

    def add_warnings(self, warnings: Sequence[str]) -> None:
        """Record concerns, and log each one as a warning."""
        for warning in warnings:
            if warning not in self.warnings:
                self.warnings.append(warning)
                logger.warning(warning)

    def to_dict(self) -> Dict[str, Any]:
        """Serializable form for the JSON report's `ancestry_diagnostics`."""
        out: Dict[str, Any] = {}
        if self.snp_overlap is not None:
            out["snp_overlap"] = self.snp_overlap.to_dict()
        if self.filled_snps_path:
            out["filled_snps_path"] = self.filled_snps_path
        if self.allele_frequency is not None:
            out["allele_frequency"] = self.allele_frequency
        if self.pc_drift is not None and not self.pc_drift.empty:
            out["pc_drift"] = self.pc_drift.to_dict(orient="list")
        if self.admixture is not None:
            out["admixture"] = self.admixture
        if self.self_test is not None:
            out["self_test"] = self.self_test
        if self.plots:
            out["plots"] = list(self.plots)
        out["warnings"] = list(self.warnings)
        return out
