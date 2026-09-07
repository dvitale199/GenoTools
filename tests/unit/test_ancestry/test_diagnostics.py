"""Unit tests for the ancestry prediction diagnostics.

Every function here is pure over frames, so these run against hand-built input
rather than a PLINK pipeline -- the `select_het_outliers` pattern. The cases are
the real failure shapes: a cohort that matched nothing, `chr1` against `1`,
hg19 against hg38, an allele-swapped match, a cohort collapsed to a point, and
an all-CAH override sitting on top of a classifier that worked.
"""

import numpy as np
import pandas as pd
import pytest

from genotools.ancestry.diagnostics import (
    AncestryDiagnostics,
    SnpOverlapReport,
    admixture_summary,
    admixture_warnings,
    allele_frequency_concordance,
    allele_frequency_warnings,
    bim_compatibility,
    format_pc_drift,
    is_low_overlap,
    pc_columns,
    pc_drift,
    pc_drift_warnings,
    snp_overlap_warnings,
)


# --- SNP overlap -----------------------------------------------------------


def test_overlap_report_counts_and_fraction():
    report = SnpOverlapReport(
        n_model_snps=1000, n_matched=750, fill_strategy="ref-mean",
        filled_snps=tuple(f"rs{i}" for i in range(250)),
    )
    assert report.n_filled == 250
    assert report.fraction_filled == pytest.approx(0.25)
    assert report.to_dict()["n_filled"] == 250
    # The ID list is a file, not a report field.
    assert "filled_snps" not in report.to_dict()
    assert "750/1000" in report.format_summary()


def test_overlap_report_survives_an_empty_model():
    report = SnpOverlapReport(n_model_snps=0, n_matched=0, fill_strategy="ref-mean")
    assert report.fraction_filled == 0.0
    assert "no SNPs" in snp_overlap_warnings(report)[0]


def test_small_overlap_loss_is_not_worth_a_warning():
    report = SnpOverlapReport(n_model_snps=1000, n_matched=990, fill_strategy="ref-mean")
    assert snp_overlap_warnings(report) == []


def test_overlap_loss_above_the_threshold_warns():
    report = SnpOverlapReport(n_model_snps=1000, n_matched=800, fill_strategy="ref-mean")
    (warning,) = snp_overlap_warnings(report)
    assert "20.00%" in warning


def test_matching_nothing_says_so_and_names_the_usual_causes():
    """The long-read failure: a cohort sharing no coordinate space at all."""
    report = SnpOverlapReport(n_model_snps=1000, n_matched=0, fill_strategy="ref-mean")
    (warning,) = snp_overlap_warnings(report)
    assert "NONE" in warning
    assert "genome build" in warning


# --- bim compatibility ------------------------------------------------------


def _bim(chroms, positions):
    return pd.DataFrame({"chr": [str(c) for c in chroms], "pos": list(positions)})


def test_matched_files_get_no_hints():
    bim = _bim(["1", "1", "2"], [100, 200, 300])
    compat = bim_compatibility(bim, bim)
    assert compat.hints == ()
    assert compat.shared_chroms == ("1", "2")


def test_chr_prefix_mismatch_is_named():
    compat = bim_compatibility(
        _bim(["1", "2"], [100, 200]), _bim(["chr1", "chr2"], [100, 200])
    )
    assert compat.shared_chroms == ()
    assert "chr1' vs '1" in compat.hints[0]
    assert "1" in compat.only_in_1 and "chr1" in compat.only_in_2


def test_disjoint_positions_on_shared_chromosomes_suggests_a_build_mismatch():
    compat = bim_compatibility(
        _bim(["1", "1"], [10_000, 20_000]),
        _bim(["1", "1"], [910_000, 920_000]),
    )
    assert compat.shared_chroms == ("1",)
    assert "hg19 vs hg38" in compat.hints[0]
    assert compat.position_ranges["1"] == (10_000, 20_000, 910_000, 920_000)


def test_unrelated_files_say_nothing_is_shared():
    compat = bim_compatibility(_bim(["1"], [100]), _bim(["X"], [100]))
    assert "no chromosome name is present in both" in compat.hints[0]


def test_report_text_lists_both_sides():
    text = bim_compatibility(
        _bim(["1"], [100]), _bim(["chr1"], [100])
    ).format_report()
    assert "variants:" in text and "hint:" in text


# --- overlap floor ----------------------------------------------------------


def test_a_panel_matched_into_a_much_larger_cohort_is_not_low_overlap():
    """The real shape: 43k of a 209k panel found in a 1.9M cohort."""
    assert not is_low_overlap(43_173, (209_517, 1_900_000))


def test_a_near_empty_match_is_low_overlap():
    assert is_low_overlap(12, (209_517, 1_900_000))


def test_zero_variants_counts_as_low_overlap():
    assert is_low_overlap(0, (0, 0))


# --- allele frequency -------------------------------------------------------


def _dosages(**columns):
    return pd.DataFrame(columns)


def test_identical_matrices_agree_perfectly():
    ref = _dosages(rs1=[0, 1, 2, 1], rs2=[2, 2, 1, 0], rs3=[0, 0, 1, 1])
    result = allele_frequency_concordance(ref, ref.copy())
    assert result["n_snps"] == 3
    assert result["correlation"] == pytest.approx(1.0)
    assert result["n_discordant"] == 0
    assert result["n_flip_signature"] == 0


def test_a_swapped_allele_is_recognised_as_a_swap_not_a_difference():
    """`2 - dosage` is the signature: frequency near `1 - ref_freq`."""
    ref = _dosages(rs1=[0, 0, 0, 1], rs2=[2, 2, 2, 1], rs3=[0, 1, 2, 1])
    cohort = ref.copy()
    cohort["rs1"] = 2 - cohort["rs1"]
    result = allele_frequency_concordance(ref, cohort)
    assert result["n_discordant"] == 1
    assert result["n_flip_signature"] == 1
    (warning,) = allele_frequency_warnings(result, min_correlation=0.0)
    assert "swap" in warning


def test_population_difference_is_discordant_without_a_swap_signature():
    ref = _dosages(rs1=[0, 0, 0, 0], rs2=[0, 1, 2, 1], rs3=[2, 1, 0, 1])
    cohort = _dosages(rs1=[1, 1, 1, 1], rs2=[0, 1, 2, 1], rs3=[2, 1, 0, 1])
    result = allele_frequency_concordance(ref, cohort)
    assert result["n_discordant"] == 1
    assert result["n_flip_signature"] == 0


def test_no_shared_snps_reports_zeros_rather_than_raising():
    result = allele_frequency_concordance(
        _dosages(rs1=[0, 1]), _dosages(rs9=[0, 1])
    )
    assert result == {
        "n_snps": 0,
        "correlation": None,
        "mean_abs_delta": None,
        "n_discordant": 0,
        "fraction_discordant": None,
        "n_flip_signature": 0,
        "fraction_flip_signature": None,
    }


def test_duplicate_cohort_probes_do_not_break_alignment():
    """A cohort can carry one site under several probe IDs."""
    ref = _dosages(rs1=[0, 1, 2], rs2=[2, 1, 0])
    cohort = pd.DataFrame(
        np.array([[0, 1, 2], [1, 1, 1], [2, 1, 0]]).T, columns=["rs1", "rs1", "rs2"]
    )
    result = allele_frequency_concordance(ref, cohort)
    assert result["n_snps"] == 2


def test_a_constant_column_yields_no_correlation_rather_than_nan():
    ref = _dosages(rs1=[1, 1], rs2=[1, 1])
    assert allele_frequency_concordance(ref, ref.copy())["correlation"] is None


def test_a_weak_correlation_warns_about_the_match_itself():
    result = {"n_snps": 500, "correlation": 0.1, "fraction_flip_signature": 0.0}
    (warning,) = allele_frequency_warnings(result)
    assert "r=0.100" in warning


# --- PC drift ---------------------------------------------------------------


def test_pc_columns_sort_numerically_and_ignore_everything_else():
    frame = pd.DataFrame(
        columns=["FID", "IID", "PC10", "PC2", "PC1", "label", "dist_to_all"]
    )
    assert pc_columns(frame) == ["PC1", "PC2", "PC10"]
    assert pc_columns(frame, 2) == ["PC1", "PC2"]


def _train_frame(n=200, seed=0):
    rng = np.random.default_rng(seed)
    return pd.DataFrame({"PC1": rng.normal(0, 1, n), "PC2": rng.normal(0, 2, n)})


def test_a_cohort_inside_the_reference_cloud_drifts_little():
    train = _train_frame()
    drift = pc_drift(train, _train_frame(seed=1))
    assert list(drift["pc"]) == ["PC1", "PC2"]
    assert abs(drift.loc[0, "shift_sd"]) < 0.5
    assert drift.loc[0, "sd_ratio"] == pytest.approx(1.0, abs=0.2)
    assert pc_drift_warnings(drift) == []


def test_a_displaced_cohort_is_measured_in_reference_sds():
    train = _train_frame()
    cohort = _train_frame(seed=1)
    cohort["PC1"] += 40
    drift = pc_drift(train, cohort)
    assert drift.loc[0, "shift_sd"] == pytest.approx(40, rel=0.1)
    (warning,) = pc_drift_warnings(drift)
    assert "too far to be ancestry" in warning


def test_a_cohort_collapsed_to_a_point_is_the_constant_fill_signature():
    train = _train_frame()
    cohort = pd.DataFrame({"PC1": [3.0] * 50, "PC2": [1.0] * 50})
    drift = pc_drift(train, cohort)
    assert drift.loc[0, "sd_ratio"] == pytest.approx(0.0)
    assert any("collapsed to a point" in w for w in pc_drift_warnings(drift))


def test_pc_drift_only_reports_pcs_present_on_both_sides():
    train = _train_frame()
    drift = pc_drift(train, pd.DataFrame({"PC1": [0.0, 1.0]}))
    assert list(drift["pc"]) == ["PC1"]


def test_formatting_handles_a_missing_value():
    drift = pd.DataFrame([{
        "pc": "PC1", "train_mean": 0.0, "train_sd": 0.0,
        "cohort_mean": 1.0, "cohort_sd": 1.0, "shift_sd": None, "sd_ratio": None,
    }])
    text = format_pc_drift(drift)
    assert "n/a" in text
    assert "no PCs in common" in format_pc_drift(pd.DataFrame())


# --- admixture override -----------------------------------------------------


def _decisions(pre, final, margin):
    return pd.DataFrame({
        "FID": range(len(pre)), "IID": range(len(pre)),
        "label_pre_admixture": pre, "label": final, "margin_to_all": margin,
    })


def test_summary_keeps_what_the_classifier_said():
    decisions = _decisions(
        pre=["EUR", "AFR", "EAS", "EUR"],
        final=["CAH", "CAH", "CAH", "CAH"],
        margin=[-2.0, -3.0, -1.0, -4.0],
    )
    summary = admixture_summary(decisions)
    assert summary["fraction_cah"] == 1.0
    assert summary["n_distinct_pre_admixture"] == 3
    assert summary["overridden_from"] == {"EUR": 2, "AFR": 1, "EAS": 1}
    assert summary["margin_to_all_quantiles"]["p50"] < 0


def test_an_all_cah_run_with_a_working_classifier_points_away_from_the_model():
    """The whole reason both labels are kept."""
    summary = admixture_summary(
        _decisions(["EUR", "AFR", "EAS"], ["CAH"] * 3, [-1.0, -1.0, -1.0])
    )
    (warning,) = admixture_warnings(summary)
    assert "3 distinct labels" in warning
    assert "not the thing that failed" in warning
    assert "--no-admixture-detection" in warning


def test_a_normal_admixture_rate_is_not_warned_about():
    summary = admixture_summary(
        _decisions(["EUR"] * 10, ["EUR"] * 9 + ["CAH"], [1.0] * 9 + [-1.0])
    )
    assert summary["fraction_cah"] == pytest.approx(0.1)
    assert admixture_warnings(summary) == []


def test_summary_of_nothing_is_empty_not_an_error():
    assert admixture_summary(pd.DataFrame())["n_samples"] == 0


# --- the collector ----------------------------------------------------------


def test_collector_deduplicates_warnings_and_serializes():
    diagnostics = AncestryDiagnostics(
        snp_overlap=SnpOverlapReport(100, 90, "ref-mean"),
        allele_frequency={"n_snps": 90},
        pc_drift=pd.DataFrame([{"pc": "PC1", "shift_sd": 0.1}]),
        admixture={"n_cah": 0},
        self_test={"balanced_accuracy": 0.9},
        plots=["/tmp/x.png"],
    )
    diagnostics.add_warnings(["same", "same", "other"])
    assert diagnostics.warnings == ["same", "other"]

    out = diagnostics.to_dict()
    assert out["snp_overlap"]["n_matched"] == 90
    assert out["pc_drift"]["pc"] == ["PC1"]
    assert out["self_test"]["balanced_accuracy"] == 0.9
    assert out["plots"] == ["/tmp/x.png"]
    assert out["warnings"] == ["same", "other"]


def test_an_empty_collector_serializes_to_just_its_warnings():
    assert AncestryDiagnostics().to_dict() == {"warnings": []}
