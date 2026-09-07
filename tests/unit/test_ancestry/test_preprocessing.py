from pathlib import Path

import pandas as pd
import pytest

from .conftest import GOLDEN


def _read_bim_snps(prefix: str) -> set:
    bim = pd.read_csv(f"{prefix}.bim", sep=r"\s+", header=None)
    return set(bim[1])


PALINDROMIC = {"AT", "TA", "CG", "GC"}


def _palindromic_snps(prefix: str) -> set:
    bim = pd.read_csv(f"{prefix}.bim", sep=r"\s+", header=None)
    return set(bim.loc[(bim[4] + bim[5]).isin(PALINDROMIC), 1])


def test_get_common_snps_matches_legacy_on_unambiguous_input(palindrome_free_bfile, tmp_path):
    """On input with nothing ambiguous to decide, the port is still exactly the
    legacy function — which is what the old version of this test was really
    checking, minus the palindromes that masked the difference."""
    from genotools.ancestry.preprocessing import get_common_snps as new_fn
    from genotools.utils import get_common_snps as legacy_fn

    g = str(palindrome_free_bfile)
    assert not _palindromic_snps(g), "fixture should have no palindromic sites"

    legacy_out, new_out = tmp_path / "legacy_common", tmp_path / "new_common"
    legacy_res = legacy_fn(g, g, str(legacy_out))
    new_res = new_fn(g, g, str(new_out))

    assert set(new_res.keys()) == set(legacy_res.keys())
    assert _read_bim_snps(str(new_out)) == _read_bim_snps(str(legacy_out))
    assert (Path(f"{new_out}.common_snps").read_text()
            == Path(f"{legacy_out}.common_snps").read_text())


def test_get_common_snps_excludes_palindromic(synth_ref_bfile, tmp_path):
    """Palindromic sites are dropped, and legacy is shown to have kept them.

    Their alleles survive strand complement unchanged, so nothing downstream
    can tell which strand they came from -- get_raw_files' allele-switch test
    can never fire for one. The legacy assertion is what makes this a
    regression test rather than a restatement of the implementation.
    """
    from genotools.ancestry.preprocessing import get_common_snps as new_fn
    from genotools.utils import get_common_snps as legacy_fn

    g = str(synth_ref_bfile)
    palindromic = _palindromic_snps(g)
    assert palindromic, "fixture must contain palindromic sites to be meaningful"

    new_res = new_fn(g, g, str(tmp_path / "new_common"))
    kept = set(pd.read_csv(new_res["common_snps"], header=None)[0])
    assert not (kept & palindromic), "palindromic sites reached the common set"

    legacy_res = legacy_fn(g, g, str(tmp_path / "legacy_common"))
    legacy_kept = set(pd.read_csv(legacy_res["common_snps"], header=None)[0])
    assert legacy_kept & palindromic, "legacy should have kept them (test is degenerate)"
    assert legacy_kept - kept == palindromic, "only palindromic sites should differ"


def test_get_common_snps_breaks_position_ties_by_missingness(dup_position_bfiles, tmp_path):
    """When one position offers several probes, the best-called one wins.

    The fixture puts the 50%-missing probe ahead of the 5%-missing one, so
    first-row-wins dedup would take `probe_hi`.
    """
    from genotools.ancestry.preprocessing import get_common_snps

    one, two = dup_position_bfiles
    res = get_common_snps(str(one), str(two), str(tmp_path / "out"))
    kept = set(pd.read_csv(res["common_snps"], header=None)[0])

    assert "probe_lo" in kept, "the better-called probe should have been chosen"
    assert "probe_hi" not in kept, "the 50%-missing probe should have lost"
    # one row per position, not one per probe
    bim = pd.read_csv(f"{tmp_path / 'out'}.bim", sep=r"\s+", header=None)
    assert not bim.duplicated(subset=[0, 3]).any(), "a position survived twice"


def _sorted_df(df):
    df = df.sort_values(["FID", "IID"]).reset_index(drop=True)
    return df.reindex(sorted(df.columns), axis=1)


def test_get_raw_files_train_golden(ref21_22_bfile, ref21_22_labels, geno21_22_pfile, tmp_path):
    """Train-mode raw_ref/raw_geno match the committed golden (== legacy at capture)."""
    from genotools.ancestry.preprocessing import get_raw_files
    new_dir = tmp_path / "new"; new_dir.mkdir()
    res = get_raw_files(
        geno_path=str(geno21_22_pfile), ref_panel=str(ref21_22_bfile),
        ref_labels=str(ref21_22_labels), out_path=str(new_dir / "out"), train=True,
    )
    pd.testing.assert_frame_equal(
        _sorted_df(res["raw_ref"]), pd.read_parquet(GOLDEN / "raw_ref_train.parquet"))
    pd.testing.assert_frame_equal(
        _sorted_df(res["raw_geno"]), pd.read_parquet(GOLDEN / "raw_geno_train.parquet"))


@pytest.fixture(scope="module")
def inference_inputs(ref21_22_bfile, ref21_22_labels, geno21_22_pfile, tmp_path_factory):
    """A model SNP list plus a cohort missing chr22, so the fill path must run.

    Module-scoped: the two PLINK stages behind it are the slow part of these
    tests, and every fill-strategy test wants the same inputs.
    """
    import shutil
    from genotools.ancestry.preprocessing import get_raw_files
    from genotools.core.executors import run_command, get_plink2

    d = tmp_path_factory.mktemp("inference_inputs")

    # common_snps from a train run on the full (chr21+22) reduced geno
    prep_dir = d / "prep"; prep_dir.mkdir()
    train_res = get_raw_files(
        geno_path=str(geno21_22_pfile), ref_panel=str(ref21_22_bfile),
        ref_labels=str(ref21_22_labels), out_path=str(prep_dir / "out"), train=True)
    common_snps = d / "model.common_snps"
    shutil.copy2(train_res["out_paths"]["common_snps"], common_snps)

    # inference geno MISSING chr22 -> ref chr22 common SNPs must be filled
    subset_geno = d / "subset_geno"
    run_command(
        f"{get_plink2()} --pfile {geno21_22_pfile} --not-chr 22 "
        f"--make-pgen psam-cols=fid,parents,sex,pheno1,phenos --out {subset_geno}",
        tool_name="plink2")
    return subset_geno, common_snps


def _run_inference(inputs, ref_bfile, ref_labels, out_dir, **kwargs):
    from genotools.ancestry.preprocessing import get_raw_files

    subset_geno, common_snps = inputs
    out_dir.mkdir(parents=True, exist_ok=True)
    return get_raw_files(
        geno_path=str(subset_geno), ref_panel=str(ref_bfile),
        ref_labels=str(ref_labels), out_path=str(out_dir / "out"), train=False,
        common_snps_file=str(common_snps), **kwargs)


def _filled_columns(golden: pd.DataFrame) -> list:
    """Columns the legacy golden filled with a constant 2 for every sample."""
    numeric = golden.drop(columns=["FID", "IID", "label"])
    return list(numeric.columns[(numeric == 2).all(axis=0)])


def test_get_raw_files_inference_constant_fill_matches_golden(
    inference_inputs, ref21_22_bfile, ref21_22_labels, tmp_path
):
    """`constant` reproduces the legacy fill exactly, int dtypes included.

    The golden was captured against 1.x's `np.repeat(2, ...)`, so this is the
    escape hatch's parity test. The overlap assertions keep it honest: without
    them it could pass on a cohort that needed no fill at all.
    """
    res = _run_inference(
        inference_inputs, ref21_22_bfile, ref21_22_labels, tmp_path / "constant",
        fill_strategy="constant",
    )
    golden = pd.read_parquet(GOLDEN / "raw_geno_inference.parquet")
    overlap = res["diagnostics"].snp_overlap
    assert overlap.n_filled == len(_filled_columns(golden)) > 0
    assert overlap.n_matched + overlap.n_filled == overlap.n_model_snps
    pd.testing.assert_frame_equal(_sorted_df(res["raw_geno"]), golden)


def test_get_raw_files_inference_defaults_to_ref_mean_fill(
    inference_inputs, ref21_22_bfile, ref21_22_labels, tmp_path
):
    """The default leaves absent SNPs NaN for the reference-fitted imputer.

    Dosage 2 is one end of the range applied identically to every sample, so a
    cohort missing much of the model's SNP list gets a large shared offset in PC
    space -- the mechanism behind an all-CAH prediction. NaN instead reaches
    `PCAReducer.transform`'s imputer, which fills the reference panel's own mean
    and so contributes nothing to the projection. Columns the cohort *did*
    carry must be untouched by the change.
    """
    res = _run_inference(
        inference_inputs, ref21_22_bfile, ref21_22_labels, tmp_path / "refmean",
    )
    raw = _sorted_df(res["raw_geno"])
    golden = pd.read_parquet(GOLDEN / "raw_geno_inference.parquet")
    filled = _filled_columns(golden)
    assert filled, "fixture no longer exercises the fill path"

    assert raw[filled].isna().all().all(), "absent SNPs should be NaN, not a dosage"
    carried = [c for c in golden.columns if c not in set(filled)]
    pd.testing.assert_frame_equal(raw[carried], golden[carried])

    overlap = res["diagnostics"].snp_overlap
    assert overlap.fill_strategy == "ref-mean"
    assert overlap.n_filled == len(filled)
    assert set(overlap.filled_snps) == set(filled)


def test_get_raw_files_writes_filled_snp_ids(
    inference_inputs, ref21_22_bfile, ref21_22_labels, tmp_path
):
    """The filled IDs are written where a user can read them."""
    keep = tmp_path / "keep"; keep.mkdir()
    res = _run_inference(
        inference_inputs, ref21_22_bfile, ref21_22_labels, tmp_path / "run",
        diagnostics_prefix=str(keep / "cohort"),
    )
    path = Path(res["diagnostics"].filled_snps_path)
    assert path == keep / "cohort_filled_snps.txt"
    written = [line.strip() for line in path.read_text().splitlines() if line.strip()]
    assert written == list(res["diagnostics"].snp_overlap.filled_snps)


def test_get_raw_files_refuses_a_mostly_filled_matrix(
    inference_inputs, ref21_22_bfile, ref21_22_labels, tmp_path
):
    """Above the limit, prediction is refused rather than reported.

    A label predicted from a feature matrix that is mostly fill describes the
    fill. The message has to name the numbers and the override, because the
    alternative -- what shipped before this round -- was a confident label.
    """
    from genotools.core.exceptions import AncestryError

    with pytest.raises(AncestryError) as excinfo:
        _run_inference(
            inference_inputs, ref21_22_bfile, ref21_22_labels, tmp_path / "strict",
            max_missing_fraction=0.01,
        )
    message = str(excinfo.value)
    assert "of the model's" in message
    assert "--ancestry-max-missing-snps" in message


def test_get_raw_files_rejects_an_unknown_fill_strategy(
    inference_inputs, ref21_22_bfile, ref21_22_labels, tmp_path
):
    from genotools.core.exceptions import AncestryError

    with pytest.raises(AncestryError, match="Unknown fill strategy"):
        _run_inference(
            inference_inputs, ref21_22_bfile, ref21_22_labels, tmp_path / "bogus",
            fill_strategy="zero",
        )


def test_get_raw_files_reports_allele_frequency_concordance(
    inference_inputs, ref21_22_bfile, ref21_22_labels, tmp_path
):
    """Matched sites are compared to the panel before anything is filled.

    Same data on both sides of the match, so agreement should be near-perfect;
    the point is that the number exists and covers the matched SNPs only.
    """
    res = _run_inference(
        inference_inputs, ref21_22_bfile, ref21_22_labels, tmp_path / "af",
    )
    concordance = res["diagnostics"].allele_frequency
    overlap = res["diagnostics"].snp_overlap
    assert concordance["n_snps"] == overlap.n_matched
    assert concordance["correlation"] > 0.99
    assert concordance["n_flip_signature"] == 0


def test_clean_up_files_removes_extensions(tmp_path):
    from genotools.ancestry.preprocessing import clean_up_files
    prefix = tmp_path / "junk"
    for ext in ("bed", "bim", "fam", "raw"):
        (tmp_path / f"junk.{ext}").write_text("x")
    clean_up_files([str(prefix)])
    for ext in ("bed", "bim", "fam", "raw"):
        assert not (tmp_path / f"junk.{ext}").exists()
