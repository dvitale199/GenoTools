from pathlib import Path

import pandas as pd

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


def test_get_raw_files_inference_golden(ref21_22_bfile, ref21_22_labels, geno21_22_pfile, tmp_path):
    """Inference-mode raw_geno matches golden; the np.repeat spy proves the
    missing-column fill loop ran (chr22 excluded from the inference geno)."""
    import shutil
    from unittest.mock import patch
    import genotools.ancestry.preprocessing as prep
    from genotools.ancestry.preprocessing import get_raw_files
    from genotools.core.executors import run_command, get_plink2

    # common_snps from a train run on the full (chr21+22) reduced geno
    prep_dir = tmp_path / "prep"; prep_dir.mkdir()
    train_res = get_raw_files(
        geno_path=str(geno21_22_pfile), ref_panel=str(ref21_22_bfile),
        ref_labels=str(ref21_22_labels), out_path=str(prep_dir / "out"), train=True)
    common_snps_src = train_res["out_paths"]["common_snps"]

    # inference geno MISSING chr22 -> ref chr22 common SNPs must be filled
    subset_geno = tmp_path / "subset_geno"
    run_command(
        f"{get_plink2()} --pfile {geno21_22_pfile} --not-chr 22 "
        f"--make-pgen psam-cols=fid,parents,sex,pheno1,phenos --out {subset_geno}",
        tool_name="plink2")

    nd = tmp_path / "new"; nd.mkdir()
    shutil.copy2(common_snps_src, nd / "model.common_snps")
    with patch.object(prep.np, "repeat", wraps=prep.np.repeat) as spy:
        res = get_raw_files(
            geno_path=str(subset_geno), ref_panel=str(ref21_22_bfile),
            ref_labels=str(ref21_22_labels), out_path=str(nd / "out"), train=False,
            common_snps_file=str(nd / "model.common_snps"))
    assert spy.called, "missing-column fill loop was not exercised (degenerate test)"
    pd.testing.assert_frame_equal(
        _sorted_df(res["raw_geno"]), pd.read_parquet(GOLDEN / "raw_geno_inference.parquet"))


def test_clean_up_files_removes_extensions(tmp_path):
    from genotools.ancestry.preprocessing import clean_up_files
    prefix = tmp_path / "junk"
    for ext in ("bed", "bim", "fam", "raw"):
        (tmp_path / f"junk.{ext}").write_text("x")
    clean_up_files([str(prefix)])
    for ext in ("bed", "bim", "fam", "raw"):
        assert not (tmp_path / f"junk.{ext}").exists()
