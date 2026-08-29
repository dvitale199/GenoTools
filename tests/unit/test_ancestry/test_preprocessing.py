from pathlib import Path

import pandas as pd

from .conftest import GOLDEN


def _read_bim_snps(prefix: str) -> set:
    bim = pd.read_csv(f"{prefix}.bim", sep=r"\s+", header=None)
    return set(bim[1])


def test_get_common_snps_matches_legacy(synth_ref_bfile, tmp_path):
    """New port yields the same extracted SNP set + common_snps file as legacy."""
    from genotools.ancestry.preprocessing import get_common_snps as new_fn
    from genotools.utils import get_common_snps as legacy_fn

    # Use the same bfile as both inputs (guarantees a nonempty common set)
    g1, g2 = str(synth_ref_bfile), str(synth_ref_bfile)

    legacy_out = tmp_path / "legacy_common"
    new_out = tmp_path / "new_common"

    legacy_res = legacy_fn(g1, g2, str(legacy_out))
    new_res = new_fn(g1, g2, str(new_out))

    assert set(new_res.keys()) == set(legacy_res.keys())
    assert _read_bim_snps(str(new_out)) == _read_bim_snps(str(legacy_out))
    assert (Path(f"{new_out}.common_snps").read_text()
            == Path(f"{legacy_out}.common_snps").read_text())


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
