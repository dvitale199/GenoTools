from pathlib import Path

import pandas as pd


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


def test_get_raw_files_train_matches_legacy(synth_ref_bfile, synth_ref_labels, synth_geno_pfile, tmp_path):
    """Train-mode raw_ref/raw_geno match legacy Ancestry.get_raw_files byte-for-byte."""
    from genotools.ancestry.preprocessing import get_raw_files as new_fn
    from genotools.ancestry.legacy import Ancestry

    ref, labels, geno = str(synth_ref_bfile), str(synth_ref_labels), str(synth_geno_pfile)

    # Legacy
    legacy_dir = tmp_path / "legacy"; legacy_dir.mkdir()
    anc = Ancestry()
    anc.geno_path = geno
    anc.ref_panel = ref
    anc.ref_labels = labels
    anc.out_path = str(legacy_dir / "out")
    anc.train = True
    anc.model_path = None
    anc.containerized = False
    legacy_res = anc.get_raw_files()

    # New
    new_dir = tmp_path / "new"; new_dir.mkdir()
    new_res = new_fn(
        geno_path=geno, ref_panel=ref, ref_labels=labels,
        out_path=str(new_dir / "out"), train=True,
    )

    pd.testing.assert_frame_equal(
        _sorted_df(legacy_res["raw_ref"]), _sorted_df(new_res["raw_ref"]),
    )
    pd.testing.assert_frame_equal(
        _sorted_df(legacy_res["raw_geno"]), _sorted_df(new_res["raw_geno"]),
    )


def test_get_raw_files_inference_matches_legacy(synth_ref_bfile, synth_ref_labels, synth_geno_pfile, tmp_path):
    """Inference-mode raw_geno matches legacy (extract-from-common-SNPs + reorder path)."""
    import shutil
    from genotools.ancestry.preprocessing import get_raw_files as new_fn
    from genotools.ancestry.legacy import Ancestry

    ref, labels, geno = str(synth_ref_bfile), str(synth_ref_labels), str(synth_geno_pfile)

    # Produce a common_snps file via a train run (new impl)
    prep = tmp_path / "prep"; prep.mkdir()
    train_res = new_fn(geno_path=geno, ref_panel=ref, ref_labels=labels,
                       out_path=str(prep / "out"), train=True)
    common_snps_src = train_res["out_paths"]["common_snps"]

    # Legacy inference derives the common_snps path from model_path: <dir>/model.common_snps
    ld = tmp_path / "legacy"; ld.mkdir()
    shutil.copy2(common_snps_src, ld / "model.common_snps")
    anc = Ancestry()
    anc.geno_path, anc.ref_panel, anc.ref_labels = geno, ref, labels
    anc.out_path = str(ld / "out")
    anc.train = False
    anc.model_path = str(ld / "model.pkl")
    anc.containerized = False
    legacy_res = anc.get_raw_files()

    # New inference takes the common_snps file directly
    nd = tmp_path / "new"; nd.mkdir()
    shutil.copy2(common_snps_src, nd / "model.common_snps")
    new_res = new_fn(geno_path=geno, ref_panel=ref, ref_labels=labels,
                     out_path=str(nd / "out"), train=False,
                     common_snps_file=str(nd / "model.common_snps"))

    pd.testing.assert_frame_equal(
        _sorted_df(legacy_res["raw_geno"]), _sorted_df(new_res["raw_geno"]),
    )


def test_clean_up_files_removes_extensions(tmp_path):
    from genotools.ancestry.preprocessing import clean_up_files
    prefix = tmp_path / "junk"
    for ext in ("bed", "bim", "fam", "raw"):
        (tmp_path / f"junk.{ext}").write_text("x")
    clean_up_files([str(prefix)])
    for ext in ("bed", "bim", "fam", "raw"):
        assert not (tmp_path / f"junk.{ext}").exists()
