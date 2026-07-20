from pathlib import Path

import pandas as pd


def test_split_cohort_matches_legacy(synth_geno_pfile, tmp_path):
    """Retained + insufficient-N pruned branches match legacy byte-for-byte."""
    from genotools.ancestry.cohort import split_cohort_by_ancestry
    from genotools.ancestry.legacy import Ancestry

    geno = str(synth_geno_pfile)
    psam = pd.read_csv(f"{geno}.psam", sep=r"\s+")
    fid = psam.get("#FID", psam.get("FID"))
    iid = psam["IID"]
    # 3 groups: a tiny "EAS" (3 samples, below min_samples) exercises the
    # insufficient-N pruned branch; EUR/AFR are retained.
    label = ["EAS" if i < 3 else ("EUR" if i % 2 == 0 else "AFR") for i in range(len(iid))]
    labels = pd.DataFrame({"FID": fid, "IID": iid, "label": label})
    labels_path = tmp_path / "pred_labels.txt"
    labels.to_csv(labels_path, sep="\t", index=False)

    ld = tmp_path / "legacy"; ld.mkdir()
    anc = Ancestry()
    anc.geno_path = geno
    anc.out_path = str(ld / "out")
    anc.subset = None
    anc.min_samples = 10
    legacy_res = anc.split_cohort_ancestry(labels_path=str(labels_path))

    nd = tmp_path / "new"; nd.mkdir()
    new_res = split_cohort_by_ancestry(
        labels_path=str(labels_path), geno_path=geno,
        out_path=str(nd / "out"), min_samples=10, subset=None,
    )

    def _sorted(d):
        return d.sort_values(["FID", "IID"]).reset_index(drop=True)

    assert sorted(new_res["labels"]) == sorted(legacy_res["labels"])
    assert "EAS" not in new_res["labels"]  # pruned for insufficient N
    assert set(new_res["pruned_samples"]["label"]) == {"EAS"}
    pd.testing.assert_frame_equal(
        _sorted(new_res["pruned_samples"]), _sorted(legacy_res["pruned_samples"]),
    )
    for lbl in new_res["labels"]:
        assert Path(f"{nd / 'out'}_{lbl}.pgen").exists()


def test_split_cohort_subset_matches_legacy(synth_geno_pfile, tmp_path):
    """subset restricts which ancestries are split (the subset branch)."""
    from genotools.ancestry.cohort import split_cohort_by_ancestry
    from genotools.ancestry.legacy import Ancestry

    geno = str(synth_geno_pfile)
    psam = pd.read_csv(f"{geno}.psam", sep=r"\s+")
    fid = psam.get("#FID", psam.get("FID"))
    iid = psam["IID"]
    labels = pd.DataFrame({
        "FID": fid, "IID": iid,
        "label": ["EUR" if i % 2 == 0 else "AFR" for i in range(len(iid))],
    })
    labels_path = tmp_path / "pred_labels.txt"
    labels.to_csv(labels_path, sep="\t", index=False)

    ld = tmp_path / "legacy"; ld.mkdir()
    anc = Ancestry()
    anc.geno_path = geno
    anc.out_path = str(ld / "out")
    anc.subset = ["EUR"]
    anc.min_samples = 10
    legacy_res = anc.split_cohort_ancestry(labels_path=str(labels_path))

    nd = tmp_path / "new"; nd.mkdir()
    new_res = split_cohort_by_ancestry(
        labels_path=str(labels_path), geno_path=geno,
        out_path=str(nd / "out"), min_samples=10, subset=["EUR"],
    )
    assert new_res["labels"] == legacy_res["labels"] == ["EUR"]
