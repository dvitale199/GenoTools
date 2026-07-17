from pathlib import Path

import pandas as pd


def test_split_cohort_matches_legacy(synth_geno_pfile, tmp_path):
    from genotools.ancestry.cohort import split_cohort_by_ancestry
    from genotools.ancestry.legacy import Ancestry

    geno = str(synth_geno_pfile)
    # Build a predicted-labels file over the real sample IDs
    psam = pd.read_csv(f"{geno}.psam", sep=r"\s+")
    fid = psam.get("#FID", psam.get("FID"))
    iid = psam["IID"]
    labels = pd.DataFrame({
        "FID": fid, "IID": iid,
        "label": ["EUR" if i % 2 == 0 else "AFR" for i in range(len(iid))],
    })
    labels_path = tmp_path / "pred_labels.txt"
    labels.to_csv(labels_path, sep="\t", index=False)

    # Legacy
    ld = tmp_path / "legacy"; ld.mkdir()
    anc = Ancestry()
    anc.geno_path = geno
    anc.out_path = str(ld / "out")
    anc.subset = None
    anc.min_samples = 10
    legacy_res = anc.split_cohort_ancestry(labels_path=str(labels_path))

    # New
    nd = tmp_path / "new"; nd.mkdir()
    new_res = split_cohort_by_ancestry(
        labels_path=str(labels_path), geno_path=geno,
        out_path=str(nd / "out"), min_samples=10, subset=None,
    )

    assert sorted(new_res["labels"]) == sorted(legacy_res["labels"])
    assert (new_res["pruned_samples"].shape == legacy_res["pruned_samples"].shape)
    # Each retained ancestry produced a pgen
    for label in new_res["labels"]:
        assert Path(f"{nd / 'out'}_{label}.pgen").exists()
