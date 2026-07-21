from pathlib import Path
import pandas as pd
from tests.regression.compare import compare_genotypes
from .conftest import GOLDEN


def _write_labels(geno: Path, out: Path):
    """Deterministic 3-group labels: tiny EAS (<min_samples) + EUR/AFR."""
    psam = pd.read_csv(f"{geno}.psam", sep=r"\s+")
    fid = psam.get("#FID", psam.get("FID")); iid = psam["IID"]
    label = ["EAS" if i < 3 else ("EUR" if i % 2 == 0 else "AFR") for i in range(len(iid))]
    df = pd.DataFrame({"FID": fid, "IID": iid, "label": label})
    df.to_csv(out, sep="\t", index=False)


def _sorted(d):
    return d.sort_values(["FID", "IID"]).reset_index(drop=True)


def test_split_cohort_golden(geno21_22_pfile, tmp_path):
    """Retained + insufficient-N pruned branches match the committed golden."""
    from genotools.ancestry.cohort import split_cohort_by_ancestry
    labels_path = tmp_path / "pred_labels.txt"
    _write_labels(geno21_22_pfile, labels_path)
    nd = tmp_path / "new"; nd.mkdir()
    res = split_cohort_by_ancestry(
        labels_path=str(labels_path), geno_path=str(geno21_22_pfile),
        out_path=str(nd / "out"), min_samples=10, subset=None)

    assert sorted(res["labels"]) == ["AFR", "EUR"]
    assert "EAS" not in res["labels"]
    assert set(res["pruned_samples"]["label"]) == {"EAS"}
    pd.testing.assert_frame_equal(
        _sorted(res["pruned_samples"]),
        pd.read_parquet(GOLDEN / "cohort_pruned_samples.parquet"))
    for lbl in ("EUR", "AFR"):
        cmp = compare_genotypes(GOLDEN / f"split_{lbl}", Path(f"{nd/'out'}_{lbl}"))
        assert cmp.equal, cmp.message


def test_split_cohort_subset_golden(geno21_22_pfile, tmp_path):
    """subset=['EUR'] restricts which ancestries are split."""
    from genotools.ancestry.cohort import split_cohort_by_ancestry
    labels_path = tmp_path / "pred_labels.txt"
    psam = pd.read_csv(f"{geno21_22_pfile}.psam", sep=r"\s+")
    fid = psam.get("#FID", psam.get("FID")); iid = psam["IID"]
    df = pd.DataFrame({"FID": fid, "IID": iid,
                       "label": ["EUR" if i % 2 == 0 else "AFR" for i in range(len(iid))]})
    df.to_csv(labels_path, sep="\t", index=False)
    nd = tmp_path / "new"; nd.mkdir()
    res = split_cohort_by_ancestry(
        labels_path=str(labels_path), geno_path=str(geno21_22_pfile),
        out_path=str(nd / "out"), min_samples=10, subset=["EUR"])
    assert res["labels"] == ["EUR"]
