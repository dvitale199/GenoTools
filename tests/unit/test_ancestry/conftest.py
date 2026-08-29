import shutil
from pathlib import Path

import pandas as pd
import pytest

from genotools.core.executors import run_command, get_plink2

SYNTH = Path("tests/data/synthetic/genotools_test")
GOLDEN = Path(__file__).parent / "golden"


@pytest.fixture(scope="session")
def synth_geno_pfile() -> Path:
    assert SYNTH.with_suffix(".pgen").exists(), "run from repo root"
    return SYNTH


@pytest.fixture(scope="session")
def synth_ref_bfile(tmp_path_factory) -> Path:
    """A bfile copy of the synthetic data, used as a reference panel."""
    out = tmp_path_factory.mktemp("refpanel") / "ref"
    run_command(
        f"{get_plink2()} --pfile {SYNTH} --make-bed --out {out}",
        tool_name="plink2",
    )
    return out


@pytest.fixture(scope="session")
def synth_ref_labels(synth_ref_bfile: Path, tmp_path_factory) -> Path:
    """FID/IID/label for the ref-panel samples (alternating EUR/AFR)."""
    fam = pd.read_csv(f"{synth_ref_bfile}.fam", sep=r"\s+", header=None)
    labels = pd.DataFrame({
        "FID": fam[0],
        "IID": fam[1],
        "label": ["EUR" if i % 2 == 0 else "AFR" for i in range(len(fam))],
    })
    path = tmp_path_factory.mktemp("reflabels") / "ref_labels.txt"
    labels.to_csv(path, sep="\t", header=False, index=False)
    return path


@pytest.fixture(scope="session")
def ref21_22_bfile(tmp_path_factory) -> Path:
    """Chr21+22 reference bfile, staged into a temp dir from the committed
    golden copy. get_common_snps/get_raw_files write side-effect files
    (.snplist, _flip.*, .log) next to this prefix, so it must never be the
    committed golden/ path itself."""
    src = GOLDEN / "ref21_22"
    assert src.with_suffix(".bed").exists(), "run the golden generator first"
    dst = tmp_path_factory.mktemp("ref21_22_src") / "ref21_22"
    for ext in (".bed", ".bim", ".fam"):
        shutil.copy2(src.with_suffix(ext), dst.with_suffix(ext))
    return dst


@pytest.fixture(scope="session")
def ref21_22_labels() -> Path:
    return GOLDEN / "ref21_22_labels.txt"


@pytest.fixture(scope="session")
def geno21_22_pfile(tmp_path_factory) -> Path:
    """Chr21+22 geno pfile (all 500 samples, small), staged into a temp dir
    from the committed golden copy for the same reason as ref21_22_bfile."""
    src = GOLDEN / "geno21_22"
    assert src.with_suffix(".pgen").exists(), "run the golden generator first"
    dst = tmp_path_factory.mktemp("geno21_22_src") / "geno21_22"
    for ext in (".pgen", ".pvar", ".psam"):
        shutil.copy2(src.with_suffix(ext), dst.with_suffix(ext))
    return dst
