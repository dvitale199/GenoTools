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


def _write_vcf(path: Path, records: list, n_samples: int) -> None:
    """Minimal VCF writer. `records` are (pos, snp_id, ref, alt, genotypes)."""
    samples = [f"S{i:03d}" for i in range(n_samples)]
    lines = [
        "##fileformat=VCFv4.2",
        "##contig=<ID=1>",
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        "#" + "\t".join(["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"] + samples),
    ]
    for pos, snp_id, ref, alt, gts in records:
        lines.append("\t".join(["1", str(pos), snp_id, ref, alt, ".", ".", ".", "GT"] + gts))
    path.write_text("\n".join(lines) + "\n")


@pytest.fixture(scope="session")
def dup_position_bfiles(tmp_path_factory) -> tuple:
    """Two bfiles sharing positions, where the *first* carries 1:1000 under two
    probe IDs with very different call rates.

    `probe_hi` (50% missing) is deliberately placed ahead of `probe_lo` (5%
    missing) in the .bim, so first-row-wins dedup selects the badly-called
    probe and a missingness-aware tie-break selects the good one. Without the
    ordering the test could pass for the wrong reason.
    """
    n = 20
    hi = ["./." if i < 10 else "0/1" for i in range(n)]          # 50% missing
    lo = ["./." if i < 1 else "0/1" for i in range(n)]           # 5% missing
    het = ["0/1"] * n
    shared = [(2000, "shared_a", "A", "G", het), (3000, "shared_b", "C", "A", het)]

    d = tmp_path_factory.mktemp("dup_pos")
    _write_vcf(d / "one.vcf", [(1000, "probe_hi", "A", "G", hi),
                               (1000, "probe_lo", "A", "G", lo)] + shared, n)
    _write_vcf(d / "two.vcf", [(1000, "ref_probe", "A", "G", het)] + shared, n)

    out = []
    for name in ("one", "two"):
        prefix = d / name
        run_command(
            f"{get_plink2()} --vcf {prefix}.vcf --make-bed --double-id --out {prefix}",
            tool_name="plink2",
        )
        out.append(prefix)
    return tuple(out)


@pytest.fixture(scope="session")
def palindrome_free_bfile(synth_ref_bfile: Path, tmp_path_factory) -> Path:
    """`synth_ref_bfile` with its palindromic (A/T, C/G) sites removed, so the
    new and legacy get_common_snps must agree exactly on it."""
    bim = pd.read_csv(f"{synth_ref_bfile}.bim", sep="\t", header=None)
    keep = bim.loc[~(bim[4] + bim[5]).isin({"AT", "TA", "CG", "GC"}), 1]
    d = tmp_path_factory.mktemp("no_palindrome")
    keep_file = d / "keep.snplist"
    keep.to_csv(keep_file, header=False, index=False)
    out = d / "ref_nopal"
    run_command(
        f"{get_plink2()} --bfile {synth_ref_bfile} --extract {keep_file} --make-bed --out {out}",
        tool_name="plink2",
    )
    return out
