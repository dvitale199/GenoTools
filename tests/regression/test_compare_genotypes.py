# Copyright 2023 The GenoTools Authors. All Rights Reserved.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
# ==============================================================================

"""Tests for genotype-content comparison used by the old-vs-new parity harness.

compare_pfiles only diffs sample/variant ID *sets*; a regression that keeps the
same IDs but corrupts genotypes or flips allele coding would pass unnoticed.
These tests cover the genotype-content comparator that closes that gap.
"""

import shutil
from pathlib import Path

import pytest

from tests.regression.compare import _compare_traw_files, compare_genotypes


def _write_traw(path: Path, samples: list[str], rows: list[list]) -> None:
    """Write a minimal PLINK2 .traw file.

    Each row is [CHR, SNP, (C)M, POS, COUNTED, ALT, geno_for_each_sample...].
    """
    header = ["CHR", "SNP", "(C)M", "POS", "COUNTED", "ALT", *samples]
    lines = ["\t".join(header)]
    for r in rows:
        lines.append("\t".join(str(x) for x in r))
    path.write_text("\n".join(lines) + "\n")


class TestCompareTrawFiles:
    """Pure genotype-matrix comparison (no PLINK required)."""

    def test_identical_matrices_are_equal(self, tmp_path: Path) -> None:
        samples = ["F1_I1", "F2_I2"]
        rows = [
            [1, "rs1", 0, 100, "A", "G", 0, 1],
            [1, "rs2", 0, 200, "C", "T", 2, 0],
        ]
        a = tmp_path / "a.traw"
        b = tmp_path / "b.traw"
        _write_traw(a, samples, rows)
        _write_traw(b, samples, rows)

        result = _compare_traw_files(a, b)

        assert result.equal is True
        assert result.sample_diff == 0
        assert result.variant_diff == 0

    def test_differing_genotype_is_detected(self, tmp_path: Path) -> None:
        samples = ["F1_I1", "F2_I2"]
        a = tmp_path / "a.traw"
        b = tmp_path / "b.traw"
        _write_traw(a, samples, [[1, "rs1", 0, 100, "A", "G", 0, 1]])
        _write_traw(b, samples, [[1, "rs1", 0, 100, "A", "G", 0, 2]])  # I2 differs

        result = _compare_traw_files(a, b)

        assert result.equal is False
        assert "rs1" in result.mismatched_variants

    def test_allele_flip_is_detected(self, tmp_path: Path) -> None:
        """Same genotype counts but swapped COUNTED/ALT alleles must not compare equal."""
        samples = ["F1_I1"]
        a = tmp_path / "a.traw"
        b = tmp_path / "b.traw"
        _write_traw(a, samples, [[1, "rs1", 0, 100, "A", "G", 0]])
        _write_traw(b, samples, [[1, "rs1", 0, 100, "G", "A", 0]])  # flipped coding

        result = _compare_traw_files(a, b)

        assert result.equal is False
        assert "rs1" in result.mismatched_variants

    def test_matching_missing_values_are_equal(self, tmp_path: Path) -> None:
        """NA in the same cell of both files counts as equal, not a mismatch."""
        samples = ["F1_I1", "F2_I2"]
        rows = [[1, "rs1", 0, 100, "A", "G", "NA", 1]]
        a = tmp_path / "a.traw"
        b = tmp_path / "b.traw"
        _write_traw(a, samples, rows)
        _write_traw(b, samples, rows)

        result = _compare_traw_files(a, b)

        assert result.equal is True

    def test_missing_variant_is_detected(self, tmp_path: Path) -> None:
        samples = ["F1_I1"]
        a = tmp_path / "a.traw"
        b = tmp_path / "b.traw"
        _write_traw(a, samples, [[1, "rs1", 0, 100, "A", "G", 0], [1, "rs2", 0, 200, "C", "T", 1]])
        _write_traw(b, samples, [[1, "rs1", 0, 100, "A", "G", 0]])  # rs2 missing

        result = _compare_traw_files(a, b)

        assert result.equal is False
        assert result.variant_diff >= 1


# PLINK2 executable for the integration-style test (skip if unavailable)
def _find_plink2():
    from_path = shutil.which("plink2")
    if from_path:
        return from_path
    candidate = Path.home() / ".genotools" / "misc" / "executables" / "plink2"
    return str(candidate) if candidate.exists() else None


class TestCompareGenotypesWithPlink:
    """PLINK-backed genotype comparison on real pfiles."""

    def test_identical_pfiles_compare_equal(self, test_geno_path: Path) -> None:
        plink2 = _find_plink2()
        if plink2 is None:
            pytest.skip("plink2 not available")

        result = compare_genotypes(test_geno_path, test_geno_path, plink2_exec=plink2)

        assert result.equal is True
        assert result.variant_diff == 0
        assert result.sample_diff == 0


_VCF_HEADER = """##fileformat=VCFv4.2
##contig=<ID=1,length=1000>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3
"""


def _vcf_to_pgen(
    tmp_path: Path,
    name: str,
    calls: list[str],
    plink2: str,
    flip_first: bool = False,
) -> Path:
    """Build a 3-sample / 3-variant pfile from explicit genotype calls.

    ``calls`` is one "GT GT GT" string per variant, so a test can change a
    single call and nothing else.
    """
    import subprocess

    variants = [("100", "v1", "A", "G"), ("200", "v2", "C", "T"), ("300", "v3", "G", "A")]
    if flip_first:
        # Swap REF/ALT on v1 only; the caller recodes its calls to match.
        pos, vid, ref, alt = variants[0]
        variants[0] = (pos, vid, alt, ref)
    rows = [
        "\t".join([c, pos, vid, ref, alt, ".", ".", ".", "GT", *gt.split()])
        for (pos, vid, ref, alt), gt, c in zip(variants, calls, "111")
    ]
    vcf = tmp_path / f"{name}.vcf"
    vcf.write_text(_VCF_HEADER + "\n".join(rows) + "\n")

    prefix = tmp_path / name
    subprocess.run(
        [plink2, "--vcf", str(vcf), "--make-pgen", "--out", str(prefix)],
        capture_output=True,
        text=True,
        check=True,
    )
    return prefix


@pytest.fixture()
def plink2_exec() -> str:
    plink2 = _find_plink2()
    if plink2 is None:
        pytest.skip("plink2 not available")
    return plink2


class TestCompareGenotypesMethodsAgree:
    """The pgen-diff and traw methods must reach the same verdict.

    compare_genotypes() used to always export both filesets to .traw text,
    which costs roughly (samples x variants x 2) bytes per side - enough to
    exhaust the disk on a real cohort. The default is now PLINK2 --pgen-diff,
    which walks the filesets in place. These tests pin the two methods
    together so the cheap path cannot quietly become the lenient one.
    """

    BASE = ["0/0 0/1 1/1", "0/1 0/1 0/0", "1/1 0/0 0/1"]

    @pytest.mark.parametrize("method", ["pgen-diff", "traw"])
    def test_identical_is_equal(self, tmp_path: Path, plink2_exec: str, method: str) -> None:
        a = _vcf_to_pgen(tmp_path, "a", self.BASE, plink2_exec)

        result = compare_genotypes(a, a, plink2_exec=plink2_exec, method=method)

        assert result.equal is True
        assert result.sample_diff == 0
        assert result.variant_diff == 0

    @pytest.mark.parametrize("method", ["pgen-diff", "traw"])
    def test_changed_call_is_detected(self, tmp_path: Path, plink2_exec: str, method: str) -> None:
        a = _vcf_to_pgen(tmp_path, "a", self.BASE, plink2_exec)
        # v2 sample S2: 0/1 -> 1/1
        changed = list(self.BASE)
        changed[1] = "0/1 1/1 0/0"
        b = _vcf_to_pgen(tmp_path, "b", changed, plink2_exec)

        result = compare_genotypes(a, b, plink2_exec=plink2_exec, method=method)

        assert result.equal is False
        assert result.mismatched_variants == ["v2"]

    @pytest.mark.parametrize("method", ["pgen-diff", "traw"])
    def test_call_becoming_missing_is_detected(
        self, tmp_path: Path, plink2_exec: str, method: str
    ) -> None:
        """A call present on one side and missing on the other is a difference."""
        a = _vcf_to_pgen(tmp_path, "a", self.BASE, plink2_exec)
        # v3 sample S3: 0/1 -> ./.
        missing = list(self.BASE)
        missing[2] = "1/1 0/0 ./."
        b = _vcf_to_pgen(tmp_path, "b", missing, plink2_exec)

        result = compare_genotypes(a, b, plink2_exec=plink2_exec, method=method)

        assert result.equal is False
        assert result.mismatched_variants == ["v3"]

    @pytest.mark.parametrize("method", ["pgen-diff", "traw"])
    def test_recoded_allele_flip_is_detected(
        self, tmp_path: Path, plink2_exec: str, method: str
    ) -> None:
        """A REF/ALT swap with genotypes recoded to match must still fail.

        The underlying calls are identical - only the representation differs -
        so nothing in the genotype content itself is wrong. It is still a
        regression worth catching, and PLINK2 refuses to run --pgen-diff at all
        when REF alleles conflict, so the .pvar allele comparison is what
        produces the verdict.

        Deliberately built from explicit VCFs rather than plink2 --maj-ref:
        that flag's effect on the written alleles differs between PLINK2
        v2.00a5.10 and v2.0.0-a.6.3, which would make this test's meaning
        depend on which plink2 happens to be resolved.
        """
        # v1 A/G with 0/0,0/1,1/1 == v1 G/A with 1/1,0/1,0/0 (same alleles).
        a = _vcf_to_pgen(tmp_path, "a", ["0/0 0/1 1/1", "0/1 0/1 0/0", "1/1 0/0 0/1"], plink2_exec)
        b = _vcf_to_pgen(
            tmp_path,
            "b",
            ["1/1 0/1 0/0", "0/1 0/1 0/0", "1/1 0/0 0/1"],
            plink2_exec,
            flip_first=True,
        )

        result = compare_genotypes(a, b, plink2_exec=plink2_exec, method=method)

        assert result.equal is False
        assert "v1" in result.mismatched_variants

    @pytest.mark.parametrize("method", ["pgen-diff", "traw"])
    def test_prefix_containing_dots(
        self, tmp_path: Path, plink2_exec: str, method: str
    ) -> None:
        """A prefix like "cohort.r12" must not lose its last dotted segment.

        Path.with_suffix() would turn "cohort.r12" into "cohort.psam"; PLINK
        appends, so the comparison has to as well.
        """
        a = _vcf_to_pgen(tmp_path, "cohort.r12", self.BASE, plink2_exec)
        assert (tmp_path / "cohort.r12.psam").exists()

        result = compare_genotypes(a, a, plink2_exec=plink2_exec, method=method)

        assert result.equal is True
