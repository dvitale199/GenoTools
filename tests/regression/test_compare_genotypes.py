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
