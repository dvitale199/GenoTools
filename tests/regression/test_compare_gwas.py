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

"""Tests for the GWAS association-table comparator used by the parity harness.

GWAS output is a PLINK2 --glm table, not a pfile, so it needs its own
comparator: align on variant ID, compare per-variant p-values and the derived
genomic-inflation lambda. These tests cover that comparator directly (no
.venv-stable / PLINK required).
"""

from pathlib import Path

from tests.regression.compare import (
    compare_gwas,
    compare_gwas_results,
    find_gwas_output,
)

_HEADER = "#CHROM\tPOS\tID\tREF\tALT\tA1\tTEST\tOBS_CT\tOR\tP\n"


def _write_glm(path: Path, pvalues: list[float]) -> None:
    """Write a minimal PLINK2 .glm.logistic.hybrid file with ADD-test rows."""
    lines = [_HEADER.rstrip("\n")]
    for i, p in enumerate(pvalues):
        lines.append(f"1\t{i + 1}\trs{i}\tA\tG\tG\tADD\t500\t1.0\t{p}")
    path.write_text("\n".join(lines) + "\n")


_PVALS = [0.5, 0.1, 0.9, 0.3, 0.7, 0.05, 0.6, 0.2, 0.8, 0.4, 0.55]


def test_identical_glm_files_are_equal(tmp_path: Path) -> None:
    a = tmp_path / "a.PHENO1.glm.logistic.hybrid"
    b = tmp_path / "b.PHENO1.glm.logistic.hybrid"
    _write_glm(a, _PVALS)
    _write_glm(b, _PVALS)

    result = compare_gwas_results(a, b)

    assert result.equal, result.message
    assert result.variant_diff == 0


def test_perturbed_pvalue_is_detected(tmp_path: Path) -> None:
    a = tmp_path / "a.PHENO1.glm.logistic.hybrid"
    b = tmp_path / "b.PHENO1.glm.logistic.hybrid"
    _write_glm(a, _PVALS)
    perturbed = list(_PVALS)
    perturbed[3] = 0.31  # beyond default p_tolerance
    _write_glm(b, perturbed)

    result = compare_gwas_results(a, b)

    assert not result.equal
    assert "rs3" in result.mismatched_variants


def test_missing_variant_is_detected(tmp_path: Path) -> None:
    a = tmp_path / "a.PHENO1.glm.logistic.hybrid"
    b = tmp_path / "b.PHENO1.glm.logistic.hybrid"
    _write_glm(a, _PVALS)
    _write_glm(b, _PVALS[:-1])  # drop one variant

    result = compare_gwas_results(a, b)

    assert not result.equal
    assert result.variant_diff >= 1


def test_lambda_divergence_flagged_even_when_ids_match(tmp_path: Path) -> None:
    """Same variant IDs but a wildly different p-distribution → lambda mismatch."""
    a = tmp_path / "a.PHENO1.glm.logistic.hybrid"
    b = tmp_path / "b.PHENO1.glm.logistic.hybrid"
    _write_glm(a, [0.5] * 20)  # lambda ~ 0
    _write_glm(b, [1e-8] * 20)  # hugely inflated
    result = compare_gwas_results(a, b, p_tolerance=1.0)  # ignore per-variant p
    assert not result.equal, result.message


def test_find_gwas_output_locates_glm(tmp_path: Path) -> None:
    prefix = tmp_path / "run"
    glm = tmp_path / "run.PHENO1.glm.logistic.hybrid"
    _write_glm(glm, _PVALS)

    found = find_gwas_output(prefix)
    assert found == glm


def test_compare_gwas_reports_missing_output(tmp_path: Path) -> None:
    a_prefix = tmp_path / "a"
    b_prefix = tmp_path / "b"
    _write_glm(tmp_path / "a.PHENO1.glm.logistic.hybrid", _PVALS)
    # b has no glm output

    result = compare_gwas(a_prefix, b_prefix)
    assert not result.equal
    assert "missing" in result.message.lower()
