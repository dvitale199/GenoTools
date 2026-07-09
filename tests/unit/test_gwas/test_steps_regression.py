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

"""Regression tests for GWAS step bugs found in the refactor hardening audit."""

from pathlib import Path

import pytest

from genotools.core.genotypes import GenotypeData
from genotools.gwas.config import GWASConfig
from genotools.gwas.steps.association import _process_logistic_results, run_gwas

SYNTHETIC = (
    Path(__file__).resolve().parents[2] / "data" / "synthetic" / "genotools_test"
)


def _plink2_available() -> bool:
    import shutil

    if shutil.which("plink2"):
        return True
    return (Path.home() / ".genotools" / "misc" / "executables" / "plink2").exists()


def _write_binary_cohort(tmp_path: Path, n: int = 10) -> GenotypeData:
    """Write a minimal .psam with both cases (PHENO1==2) and controls (PHENO1==1).

    Only the .psam is needed: _process_logistic_results reads case/control
    counts from data.path.psam and never opens the .pgen/.pvar.
    """
    prefix = tmp_path / "cohort"
    rows = ["#FID\tIID\tSEX\tPHENO1"]
    for i in range(n):
        pheno = "1" if i % 2 == 0 else "2"  # alternate control/case
        rows.append(f"F{i}\tI{i}\t1\t{pheno}")
    prefix.with_suffix(".psam").write_text("\n".join(rows) + "\n")
    return GenotypeData(path=prefix, format="pfile", sample_count=n, variant_count=50)


def _write_glm_logistic(tmp_path: Path, n_variants: int = 50) -> Path:
    """Write a minimal PLINK2 .glm.logistic.hybrid file with ADD-test rows."""
    glm = tmp_path / "cohort.PHENO1.glm.logistic.hybrid"
    rows = ["#CHROM\tPOS\tID\tTEST\tP"]
    for i in range(n_variants):
        rows.append(f"1\t{i + 1}\trs{i}\tADD\t0.5")
    glm.write_text("\n".join(rows) + "\n")
    return glm


class TestLogisticResultSummaryFormatting:
    """Regression: association.py summary log used an invalid f-string format spec.

    `f"lambda1000={inflation.lambda_1000:.4f if inflation.lambda_1000 else 'N/A'}"`
    treats everything after the ':' as a format spec, raising ValueError (float) or
    TypeError (None) on EVERY binary-phenotype GWAS, right before returning the result.
    """

    def test_binary_gwas_summary_does_not_crash(self, tmp_path: Path) -> None:
        """_process_logistic_results completes and returns a successful result."""
        data = _write_binary_cohort(tmp_path)
        glm = _write_glm_logistic(tmp_path)

        result = _process_logistic_results(
            output_path=glm,
            data=data,
            config=GWASConfig(),
            log="",
        )

        assert result.success is True
        assert result.phenotype_type == "binary"
        assert result.inflation is not None
        assert result.inflation.lambda_1000 is not None  # cases+controls present


class TestGwasGlmOptionsAreTokenized:
    """Regression: run_gwas passed the multi-modifier --glm options as a SINGLE
    argv token. shell_do (old code) split the command on whitespace; run_command
    (new code) treats each list element as one argument, so PLINK2 received
    'hide-covar firth-fallback no-x-sex cols=...' as one token and rejected it
    ('Invalid --glm argument ...'), producing no GWAS output. This broke every
    GWAS run in the new CLI, silently under --warn.
    """

    @pytest.mark.skipif(not _plink2_available(), reason="plink2 not available")
    def test_run_gwas_produces_output(self, tmp_path: Path) -> None:
        """A real PLINK2 --glm run must succeed and write a results file."""
        if not SYNTHETIC.with_suffix(".pgen").exists():
            pytest.skip("synthetic test data not found")

        data = GenotypeData.from_path(SYNTHETIC)
        out = tmp_path / "gwas_out"

        # No covariates -> exercises the allow-no-covars branch of the command.
        result = run_gwas(data, GWASConfig(), out)

        assert result.success is True
        assert result.output_path is not None and result.output_path.exists()
        assert result.n_variants_tested > 0
