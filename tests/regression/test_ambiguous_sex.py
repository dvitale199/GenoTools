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

"""QC steps must work on cohorts holding ambiguous-sex samples.

Real cohorts contain samples with SEX=0 and a recorded phenotype. PLINK 1.9
refuses to combine a genotype write with any other command for such data:

    Error: When ambiguous-sex samples with phenotype data are present,
    --make-bed/--make-just-fam/--recode/--write-covar usually cannot be
    combined with other commands.

Every QC step that shells out to PLINK 1.9 - sex, case_control, haplotype, hwe
- tripped over this once run_plink() started appending --make-bed, and with
--warn the pipeline reported the steps as failed-but-continuing and produced
output that had silently skipped them. On a 10k GP2 r12 subset that dropped 4
QC steps for 3 of 11 ancestry groups: 2 ambiguous-sex samples out of 6,802 were
enough to skip 88 sex outliers and ~98k variants in EUR alone.

The synthetic fixture has no ambiguous-sex samples at all (every sample is
SEX=1 or 2 with a phenotype), so nothing in the suite exercised this. These
tests derive an ambiguous-sex fixture at run time - only the .psam is rewritten,
so no golden file is affected.
"""

from pathlib import Path

import pandas as pd
import pytest


@pytest.fixture()
def ambiguous_sex_geno(test_geno_path: Path, tmp_path: Path) -> Path:
    """Copy the synthetic pfile, setting SEX=0 on samples that have phenotypes.

    Sex lives in the .psam only, so the .pgen/.pvar are copied untouched.
    """
    import shutil

    prefix = tmp_path / "ambiguous_sex"
    for ext in (".pgen", ".pvar"):
        shutil.copy2(test_geno_path.with_suffix(ext), prefix.with_suffix(ext))

    psam = pd.read_csv(test_geno_path.with_suffix(".psam"), sep="\t", dtype=str)
    sex_col = "SEX"
    pheno_col = "PHENO1"
    assert sex_col in psam.columns and pheno_col in psam.columns, (
        f"unexpected psam columns: {list(psam.columns)}"
    )
    # Pick samples that carry a real phenotype - that is the combination PLINK
    # objects to. An ambiguous-sex sample with a missing phenotype does not
    # trigger it.
    with_pheno = psam.index[psam[pheno_col].isin(["1", "2"])]
    assert len(with_pheno) >= 3, "fixture needs phenotyped samples"
    psam.loc[with_pheno[:3], sex_col] = "0"
    psam.to_csv(prefix.with_suffix(".psam"), sep="\t", index=False)
    return prefix


def _assert_ambiguous(prefix: Path) -> None:
    """Guard the fixture itself: the risky combination must really be present."""
    psam = pd.read_csv(prefix.with_suffix(".psam"), sep="\t", dtype=str)
    risky = psam[(psam["SEX"] == "0") & (psam["PHENO1"].isin(["1", "2"]))]
    assert len(risky) == 3, (
        "fixture no longer contains ambiguous-sex samples with phenotype data, "
        "so these tests would pass vacuously"
    )


class TestAmbiguousSexQCSteps:
    """Each PLINK 1.9-backed step must complete, not fail-and-skip."""

    def test_fixture_has_ambiguous_sex_with_phenotypes(
        self, ambiguous_sex_geno: Path
    ) -> None:
        _assert_ambiguous(ambiguous_sex_geno)

    def test_sex_prune_runs(self, ambiguous_sex_geno: Path, tmp_path: Path) -> None:
        from genotools.core.genotypes import GenotypeData
        from genotools.qc.config import SexConfig
        from genotools.qc.steps.sex import filter_sex

        _assert_ambiguous(ambiguous_sex_geno)
        data = GenotypeData.from_path(ambiguous_sex_geno)

        result = filter_sex(data, SexConfig(), tmp_path / "sex_out")

        assert result.output.path.with_suffix(".pgen").exists()

    def test_case_control_prune_runs(
        self, ambiguous_sex_geno: Path, tmp_path: Path
    ) -> None:
        from genotools.core.genotypes import GenotypeData
        from genotools.qc.config import CaseControlConfig
        from genotools.qc.steps.case_control import filter_case_control

        _assert_ambiguous(ambiguous_sex_geno)
        data = GenotypeData.from_path(ambiguous_sex_geno)

        result = filter_case_control(
            data, CaseControlConfig(), tmp_path / "cc_out"
        )

        assert result.output.path.with_suffix(".pgen").exists()

    def test_haplotype_prune_runs(
        self, ambiguous_sex_geno: Path, tmp_path: Path
    ) -> None:
        from genotools.core.genotypes import GenotypeData
        from genotools.qc.config import HaplotypeConfig
        from genotools.qc.steps.haplotype import filter_haplotype

        _assert_ambiguous(ambiguous_sex_geno)
        data = GenotypeData.from_path(ambiguous_sex_geno)

        result = filter_haplotype(
            data, HaplotypeConfig(), tmp_path / "hap_out"
        )

        assert result.output.path.with_suffix(".pgen").exists()

    def test_hwe_prune_runs(self, ambiguous_sex_geno: Path, tmp_path: Path) -> None:
        from genotools.core.genotypes import GenotypeData
        from genotools.qc.config import HWEConfig
        from genotools.qc.steps.hwe import filter_hwe

        _assert_ambiguous(ambiguous_sex_geno)
        data = GenotypeData.from_path(ambiguous_sex_geno)

        result = filter_hwe(data, HWEConfig(), tmp_path / "hwe_out")

        assert result.output.path.with_suffix(".pgen").exists()
