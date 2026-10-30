"""The haplotype step reports the p-value threshold it actually applied."""

import dataclasses
import shutil
from pathlib import Path

import pytest

from genotools.core.genotypes import GenotypeData
from genotools.qc.config import HaplotypeConfig
from genotools.qc.steps.haplotype import effective_p_threshold, filter_haplotype

SYNTHETIC = Path(__file__).resolve().parents[2] / "data" / "synthetic" / "genotools_test"
_EXE = Path.home() / ".genotools" / "misc" / "executables"
needs_plink = pytest.mark.skipif(
    not ((_EXE / "plink").exists() or shutil.which("plink")),
    reason="PLINK 1.9 not available",
)


class TestEffectiveThreshold:
    def test_small_cohort_keeps_the_configured_value(self) -> None:
        assert effective_p_threshold(HaplotypeConfig(p_threshold=1e-4), 10000) == 1e-4

    def test_large_cohort_is_scaled_by_sample_count(self) -> None:
        assert effective_p_threshold(HaplotypeConfig(p_threshold=1e-4), 20000) == pytest.approx(2.5e-6)

    def test_adjustment_can_be_turned_off(self) -> None:
        config = HaplotypeConfig(p_threshold=1e-4, adjust_for_sample_size=False)
        assert effective_p_threshold(config, 50000) == 1e-4


@needs_plink
class TestStepRecordsAppliedThreshold:
    def test_unadjusted_value_is_recorded(self, tmp_path: Path) -> None:
        data = GenotypeData.from_path(SYNTHETIC)
        result = filter_haplotype(data, HaplotypeConfig(p_threshold=1e-4), tmp_path / "out")
        assert result.parameters == {"p_threshold": 1e-4}
        assert "p_threshold" not in result.metrics  # counts only

    def test_adjusted_value_is_recorded(self, tmp_path: Path) -> None:
        data = dataclasses.replace(GenotypeData.from_path(SYNTHETIC), sample_count=20000)
        result = filter_haplotype(data, HaplotypeConfig(p_threshold=1e-4), tmp_path / "out")
        assert result.parameters == {"p_threshold": pytest.approx(0.05 / 20000)}
        assert result.to_dict()["parameters"]["p_threshold"] == pytest.approx(2.5e-6)
