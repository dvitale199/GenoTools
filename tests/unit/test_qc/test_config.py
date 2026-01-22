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

"""Tests for QC configuration classes."""

import pytest

from genotools.qc.config import (
    CallrateConfig,
    CaseControlConfig,
    GenoConfig,
    HaplotypeConfig,
    HetConfig,
    HWEConfig,
    LDConfig,
    RelatedConfig,
    SexConfig,
)


class TestCallrateConfig:
    """Tests for CallrateConfig."""

    def test_default_values(self) -> None:
        """Default values are set correctly."""
        config = CallrateConfig()
        assert config.mind == 0.05

    def test_custom_value(self) -> None:
        """Custom values are accepted."""
        config = CallrateConfig(mind=0.02)
        assert config.mind == 0.02

    def test_invalid_mind_negative(self) -> None:
        """Negative mind raises ValueError."""
        with pytest.raises(ValueError, match="mind"):
            CallrateConfig(mind=-0.1)

    def test_invalid_mind_too_large(self) -> None:
        """Mind > 1 raises ValueError."""
        with pytest.raises(ValueError, match="mind"):
            CallrateConfig(mind=1.5)

    def test_frozen(self) -> None:
        """Config is immutable."""
        config = CallrateConfig()
        with pytest.raises(AttributeError):
            config.mind = 0.1  # type: ignore[misc]


class TestSexConfig:
    """Tests for SexConfig."""

    def test_default_values(self) -> None:
        """Default values are set correctly."""
        config = SexConfig()
        assert config.female_max_f == 0.25
        assert config.male_min_f == 0.75

    def test_custom_values(self) -> None:
        """Custom values are accepted."""
        config = SexConfig(female_max_f=0.2, male_min_f=0.8)
        assert config.female_max_f == 0.2
        assert config.male_min_f == 0.8

    def test_invalid_female_max_f_negative(self) -> None:
        """Negative female_max_f raises ValueError."""
        with pytest.raises(ValueError, match="female_max_f"):
            SexConfig(female_max_f=-0.1)

    def test_invalid_male_min_f_too_large(self) -> None:
        """male_min_f > 1 raises ValueError."""
        with pytest.raises(ValueError, match="male_min_f"):
            SexConfig(male_min_f=1.5)

    def test_invalid_thresholds_overlap(self) -> None:
        """female_max_f >= male_min_f raises ValueError."""
        with pytest.raises(ValueError, match="female_max_f.*must be less than"):
            SexConfig(female_max_f=0.6, male_min_f=0.5)

    def test_invalid_thresholds_equal(self) -> None:
        """female_max_f == male_min_f raises ValueError."""
        with pytest.raises(ValueError, match="female_max_f.*must be less than"):
            SexConfig(female_max_f=0.5, male_min_f=0.5)


class TestHetConfig:
    """Tests for HetConfig."""

    def test_default_values(self) -> None:
        """Default values are set correctly."""
        config = HetConfig()
        assert config.f_lower == -0.15
        assert config.f_upper == 0.15
        assert config.auto_detect is False

    def test_auto_detect_mode(self) -> None:
        """Auto-detect mode can be enabled."""
        config = HetConfig(auto_detect=True)
        assert config.auto_detect is True

    def test_auto_detect_via_sentinel(self) -> None:
        """Setting both bounds to -1 enables auto-detect mode."""
        config = HetConfig(f_lower=-1.0, f_upper=-1.0)
        # Should not raise validation error
        assert config.f_lower == -1.0
        assert config.f_upper == -1.0

    def test_invalid_bounds_order(self) -> None:
        """f_lower >= f_upper raises ValueError."""
        with pytest.raises(ValueError, match="f_lower.*must be less than"):
            HetConfig(f_lower=0.2, f_upper=0.1)

    def test_invalid_bounds_equal(self) -> None:
        """f_lower == f_upper raises ValueError."""
        with pytest.raises(ValueError, match="f_lower.*must be less than"):
            HetConfig(f_lower=0.15, f_upper=0.15)


class TestRelatedConfig:
    """Tests for RelatedConfig."""

    def test_default_values(self) -> None:
        """Default values are set correctly."""
        config = RelatedConfig()
        assert config.related_cutoff == 0.0884
        assert config.duplicated_cutoff == 0.354
        assert config.prune_related is True
        assert config.prune_duplicated is True

    def test_custom_values(self) -> None:
        """Custom values are accepted."""
        config = RelatedConfig(
            related_cutoff=0.1,
            duplicated_cutoff=0.4,
            prune_related=False,
            prune_duplicated=False,
        )
        assert config.related_cutoff == 0.1
        assert config.duplicated_cutoff == 0.4
        assert config.prune_related is False
        assert config.prune_duplicated is False

    def test_invalid_related_without_duplicated(self) -> None:
        """Cannot prune related without also pruning duplicated."""
        with pytest.raises(ValueError, match="Cannot prune related.*without"):
            RelatedConfig(prune_related=True, prune_duplicated=False)

    def test_invalid_cutoff_negative(self) -> None:
        """Negative cutoff raises ValueError."""
        with pytest.raises(ValueError, match="related_cutoff"):
            RelatedConfig(related_cutoff=-0.1)

    def test_invalid_cutoff_too_large(self) -> None:
        """Cutoff > 0.5 raises ValueError."""
        with pytest.raises(ValueError, match="duplicated_cutoff"):
            RelatedConfig(duplicated_cutoff=0.6)


class TestGenoConfig:
    """Tests for GenoConfig."""

    def test_default_values(self) -> None:
        """Default values are set correctly."""
        config = GenoConfig()
        assert config.geno == 0.05

    def test_custom_value(self) -> None:
        """Custom values are accepted."""
        config = GenoConfig(geno=0.02)
        assert config.geno == 0.02

    def test_invalid_geno_negative(self) -> None:
        """Negative geno raises ValueError."""
        with pytest.raises(ValueError, match="geno"):
            GenoConfig(geno=-0.1)


class TestCaseControlConfig:
    """Tests for CaseControlConfig."""

    def test_default_values(self) -> None:
        """Default values are set correctly."""
        config = CaseControlConfig()
        assert config.p_threshold == 1e-4

    def test_custom_value(self) -> None:
        """Custom values are accepted."""
        config = CaseControlConfig(p_threshold=1e-6)
        assert config.p_threshold == 1e-6

    def test_invalid_threshold_negative(self) -> None:
        """Negative threshold raises ValueError."""
        with pytest.raises(ValueError, match="p_threshold"):
            CaseControlConfig(p_threshold=-0.1)


class TestHaplotypeConfig:
    """Tests for HaplotypeConfig."""

    def test_default_values(self) -> None:
        """Default values are set correctly."""
        config = HaplotypeConfig()
        assert config.p_threshold == 1e-4
        assert config.adjust_for_sample_size is True

    def test_disable_sample_size_adjustment(self) -> None:
        """Sample size adjustment can be disabled."""
        config = HaplotypeConfig(adjust_for_sample_size=False)
        assert config.adjust_for_sample_size is False


class TestHWEConfig:
    """Tests for HWEConfig."""

    def test_default_values(self) -> None:
        """Default values are set correctly."""
        config = HWEConfig()
        assert config.hwe_threshold == 1e-4
        assert config.filter_controls is False

    def test_filter_controls_enabled(self) -> None:
        """Filter controls can be enabled."""
        config = HWEConfig(filter_controls=True)
        assert config.filter_controls is True


class TestLDConfig:
    """Tests for LDConfig."""

    def test_default_values(self) -> None:
        """Default values are set correctly."""
        config = LDConfig()
        assert config.window_size == 50
        assert config.step_size == 5
        assert config.r2_threshold == 0.5

    def test_custom_values(self) -> None:
        """Custom values are accepted."""
        config = LDConfig(window_size=100, step_size=10, r2_threshold=0.2)
        assert config.window_size == 100
        assert config.step_size == 10
        assert config.r2_threshold == 0.2

    def test_invalid_r2_too_large(self) -> None:
        """r2 > 1 raises ValueError."""
        with pytest.raises(ValueError, match="r2_threshold"):
            LDConfig(r2_threshold=1.5)

    def test_invalid_window_size_negative(self) -> None:
        """Negative window_size raises ValueError."""
        with pytest.raises(ValueError, match="window_size"):
            LDConfig(window_size=-10)


class TestConfigImmutability:
    """Tests that all configs are immutable (frozen)."""

    def test_all_configs_frozen(self) -> None:
        """All config classes are frozen dataclasses."""
        configs = [
            CallrateConfig(),
            SexConfig(),
            HetConfig(),
            RelatedConfig(),
            GenoConfig(),
            CaseControlConfig(),
            HaplotypeConfig(),
            HWEConfig(),
            LDConfig(),
        ]

        for config in configs:
            # Try to modify an attribute - should raise
            attr = list(config.__dataclass_fields__.keys())[0]
            with pytest.raises(AttributeError):
                setattr(config, attr, "new_value")
