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

"""Tests for heterozygosity outlier selection.

``select_het_outliers`` is split out of ``filter_heterozygosity`` so the bound
arithmetic can be exercised against a frame with a known mean and sigma,
without PLINK in the loop. These pin the ``sd`` multiplier to real behaviour
rather than to config plumbing.
"""

from typing import List

import pandas as pd
import pytest

from genotools.qc.config import HetConfig
from genotools.qc.steps.heterozygosity import select_het_outliers


def _het_frame(f_values: List[float], obs_hom: List[int] = None) -> pd.DataFrame:
    """Build a frame shaped like PLINK's .het output.

    F and the derived rate are independent here, which is the point: the two
    modes threshold different statistics, so a frame can separate them.
    """
    n = len(f_values)
    hom = obs_hom if obs_hom is not None else [5000] * n
    return pd.DataFrame(
        {
            "#FID": [f"F{i}" for i in range(n)],
            "IID": [f"I{i}" for i in range(n)],
            "O(HOM)": hom,
            "E(HOM)": [5000.0] * n,
            "OBS_CT": [10000] * n,
            "F": f_values,
        }
    )


class TestFixedBounds:
    """The default path, unchanged."""

    def test_bounds_are_inclusive(self) -> None:
        het = _het_frame([-0.2, -0.15, 0.0, 0.15, 0.2])
        outliers, metrics = select_het_outliers(het, HetConfig())
        assert sorted(outliers["F"]) == [-0.2, -0.15, 0.15, 0.2]
        assert metrics["het_mode"] == "fixed"
        assert metrics["het_statistic"] == "F"

    def test_custom_bounds(self) -> None:
        het = _het_frame([-0.25, 0.0, 0.25])
        outliers, _ = select_het_outliers(
            het, HetConfig(f_lower=-0.3, f_upper=0.3)
        )
        assert outliers.empty


class TestAutoDetectThresholdsF:
    """``auto_detect=True`` bounds F, so both CLI spellings of --het bound the
    same statistic. Selected without the legacy [-1, -1] sentinel."""

    def test_auto_detect_without_the_sentinel(self) -> None:
        het = _het_frame([0.0] * 20 + [1.0])
        outliers, metrics = select_het_outliers(het, HetConfig(auto_detect=True))
        assert metrics["het_mode"] == "sd"
        assert metrics["het_statistic"] == "F"
        assert metrics["het_sd"] == 3.0
        assert not outliers.empty

    def test_auto_detect_ignores_the_rate(self) -> None:
        """F constant, rate wildly variable: the F branch must prune nothing.

        If auto_detect ever fell through to the rate branch this would flag the
        rate outlier, silently swapping which statistic is bounded.
        """
        het = _het_frame([0.01] * 20, obs_hom=[5000] * 19 + [100])
        outliers, metrics = select_het_outliers(het, HetConfig(auto_detect=True))
        assert outliers.empty
        assert metrics["het_statistic"] == "F"

    def test_fixed_bounds_are_ignored_in_auto_mode(self) -> None:
        """f_lower/f_upper ride along unused rather than being consulted."""
        het = _het_frame([0.0] * 20 + [1.0])
        wide, _ = select_het_outliers(
            het, HetConfig(f_lower=-10.0, f_upper=10.0, auto_detect=True)
        )
        narrow, _ = select_het_outliers(
            het, HetConfig(f_lower=-0.001, f_upper=0.001, auto_detect=True)
        )
        assert len(wide) == len(narrow)


class TestAutoSdMovesTheCutoff:
    """The multiplier has to move real behaviour, not just ride in the config."""

    @staticmethod
    def _frame_with_probe_between_2_and_3_sd():
        """A frame whose last sample sits between 2 and 3 sd of the whole set.

        Built by construction rather than assumed: mean and sd are computed
        over every sample, the probe included, so the probe widens the very
        bounds that judge it.
        """
        values = [-0.01, 0.01] * 50 + [0.025]
        het = _het_frame(values)
        mean = het["F"].mean()
        std = het["F"].std()
        z = (0.025 - mean) / std
        assert 2.0 < z < 3.0, f"probe must sit between 2 and 3 sd, got {z}"
        return het

    def test_probe_kept_at_sd_3_pruned_at_sd_2(self) -> None:
        het = self._frame_with_probe_between_2_and_3_sd()

        kept, _ = select_het_outliers(
            het, HetConfig(auto_detect=True, auto_sd=3.0)
        )
        pruned, _ = select_het_outliers(
            het, HetConfig(auto_detect=True, auto_sd=2.0)
        )

        assert kept.empty, "a 2.5 sd sample survives a 3 sd cut"
        assert 0.025 in list(pruned["F"]), "and is pruned by a 2 sd cut"

    def test_tightening_the_multiplier_can_only_prune_more(self) -> None:
        het = _het_frame([-0.03, -0.01, 0.0, 0.01, 0.02, 0.03] * 10 + [0.09])
        counts = [
            len(select_het_outliers(het, HetConfig(auto_detect=True, auto_sd=sd))[0])
            for sd in (1.0, 1.5, 2.0, 3.0, 4.0)
        ]
        assert counts == sorted(counts, reverse=True), counts

    def test_reported_bounds_match_the_multiplier(self) -> None:
        het = _het_frame([-0.01, 0.01] * 50)
        _, metrics = select_het_outliers(
            het, HetConfig(auto_detect=True, auto_sd=2.0)
        )
        mean, std = het["F"].mean(), het["F"].std()
        assert metrics["het_lower"] == pytest.approx(mean - 2.0 * std)
        assert metrics["het_upper"] == pytest.approx(mean + 2.0 * std)


class TestLegacySentinel:
    """``f_lower == f_upper == -1.0`` is the 1.x way to ask for derived bounds.

    It thresholds the derived heterozygosity *rate* at 3 sd. Preserved verbatim
    for Python-API callers that already pass it; the CLI no longer produces it.
    """

    def test_sentinel_selects_the_rate_branch(self) -> None:
        het = _het_frame([0.01] * 20, obs_hom=[5000] * 19 + [100])
        outliers, metrics = select_het_outliers(
            het, HetConfig(f_lower=-1.0, f_upper=-1.0)
        )
        assert metrics["het_statistic"] == "rate"
        assert metrics["het_sd"] == 3.0
        assert len(outliers) == 1

    def test_sentinel_keeps_the_HET_column(self) -> None:
        """Existing callers read it out of the .outliers file."""
        het = _het_frame([0.01] * 20, obs_hom=[5000] * 19 + [100])
        outliers, _ = select_het_outliers(
            het, HetConfig(f_lower=-1.0, f_upper=-1.0)
        )
        assert "HET" in outliers.columns

    def test_sentinel_ignores_auto_sd(self) -> None:
        """The legacy path is fixed at 3 sd; the knob is a CLI-era addition."""
        het = _het_frame([0.01] * 20, obs_hom=[5000] * 19 + [100])
        _, metrics = select_het_outliers(
            het, HetConfig(f_lower=-1.0, f_upper=-1.0, auto_sd=1.0)
        )
        assert metrics["het_sd"] == 3.0

    def test_auto_detect_wins_over_the_sentinel(self) -> None:
        """Precedence has to be explicit, or a config setting both would pick
        its statistic by accident of branch order."""
        het = _het_frame([0.01] * 20, obs_hom=[5000] * 19 + [100])
        _, metrics = select_het_outliers(
            het, HetConfig(f_lower=-1.0, f_upper=-1.0, auto_detect=True)
        )
        assert metrics["het_statistic"] == "F"


class TestAutoSdValidation:
    """HetConfig rejects a non-positive multiplier at construction."""

    @pytest.mark.parametrize("bad", [0.0, -1.0])
    def test_non_positive_rejected(self, bad: float) -> None:
        with pytest.raises(ValueError, match="auto_sd"):
            HetConfig(auto_detect=True, auto_sd=bad)

    def test_default_is_three(self) -> None:
        assert HetConfig().auto_sd == 3.0

    def test_validated_even_under_the_sentinel(self) -> None:
        with pytest.raises(ValueError, match="auto_sd"):
            HetConfig(f_lower=-1.0, f_upper=-1.0, auto_sd=0.0)
