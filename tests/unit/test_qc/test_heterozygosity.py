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

from pathlib import Path
from typing import Any, List

import pandas as pd
import pytest

from genotools.qc import steps
from genotools.qc.config import HetConfig
from genotools.qc.steps.heterozygosity import (
    filter_heterozygosity,
    select_het_outliers,
)


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


class _FakeGenotypeData:
    """Enough of GenotypeData for filter_heterozygosity."""

    format = "pfile"

    def __init__(self, path: Path, sample_count: int) -> None:
        self.path = path
        self.sample_count = sample_count

    @classmethod
    def from_path(cls, path: Path) -> "_FakeGenotypeData":
        return cls(Path(path), 92)


class TestDerivedBoundsAreReported:
    """The derived bounds have to leave the step, not just the log.

    ``select_het_outliers`` works them out and hands them back; before this the
    step logged them and dropped them, so a report could say 8 samples were
    pruned and nothing about the rule that pruned them. With per-group bounds
    that makes two rows of the same report incomparable.
    """

    @staticmethod
    def _stub_plink(monkeypatch: pytest.MonkeyPatch, tmp_path: Path, f_values: List[float]) -> None:
        """Stand in for the three PLINK calls, writing what each produces."""
        het_mod = steps.heterozygosity

        def fake_run_plink2(
            input_path: Any, output_path: Any, extra_args: List[str] = None, **kwargs: Any
        ) -> None:
            extra_args = extra_args or []
            out = Path(output_path)
            if "--indep-pairwise" in extra_args:
                out.with_suffix(".prune.in").write_text("rs1\nrs2\n")
            elif "--het" in extra_args:
                _het_frame(f_values).to_csv(
                    out.with_suffix(".het"), sep="\t", index=False
                )
            else:
                for ext in (".pgen", ".pvar", ".psam"):
                    out.with_suffix(ext).write_text("")

        monkeypatch.setattr(het_mod, "run_plink2", fake_run_plink2)
        monkeypatch.setattr(het_mod, "GenotypeData", _FakeGenotypeData)

    def _run(self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch, config: HetConfig):
        geno = tmp_path / "geno"
        for ext in (".pgen", ".pvar", ".psam"):
            geno.with_suffix(ext).write_text("")
        # 100 samples: two at +/-1.0 sit far outside any of these bounds.
        f_values = [0.0] * 98 + [-1.0, 1.0]
        self._stub_plink(monkeypatch, tmp_path, f_values)
        data = _FakeGenotypeData(geno, 100)
        return filter_heterozygosity(data, config, tmp_path / "out")

    def test_sd_mode_reports_the_bounds_it_derived(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """Revert check: drop ``parameters=het_metrics`` from the FilterResult
        and this fails with an empty dict."""
        result = self._run(
            tmp_path, monkeypatch, HetConfig(auto_detect=True, auto_sd=2.0)
        )

        assert result.parameters["het_mode"] == "sd"
        assert result.parameters["het_statistic"] == "F"
        assert result.parameters["het_sd"] == 2.0
        # Derived from the frame, so not predictable from the config alone -
        # which is exactly why the step has to report them.
        assert result.parameters["het_lower"] < 0 < result.parameters["het_upper"]

    def test_fixed_mode_reports_the_bounds_it_was_given(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """Fixed bounds are reported too: without them a mixed run cannot say
        which groups were cut by which rule."""
        result = self._run(
            tmp_path, monkeypatch, HetConfig(f_lower=-0.2, f_upper=0.2)
        )

        assert result.parameters["het_mode"] == "fixed"
        assert result.parameters["het_sd"] is None
        assert (result.parameters["het_lower"], result.parameters["het_upper"]) == (
            -0.2,
            0.2,
        )

    def test_metrics_stay_counts_only(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """The reason parameters is a separate channel: every metric becomes a
        row under a column meaning "how many were pruned"."""
        result = self._run(
            tmp_path, monkeypatch, HetConfig(auto_detect=True, auto_sd=2.0)
        )

        assert set(result.metrics) == {"outlier_count"}
        assert isinstance(result.metrics["outlier_count"], int)
