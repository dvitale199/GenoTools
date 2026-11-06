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

"""#323: a step decides "this removed everything" from its own count.

PLINK2's "No samples remaining" message is matched in ``run_plink2`` too, but
only as a backstop. These pin the primary signal: a failed filtering call is
classified by what the step had already worked out it was removing, so a
reworded PLINK message cannot quietly bring back the bug where later steps QC
the unfiltered input.
"""

from pathlib import Path
from typing import List

import pandas as pd
import pytest

from genotools.core.exceptions import ExternalToolError, NoDataRemainingError
from genotools.core.executors import emptying_filter
from genotools.core.genotypes import GenotypeData
from genotools.qc.config import CallrateConfig, GenoConfig
from genotools.qc.steps import callrate as callrate_mod
from genotools.qc.steps import variant_missingness as geno_mod
from genotools.qc.steps.callrate import filter_callrate
from genotools.qc.steps.variant_missingness import filter_variant_missingness

# Deliberately not PLINK2's wording, so only the step's own count can match.
REWORDED = "Error: The filters left an empty dataset."


def _tool_error() -> ExternalToolError:
    return ExternalToolError("plink2", "plink2 --mind 0", 13, REWORDED)


class TestEmptyingFilter:
    def test_failure_is_classified_by_the_count(self, tmp_path: Path) -> None:
        out = tmp_path / "cohort.v2_callrate"
        with pytest.raises(NoDataRemainingError) as exc:
            with emptying_filter(True, "samples", out):
                raise _tool_error()

        assert exc.value.remaining == "samples"
        assert exc.value.output_path == str(out)
        assert exc.value.returncode == 13
        assert REWORDED in exc.value.stderr

    def test_failure_with_data_left_stays_plain(self, tmp_path: Path) -> None:
        with pytest.raises(ExternalToolError) as exc:
            with emptying_filter(False, "samples", tmp_path / "out"):
                raise _tool_error()

        assert not isinstance(exc.value, NoDataRemainingError)

    def test_already_classified_error_passes_through(self, tmp_path: Path) -> None:
        original = NoDataRemainingError(
            "plink2", "cmd", 13, "No variants remaining", "variants", str(tmp_path / "a")
        )
        with pytest.raises(NoDataRemainingError) as exc:
            with emptying_filter(False, "samples", tmp_path / "b"):
                raise original

        assert exc.value is original

    def test_success_is_untouched(self, tmp_path: Path) -> None:
        with emptying_filter(True, "samples", tmp_path / "out"):
            pass


def _input(tmp_path: Path, samples: int, variants: int) -> GenotypeData:
    prefix = tmp_path / "cohort.v2"
    Path(f"{prefix}.pgen").touch()
    return GenotypeData(
        path=prefix, format="pfile", sample_count=samples, variant_count=variants
    )


def _stub_plink2(monkeypatch, module, report_ext: str, report: pd.DataFrame) -> None:
    """The --missing call writes `report`; the filtering call fails, reworded."""

    def fake_run_plink2(*, output_path: Path, extra_args: List[str], **_) -> None:
        if "--missing" in extra_args:
            report.to_csv(f"{output_path}.{report_ext}", sep="\t", index=False)
            return
        raise _tool_error()

    monkeypatch.setattr(module, "run_plink2", fake_run_plink2)


def _smiss(f_miss: List[float]) -> pd.DataFrame:
    n = len(f_miss)
    return pd.DataFrame(
        {"#FID": ["0"] * n, "IID": [f"s{i}" for i in range(n)], "F_MISS": f_miss}
    )


def _vmiss(f_miss: List[float]) -> pd.DataFrame:
    return pd.DataFrame({"ID": [f"v{i}" for i in range(len(f_miss))], "F_MISS": f_miss})


class TestCallrateCountsItsOwnRemovals:
    def test_every_sample_over_mind_is_no_samples_remaining(
        self, monkeypatch, tmp_path: Path
    ) -> None:
        _stub_plink2(monkeypatch, callrate_mod, "smiss", _smiss([0.01, 0.2, 0.5]))
        out = tmp_path / "cohort.v2_callrate"

        with pytest.raises(NoDataRemainingError) as exc:
            filter_callrate(_input(tmp_path, 3, 10), CallrateConfig(mind=0.0), out)

        assert exc.value.remaining == "samples"
        assert exc.value.output_path == str(out)

    def test_a_sample_at_mind_survives_so_the_failure_is_ordinary(
        self, monkeypatch, tmp_path: Path
    ) -> None:
        # --mind removes rates *above* the threshold: 0.0 at mind 0 is kept.
        _stub_plink2(monkeypatch, callrate_mod, "smiss", _smiss([0.0, 0.2, 0.5]))

        with pytest.raises(ExternalToolError) as exc:
            filter_callrate(
                _input(tmp_path, 3, 10), CallrateConfig(mind=0.0), tmp_path / "out"
            )

        assert not isinstance(exc.value, NoDataRemainingError)


class TestGenoCountsItsOwnRemovals:
    def test_every_variant_over_geno_is_no_variants_remaining(
        self, monkeypatch, tmp_path: Path
    ) -> None:
        _stub_plink2(monkeypatch, geno_mod, "vmiss", _vmiss([0.1, 0.2, 0.3]))
        out = tmp_path / "cohort.v2_geno"

        with pytest.raises(NoDataRemainingError) as exc:
            filter_variant_missingness(
                _input(tmp_path, 10, 3), GenoConfig(geno=0.05), out
            )

        assert exc.value.remaining == "variants"
        assert exc.value.output_path == str(out)

    def test_a_variant_at_geno_survives_so_the_failure_is_ordinary(
        self, monkeypatch, tmp_path: Path
    ) -> None:
        _stub_plink2(monkeypatch, geno_mod, "vmiss", _vmiss([0.05, 0.2, 0.3]))

        with pytest.raises(ExternalToolError) as exc:
            filter_variant_missingness(
                _input(tmp_path, 10, 3), GenoConfig(geno=0.05), tmp_path / "out"
            )

        assert not isinstance(exc.value, NoDataRemainingError)
