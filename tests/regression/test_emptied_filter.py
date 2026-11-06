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

"""#323: a filter that removes every sample ends the run, through the real CLI.

PLINK2 refuses to write an empty fileset ("No samples remaining after main
filters", exit 13). Before the fix that was an ordinary step failure, so under
the default warn-and-continue the next steps QC'd the unfiltered input and
reported PASS. Unit tests cover the runner with the error faked; this covers
PLINK2's real wording reaching it, and everything written on the way out.

Two signals detect the emptying, and each has a test here: the step's own count
of what it removes (primary), and PLINK2's message (backstop). Each is checked
with the other out of the way, so a PLINK upgrade that rewords the message fails
the backstop test by name instead of silently weakening the fix.

Requires the synthetic test data and real PLINK; skips cleanly otherwise.
"""

import json
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

from genotools.core import executors
from genotools.core.exceptions import NoDataRemainingError
from genotools.core.genotypes import GenotypeData
from genotools.qc.config import CallrateConfig
from genotools.qc.steps.callrate import filter_callrate

REPO_ROOT = Path(__file__).resolve().parents[2]


def _run_cli(args: list[str]) -> subprocess.CompletedProcess:
    return subprocess.run(
        [sys.executable, "-m", "genotools", *args],
        cwd=str(REPO_ROOT),
        capture_output=True,
        text=True,
        timeout=300,
    )


@pytest.mark.parametrize("no_warn", [False, True])
def test_mind_removing_every_sample_stops_the_run(
    test_geno_path: Path, tmp_path: Path, no_warn: bool
) -> None:
    out = tmp_path / f"emptied{int(no_warn)}"
    args = [
        "--pfile", str(test_geno_path), "--out", str(out),
        "--callrate", "0", "--geno", "0.05", "--hwe", "1e-4", "--plots",
    ]
    if no_warn:
        args.append("--no-warn")

    res = _run_cli(args)

    assert res.returncode == 1, f"{res.stdout}\n{res.stderr}"
    assert "ERROR: No data remain after QC: callrate removed all 500 samples" in res.stderr
    assert "Traceback" not in res.stderr, res.stderr

    # The report is written before the run fails, and says what happened.
    report = json.loads(Path(f"{out}.json").read_text())
    qc = pd.DataFrame(report["QC"])
    rows = qc.set_index("step")
    assert rows.loc["callrate_prune", "outcome"] == "pass"
    assert rows.loc["callrate_prune", "pruned_count"] == 500
    for step in ("geno_prune", "hwe_prune"):
        assert rows.loc[step, "outcome"] == "skipped"
        assert rows.loc[step, "reason"] == "no samples remain after callrate"
    params = pd.DataFrame(report["parameters"])
    assert params.loc[params["parameter"] == "mind", "value"].item() == 0.0

    # No fileset: nothing passed QC, and the unfiltered input must not stand in.
    assert not Path(f"{out}.pgen").exists()

    # Every sample is accounted for as removed by callrate, with the
    # missingness that removed it.
    fate = pd.read_csv(f"{out}_sample_fate.tsv", sep="\t")
    assert len(fate) == 500
    assert set(fate["removed_by"]) == {"callrate"}
    assert fate["F_MISS"].notna().all()

    assert (Path(f"{out}_plots") / "static" / "sample_missingness.png").is_file()


def test_callrate_detects_it_without_plink_wording(
    test_geno_path: Path, tmp_path: Path, monkeypatch
) -> None:
    # As if PLINK2 reworded its message: only callrate's own count can tell.
    monkeypatch.setattr(executors, "_no_data_remaining", lambda text: None)
    out = tmp_path / "emptied_callrate"

    with pytest.raises(NoDataRemainingError) as exc:
        filter_callrate(
            GenotypeData.from_path(test_geno_path), CallrateConfig(mind=0.0), out
        )

    assert exc.value.remaining == "samples"
    assert exc.value.output_path == str(out)


def test_plink2_still_words_it_the_way_the_backstop_matches(
    test_geno_path: Path, tmp_path: Path
) -> None:
    # The pinned PLINK2's own message, with no step count involved. If a PLINK2
    # upgrade fails this, update _no_data_remaining in core/executors.py.
    out = tmp_path / "emptied_mind"

    with pytest.raises(NoDataRemainingError) as exc:
        executors.run_plink2(
            input_path=test_geno_path, output_path=out, extra_args=["--mind", "0"]
        )

    assert exc.value.remaining == "samples"
