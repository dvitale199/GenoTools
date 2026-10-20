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

"""Input and output prefixes whose names contain periods, through the real CLI.

--pfile and --out are PLINK prefixes: GenoTools appends .pgen, .json, ... to
them. It used to build those names with ``Path.with_suffix``, which replaces
whatever follows the last dot -- ``genotools.test.v2`` became
``genotools.test.pgen``. Every step then failed to find its input, and the run
exited 1 with the raw input copied out as ``{out}.pgen``.

The exit code catches a step that fails outright, but not a path that is wrong
without failing -- the report written to ``clean.json``, or a step reading a
different file that happens to exist. So this runs the same QC twice, under
plain names and under dotted ones, and requires the same per-step outcomes and
the same surviving samples and variants.

Revert check: restoring ``data.path.with_suffix(...)`` for ``input_file`` in
qc/steps/callrate.py, or ``out_path.with_suffix(".json")`` in cli/__init__.py,
fails this test.

Requires the synthetic test data and real PLINK.
"""

import json
import shutil
import subprocess
import sys
from pathlib import Path

from tests.regression.compare import compare_pfiles

REPO_ROOT = Path(__file__).resolve().parents[2]


def _run_qc(geno: Path, out: Path) -> dict:
    """Run --all-sample --all-variant; return the report's per-step pass_fail."""
    out.parent.mkdir(parents=True)
    result = subprocess.run(
        [sys.executable, "-m", "genotools", "--pfile", str(geno), "--out", str(out),
         "--all-sample", "--all-variant"],
        cwd=str(REPO_ROOT),
        capture_output=True,
        text=True,
        timeout=300,
    )
    assert result.returncode == 0, (
        f"genotools failed:\nstdout: {result.stdout}\nstderr: {result.stderr}"
    )
    report = Path(f"{out}.json")
    assert report.exists(), f"report not written at {report}"
    return json.loads(report.read_text())["pass_fail"]


def _outcomes(pass_fail: dict) -> dict:
    return {step: info["outcome"] for step, info in pass_fail.items()}


def test_dotted_names_give_the_same_qc_as_plain_names(
    test_geno_path: Path, tmp_path: Path
) -> None:
    dotted_in = tmp_path / "in" / "genotools.test.v2"
    dotted_in.parent.mkdir()
    for ext in (".pgen", ".pvar", ".psam"):
        shutil.copy2(f"{test_geno_path}{ext}", f"{dotted_in}{ext}")

    plain_out = tmp_path / "plain" / "clean"
    dotted_out = tmp_path / "dotted" / "clean.v2"

    plain = _run_qc(test_geno_path, plain_out)
    dotted = _run_qc(dotted_in, dotted_out)

    # with_suffix wrote the report to clean.json
    assert not (dotted_out.parent / "clean.json").exists()

    # Guard against a vacuous match where both runs fail the same way.
    assert _outcomes(plain)["callrate"] == "pass"
    assert _outcomes(dotted) == _outcomes(plain)

    comparison = compare_pfiles(plain_out, dotted_out)
    assert comparison.equal, comparison.message
