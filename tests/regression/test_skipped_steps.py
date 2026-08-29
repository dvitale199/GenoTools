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

"""End-to-end behavior of a data-driven step skip, through the real CLI.

A step the data rules out is reported, not failed: the run exits 0, the JSON
carries the step with ``outcome="skipped"`` and a reason, and the pfile chain
still reaches ``{out}``. Unit tests cover the decision; this covers the whole
path, which is where both skip bugs so far were actually found.

Requires the synthetic test data and real PLINK; skips cleanly otherwise.
"""

import json
import subprocess
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]


def _run_cli(args: list[str]) -> subprocess.CompletedProcess:
    return subprocess.run(
        [sys.executable, "-m", "genotools", *args],
        cwd=str(REPO_ROOT),
        capture_output=True,
        text=True,
        timeout=300,
    )


def _pfile_with_column(src: Path, dst: Path, column: str, value_of) -> Path:
    """Copy src's pfiles to dst, rewriting one psam column as value_of(row)."""
    dst.with_suffix(".pgen").write_bytes(src.with_suffix(".pgen").read_bytes())
    dst.with_suffix(".pvar").write_bytes(src.with_suffix(".pvar").read_bytes())

    lines = src.with_suffix(".psam").read_text().splitlines()
    col = lines[0].split().index(column)
    rows = [lines[0]]
    for i, line in enumerate(lines[1:]):
        fields = line.split()
        fields[col] = str(value_of(i))
        rows.append("\t".join(fields))
    dst.with_suffix(".psam").write_text("\n".join(rows) + "\n")
    return dst


def _qc_rows(out: Path) -> list[dict]:
    """QC rows from the results JSON, which stores a column-oriented frame."""
    qc = json.loads(Path(f"{out}.json").read_text())["QC"]
    return [
        {col: values[key] for col, values in qc.items()}
        for key in qc["step"]
    ]


@pytest.mark.parametrize("skip_fails", [False, True])
def test_case_control_without_both_phenotypes_is_skipped(
    test_geno_path: Path, tmp_path: Path, skip_fails: bool
) -> None:
    """One decision, one behavior. Case-control missingness needs both cases
    and controls; with only cases it must be reported as skipped rather than
    raising - and parametrized over --skip-fails because that flag suppresses
    the cohort-level decision, leaving only the per-dataset one.
    """
    cases_only = _pfile_with_column(
        test_geno_path, tmp_path / "cases_only", "PHENO1", lambda i: 2
    )
    out = tmp_path / f"skipcc{int(skip_fails)}"
    args = ["--pfile", str(cases_only), "--out", str(out), "--geno", "--case-control"]
    if skip_fails:
        args.append("--skip-fails")

    res = _run_cli(args)

    assert res.returncode == 0, (
        f"a skipped step must not fail the run:\n{res.stdout}\n{res.stderr}"
    )
    assert "Traceback" not in res.stderr, res.stderr

    rows = {row["step"]: row for row in _qc_rows(out)}
    skipped = rows["case_control_missingness_prune"]
    assert skipped["outcome"] == "skipped"
    assert skipped["pruned_count"] == 0
    assert "not both" in skipped["reason"]
    assert rows["geno_prune"]["outcome"] == "pass"

    # The trailing skip writes nothing, so the last step that ran owns {out}.
    assert Path(f"{out}.pgen").exists(), "trailing skip swallowed the output"


def test_case_control_runs_when_both_phenotypes_present(
    test_geno_path: Path, tmp_path: Path
) -> None:
    """The negative control: the skip must not fire on data that supports the
    step, or it would silently disable case-control QC for every real cohort.
    """
    both = _pfile_with_column(
        test_geno_path, tmp_path / "both_phenos", "PHENO1", lambda i: 1 + (i % 2)
    )
    out = tmp_path / "runcc"

    res = _run_cli(
        ["--pfile", str(both), "--out", str(out), "--geno", "--case-control"]
    )

    assert res.returncode == 0, f"{res.stdout}\n{res.stderr}"
    rows = {row["step"]: row for row in _qc_rows(out)}
    assert rows["case_control_missingness_prune"]["outcome"] == "pass"
    assert rows["case_control_missingness_prune"]["reason"] is None


@pytest.mark.parametrize("skip_fails", [False, True])
def test_sex_without_recorded_sex_is_skipped(
    test_geno_path: Path, tmp_path: Path, skip_fails: bool
) -> None:
    """--check-sex has nothing to compare its calls against when no sample sex
    is recorded, so the step is reported as skipped rather than left to fail.
    """
    unsexed = _pfile_with_column(
        test_geno_path, tmp_path / "unsexed", "SEX", lambda i: 0
    )
    out = tmp_path / f"skipsex{int(skip_fails)}"
    args = ["--pfile", str(unsexed), "--out", str(out), "--geno", "--sex"]
    if skip_fails:
        args.append("--skip-fails")

    res = _run_cli(args)

    assert res.returncode == 0, (
        f"a skipped step must not fail the run:\n{res.stdout}\n{res.stderr}"
    )
    assert "Traceback" not in res.stderr, res.stderr

    rows = {row["step"]: row for row in _qc_rows(out)}
    skipped = rows["sex_prune"]
    assert skipped["outcome"] == "skipped"
    assert skipped["pruned_count"] == 0
    assert skipped["reason"] == "no sample sex data is available"
    assert Path(f"{out}.pgen").exists(), "trailing skip swallowed the output"


def test_sex_runs_when_sex_is_recorded(
    test_geno_path: Path, tmp_path: Path
) -> None:
    """The negative control: the skip must not fire on data that supports the
    step, or sex QC would silently stop running for every real cohort.
    """
    out = tmp_path / "runsex"

    res = _run_cli(
        ["--pfile", str(test_geno_path), "--out", str(out), "--geno", "--sex"]
    )

    assert res.returncode == 0, f"{res.stdout}\n{res.stderr}"
    rows = {row["step"]: row for row in _qc_rows(out)}
    assert rows["sex_prune"]["outcome"] == "pass"
    assert rows["sex_prune"]["reason"] is None
