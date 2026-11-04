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

"""PCA leaves no intermediates behind when PLINK fails (#329).

On a real ancestry run, `--indep-pairwise` exited 13 on every group under 50
samples. The error raised out of `run_pca_pruning` before either cleanup path,
leaving `{out}_pca_pruned_tmp.*` and `{out}_pca_pruned_hg38.exclusion` beside
the run's real outputs. PLINK is stubbed here: each command writes what PLINK
would, and the one under test fails.
"""

from pathlib import Path
from typing import List

import pytest

from genotools.core.exceptions import ExternalToolError
from genotools.core.genotypes import GenotypeData
from genotools.gwas.config import PCAConfig, PCAPruningConfig
from genotools.gwas.results import PruningResult
from genotools.gwas.steps import pca as pca_module


def _pfiles(prefix: str) -> None:
    for ext in ("pgen", "pvar", "psam", "log"):
        Path(f"{prefix}.{ext}").write_text("x\n")


def _cohort(tmp_path: Path) -> GenotypeData:
    prefix = tmp_path / "cohort"
    _pfiles(str(prefix))
    return GenotypeData(path=prefix, format="pfile", sample_count=11, variant_count=100)


def _fake_plink(fail_on: str):
    """A run_command that writes PLINK's outputs and fails on one command."""

    def run(cmd: List[str], tool_name: str = "plink2", **_):
        out = cmd[cmd.index("--out") + 1]
        if fail_on in cmd:
            Path(f"{out}.log").write_text("Error\n")
            raise ExternalToolError(tool_name, " ".join(cmd), 13, "fewer than 50 samples")
        if "--indep-pairwise" in cmd:
            for ext in ("prune.in", "prune.out", "log"):
                Path(f"{out}.{ext}").write_text("rs1\n")
        else:
            _pfiles(out)

    return run


@pytest.fixture
def plink(monkeypatch):
    monkeypatch.setattr(pca_module, "get_plink2", lambda: "plink2")

    def install(fail_on: str) -> None:
        monkeypatch.setattr(pca_module, "run_command", _fake_plink(fail_on))

    return install


def _left_behind(tmp_path: Path) -> List[str]:
    return sorted(p.name for p in tmp_path.iterdir() if not p.name.startswith("cohort."))


@pytest.mark.parametrize("fail_on", ["--maf", "--indep-pairwise", "--extract"])
def test_pruning_removes_its_intermediates_when_plink_fails(
    tmp_path: Path, plink, fail_on: str
) -> None:
    plink(fail_on)
    out = tmp_path / "r10_EUR_pca_pruned"
    with pytest.raises(ExternalToolError):
        pca_module.run_pca_pruning(_cohort(tmp_path), PCAPruningConfig(), out)
    leftovers = [n for n in _left_behind(tmp_path) if not n.startswith("r10_EUR_pca_pruned.")]
    assert leftovers == []


def test_pruning_success_keeps_only_its_output(tmp_path: Path, plink, monkeypatch) -> None:
    plink("never")
    monkeypatch.setattr(
        pca_module.GenotypeData, "from_path",
        classmethod(lambda cls, p: GenotypeData(path=p, format="pfile", sample_count=60, variant_count=40)),
    )
    out = tmp_path / "r10_AMR_pca_pruned"
    result = pca_module.run_pca_pruning(_cohort(tmp_path), PCAPruningConfig(), out)
    assert result.success
    assert _left_behind(tmp_path) == [f"r10_AMR_pca_pruned.{e}" for e in ("log", "pgen", "psam", "pvar")]


def test_pca_removes_the_pruned_pfiles_when_plink_fails(
    tmp_path: Path, plink, monkeypatch
) -> None:
    plink("--pca")
    out = tmp_path / "r10_EUR"

    def pruned(data, config, out_path, build):
        _pfiles(str(out_path))
        return PruningResult(output_path=out_path, variants_before=100, variants_after=40,
                             exclusion_file=None, success=True, log="")

    monkeypatch.setattr(pca_module, "run_pca_pruning", pruned)
    with pytest.raises(ExternalToolError):
        pca_module.run_pca(_cohort(tmp_path), PCAConfig(n_pcs=10), out)
    assert [n for n in _left_behind(tmp_path) if "pca_pruned" in n] == []
