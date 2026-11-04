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

"""#303: the per-sample fate table, assembled from hand-built frames."""

from pathlib import Path

import numpy as np
import pandas as pd

from genotools.cli.output import PipelineOutput
from genotools.cli.sample_fate import (
    FATE_COLUMNS,
    assemble_sample_fate,
    collect_sample_stats,
    kinship_max,
    read_plot_tables,
    write_sample_fate,
)


def _psam(n: int = 6) -> pd.DataFrame:
    return pd.DataFrame({"#FID": [f"f{i}" for i in range(n)], "IID": [f"i{i}" for i in range(n)]})


def _pruned(rows) -> pd.DataFrame:
    return pd.DataFrame(rows, columns=["FID", "IID", "step"])


def _removed(fate: pd.DataFrame) -> dict:
    return dict(zip(fate["IID"], fate["removed_by"]))


class TestAssemble:
    def test_every_input_sample_exactly_once_in_input_order(self) -> None:
        fate = assemble_sample_fate(_psam(), _pruned([("f2", "i2", "het")]))
        assert list(fate.columns) == FATE_COLUMNS
        assert list(fate["IID"]) == [f"i{i}" for i in range(6)]
        assert fate.duplicated(["FID", "IID"]).sum() == 0

    def test_removed_by_matches_pruned_samples(self) -> None:
        pruned = _pruned([("f1", "i1", "callrate"), ("f4", "i4", "het")])
        fate = assemble_sample_fate(_psam(), pruned)
        assert _removed(fate) == {
            "i0": "kept", "i1": "callrate", "i2": "kept",
            "i3": "kept", "i4": "het", "i5": "kept",
        }

    def test_first_step_wins_when_listed_twice(self) -> None:
        pruned = _pruned([("f1", "i1", "callrate"), ("f1", "i1", "het")])
        assert _removed(assemble_sample_fate(_psam(), pruned))["i1"] == "callrate"

    def test_min_samples_drops_reported(self) -> None:
        pruned = pd.DataFrame(
            {
                "FID": ["f0", "f3"],
                "IID": ["i0", "i3"],
                "step": "insufficient_ancestry_sample_n",
                "label": "SAS",
            }
        )
        labels = pd.DataFrame(
            {"FID": [f"f{i}" for i in range(6)], "IID": [f"i{i}" for i in range(6)],
             "label": ["SAS", "EUR", "EUR", "SAS", "EUR", "EUR"]}
        )
        fate = assemble_sample_fate(_psam(), pruned, labels)
        assert _removed(fate)["i0"] == "insufficient_ancestry_sample_n"
        assert _removed(fate)["i3"] == "insufficient_ancestry_sample_n"
        assert dict(zip(fate["IID"], fate["ancestry"]))["i0"] == "SAS"

    def test_flat_run_ancestry_is_all(self) -> None:
        fate = assemble_sample_fate(_psam(), None, None)
        assert set(fate["ancestry"]) == {"all"}
        assert set(fate["removed_by"]) == {"kept"}

    def test_ancestry_run_uses_predicted_label(self) -> None:
        labels = pd.DataFrame(
            {"FID": ["f0", "f1"], "IID": ["i0", "i1"], "label": ["EUR", "AFR"]}
        )
        fate = assemble_sample_fate(_psam(2), None, labels)
        assert list(fate["ancestry"]) == ["EUR", "AFR"]

    def test_numeric_ids_join_across_dtypes(self) -> None:
        psam = pd.DataFrame({"#FID": [1, 2], "IID": [10, 20]})
        pruned = pd.DataFrame({"FID": ["2"], "IID": ["20"], "step": ["sex"]})
        assert _removed(assemble_sample_fate(psam, pruned)) == {"10": "kept", "20": "sex"}

    def test_psam_without_fid(self) -> None:
        psam = pd.DataFrame({"#IID": ["a", "b"]})
        pruned = pd.DataFrame({"FID": ["0"], "IID": ["b"], "step": ["het"]})
        fate = assemble_sample_fate(psam, pruned)
        assert _removed(fate) == {"a": "kept", "b": "het"}

    def test_missing_stats_are_nan_not_a_crash(self) -> None:
        fate = assemble_sample_fate(_psam(), None, None, None)
        assert fate[["F_MISS", "sex_F", "het_F", "max_kinship"]].isna().all().all()

    def test_stats_attach_by_sample(self) -> None:
        tables = [{
            "callrate_sample_missingness": pd.DataFrame(
                {"#FID": ["f0", "f1"], "IID": ["i0", "i1"], "F_MISS": [0.01, 0.2]}
            ),
            "het_f": pd.DataFrame({"#FID": ["f1"], "IID": ["i1"], "F": [-0.3]}),
        }]
        fate = assemble_sample_fate(_psam(3), None, None, collect_sample_stats(tables))
        by = fate.set_index("IID")
        assert by.loc["i1", "F_MISS"] == 0.2
        assert by.loc["i1", "het_F"] == -0.3
        assert np.isnan(by.loc["i0", "het_F"]) and np.isnan(by.loc["i2", "F_MISS"])


class TestStats:
    def test_kinship_max_over_both_sides_of_a_pair(self) -> None:
        kin = pd.DataFrame(
            {"#FID1": ["a", "a"], "IID1": ["1", "1"], "FID2": ["b", "c"],
             "IID2": ["2", "3"], "IBS0": [0, 0], "KINSHIP": [0.3, 0.1]}
        )
        got = kinship_max(kin).set_index("IID")["max_kinship"]
        assert got.to_dict() == {"1": 0.3, "2": 0.3, "3": 0.1}

    def test_sex_uses_all_x_pass(self) -> None:
        sex = pd.DataFrame(
            {"#FID": ["f0", "f0"], "IID": ["i0", "i0"], "F": [0.9, 0.1],
             "check": ["all_x", "filtered_x"]}
        )
        stats = collect_sample_stats([{"sex_check": sex}])
        assert stats.loc[0, "sex_F"] == 0.9

    def test_groups_are_concatenated(self) -> None:
        a = {"het_f": pd.DataFrame({"#FID": ["f0"], "IID": ["i0"], "F": [0.1]})}
        b = {"het_f": pd.DataFrame({"#FID": ["f1"], "IID": ["i1"], "F": [0.2]})}
        stats = collect_sample_stats([a, b])
        assert sorted(stats["IID"]) == ["i0", "i1"]

    def test_no_tables(self) -> None:
        assert collect_sample_stats([]).empty

    def test_read_plot_tables_skips_missing_and_ancestry_dir(self, tmp_path: Path) -> None:
        assert read_plot_tables(None) == []
        assert read_plot_tables(tmp_path / "nope") == []
        (tmp_path / "ancestry").mkdir()
        (tmp_path / "ancestry" / "het_f.tsv").write_text("#FID\tIID\tF\nx\ty\t1\n")
        (tmp_path / "EUR").mkdir()
        (tmp_path / "EUR" / "het_f.tsv").write_text("#FID\tIID\tF\nf0\ti0\t0.5\n")
        groups = read_plot_tables(tmp_path)
        assert len(groups) == 1 and groups[0]["het_f"]["F"].iloc[0] == 0.5


class TestWrite:
    def test_dotted_prefix_not_rewritten(self, tmp_path: Path) -> None:
        out = tmp_path / "cohort.v2"
        path = write_sample_fate(out, _psam(2), None, None, None)
        assert path == Path(f"{out}_sample_fate.tsv")
        assert path.name == "cohort.v2_sample_fate.tsv"
        table = pd.read_csv(path, sep="\t")
        assert list(table.columns) == FATE_COLUMNS and len(table) == 2

    def test_no_input_samples_writes_nothing(self, tmp_path: Path) -> None:
        assert write_sample_fate(tmp_path / "o", None, None, None, None) is None

    def test_unwritable_location_warns_not_raises(self, tmp_path: Path) -> None:
        assert write_sample_fate(tmp_path / "nodir" / "o", _psam(2), None, None, None) is None


class TestReportKey:
    def test_absent_by_default(self) -> None:
        assert "sample_fate" not in PipelineOutput().to_dict()

    def test_recorded(self) -> None:
        assert PipelineOutput(sample_fate="x_sample_fate.tsv").to_dict()["sample_fate"] == "x_sample_fate.tsv"
