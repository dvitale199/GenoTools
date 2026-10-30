"""Plot-ready tables the QC steps keep, tested on hand-built frames."""

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from genotools.qc.plot_data import (
    MISSINGNESS_EDGES,
    NEG_LOG10_P_EDGES,
    failing_rows,
    histogram,
    plot_data_path,
    pvalue_histogram,
    write_plot_data,
)
from genotools.qc.steps.hwe import hwe_plot_tables


class TestHistogram:
    def test_every_bin_present_and_nans_dropped(self) -> None:
        frame = histogram(pd.Series([0.0, 0.001, 0.012, np.nan, 1.0]), MISSINGNESS_EDGES)
        assert len(frame) == len(MISSINGNESS_EDGES) - 1
        assert frame["count"].sum() == 4
        assert frame.loc[0, "count"] == 2  # [0, 0.005)
        assert frame["count"].iloc[-1] == 1  # 1.0 lands in the closed last bin

    def test_size_does_not_grow_with_input(self) -> None:
        rng = np.random.default_rng(0)
        small = histogram(pd.Series(rng.random(10)), MISSINGNESS_EDGES)
        large = histogram(pd.Series(rng.random(1_000_000)), MISSINGNESS_EDGES)
        assert len(small) == len(large)


class TestPvalueHistogram:
    def test_zero_p_lands_in_open_last_bin(self) -> None:
        frame = pvalue_histogram(pd.Series([1.0, 0.5, 1e-300, 0.0]))
        assert len(frame) == len(NEG_LOG10_P_EDGES) - 1
        assert frame["count"].sum() == 4
        assert frame["bin_upper"].iloc[-1] == np.inf
        assert frame["count"].iloc[-1] == 2  # 1e-300 and 0 are both past 50

    def test_one_histogram_per_group(self) -> None:
        p = pd.Series([0.5, 0.5, 1e-8])
        by = pd.Series(["ALL", "UNAFF", "UNAFF"], name="test")
        frame = pvalue_histogram(p, by=by)
        assert list(frame.columns) == ["test", "bin_lower", "bin_upper", "count"]
        assert frame.groupby("test")["count"].sum().to_dict() == {"ALL": 1, "UNAFF": 2}

    def test_unparseable_p_is_dropped(self) -> None:
        assert pvalue_histogram(pd.Series(["NA", 0.1]))["count"].sum() == 1


def test_failing_rows_most_extreme_first_and_capped() -> None:
    frame = pd.DataFrame({"ID": list("abcde"), "F_MISS": [0.01, 0.3, 0.2, 0.9, 0.06]})
    out = failing_rows(frame, frame["F_MISS"] > 0.05, "F_MISS", ascending=False, limit=3)
    assert list(out["ID"]) == ["d", "b", "c"]


def test_written_beside_a_dotted_prefix(tmp_path: Path) -> None:
    """#296: a prefix is not a filename."""
    out = tmp_path / "cohort.v2"
    written = write_plot_data({"f": pd.DataFrame({"x": [1]})}, out)
    assert written == {"f": tmp_path / "cohort.v2_plotdata_f.tsv"}
    assert plot_data_path(out, "f").read_text() == "x\n1\n"


class TestHwePlotTables:
    @staticmethod
    def _hardy(rows) -> pd.DataFrame:
        return pd.DataFrame(rows, columns=["SNP", "TEST", "P"])

    def test_case_control_failures_come_from_unaff(self) -> None:
        """PLINK 1.9 --hwe tests controls only on case/control data."""
        hardy = self._hardy([
            ("rs1", "ALL", 1e-9), ("rs1", "AFF", 1e-9), ("rs1", "UNAFF", 0.4),
            ("rs2", "ALL", 0.3), ("rs2", "AFF", 0.3), ("rs2", "UNAFF", 1e-7),
        ])
        tables = hwe_plot_tables(hardy, kept={"rs1"})
        failing = tables["failing"]
        assert list(failing["SNP"]) == ["rs2"]
        assert failing["test"].item() == "UNAFF"
        assert failing["P"].item() == 1e-7
        assert set(tables["pvalue_bins"]["test"]) == {"ALL", "AFF", "UNAFF"}

    def test_failures_follow_the_filter_not_the_threshold(self) -> None:
        """Whatever PLINK removed is failing, read from its snplist."""
        hardy = self._hardy([("rs1", "ALL", 0.9), ("rs2", "ALL", 0.8)])
        tables = hwe_plot_tables(hardy, kept={"rs1"})
        assert list(tables["failing"]["SNP"]) == ["rs2"]

    def test_no_phenotype_uses_the_only_test(self) -> None:
        hardy = self._hardy([("rs1", "ALL(NP)", 1e-9), ("rs2", "ALL(NP)", 0.5)])
        tables = hwe_plot_tables(hardy, kept={"rs2"})
        assert tables["failing"]["test"].item() == "ALL(NP)"


def test_related_file_keeps_only_pairs_at_the_cutoff() -> None:
    """The KING table reaches 3rd degree for the plot; .related must not."""
    from genotools.qc.steps.relatedness import split_kinship_table

    kin0 = pd.DataFrame({
        "#FID1": ["a", "b", "c"], "IID1": ["a", "b", "c"],
        "FID2": ["x", "y", "z"], "IID2": ["x", "y", "z"],
        "NSNP": [100] * 3, "HETHET": [0.1] * 3,
        "IBS0": [0.001, 0.01, 0.02],
        "KINSHIP": [0.49, 0.0884, 0.05],  # duplicate, at the cutoff, 3rd degree
    })
    plot_table, reported = split_kinship_table(kin0, related_cutoff=0.0884)
    assert list(plot_table["#FID1"]) == ["a", "b", "c"]
    assert list(plot_table.columns) == ["#FID1", "IID1", "FID2", "IID2", "IBS0", "KINSHIP"]
    assert list(reported["#FID1"]) == ["a", "b"]
    assert "NSNP" in reported.columns  # .related keeps the full kin0 shape
