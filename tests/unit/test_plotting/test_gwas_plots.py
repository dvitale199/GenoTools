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

"""GWAS plots: .glm parsing (logistic and linear, covariate rows), the Manhattan,
QQ, PC and scree frames on hand-built files, availability, and render smoke tests."""

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import pytest  # noqa: E402

from genotools.plotting.plots import gwas  # noqa: E402
from genotools.plotting.report import load_report  # noqa: E402

NAMES = ["gwas_manhattan", "gwas_qq", "gwas_pcs_by_phenotype", "gwas_pca_scree", "ancestry_pca_scree"]

LOGISTIC_HEADER = "#CHROM\tPOS\tID\tREF\tALT\tA1\tTEST\tOBS_CT\tOR\tP\tERRCODE"
LINEAR_HEADER = "#CHROM\tPOS\tID\tREF\tALT\tA1\tTEST\tOBS_CT\tBETA\tSE\tT_STAT\tP\tERRCODE"


def write_logistic(path: Path, n: int = 40) -> Path:
    rows = [LOGISTIC_HEADER]
    rng = np.random.default_rng(0)
    for i in range(n):
        chrom = "X" if i >= n - 5 else str(1 + i % 3)
        p = rng.uniform(1e-6, 1)
        rows.append(f"{chrom}\t{1000 + i * 10}\trs{i}\tA\tG\tG\tADD\t100\t1.1\t{p}\t.")
        rows.append(f"{chrom}\t{1000 + i * 10}\trs{i}\tA\tG\tG\tPC1\t100\t1.1\t0.5\t.")
    path.write_text("\n".join(rows) + "\n")
    return path


def write_linear(path: Path, n: int = 40) -> Path:
    rows = [LINEAR_HEADER]
    rng = np.random.default_rng(1)
    for i in range(n):
        p = rng.uniform(1e-6, 1)
        rows.append(f"{1 + i % 3}\t{1000 + i * 10}\trs{i}\tA\tG\tG\tADD\t100\t0.1\t0.2\t0.3\t{p}\t.")
        rows.append(f"{1 + i % 3}\t{1000 + i * 10}\trs{i}\tA\tG\tG\tCOV1\t100\t0.1\t0.2\t0.3\t0.001\t.")
    path.write_text("\n".join(rows) + "\n")
    return path


def make_report(tmp_path: Path, glm: Path = None, with_pcs: bool = True, pheno=(1, 2), ancestry_eig: bool = False,
                gwas_block: bool = True, quantitative: bool = False) -> Path:
    """A report whose assoc parameters point at files under tmp_path."""
    ids = [f"s{i}" for i in range(20)]
    vec = pd.DataFrame({"#FID": ids, "IID": ids, "PC1": np.linspace(-1, 1, 20), "PC2": np.linspace(1, -1, 20)})
    eigenvec = tmp_path / "run.eigenvec"
    vec.to_csv(eigenvec, sep="\t", index=False)
    eigenval = tmp_path / "run.eigenval"
    eigenval.write_text("3.0\n2.0\n1.5\n")
    if quantitative:
        values = list(np.linspace(0.5, 9.5, 20))
    else:
        values = [pheno[i % 2] for i in range(20)]
        values[0] = -9
    params = [
        ("assoc", "gwas.pheno_name", "PHENO1", "requested"),
        ("assoc", "eigenvec_path", str(eigenvec) if with_pcs else None, "resolved"),
        ("assoc", "eigenval_path", str(eigenval), "resolved"),
        ("assoc", "gwas_output_path", str(glm) if glm else None, "resolved"),
    ]
    raw = {
        "QC": pd.DataFrame([{"step": "assoc", "pruned_count": 0, "metric": "x", "ancestry": "all",
                             "level": "sample", "pass": True, "outcome": "pass", "reason": None}]).to_dict(),
        "input_samples": pd.DataFrame({"#FID": ids, "IID": ids, "PHENO1": values}).to_dict(),
        "parameters": pd.DataFrame(params, columns=["step", "parameter", "value", "source"])
        .assign(ancestry="all").to_dict(),
        "plot_data": str(tmp_path / "run_plot_data"),
    }
    if gwas_block:
        raw["GWAS"] = pd.DataFrame({
            "value": [1.2345, 1.0987, 1.5, 1.4], "ancestry": "all",
            "metric": ["lambda", "lambda1000", "lambda_maf", "lambda1000_maf"],
        }).to_dict()
    if ancestry_eig:
        d = tmp_path / "run_plot_data" / "ancestry"
        d.mkdir(parents=True)
        pd.DataFrame({"PC": ["PC1", "PC2", "PC3"], "eigenvalue": [5.0, 2.0, 1.0],
                      "explained_variance_ratio": [0.5, 0.3, 0.1]}).to_csv(
            d / "pca_eigenvalues.tsv", sep="\t", index=False)
    path = tmp_path / "run.json"
    path.write_text(json.dumps(raw))
    return path


def by_name(name: str):
    return next(p for p in gwas.PLOTS if p.name == name)


class TestReadGlm:
    def test_logistic_keeps_add_rows(self, tmp_path) -> None:
        glm = gwas.read_glm(write_logistic(tmp_path / "a.glm.logistic.hybrid"))
        assert len(glm) == 40
        assert set(glm.columns) == {"chrom", "pos", "id", "p"}

    def test_linear_drops_covariate_rows(self, tmp_path) -> None:
        glm = gwas.read_glm(write_linear(tmp_path / "a.PHENO1.glm.linear"))
        assert len(glm) == 40
        assert (glm["p"] != 0.001).all()  # the COV1 rows

    def test_bad_p_dropped_and_zero_floored(self, tmp_path) -> None:
        path = tmp_path / "x.glm.linear"
        path.write_text(
            LINEAR_HEADER + "\n"
            "1\t1\ta\tA\tG\tG\tADD\t9\t.\t.\t.\t0\t.\n"
            "1\t2\tb\tA\tG\tG\tADD\t9\t.\t.\t.\tNA\t.\n"
            "1\t3\tc\tA\tG\tG\tADD\t9\t.\t.\t.\t0.5\t.\n"
        )
        glm = gwas.read_glm(path)
        assert list(glm["id"]) == ["a", "c"]
        assert glm["p"].iloc[0] == gwas.P_FLOOR


class TestManhattan:
    def test_positions_increase_across_chromosomes(self, tmp_path) -> None:
        points, ticks = gwas.manhattan_frame(gwas.read_glm(write_logistic(tmp_path / "a.glm")))
        assert list(ticks["label"]) == ["1", "2", "3", "X"]
        assert points["x"].is_monotonic_increasing
        assert set(points["band"]) == {0, 1}
        assert np.allclose(points["neglogp"], -np.log10(points["p"]))
        # Each tick sits inside its own chromosome's span.
        for _, t in ticks.iterrows():
            xs = points.loc[points["chrom"] == t["chrom"], "x"]
            assert xs.min() <= t["x"] <= xs.max()

    def test_unplaceable_chromosomes_dropped(self) -> None:
        glm = pd.DataFrame({"chrom": ["1", "chrUn", "chr2"], "pos": [1, 2, 3], "id": list("abc"), "p": [0.1] * 3})
        points, _ = gwas.manhattan_frame(glm)
        assert list(points["id"]) == ["a", "c"]

    def test_empty(self) -> None:
        points, ticks = gwas.manhattan_frame(pd.DataFrame({"chrom": [], "pos": [], "id": [], "p": []}))
        assert points.empty and ticks.empty


class TestQQ:
    def test_expected_and_order(self) -> None:
        frame = gwas.qq_frame(pd.DataFrame({"p": [0.5, 0.01, 0.2, 0.9]}))
        assert len(frame) == 4
        # Most significant last: p = 0.01 against the top expected value.
        assert frame["observed"].iloc[-1] == pytest.approx(2.0)
        assert frame["expected"].iloc[-1] == pytest.approx(-np.log10(0.5 / 4))
        assert frame["observed"].iloc[0] == pytest.approx(-np.log10(0.9))

    def test_thinning_keeps_tail(self) -> None:
        p = np.random.default_rng(0).uniform(size=50_000)
        full = gwas.qq_frame(pd.DataFrame({"p": p}))
        thin = gwas.qq_frame(pd.DataFrame({"p": p}), max_points=1000)
        assert len(thin) < len(full)
        assert (thin["expected"] >= 3.1).sum() == (full["expected"] >= 3.1).sum()

    def test_lambda_values_all_present(self, tmp_path) -> None:
        report = load_report(make_report(tmp_path))
        values = dict(gwas.lambda_values(report))
        assert values["λ GC"] == pytest.approx(1.2345)
        assert len(values) == 4
        assert "1.234" in gwas._lambda_text(report) or "1.235" in gwas._lambda_text(report)

    def test_lambda_absent(self, tmp_path) -> None:
        assert gwas.lambda_values(load_report(make_report(tmp_path, gwas_block=False))) == []


class TestPcs:
    def test_case_control_join(self, tmp_path) -> None:
        frame = gwas.pc_pheno_frame(load_report(make_report(tmp_path)))
        assert len(frame) == 20
        counts = frame["group"].value_counts()
        assert counts["missing"] == 1 and counts["control"] + counts["case"] == 19

    def test_quantitative(self, tmp_path) -> None:
        frame = gwas.pc_pheno_frame(load_report(make_report(tmp_path, quantitative=True)))
        assert (frame["group"] == "quantitative").all()


class TestScree:
    def test_eigenval_file(self, tmp_path) -> None:
        path = tmp_path / "x.eigenval"
        path.write_text("3.0\n2.0\n")
        frame = gwas.scree_frame(path)
        assert list(frame["eigenvalue"]) == [3.0, 2.0]
        assert list(frame["label"]) == ["PC1", "PC2"]
        assert frame["ratio"].isna().all()

    def test_ancestry_tsv(self, tmp_path) -> None:
        make_report(tmp_path, ancestry_eig=True)
        frame = gwas.scree_frame(tmp_path / "run_plot_data" / "ancestry" / "pca_eigenvalues.tsv")
        assert list(frame["ratio"]) == [0.5, 0.3, 0.1]


class TestAvailability:
    def test_missing_files_are_unavailable_not_errors(self, tmp_path) -> None:
        report = load_report(make_report(tmp_path, glm=tmp_path / "nope.glm.linear"))
        assert not by_name("gwas_manhattan").available(report)
        assert not by_name("gwas_qq").available(report)
        assert not by_name("ancestry_pca_scree").available(report)
        assert by_name("gwas_pca_scree").available(report)
        assert by_name("gwas_pcs_by_phenotype").available(report)

    def test_report_without_parameters(self, tmp_path) -> None:
        path = tmp_path / "bare.json"
        path.write_text(json.dumps({}))
        report = load_report(path)
        assert [p.name for p in gwas.PLOTS if p.available(report)] == []

    def test_no_eigenvec_path(self, tmp_path) -> None:
        report = load_report(make_report(tmp_path, with_pcs=False))
        assert not by_name("gwas_pcs_by_phenotype").available(report)

    def test_ancestry_scree_available_with_file(self, tmp_path) -> None:
        report = load_report(make_report(tmp_path, ancestry_eig=True))
        assert by_name("ancestry_pca_scree").available(report)

    def test_names(self) -> None:
        assert [p.name for p in gwas.PLOTS] == NAMES


@pytest.fixture(params=["logistic", "linear"])
def full_report(request, tmp_path) -> Path:
    glm = (write_logistic if request.param == "logistic" else write_linear)(tmp_path / f"r.glm.{request.param}")
    return make_report(tmp_path, glm=glm, ancestry_eig=True)


class TestRender:
    def test_every_gwas_plot_renders_png_and_plotly(self, full_report, tmp_path) -> None:
        report = load_report(full_report)
        for plot in gwas.PLOTS:
            assert plot.available(report), plot.name
            out = plot.png(report, plt, tmp_path / f"{plot.name}.png")
            assert out.stat().st_size > 1000
            fig = plot.plotly(report)
            assert len(fig.data) >= 1
            fig.to_json()

    def test_quantitative_pcs_render(self, tmp_path) -> None:
        report = load_report(make_report(tmp_path, quantitative=True))
        plot = by_name("gwas_pcs_by_phenotype")
        assert plot.png(report, plt, tmp_path / "q.png").exists()
        plot.plotly(report).to_json()

    def test_qq_annotation_carries_lambda(self, full_report) -> None:
        fig = by_name("gwas_qq").plotly(load_report(full_report))
        assert "1.234" in fig.layout.annotations[0].text or "1.235" in fig.layout.annotations[0].text
