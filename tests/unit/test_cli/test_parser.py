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

"""Tests for CLI argument parsing."""

from pathlib import Path
from typing import List

import pytest

from genotools.cli.parser import (
    InputArgs,
    SampleQCArgs,
    VariantQCArgs,
    AncestryArgs,
    GWASArgs,
    OutputArgs,
    PipelineArgs,
    create_parser,
    parse_args,
    _parse_sex_args,
    _parse_het_args,
    _parse_ld_args,
)
from genotools.ancestry.config import InferenceMode


class TestInputArgs:
    """Tests for InputArgs dataclass."""

    def test_pfile_input(self) -> None:
        """InputArgs accepts pfile input."""
        args = InputArgs(pfile=Path("/data/test"))
        assert args.geno_path == Path("/data/test")
        assert args.input_format == "pfile"

    def test_bfile_input(self) -> None:
        """InputArgs accepts bfile input."""
        args = InputArgs(bfile=Path("/data/test"))
        assert args.geno_path == Path("/data/test")
        assert args.input_format == "bfile"

    def test_vcf_input(self) -> None:
        """InputArgs accepts vcf input."""
        args = InputArgs(vcf=Path("/data/test.vcf"))
        assert args.geno_path == Path("/data/test")
        assert args.input_format == "vcf"

    def test_vcf_gz_input(self) -> None:
        """InputArgs handles .vcf.gz input."""
        args = InputArgs(vcf=Path("/data/test.vcf.gz"))
        assert args.geno_path == Path("/data/test")

    def test_no_input_raises(self) -> None:
        """InputArgs raises without any input."""
        with pytest.raises(ValueError, match="At least one input"):
            InputArgs()

    def test_pfile_takes_precedence(self) -> None:
        """pfile takes precedence over bfile."""
        args = InputArgs(pfile=Path("/data/pfile"), bfile=Path("/data/bfile"))
        assert args.geno_path == Path("/data/pfile")
        assert args.input_format == "pfile"


class TestSampleQCArgs:
    """Tests for SampleQCArgs dataclass."""

    def test_default_values(self) -> None:
        """SampleQCArgs has correct defaults."""
        args = SampleQCArgs()
        assert args.run_callrate is False
        assert args.callrate_threshold == 0.02
        assert args.run_sex is False
        assert args.female_max_f == 0.25
        assert args.male_min_f == 0.75
        assert args.run_het is False
        assert args.het_lower == -0.15
        assert args.het_upper == 0.15
        assert args.run_related is False
        assert args.related_cutoff == 0.0884
        assert args.duplicated_cutoff == 0.354

    def test_to_callrate_config(self) -> None:
        """to_callrate_config creates correct config."""
        args = SampleQCArgs(callrate_threshold=0.05)
        config = args.to_callrate_config()
        assert config.mind == 0.05

    def test_to_sex_config(self) -> None:
        """to_sex_config creates correct config."""
        args = SampleQCArgs(female_max_f=0.2, male_min_f=0.8)
        config = args.to_sex_config()
        assert config.female_max_f == 0.2
        assert config.male_min_f == 0.8

    def test_to_het_config(self) -> None:
        """to_het_config creates correct config."""
        args = SampleQCArgs(het_lower=-0.2, het_upper=0.2)
        config = args.to_het_config()
        assert config.f_lower == -0.2
        assert config.f_upper == 0.2

    def test_to_related_config(self) -> None:
        """to_related_config creates correct config."""
        args = SampleQCArgs(
            related_cutoff=0.1,
            duplicated_cutoff=0.4,
            prune_related=True,
            prune_duplicated=True,
        )
        config = args.to_related_config()
        assert config.related_cutoff == 0.1
        assert config.duplicated_cutoff == 0.4
        assert config.prune_related is True
        assert config.prune_duplicated is True


class TestVariantQCArgs:
    """Tests for VariantQCArgs dataclass."""

    def test_default_values(self) -> None:
        """VariantQCArgs has correct defaults."""
        args = VariantQCArgs()
        assert args.run_geno is False
        assert args.geno_threshold == 0.05
        assert args.run_case_control is False
        assert args.run_haplotype is False
        assert args.run_hwe is False
        assert args.filter_controls is False
        assert args.run_ld is False
        assert args.ld_window_size == 50
        assert args.ld_step_size == 5
        assert args.ld_r2_threshold == 0.5

    def test_to_geno_config(self) -> None:
        """to_geno_config creates correct config."""
        args = VariantQCArgs(geno_threshold=0.02)
        config = args.to_geno_config()
        assert config.geno == 0.02

    def test_to_ld_config(self) -> None:
        """to_ld_config creates correct config."""
        args = VariantQCArgs(ld_window_size=100, ld_step_size=10, ld_r2_threshold=0.2)
        config = args.to_ld_config()
        assert config.window_size == 100
        assert config.step_size == 10
        assert config.r2_threshold == 0.2


class TestAncestryArgs:
    """Tests for AncestryArgs dataclass."""

    def test_default_values(self) -> None:
        """AncestryArgs has correct defaults."""
        args = AncestryArgs()
        assert args.run_ancestry is False
        assert args.ref_panel is None
        assert args.use_container is False
        assert args.use_singularity is False
        assert args.use_cloud is False
        assert args.inference_mode == InferenceMode.LOCAL

    def test_container_mode(self) -> None:
        """Container mode is detected correctly."""
        args = AncestryArgs(use_container=True)
        assert args.inference_mode == InferenceMode.CONTAINER

    def test_singularity_mode(self) -> None:
        """Singularity mode is detected correctly."""
        args = AncestryArgs(use_singularity=True)
        assert args.inference_mode == InferenceMode.SINGULARITY

    def test_cloud_mode(self) -> None:
        """Cloud mode is detected correctly."""
        args = AncestryArgs(use_cloud=True)
        assert args.inference_mode == InferenceMode.CLOUD

    def test_container_with_model_raises(self) -> None:
        """Using container with model path raises error."""
        with pytest.raises(ValueError, match="Cannot use both"):
            AncestryArgs(use_container=True, model_path=Path("/path/to/model.pkl"))


class TestGWASArgs:
    """Tests for GWASArgs dataclass."""

    def test_default_values(self) -> None:
        """GWASArgs has correct defaults."""
        args = GWASArgs()
        assert args.run_pca is False
        assert args.n_pcs == 10
        assert args.build == "hg38"
        assert args.run_gwas is False


class TestOutputArgs:
    """Tests for OutputArgs dataclass."""

    def test_default_values(self) -> None:
        """OutputArgs has correct defaults."""
        args = OutputArgs()
        assert args.full_output is False
        assert args.skip_fails is False
        assert args.warn_only is True


class TestPipelineArgs:
    """Tests for PipelineArgs dataclass."""

    @pytest.fixture
    def basic_args(self) -> PipelineArgs:
        """Create basic pipeline args."""
        return PipelineArgs(
            input=InputArgs(pfile=Path("/data/test")),
            output=OutputArgs(out_path=Path("/output/test")),
        )

    def test_geno_path_property(self, basic_args: PipelineArgs) -> None:
        """geno_path property returns correct path."""
        assert basic_args.geno_path == Path("/data/test")

    def test_out_path_property(self, basic_args: PipelineArgs) -> None:
        """out_path property returns correct path."""
        assert basic_args.out_path == Path("/output/test")

    def test_get_enabled_sample_steps(self) -> None:
        """get_enabled_sample_steps returns enabled steps."""
        args = PipelineArgs(
            input=InputArgs(pfile=Path("/data/test")),
            output=OutputArgs(out_path=Path("/output/test")),
            sample_qc=SampleQCArgs(run_callrate=True, run_sex=True),
        )
        steps = args.get_enabled_sample_steps()
        assert steps == ["callrate", "sex"]

    def test_get_enabled_variant_steps(self) -> None:
        """get_enabled_variant_steps returns enabled steps."""
        args = PipelineArgs(
            input=InputArgs(pfile=Path("/data/test")),
            output=OutputArgs(out_path=Path("/output/test")),
            variant_qc=VariantQCArgs(run_geno=True, run_hwe=True),
        )
        steps = args.get_enabled_variant_steps()
        assert steps == ["hwe", "geno"]

    def test_get_all_enabled_steps_with_assoc(self) -> None:
        """get_all_enabled_steps includes assoc when PCA/GWAS enabled."""
        args = PipelineArgs(
            input=InputArgs(pfile=Path("/data/test")),
            output=OutputArgs(out_path=Path("/output/test")),
            gwas=GWASArgs(run_pca=True),
        )
        steps = args.get_all_enabled_steps()
        assert "assoc" in steps

    def test_has_any_qc_steps(self) -> None:
        """has_any_qc_steps returns correct value."""
        args_no_qc = PipelineArgs(
            input=InputArgs(pfile=Path("/data/test")),
            output=OutputArgs(out_path=Path("/output/test")),
        )
        assert args_no_qc.has_any_qc_steps() is False

        args_with_qc = PipelineArgs(
            input=InputArgs(pfile=Path("/data/test")),
            output=OutputArgs(out_path=Path("/output/test")),
            sample_qc=SampleQCArgs(run_callrate=True),
        )
        assert args_with_qc.has_any_qc_steps() is True

    def test_to_legacy_dict(self, basic_args: PipelineArgs) -> None:
        """to_legacy_dict produces correct format."""
        legacy = basic_args.to_legacy_dict()
        assert legacy["pfile"] == "/data/test"
        assert legacy["out"] == "/output/test"
        assert legacy["warn"] is True
        assert legacy["callrate"] is None  # Not enabled
        assert legacy["ancestry"] is False


class TestParseSexArgs:
    """Tests for _parse_sex_args helper."""

    def test_none_input(self) -> None:
        """None input returns disabled state."""
        run, female, male = _parse_sex_args(None)
        assert run is False
        assert female == 0.25
        assert male == 0.75

    def test_empty_list(self) -> None:
        """Empty list enables with defaults."""
        run, female, male = _parse_sex_args([])
        assert run is True
        assert female == 0.25
        assert male == 0.75

    def test_two_values(self) -> None:
        """Two values are parsed correctly."""
        run, female, male = _parse_sex_args([0.2, 0.8])
        assert run is True
        assert female == 0.2
        assert male == 0.8

    def test_wrong_count_raises(self) -> None:
        """Wrong value count raises error."""
        with pytest.raises(ValueError, match="requires 0 or 2 values"):
            _parse_sex_args([0.2])


class TestParseHetArgs:
    """Tests for _parse_het_args helper."""

    def test_none_input(self) -> None:
        """None input returns disabled state."""
        run, lower, upper = _parse_het_args(None)
        assert run is False
        assert lower == -0.15
        assert upper == 0.15

    def test_empty_list(self) -> None:
        """Empty list enables with defaults."""
        run, lower, upper = _parse_het_args([])
        assert run is True
        assert lower == -0.15
        assert upper == 0.15

    def test_two_values(self) -> None:
        """Two values are parsed correctly."""
        run, lower, upper = _parse_het_args([-0.2, 0.2])
        assert run is True
        assert lower == -0.2
        assert upper == 0.2


class TestParseLdArgs:
    """Tests for _parse_ld_args helper."""

    def test_none_input(self) -> None:
        """None input returns disabled state."""
        run, window, step, r2 = _parse_ld_args(None)
        assert run is False
        assert window == 50
        assert step == 5
        assert r2 == 0.5

    def test_empty_list(self) -> None:
        """Empty list enables with defaults."""
        run, window, step, r2 = _parse_ld_args([])
        assert run is True
        assert window == 50
        assert step == 5
        assert r2 == 0.5

    def test_three_values(self) -> None:
        """Three values are parsed correctly."""
        run, window, step, r2 = _parse_ld_args([100.0, 10.0, 0.2])
        assert run is True
        assert window == 100
        assert step == 10
        assert r2 == 0.2

    def test_wrong_count_raises(self) -> None:
        """Wrong value count raises error."""
        with pytest.raises(ValueError, match="requires 0 or 3 values"):
            _parse_ld_args([50.0, 5.0])


class TestCreateParser:
    """Tests for create_parser function."""

    def test_parser_created(self) -> None:
        """Parser is created successfully."""
        parser = create_parser()
        assert parser.prog == "genotools"

    def test_help_shows_groups(self) -> None:
        """Parser has expected argument groups."""
        parser = create_parser()
        group_names = [g.title for g in parser._action_groups]
        assert "Input files" in group_names
        assert "Output configuration" in group_names
        assert "Sample-level QC" in group_names
        assert "Variant-level QC" in group_names
        assert "Ancestry prediction" in group_names
        assert "GWAS and PCA" in group_names


class TestParseArgs:
    """Tests for parse_args function."""

    def test_minimal_args(self) -> None:
        """Minimal args are parsed correctly."""
        args = parse_args(["--pfile", "/data/test", "--out", "/output/test"])
        assert args.geno_path == Path("/data/test")
        assert args.out_path == Path("/output/test")

    def test_no_input_raises(self) -> None:
        """Missing input raises error."""
        with pytest.raises(SystemExit):
            parse_args(["--out", "/output/test"])

    def test_callrate_flag_only(self) -> None:
        """--callrate without value uses default."""
        args = parse_args([
            "--pfile", "/data/test",
            "--out", "/output/test",
            "--callrate",
        ])
        assert args.sample_qc.run_callrate is True
        assert args.sample_qc.callrate_threshold == 0.02

    def test_callrate_with_value(self) -> None:
        """--callrate with value uses provided value."""
        args = parse_args([
            "--pfile", "/data/test",
            "--out", "/output/test",
            "--callrate", "0.05",
        ])
        assert args.sample_qc.run_callrate is True
        assert args.sample_qc.callrate_threshold == 0.05

    def test_all_sample_flag(self) -> None:
        """--all-sample enables all sample QC."""
        args = parse_args([
            "--pfile", "/data/test",
            "--out", "/output/test",
            "--all-sample",
        ])
        assert args.sample_qc.run_callrate is True
        assert args.sample_qc.run_sex is True
        assert args.sample_qc.run_het is True
        assert args.sample_qc.run_related is True
        # all_sample uses different default
        assert args.sample_qc.callrate_threshold == 0.05

    def test_all_variant_flag(self) -> None:
        """--all-variant enables variant QC (except LD)."""
        args = parse_args([
            "--pfile", "/data/test",
            "--out", "/output/test",
            "--all-variant",
        ])
        assert args.variant_qc.run_geno is True
        assert args.variant_qc.run_case_control is True
        assert args.variant_qc.run_haplotype is True
        assert args.variant_qc.run_hwe is True
        assert args.variant_qc.run_ld is False  # LD not included

    def test_sex_with_values(self) -> None:
        """--sex with values parses correctly."""
        args = parse_args([
            "--pfile", "/data/test",
            "--out", "/output/test",
            "--sex", "0.2", "0.8",
        ])
        assert args.sample_qc.run_sex is True
        assert args.sample_qc.female_max_f == 0.2
        assert args.sample_qc.male_min_f == 0.8

    def test_ld_with_values(self) -> None:
        """--ld with values parses correctly."""
        args = parse_args([
            "--pfile", "/data/test",
            "--out", "/output/test",
            "--ld", "100", "10", "0.2",
        ])
        assert args.variant_qc.run_ld is True
        assert args.variant_qc.ld_window_size == 100
        assert args.variant_qc.ld_step_size == 10
        assert args.variant_qc.ld_r2_threshold == 0.2

    def test_ancestry_flag(self) -> None:
        """--ancestry enables ancestry prediction."""
        args = parse_args([
            "--pfile", "/data/test",
            "--out", "/output/test",
            "--ancestry",
        ])
        assert args.ancestry.run_ancestry is True

    def test_container_flag(self) -> None:
        """--container sets container mode."""
        args = parse_args([
            "--pfile", "/data/test",
            "--out", "/output/test",
            "--ancestry",
            "--container",
        ])
        assert args.ancestry.use_container is True
        assert args.ancestry.inference_mode == InferenceMode.CONTAINER

    def test_warn_default(self) -> None:
        """warn_only defaults to True."""
        args = parse_args([
            "--pfile", "/data/test",
            "--out", "/output/test",
        ])
        assert args.warn_only is True

    def test_no_warn_flag(self) -> None:
        """--no-warn disables warn mode."""
        args = parse_args([
            "--pfile", "/data/test",
            "--out", "/output/test",
            "--no-warn",
        ])
        assert args.warn_only is False

    def test_pca_flag(self) -> None:
        """--pca enables PCA with default components."""
        args = parse_args([
            "--pfile", "/data/test",
            "--out", "/output/test",
            "--pca",
        ])
        assert args.gwas.run_pca is True
        assert args.gwas.n_pcs == 10

    def test_pca_with_value(self) -> None:
        """--pca with value sets custom components."""
        args = parse_args([
            "--pfile", "/data/test",
            "--out", "/output/test",
            "--pca", "20",
        ])
        assert args.gwas.run_pca is True
        assert args.gwas.n_pcs == 20
