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

import logging
from pathlib import Path
from typing import List

import pytest

from genotools.cli.parser import (
    HetSpec,
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
    _parse_het_ancestry_args,
    _parse_het_spec,
    _parse_ld_args,
)


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
        assert args.het_auto is False
        assert args.het_auto_sd == 3.0
        assert args.het_by_ancestry == {}
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

    def test_het_config_for_no_label_uses_base(self) -> None:
        """A flat run has no label, so only the base can apply."""
        args = SampleQCArgs(het_lower=-0.2, het_upper=0.2)
        config = args.het_config_for(None)
        assert config.auto_detect is False
        assert (config.f_lower, config.f_upper) == (-0.2, 0.2)

    def test_het_config_for_adaptive_base(self) -> None:
        args = SampleQCArgs(het_auto=True, het_auto_sd=2.0)
        config = args.het_config_for(None)
        assert config.auto_detect is True
        assert config.auto_sd == 2.0

    def test_adaptive_base_avoids_the_legacy_sentinel(self) -> None:
        """auto_detect must not arrive as [-1, -1].

        HetConfig reads that sentinel as the legacy *rate*-based mode; the CLI's
        sd mode thresholds F. If the spec ever forced the sentinel the two would
        silently swap statistics.
        """
        config = SampleQCArgs(het_auto=True).het_config_for(None)
        assert (config.f_lower, config.f_upper) != (-1.0, -1.0)

    def test_het_config_for_unknown_label_falls_back_to_base(self) -> None:
        args = SampleQCArgs(
            het_lower=-0.2,
            het_upper=0.2,
            het_by_ancestry={"AMR": HetSpec(auto=True)},
        )
        config = args.het_config_for("EUR")
        assert config.auto_detect is False
        assert (config.f_lower, config.f_upper) == (-0.2, 0.2)

    def test_override_beats_fixed_base(self) -> None:
        args = SampleQCArgs(
            het_lower=-0.2,
            het_upper=0.2,
            het_by_ancestry={"AMR": HetSpec(auto=True, sd=3.0)},
        )
        assert args.het_config_for("AMR").auto_detect is True
        assert args.het_config_for("EUR").auto_detect is False

    def test_fixed_override_beats_adaptive_base(self) -> None:
        """The reverse direction: an adaptive base pinned back for one group."""
        args = SampleQCArgs(
            het_auto=True,
            het_auto_sd=2.0,
            het_by_ancestry={"CAH": HetSpec(lower=-0.3, upper=0.3)},
        )
        cah = args.het_config_for("CAH")
        assert cah.auto_detect is False
        assert (cah.f_lower, cah.f_upper) == (-0.3, 0.3)
        assert args.het_config_for("EUR").auto_detect is True

    def test_per_group_multiplier_differs_from_base(self) -> None:
        args = SampleQCArgs(
            het_auto=True,
            het_auto_sd=2.0,
            het_by_ancestry={"AMR": HetSpec(auto=True, sd=1.5)},
        )
        assert args.het_config_for("AMR").auto_sd == 1.5
        assert args.het_config_for("EUR").auto_sd == 2.0

    def test_to_het_config_matches_the_base(self) -> None:
        """to_het_config is het_config_for(None); they cannot drift."""
        args = SampleQCArgs(
            het_auto=True,
            het_auto_sd=2.5,
            het_by_ancestry={"AMR": HetSpec(auto=True, sd=1.0)},
        )
        assert args.to_het_config() == args.het_config_for(None)

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

    @pytest.mark.parametrize(
        "attr,flag",
        [
            ("use_container", "--container"),
            ("use_singularity", "--singularity"),
            ("use_cloud", "--cloud"),
        ],
    )
    def test_remote_inference_flags_rejected(self, attr: str, flag: str) -> None:
        """The three remote-execution flags are refused, not silently ignored.

        No execution path reads them, so accepting one would run local
        prediction while the user believed it ran in a container or on cloud.
        """
        with pytest.raises(ValueError, match=f"{flag} is not supported"):
            AncestryArgs(**{attr: True})


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
        # Verbosity: curated console on, INFO level.
        assert args.quiet is False
        assert args.debug is False

    def test_verbosity_flags_parse(self) -> None:
        """--quiet / --debug reach OutputArgs (defaults are off)."""
        default = parse_args(["--pfile", "/data/test", "--out", "/output/test"])
        assert default.output.quiet is False
        assert default.output.debug is False

        verbose = parse_args(
            ["--pfile", "/data/test", "--out", "/output/test", "--quiet", "--debug"]
        )
        assert verbose.output.quiet is True
        assert verbose.output.debug is True


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
        assert steps == ["geno", "hwe"]

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


class TestParseHetSpec:
    """Tests for _parse_het_spec, the shared grammar for both het flags.

    One resolver backs --het and --het-ancestry, so these rows are the whole
    grammar surface for both.
    """

    def test_no_tokens_is_fixed_defaults(self) -> None:
        spec = _parse_het_spec([], "--het")
        assert spec == HetSpec(auto=False, lower=-0.15, upper=0.15)

    def test_two_numbers_are_fixed_bounds(self) -> None:
        spec = _parse_het_spec(["-0.2", "0.2"], "--het")
        assert spec.auto is False
        assert (spec.lower, spec.upper) == (-0.2, 0.2)

    def test_bare_sd_uses_default_multiplier(self) -> None:
        spec = _parse_het_spec(["sd"], "--het")
        assert spec.auto is True
        assert spec.sd == 3.0

    def test_sd_with_multiplier(self) -> None:
        spec = _parse_het_spec(["sd", "2.5"], "--het")
        assert spec.auto is True
        assert spec.sd == 2.5

    def test_sd_is_case_insensitive(self) -> None:
        assert _parse_het_spec(["SD"], "--het").auto is True

    def test_legacy_sentinel_maps_to_sd(self, caplog) -> None:
        """--het -1 -1 was the 1.x way to ask for derived bounds.

        Still accepted, since old command lines carry it, but redirected to the
        keyword: a silent second spelling would be worse than a deprecation.
        """
        with caplog.at_level(logging.WARNING, logger="genotools"):
            spec = _parse_het_spec(["-1", "-1"], "--het")
        assert spec == HetSpec(auto=True, sd=3.0)
        assert "deprecated" in caplog.text
        assert "sd" in caplog.text

    def test_sd_rejects_two_multipliers(self) -> None:
        with pytest.raises(ValueError, match="at most one multiplier"):
            _parse_het_spec(["sd", "2", "3"], "--het")

    def test_sd_rejects_non_positive_multiplier(self) -> None:
        with pytest.raises(ValueError, match="greater than 0"):
            _parse_het_spec(["sd", "0"], "--het")

    def test_sd_rejects_negative_multiplier(self) -> None:
        with pytest.raises(ValueError, match="greater than 0"):
            _parse_het_spec(["sd", "-2"], "--het")

    def test_sd_rejects_unparseable_multiplier(self) -> None:
        with pytest.raises(ValueError, match="expects numbers or 'sd'"):
            _parse_het_spec(["sd", "abc"], "--het")

    def test_lone_number_is_an_arity_error(self) -> None:
        with pytest.raises(ValueError, match="requires 0 or 2 values"):
            _parse_het_spec(["0.1"], "--het")

    def test_three_numbers_is_an_arity_error(self) -> None:
        with pytest.raises(ValueError, match="requires 0 or 2 values"):
            _parse_het_spec(["0.1", "0.2", "0.3"], "--het")

    def test_unparseable_token_reports_type_not_arity(self) -> None:
        """--het abc lost argparse's type check when it became type=str.

        The helper has to supply it, and the type error has to win over the
        arity one: 'abc' wants explaining, not counting.
        """
        with pytest.raises(ValueError, match="expects numbers or 'sd'"):
            _parse_het_spec(["abc"], "--het")

    def test_unordered_bounds_rejected(self) -> None:
        with pytest.raises(ValueError, match="must be less than"):
            _parse_het_spec(["0.3", "0.1"], "--het")

    def test_error_names_the_calling_flag(self) -> None:
        """Both flags share the resolver, so errors must name the right one."""
        with pytest.raises(ValueError, match="--het-ancestry AMR"):
            _parse_het_spec(["0.1"], "--het-ancestry AMR")


class TestParseHetArgs:
    """Tests for _parse_het_args helper."""

    def test_none_input(self) -> None:
        """None input returns disabled state."""
        run, spec = _parse_het_args(None)
        assert run is False
        assert spec == HetSpec()

    def test_empty_list(self) -> None:
        """Empty list enables with defaults."""
        run, spec = _parse_het_args([])
        assert run is True
        assert (spec.auto, spec.lower, spec.upper) == (False, -0.15, 0.15)

    def test_two_values(self) -> None:
        """Two values are parsed correctly."""
        run, spec = _parse_het_args(["-0.2", "0.2"])
        assert run is True
        assert (spec.lower, spec.upper) == (-0.2, 0.2)

    def test_sd(self) -> None:
        run, spec = _parse_het_args(["sd"])
        assert run is True
        assert (spec.auto, spec.sd) == (True, 3.0)

    def test_sd_with_multiplier(self) -> None:
        run, spec = _parse_het_args(["sd", "2"])
        assert run is True
        assert (spec.auto, spec.sd) == (True, 2.0)

    def test_wrong_count_raises(self) -> None:
        """Wrong value count raises error.

        TestParseSexArgs and TestParseLdArgs both have this; --het's was
        missing.
        """
        with pytest.raises(ValueError, match="requires 0 or 2 values"):
            _parse_het_args(["0.2"])


class TestParseHetAncestryArgs:
    """Tests for _parse_het_ancestry_args helper."""

    def test_none_is_empty(self) -> None:
        assert _parse_het_ancestry_args(None) == {}

    def test_single_sd_override(self) -> None:
        specs = _parse_het_ancestry_args([["AMR", "sd"]])
        assert specs == {"AMR": HetSpec(auto=True, sd=3.0)}

    def test_sd_with_multiplier(self) -> None:
        specs = _parse_het_ancestry_args([["AMR", "sd", "1.5"]])
        assert specs["AMR"].sd == 1.5

    def test_fixed_range_override(self) -> None:
        specs = _parse_het_ancestry_args([["CAH", "-0.3", "0.3"]])
        assert specs["CAH"] == HetSpec(auto=False, lower=-0.3, upper=0.3)

    def test_repeats_accumulate(self) -> None:
        specs = _parse_het_ancestry_args([["AMR", "sd"], ["CAH", "-0.3", "0.3"]])
        assert set(specs) == {"AMR", "CAH"}
        assert specs["AMR"].auto is True
        assert specs["CAH"].auto is False

    def test_duplicate_label_rejected(self) -> None:
        with pytest.raises(ValueError, match="given more than once"):
            _parse_het_ancestry_args([["AMR", "sd"], ["AMR", "-0.3", "0.3"]])

    def test_label_without_spec_rejected(self) -> None:
        """A bare label would silently mean "fixed defaults for this group" -
        the class of quiet no-op this flag exists to remove."""
        with pytest.raises(ValueError, match="needs a spec"):
            _parse_het_ancestry_args([["AMR"]])

    def test_missing_label_rejected(self) -> None:
        """--het-ancestry sd 2 names no group; catch the transposition."""
        with pytest.raises(ValueError, match="expects an ancestry label first"):
            _parse_het_ancestry_args([["sd", "2"]])

    def test_numeric_label_rejected(self) -> None:
        with pytest.raises(ValueError, match="expects an ancestry label first"):
            _parse_het_ancestry_args([["-0.3", "0.3"]])


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

    @pytest.mark.parametrize("flag", ["--container", "--singularity", "--cloud"])
    def test_remote_inference_flags_rejected(self, flag: str) -> None:
        """Passing a remote-execution flag on the command line is an error."""
        with pytest.raises(ValueError, match=f"{flag} is not supported"):
            parse_args([
                "--pfile", "/data/test",
                "--out", "/output/test",
                "--ancestry",
                flag,
            ])

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


class TestHetGrammarEndToEnd:
    """Every row of the --het / --het-ancestry semantics table, through
    parse_args - the layer users actually reach."""

    BASE = ["--pfile", "/data/test", "--out", "/tmp/out"]

    def _sample_qc(self, extra: List[str]):
        return parse_args(self.BASE + extra).sample_qc

    def test_bare_het_is_unchanged(self) -> None:
        sq = self._sample_qc(["--het"])
        assert sq.run_het is True
        assert (sq.het_lower, sq.het_upper) == (-0.15, 0.15)
        assert sq.het_auto is False

    def test_fixed_bounds(self) -> None:
        sq = self._sample_qc(["--het", "-0.2", "0.2"])
        assert (sq.het_lower, sq.het_upper) == (-0.2, 0.2)

    def test_negative_bounds_survive_type_str(self) -> None:
        """--het is type=str now; argparse's negative-number heuristic keys off
        the parser's option strings, not the argument type, so a leading
        negative must still parse rather than read as a flag."""
        sq = self._sample_qc(["--het", "-0.5", "-0.1"])
        assert (sq.het_lower, sq.het_upper) == (-0.5, -0.1)

    def test_sd_base(self) -> None:
        sq = self._sample_qc(["--het", "sd"])
        assert sq.het_auto is True
        assert sq.het_auto_sd == 3.0

    def test_sd_base_with_multiplier(self) -> None:
        sq = self._sample_qc(["--het", "sd", "2.5"])
        assert (sq.het_auto, sq.het_auto_sd) == (True, 2.5)

    def test_amr_override_on_fixed_base(self) -> None:
        """Today's --amr-het, generalized."""
        sq = self._sample_qc(
            ["--ancestry", "--het", "-0.2", "0.2", "--het-ancestry", "AMR", "sd"]
        )
        assert sq.het_config_for("AMR").auto_detect is True
        assert sq.het_config_for("EUR").f_lower == -0.2

    def test_per_group_multiplier(self) -> None:
        sq = self._sample_qc(
            ["--ancestry", "--het", "sd", "2", "--het-ancestry", "AMR", "sd", "1.5"]
        )
        assert sq.het_config_for("AMR").auto_sd == 1.5
        assert sq.het_config_for("EUR").auto_sd == 2.0

    def test_mixed_adaptive_and_fixed_overrides(self) -> None:
        sq = self._sample_qc([
            "--ancestry",
            "--het-ancestry", "AMR", "sd",
            "--het-ancestry", "CAH", "-0.3", "0.3",
        ])
        assert sq.het_config_for("AMR").auto_detect is True
        cah = sq.het_config_for("CAH")
        assert (cah.f_lower, cah.f_upper) == (-0.3, 0.3)

    def test_override_implies_the_step_runs(self) -> None:
        """No --het, no --all-sample: the override still turns het on.

        An override implies the system it overrides is on. --amr-het did not do
        this, so --ancestry --amr-het alone was a second silent no-op.
        """
        sq = self._sample_qc(["--ancestry", "--het-ancestry", "AMR", "sd"])
        assert sq.run_het is True
        assert "het" in parse_args(
            self.BASE + ["--ancestry", "--het-ancestry", "AMR", "sd"]
        ).get_enabled_sample_steps()

    def test_all_sample_with_sd_base(self) -> None:
        """The production per-ancestry job."""
        sq = self._sample_qc(["--all-sample", "--all-variant", "--het", "sd"])
        assert sq.run_het is True
        assert sq.het_auto is True

    def test_het_ancestry_without_ancestry_errors(self) -> None:
        """The whole point of the change: a per-label override in a run with no
        labels is impossible to honour, so it must not be accepted silently."""
        with pytest.raises(ValueError, match="requires --ancestry"):
            parse_args(self.BASE + ["--het-ancestry", "AMR", "sd"])

    def test_het_ancestry_error_names_the_label_and_the_fix(self) -> None:
        with pytest.raises(ValueError, match="AMR") as excinfo:
            parse_args(self.BASE + ["--het-ancestry", "AMR", "sd"])
        assert "--het" in str(excinfo.value)

    def test_no_underscore_alias(self) -> None:
        """Underscore spellings are deprecated 1.x aliases only; a flag added
        in 2.0.1 gets just the hyphen form, like every other 2.0 addition."""
        with pytest.raises(SystemExit):
            parse_args(self.BASE + ["--ancestry", "--het_ancestry", "AMR", "sd"])

    def test_legacy_sentinel_still_parses_and_warns(self, caplog) -> None:
        with caplog.at_level(logging.WARNING, logger="genotools"):
            sq = self._sample_qc(["--het", "-1", "-1"])
        assert sq.het_auto is True
        assert "deprecated" in caplog.text

    def test_legacy_dict_carries_the_new_keys(self) -> None:
        legacy = parse_args(
            self.BASE + ["--ancestry", "--het", "sd", "2", "--het-ancestry", "AMR", "sd"]
        ).to_legacy_dict()
        assert legacy["het_auto"] is True
        assert legacy["het_auto_sd"] == 2.0
        assert set(legacy["het_by_ancestry"]) == {"AMR"}
        assert "amr_het" not in legacy


class TestAmrHetRemoved:
    """--amr-het was removed, not aliased.

    It was read in exactly one place, inside the --ancestry branch, so it did
    nothing in a flat run - which is how the per-ancestry production workflow
    runs QC. A command line still carrying it gets a targeted error naming the
    replacement, rather than argparse's bare "unrecognized arguments".
    """

    BASE = ["--pfile", "/data/test", "--out", "/tmp/out"]

    @pytest.mark.parametrize("flag", ["--amr-het", "--amr_het"])
    def test_rejected_with_a_pointer_to_the_replacement(self, flag: str) -> None:
        with pytest.raises(ValueError) as excinfo:
            parse_args(self.BASE + ["--all-sample", flag])
        message = str(excinfo.value)
        assert "removed" in message
        assert "--het-ancestry AMR sd" in message

    def test_rejected_even_with_a_trailing_value(self) -> None:
        """--amr-het False must report the removal, not the presence-flag rule.

        The removal check runs before _reject_boolean_values for this reason.
        """
        with pytest.raises(ValueError, match="removed"):
            parse_args(self.BASE + ["--amr-het", "False"])

    def test_not_in_the_spelling_rename_table(self) -> None:
        """_DEPRECATED_SPELLINGS models renames only - a removed flag there
        would claim the underscore form still works."""
        from genotools.cli.parser import _DEPRECATED_SPELLINGS

        assert "--amr_het" not in _DEPRECATED_SPELLINGS

    def test_parser_does_not_define_it(self) -> None:
        from genotools.cli.parser import create_parser

        options = {
            opt for action in create_parser()._actions for opt in action.option_strings
        }
        assert "--amr-het" not in options
        assert "--amr_het" not in options


class TestDeprecatedFlagSpellings:
    """Pre-refactor underscore flags must keep working, with a warning.

    The refactor renamed 18 flags from underscore to hyphen and dropped two
    others, so every existing user invocation failed with "unrecognized
    arguments". The old spellings are accepted as aliases so scripts keep
    running through the 2.0 transition.
    """

    BASE = ["--pfile", "/data/test", "--out", "/tmp/out"]

    def test_every_renamed_flag_is_still_accepted(self) -> None:
        """No old spelling may raise SystemExit from argparse."""
        from genotools.cli.parser import _DEPRECATED_SPELLINGS, create_parser

        parser = create_parser()
        takes_value = {
            "--ref_panel": "/tmp/ref",
            "--ref_labels": "/tmp/lab",
            "--related_cutoff": "0.0884",
            "--duplicated_cutoff": "0.354",
            "--min_samples": "50",
            "--covar_names": "AGE,SEX",
            "--subset_ancestry": "EUR",
        }
        for old in _DEPRECATED_SPELLINGS:
            argv = list(self.BASE) + [old]
            if old in takes_value:
                argv.append(takes_value[old])
            # Must not raise; argparse exits on an unrecognized flag.
            parser.parse_args(argv)

    def test_alias_sets_the_same_destination(self) -> None:
        """Old and new spellings must land on the same parsed value."""
        from genotools.cli.parser import create_parser

        parser = create_parser()
        old = parser.parse_args(self.BASE + ["--ref_panel", "/tmp/ref", "--full_output"])
        new = parser.parse_args(self.BASE + ["--ref-panel", "/tmp/ref", "--full-output"])

        assert old.ref_panel == new.ref_panel
        assert old.full_output == new.full_output is True

    def test_positive_form_is_accepted_as_a_noop(self) -> None:
        """--warn / --prune_duplicated asked for what is now the default."""
        from genotools.cli.parser import create_parser

        parser = create_parser()
        for argv in (["--warn"], ["--prune_duplicated"]):
            ns = parser.parse_args(self.BASE + argv)
            assert ns.no_warn is False
            assert ns.no_prune_duplicated is False

    @pytest.mark.parametrize(
        "argv",
        [
            ["--all_sample", "False"],
            ["--all-sample", "False"],
            ["--ancestry", "True"],
            ["--warn", "False"],
            ["--prune_duplicated", "False"],
        ],
    )
    def test_boolean_value_form_is_rejected_with_guidance(self, argv) -> None:
        """1.x accepted `--flag True/False`; 2.0 flags are presence-only.

        1.x declared every boolean as `type=str, nargs='?', const='True'`,
        copying the threshold flags (--callrate and friends) where an optional
        value is actually wanted. That is why `--all_sample False` parsed at all.
        Rejecting it is correct, but the error has to say what to do instead.
        """
        from genotools.cli.parser import _reject_boolean_values, create_parser

        parser = create_parser()
        with pytest.raises(SystemExit):
            _reject_boolean_values(parser, self.BASE + argv)

    def test_value_taking_flags_still_accept_bool_like_values(self) -> None:
        """The rejection must only fire on presence-only flags.

        --out is a path; a path that happens to be named "False" is legal.
        """
        from genotools.cli.parser import _reject_boolean_values, create_parser

        parser = create_parser()
        # Must not raise.
        _reject_boolean_values(parser, ["--pfile", "/data/x", "--out", "False"])
        _reject_boolean_values(parser, ["--pfile", "/data/x", "--out", "/tmp/o",
                                        "--callrate", "0.05"])

    def test_prune_defaults_match_1x(self) -> None:
        """2.0 defaults are the 1.x defaults: related off, duplicated on.

        1.x declared --prune_related default='False' and --prune_duplicated
        default='True' as strings, then converted 'True'/'False' to real bools,
        so under --all_sample related samples were reported but not pruned.
        """
        args = parse_args(self.BASE + ["--all_sample"])

        assert args.sample_qc.prune_related is False
        assert args.sample_qc.prune_duplicated is True

    def test_deprecated_use_is_warned(self, caplog) -> None:
        import logging

        from genotools.cli.parser import _warn_deprecated_flags

        with caplog.at_level(logging.WARNING, logger="genotools"):
            _warn_deprecated_flags(["--ref_panel", "/tmp/ref", "--warn"])

        text = " ".join(r.message for r in caplog.records)
        assert "--ref_panel is deprecated" in text
        assert "--ref-panel" in text
        assert "--warn is deprecated" in text

    def test_new_spellings_are_not_warned(self, caplog) -> None:
        import logging

        from genotools.cli.parser import _warn_deprecated_flags

        with caplog.at_level(logging.WARNING, logger="genotools"):
            _warn_deprecated_flags(["--ref-panel", "/tmp/ref", "--no-warn"])

        assert [r for r in caplog.records if "deprecated" in r.message] == []

    def test_equals_form_is_detected(self) -> None:
        """--ref_panel=/x must warn just like --ref_panel /x."""
        import logging

        from genotools.cli.parser import _warn_deprecated_flags

        logger = logging.getLogger("genotools.cli.parser")
        records: list = []
        handler = logging.Handler()
        handler.emit = records.append  # type: ignore[method-assign]
        logger.addHandler(handler)
        try:
            _warn_deprecated_flags(["--ref_panel=/tmp/ref"])
        finally:
            logger.removeHandler(handler)

        assert any("--ref_panel is deprecated" in r.getMessage() for r in records)
