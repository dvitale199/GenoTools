"""Unit tests for genotools.core module."""

import logging
import tempfile
from pathlib import Path

import pytest

from genotools.core import (
    GenotypeData,
    GenoToolsError,
    QCError,
    AncestryError,
    ExternalToolError,
    FileFormatError,
    ValidationError,
    DependencyError,
    setup_logging,
    get_logger,
    step_context,
    current_step,
    BaseConfig,
    ThresholdConfig,
)
from genotools.core.logging import RunLog, install_run_logging, raw_sink


class TestExceptions:
    """Tests for exception hierarchy."""

    def test_genotools_error_is_base(self):
        """All custom exceptions inherit from GenoToolsError."""
        assert issubclass(QCError, GenoToolsError)
        assert issubclass(AncestryError, GenoToolsError)
        assert issubclass(ExternalToolError, GenoToolsError)
        assert issubclass(FileFormatError, GenoToolsError)
        assert issubclass(ValidationError, GenoToolsError)
        assert issubclass(DependencyError, GenoToolsError)

    def test_external_tool_error_attributes(self):
        """ExternalToolError stores command details."""
        err = ExternalToolError(
            tool="plink2",
            command="plink2 --version",
            returncode=1,
            stderr="error message",
            stdout="output"
        )
        assert err.tool == "plink2"
        assert err.command == "plink2 --version"
        assert err.returncode == 1
        assert err.stderr == "error message"
        assert err.stdout == "output"
        assert "plink2" in str(err)
        assert "exit code 1" in str(err)

    def test_external_tool_error_truncates_long_stderr(self):
        """Long stderr is truncated in error message."""
        long_stderr = "x" * 1000
        err = ExternalToolError(
            tool="plink2",
            command="cmd",
            returncode=1,
            stderr=long_stderr
        )
        # Full stderr is preserved in attribute
        assert len(err.stderr) == 1000
        # But message is truncated
        assert len(str(err)) < 1000


class TestLogging:
    """Tests for logging utilities."""

    def test_step_context_sets_and_clears(self):
        """step_context sets current_step and clears on exit."""
        assert current_step.get() == ""

        with step_context("test_step"):
            assert current_step.get() == "test_step"

        assert current_step.get() == ""

    def test_step_context_clears_on_exception(self):
        """step_context clears even if exception occurs."""
        try:
            with step_context("failing_step"):
                assert current_step.get() == "failing_step"
                raise ValueError("test error")
        except ValueError:
            pass

        assert current_step.get() == ""

    def test_get_logger_returns_logger(self):
        """get_logger returns a logger instance."""
        logger = get_logger(__name__)
        assert logger is not None
        assert "genotools" in logger.name


def _info_record(msg: str, step: str = "") -> logging.LogRecord:
    """Build an INFO LogRecord with an optional bracketed step tag."""
    rec = logging.LogRecord(
        name="genotools.test", level=logging.INFO, pathname=__file__,
        lineno=1, msg=msg, args=(), exc_info=None,
    )
    rec.step = f"[{step}]" if step else ""
    return rec


class TestRunLog:
    """Tests for the consolidated-log RunLog writer (round 7)."""

    def test_banner_header_structured_then_raw_order(self, tmp_path: Path):
        """Within a section: header, then structured lines, then buffered raw."""
        log_path = tmp_path / "out_all_logs.log"
        rl = RunLog(log_path)
        rl.write_banner()
        rl.begin_section("callrate_prune")
        rl.write_record(_info_record("Filtering samples (mind=0.02)", "callrate_prune"))
        rl.write_record(_info_record("filtering complete: 5 removed", "callrate_prune"))
        rl.append_raw("plink2 --geno --out x", "PLINK log line: 5 samples removed.")
        rl.end_section(None)
        rl.close()

        content = log_path.read_text()
        # Banner present.
        assert "GENOTOOLS" in content or "██" in content
        # Section header present.
        i_header = content.index("===== callrate_prune =====")
        # Structured lines present and after header.
        i_struct = content.index("filtering complete: 5 removed")
        # Raw present and AFTER the structured summary.
        i_raw = content.index("PLINK log line: 5 samples removed.")
        assert i_header < i_struct < i_raw
        # Structured line carries the step tag.
        assert "[callrate_prune]" in content

    def test_end_section_writes_per_step_raw_file(self, tmp_path: Path):
        """end_section(raw_path) writes the buffered raw to a per-step file too."""
        log_path = tmp_path / "out_all_logs.log"
        raw_path = tmp_path / "out_callrate.log"
        rl = RunLog(log_path)
        rl.begin_section("callrate_prune")
        rl.append_raw("plink2 --out x", "raw plink content here")
        rl.end_section(raw_path)
        rl.close()

        assert raw_path.exists()
        assert "raw plink content here" in raw_path.read_text()
        # Consolidated also contains it.
        assert "raw plink content here" in log_path.read_text()

    def test_append_raw_buffers_until_end_section(self, tmp_path: Path):
        """append_raw does not hit disk until end_section flushes it."""
        log_path = tmp_path / "out_all_logs.log"
        rl = RunLog(log_path)
        rl.begin_section("geno_prune")
        rl.append_raw("plink2 --out x", "DEFERRED_RAW_MARKER")
        # Not flushed yet.
        assert "DEFERRED_RAW_MARKER" not in log_path.read_text()
        rl.end_section(None)
        assert "DEFERRED_RAW_MARKER" in log_path.read_text()
        rl.close()

    def test_write_summary_appends_block(self, tmp_path: Path):
        """write_summary appends a recognizable summary block at the tail."""
        log_path = tmp_path / "out_all_logs.log"
        rl = RunLog(log_path)
        rl.write_banner()
        rl.write_summary([("callrate_prune", 5, True), ("geno_prune", 112, True)])
        rl.close()
        content = log_path.read_text()
        assert "summary" in content.lower()
        assert "callrate_prune" in content
        assert "112" in content

    def test_close_is_idempotent(self, tmp_path: Path):
        """Double close does not raise."""
        rl = RunLog(tmp_path / "out_all_logs.log")
        rl.close()
        rl.close()

    def test_install_run_logging_sets_sink_and_banner(self, tmp_path: Path):
        """install_run_logging writes the banner and registers the raw sink."""
        out = tmp_path / "out"
        rl = install_run_logging(str(out), console=False)
        try:
            assert raw_sink.get() is rl
            assert Path(f"{out}_all_logs.log").exists()
            # A genotools logger record now flows into the consolidated log.
            get_logger("genotools.test").info("hello from a step")
            for h in logging.getLogger("genotools").handlers:
                h.flush()
            assert "hello from a step" in Path(f"{out}_all_logs.log").read_text()
        finally:
            rl.close()
            raw_sink.set(None)


class TestGenotypeData:
    """Tests for GenotypeData class."""

    @pytest.fixture
    def test_data_path(self) -> Path:
        """Path to test genotype data."""
        path = Path("tests/data/synthetic/genotools_test")
        if not path.with_suffix(".pgen").exists():
            pytest.skip("Test data not found")
        return path

    def test_from_path_detects_pfile(self, test_data_path: Path):
        """from_path correctly detects pfile format."""
        data = GenotypeData.from_path(test_data_path)
        assert data.format == "pfile"
        assert data.path == test_data_path

    def test_from_path_counts_samples_variants(self, test_data_path: Path):
        """from_path correctly counts samples and variants."""
        data = GenotypeData.from_path(test_data_path)
        assert data.sample_count == 500
        assert data.variant_count == 40500

    def test_from_path_accepts_string(self, test_data_path: Path):
        """from_path accepts string paths."""
        data = GenotypeData.from_path(str(test_data_path))
        assert data.path == test_data_path

    def test_from_path_raises_on_missing_files(self):
        """from_path raises FileFormatError for missing files."""
        with pytest.raises(FileFormatError):
            GenotypeData.from_path(Path("/nonexistent/path"))

    def test_is_immutable(self, test_data_path: Path):
        """GenotypeData is immutable (frozen)."""
        data = GenotypeData.from_path(test_data_path)
        with pytest.raises(AttributeError):
            data.sample_count = 999  # type: ignore

    def test_pfile_properties(self, test_data_path: Path):
        """pfile properties return correct paths."""
        data = GenotypeData.from_path(test_data_path)
        assert data.pgen == test_data_path.with_suffix(".pgen")
        assert data.pvar == test_data_path.with_suffix(".pvar")
        assert data.psam == test_data_path.with_suffix(".psam")
        assert data.bed is None
        assert data.bim is None
        assert data.fam is None

    def test_to_bfile_creates_files(self, test_data_path: Path):
        """to_bfile creates bed/bim/fam files."""
        data = GenotypeData.from_path(test_data_path)

        with tempfile.TemporaryDirectory() as tmpdir:
            out_path = Path(tmpdir) / "converted"
            bfile_data = data.to_bfile(out_path)

            assert bfile_data.format == "bfile"
            assert bfile_data.bed is not None
            assert bfile_data.bed.exists()
            assert bfile_data.bim is not None
            assert bfile_data.bim.exists()
            assert bfile_data.fam is not None
            assert bfile_data.fam.exists()

    def test_to_bfile_preserves_counts(self, test_data_path: Path):
        """to_bfile preserves sample and variant counts."""
        data = GenotypeData.from_path(test_data_path)

        with tempfile.TemporaryDirectory() as tmpdir:
            out_path = Path(tmpdir) / "converted"
            bfile_data = data.to_bfile(out_path)

            assert bfile_data.sample_count == data.sample_count
            assert bfile_data.variant_count == data.variant_count

    def test_ensure_pfile_returns_self_if_already_pfile(self, test_data_path: Path):
        """ensure_pfile returns self if already pfile format."""
        data = GenotypeData.from_path(test_data_path)

        with tempfile.TemporaryDirectory() as tmpdir:
            result = data.ensure_pfile(Path(tmpdir))
            assert result is data  # Same object

    def test_validate_required_columns_passes(self, test_data_path: Path):
        """validate_required_columns passes for valid data."""
        data = GenotypeData.from_path(test_data_path)
        # Should not raise
        data.validate_required_columns()

    def test_str_representation(self, test_data_path: Path):
        """__str__ returns informative representation."""
        data = GenotypeData.from_path(test_data_path)
        s = str(data)
        assert "GenotypeData" in s
        assert "pfile" in s
        assert "500" in s  # sample count

    def test_from_vcf_creates_pfile(self, test_data_path: Path, tmp_path: Path):
        """from_vcf converts a VCF to pfile and cleans up intermediates."""
        from genotools.core import GenotypeData
        from genotools.core.executors import get_plink2
        # Make a VCF from the synthetic pfile
        vcf_prefix = tmp_path / "from_vcf_input"
        import subprocess
        subprocess.run(
            [str(get_plink2()), "--pfile", str(test_data_path),
             "--autosome", "--export", "vcf", "--out", str(vcf_prefix)],
            check=True, capture_output=True,
        )
        vcf_file = vcf_prefix.with_suffix(".vcf")
        assert vcf_file.exists()

        data = GenotypeData.from_vcf(vcf_file)

        assert data.format == "pfile"
        assert (vcf_prefix.with_suffix(".pgen")).exists()
        # intermediate bfile removed
        assert not (vcf_prefix.with_suffix(".bed")).exists()

    def test_from_vcf_raises_on_missing(self, tmp_path: Path):
        from genotools.core import GenotypeData
        with pytest.raises(FileNotFoundError):
            GenotypeData.from_vcf(tmp_path / "nope.vcf")

    def test_to_pfile_from_bfile_in_place(self, test_data_path: Path, tmp_path: Path):
        """bfile -> pfile at the same prefix (the runner's bfile input branch)."""
        from genotools.core import GenotypeData
        bfile_prefix = tmp_path / "bf_input"
        src = GenotypeData.from_path(test_data_path)
        src.to_bfile(bfile_prefix)
        assert (bfile_prefix.with_suffix(".bed")).exists()

        bf = GenotypeData.from_path(bfile_prefix)
        assert bf.format == "bfile"
        result = bf.to_pfile(bfile_prefix)
        assert result.format == "pfile"
        assert (bfile_prefix.with_suffix(".pgen")).exists()
        assert result.sample_count == src.sample_count
        assert result.variant_count == src.variant_count


class TestConfig:
    """Tests for configuration classes."""

    def test_base_config_to_dict(self):
        """BaseConfig.to_dict serializes correctly."""
        from dataclasses import dataclass

        @dataclass(frozen=True)
        class TestConfig(BaseConfig):
            value: int = 10
            name: str = "test"

        config = TestConfig()
        d = config.to_dict()
        assert d == {"value": 10, "name": "test"}

    def test_base_config_from_dict(self):
        """BaseConfig.from_dict deserializes correctly."""
        from dataclasses import dataclass

        @dataclass(frozen=True)
        class TestConfig(BaseConfig):
            value: int = 10
            name: str = "test"

        config = TestConfig.from_dict({"value": 20, "name": "custom"})
        assert config.value == 20
        assert config.name == "custom"

    def test_base_config_from_dict_ignores_extras(self):
        """from_dict ignores unknown keys."""
        from dataclasses import dataclass

        @dataclass(frozen=True)
        class TestConfig(BaseConfig):
            value: int = 10

        config = TestConfig.from_dict({"value": 20, "unknown": "ignored"})
        assert config.value == 20
        assert not hasattr(config, "unknown")

    def test_base_config_replace(self):
        """replace creates modified copy."""
        from dataclasses import dataclass

        @dataclass(frozen=True)
        class TestConfig(BaseConfig):
            value: int = 10
            name: str = "test"

        config = TestConfig()
        modified = config.replace(value=99)
        assert modified.value == 99
        assert modified.name == "test"
        assert config.value == 10  # Original unchanged

    def test_threshold_config_validation(self):
        """ThresholdConfig validation methods work."""
        from dataclasses import dataclass

        @dataclass(frozen=True)
        class TestConfig(ThresholdConfig):
            proportion: float = 0.5

            def validate(self) -> None:
                self._validate_proportion(self.proportion, "proportion")

        # Valid config
        config = TestConfig(proportion=0.5)
        config.validate()  # Should not raise

        # Invalid config
        config_invalid = TestConfig(proportion=1.5)
        with pytest.raises(ValueError, match="proportion must be between"):
            config_invalid.validate()


class TestBanner:
    """Tests for the ASCII banner."""

    def test_banner_returns_ascii_header(self):
        from genotools.core.logging import banner
        b = banner()
        assert isinstance(b, str)
        assert "█" in b            # box-drawing art present
        assert b.count("\n") >= 5  # multi-line

    def test_banner_matches_legacy(self):
        # Faithful copy of legacy gt_header (differential; drop when utils.py is removed)
        from genotools.core.logging import banner
        from genotools.utils import gt_header
        assert banner() == gt_header()


class TestImportPerformance:
    """Tests for import performance."""

    def test_import_time_under_threshold(self):
        """Core module imports in under 100ms."""
        import time
        import importlib
        import sys

        # Remove cached imports
        modules_to_remove = [k for k in sys.modules if k.startswith("genotools.core")]
        for mod in modules_to_remove:
            del sys.modules[mod]

        start = time.time()
        importlib.import_module("genotools.core")
        elapsed = (time.time() - start) * 1000

        assert elapsed < 100, f"Import took {elapsed:.2f}ms, expected < 100ms"

    def test_no_side_effects_on_import(self):
        """Importing does not trigger downloads or other side effects."""
        import sys
        import io
        import importlib

        # Remove cached imports
        modules_to_remove = [k for k in sys.modules if k.startswith("genotools.core")]
        for mod in modules_to_remove:
            del sys.modules[mod]

        # Capture stderr
        old_stderr = sys.stderr
        sys.stderr = io.StringIO()

        try:
            importlib.import_module("genotools.core")
            captured = sys.stderr.getvalue()
        finally:
            sys.stderr = old_stderr

        assert captured == "", f"Import produced unexpected output: {captured}"


class TestValidation:
    """Tests for core.validation.validate_input."""

    @pytest.fixture
    def geno(self) -> Path:
        return Path("tests/data/synthetic/genotools_test")

    def test_passes_on_valid_input(self, geno: Path, tmp_path: Path, capsys):
        from genotools.core.validation import validate_input
        validate_input(geno, tmp_path / "out", skip_fails=False)  # no raise
        out = capsys.readouterr().out
        assert "breakdown" in out.lower()

    def test_raises_on_missing_pgen(self, tmp_path: Path):
        from genotools.core.validation import validate_input
        with pytest.raises(FileNotFoundError):
            validate_input(tmp_path / "nope", tmp_path / "out", skip_fails=False)

    def test_raises_when_log_exists(self, geno: Path, tmp_path: Path):
        from genotools.core.validation import validate_input
        from genotools.core.exceptions import ValidationError
        log = tmp_path / "out_all_logs.log"
        log.write_text("prior run\n")
        with pytest.raises(ValidationError):
            validate_input(geno, tmp_path / "out", skip_fails=False)

    def test_skip_fails_ignores_existing_log(self, geno: Path, tmp_path: Path):
        from genotools.core.validation import validate_input
        (tmp_path / "out_all_logs.log").write_text("prior\n")
        validate_input(geno, tmp_path / "out", skip_fails=True)  # no raise

    def test_raises_on_missing_sex_column(self, geno: Path, tmp_path: Path):
        from genotools.core.validation import validate_input
        from genotools.core.exceptions import ValidationError
        import pandas as pd
        # Build a pfile-shaped fixture whose psam lacks SEX
        bad = tmp_path / "bad"
        bad.with_suffix(".pgen").write_bytes((geno.with_suffix(".pgen")).read_bytes())
        bad.with_suffix(".pvar").write_bytes((geno.with_suffix(".pvar")).read_bytes())
        psam = pd.read_csv(geno.with_suffix(".psam"), sep=r"\s+")
        psam.drop(columns=["SEX"]).to_csv(bad.with_suffix(".psam"), sep="\t", index=False)
        with pytest.raises(ValidationError):
            validate_input(bad, tmp_path / "out2", skip_fails=False)

    def test_raises_on_missing_pheno_column(self, geno: Path, tmp_path: Path):
        from genotools.core.validation import validate_input
        from genotools.core.exceptions import ValidationError
        import pandas as pd
        bad = tmp_path / "bad_pheno"
        bad.with_suffix(".pgen").write_bytes((geno.with_suffix(".pgen")).read_bytes())
        bad.with_suffix(".pvar").write_bytes((geno.with_suffix(".pvar")).read_bytes())
        psam = pd.read_csv(geno.with_suffix(".psam"), sep=r"\s+")
        psam.drop(columns=["PHENO1"]).to_csv(bad.with_suffix(".psam"), sep="\t", index=False)
        with pytest.raises(ValidationError):
            validate_input(bad, tmp_path / "outp", skip_fails=False)

    # --- data-driven step-skip decisions (ported from legacy upfront_check) ---

    def _pfile_with(self, geno: Path, dst: Path, *, psam=None, pvar_rows=None):
        """Copy geno's pgen; write a crafted psam (DataFrame) and/or truncate pvar."""
        import pandas as pd
        dst.with_suffix(".pgen").write_bytes(geno.with_suffix(".pgen").read_bytes())
        if psam is None:
            dst.with_suffix(".psam").write_bytes(geno.with_suffix(".psam").read_bytes())
        else:
            psam.to_csv(dst.with_suffix(".psam"), sep="\t", index=False)
        src_pvar = geno.with_suffix(".pvar").read_text().splitlines()
        if pvar_rows is None:
            dst.with_suffix(".pvar").write_text("\n".join(src_pvar) + "\n")
        else:
            header = [ln for ln in src_pvar if ln.startswith("#")]
            body = [ln for ln in src_pvar if not ln.startswith("#")][:pvar_rows]
            dst.with_suffix(".pvar").write_text("\n".join(header + body) + "\n")

    def test_skip_sex_when_no_sex_data(self, geno: Path, tmp_path: Path):
        import pandas as pd
        from genotools.core.validation import validate_input
        psam = pd.read_csv(geno.with_suffix(".psam"), sep=r"\s+")
        psam["SEX"] = 0
        bad = tmp_path / "nosex"
        self._pfile_with(geno, bad, psam=psam)
        d = validate_input(bad, tmp_path / "o1", sex_requested=True)
        assert d.skip_sex is True

    def test_skip_sex_when_no_x_chrom(self, geno: Path, tmp_path: Path):
        from genotools.core.validation import validate_input
        # Keep only chr-1 variants (no X) by truncating to the first block.
        bad = tmp_path / "nox"
        self._pfile_with(geno, bad, pvar_rows=100)
        d = validate_input(bad, tmp_path / "o2", sex_requested=True)
        assert d.skip_sex is True

    def test_disable_filter_controls_when_no_controls(self, geno: Path, tmp_path: Path):
        import pandas as pd
        from genotools.core.validation import validate_input
        psam = pd.read_csv(geno.with_suffix(".psam"), sep=r"\s+")
        psam["PHENO1"] = 2  # all cases, no controls
        bad = tmp_path / "nocontrol"
        self._pfile_with(geno, bad, psam=psam)
        d = validate_input(bad, tmp_path / "o3", hwe_requested=True, filter_controls=True)
        assert d.disable_filter_controls is True

    def test_skip_case_control_when_only_cases(self, geno: Path, tmp_path: Path):
        import pandas as pd
        from genotools.core.validation import validate_input
        psam = pd.read_csv(geno.with_suffix(".psam"), sep=r"\s+")
        psam["PHENO1"] = 2
        bad = tmp_path / "onlycases"
        self._pfile_with(geno, bad, psam=psam)
        d = validate_input(bad, tmp_path / "o4", case_control_requested=True)
        assert d.skip_case_control is True

    def test_skip_het_when_few_variants(self, geno: Path, tmp_path: Path):
        from genotools.core.validation import validate_input
        bad = tmp_path / "fewvar"
        self._pfile_with(geno, bad, pvar_rows=40)  # < 50 variants
        d = validate_input(bad, tmp_path / "o5", het_requested=True)
        assert d.skip_het is True

    def test_no_decisions_when_skip_fails(self, geno: Path, tmp_path: Path):
        import pandas as pd
        from genotools.core.validation import validate_input
        psam = pd.read_csv(geno.with_suffix(".psam"), sep=r"\s+")
        psam["SEX"] = 0
        bad = tmp_path / "sf"
        self._pfile_with(geno, bad, psam=psam)
        d = validate_input(bad, tmp_path / "o6", skip_fails=True, sex_requested=True)
        assert d.skip_sex is False

    def test_no_decisions_when_not_requested(self, geno: Path, tmp_path: Path):
        from genotools.core.validation import validate_input
        d = validate_input(geno, tmp_path / "o7")  # nothing requested
        assert (d.skip_sex, d.skip_case_control, d.skip_het, d.disable_filter_controls) \
            == (False, False, False, False)
