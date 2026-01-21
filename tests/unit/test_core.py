"""Unit tests for genotools.core module."""

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
