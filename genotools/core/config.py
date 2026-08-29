"""Base configuration classes for GenoTools.

This module provides base classes and utilities for configuration
management across the codebase. Domain-specific configurations
(QC, ancestry, GWAS) inherit from these base classes.

Usage:
    from genotools.core.config import BaseConfig

    @dataclass(frozen=True)
    class MyStepConfig(BaseConfig):
        threshold: float = 0.05
        min_samples: int = 10
"""

from dataclasses import dataclass, fields, asdict
from pathlib import Path
from typing import Any, Dict, Optional, TypeVar, Type
import json

T = TypeVar("T", bound="BaseConfig")


@dataclass(frozen=True)
class BaseConfig:
    """Base class for all configuration dataclasses.

    Provides common functionality for configuration objects:
    - Immutability (frozen=True)
    - Serialization to/from dict and JSON
    - Validation hooks
    - Copy with modifications

    Subclasses should also be frozen dataclasses with typed fields
    and sensible defaults.
    """

    def to_dict(self) -> Dict[str, Any]:
        """Convert configuration to dictionary.

        Returns:
            Dictionary representation of the config.
        """
        return asdict(self)

    def to_json(self, indent: int = 2) -> str:
        """Convert configuration to JSON string.

        Args:
            indent: Indentation level for pretty printing.

        Returns:
            JSON string representation.
        """
        return json.dumps(self.to_dict(), indent=indent, default=str)

    @classmethod
    def from_dict(cls: Type[T], data: Dict[str, Any]) -> T:
        """Create configuration from dictionary.

        Only uses keys that match field names, ignoring extras.

        Args:
            data: Dictionary with configuration values.

        Returns:
            New configuration instance.
        """
        field_names = {f.name for f in fields(cls)}
        filtered_data = {k: v for k, v in data.items() if k in field_names}
        return cls(**filtered_data)

    @classmethod
    def from_json(cls: Type[T], json_str: str) -> T:
        """Create configuration from JSON string.

        Args:
            json_str: JSON string with configuration values.

        Returns:
            New configuration instance.
        """
        data = json.loads(json_str)
        return cls.from_dict(data)

    def replace(self: T, **changes: Any) -> T:
        """Create a copy with some fields replaced.

        Args:
            **changes: Field names and new values.

        Returns:
            New configuration instance with changes applied.

        Example:
            new_config = config.replace(threshold=0.1)
        """
        from dataclasses import replace as dc_replace
        return dc_replace(self, **changes)

    def validate(self) -> None:
        """Validate configuration values.

        Override in subclasses to add validation logic.
        Called automatically by __post_init__ if defined.

        Raises:
            ValueError: If configuration is invalid.
        """
        pass


@dataclass(frozen=True)
class PathConfig(BaseConfig):
    """Configuration for file paths.

    Provides a standardized way to configure input/output paths
    with optional working directory for intermediate files.
    """

    input_path: Optional[Path] = None
    output_path: Optional[Path] = None
    work_dir: Optional[Path] = None

    def __post_init__(self) -> None:
        """Convert string paths to Path objects."""
        # Since frozen, we need to use object.__setattr__
        if isinstance(self.input_path, str):
            object.__setattr__(self, "input_path", Path(self.input_path))
        if isinstance(self.output_path, str):
            object.__setattr__(self, "output_path", Path(self.output_path))
        if isinstance(self.work_dir, str):
            object.__setattr__(self, "work_dir", Path(self.work_dir))


@dataclass(frozen=True)
class ThresholdConfig(BaseConfig):
    """Base class for configurations with numeric thresholds.

    Provides validation for threshold values that must be
    within a specific range (typically 0-1 for proportions).
    """

    def _validate_proportion(
        self,
        value: float,
        name: str,
        min_val: float = 0.0,
        max_val: float = 1.0
    ) -> None:
        """Validate a proportion is within bounds.

        Args:
            value: The value to validate.
            name: Parameter name for error message.
            min_val: Minimum allowed value (inclusive).
            max_val: Maximum allowed value (inclusive).

        Raises:
            ValueError: If value is out of bounds.
        """
        if not min_val <= value <= max_val:
            raise ValueError(
                f"{name} must be between {min_val} and {max_val}, got {value}"
            )

    def _validate_positive(self, value: float, name: str) -> None:
        """Validate a value is positive.

        Args:
            value: The value to validate.
            name: Parameter name for error message.

        Raises:
            ValueError: If value is not positive.
        """
        if value <= 0:
            raise ValueError(f"{name} must be positive, got {value}")

    def _validate_non_negative(self, value: float, name: str) -> None:
        """Validate a value is non-negative.

        Args:
            value: The value to validate.
            name: Parameter name for error message.

        Raises:
            ValueError: If value is negative.
        """
        if value < 0:
            raise ValueError(f"{name} must be non-negative, got {value}")


@dataclass(frozen=True)
class PipelineConfig(BaseConfig):
    """Base configuration for pipeline execution.

    Controls overall pipeline behavior like error handling
    and output verbosity.
    """

    warn_only: bool = False
    """If True, continue on step failure instead of stopping."""

    verbose: bool = False
    """If True, enable verbose output."""

    full_output: bool = False
    """If True, preserve all intermediate files."""

    def should_continue_on_error(self) -> bool:
        """Check if pipeline should continue after an error.

        Returns:
            True if warn_only is enabled.
        """
        return self.warn_only
