"""Core infrastructure for GenoTools.

This module provides the foundational components used across GenoTools:

- **GenotypeData**: Immutable representation of genotype files
- **Exceptions**: Structured exception hierarchy
- **Logging**: Context-aware structured logging
- **Executors**: Safe external tool execution (PLINK, KING)
- **Config**: Base configuration classes

Usage:
    from genotools.core import GenotypeData, setup_logging, get_logger
    from genotools.core import GenoToolsError, QCError, ExternalToolError

    # Setup logging
    setup_logging(level="INFO")
    logger = get_logger(__name__)

    # Load genotype data
    data = GenotypeData.from_path("/path/to/genotypes")
    print(f"Loaded {data.sample_count} samples")
"""

# Genotype data representation
from .genotypes import GenotypeData

# Exception hierarchy
from .exceptions import (
    GenoToolsError,
    ValidationError,
    QCError,
    AncestryError,
    GWASError,
    ExternalToolError,
    FileFormatError,
    DependencyError,
)

# Logging utilities
from .logging import (
    setup_logging,
    get_logger,
    set_step,
    clear_step,
    step_context,
    current_step,
    banner,
)

# External tool execution
from .executors import (
    get_plink,
    get_plink2,
    get_king,
    run_command,
    run_plink,
    run_plink2,
    run_king,
    temp_plink_files,
    CommandResult,
    PLINK2_PSAM_COLS,
)

# Configuration base classes
from .config import (
    BaseConfig,
    PathConfig,
    ThresholdConfig,
    PipelineConfig,
)

__all__ = [
    # Genotype data
    "GenotypeData",
    # Exceptions
    "GenoToolsError",
    "ValidationError",
    "QCError",
    "AncestryError",
    "GWASError",
    "ExternalToolError",
    "FileFormatError",
    "DependencyError",
    # Logging
    "setup_logging",
    "get_logger",
    "set_step",
    "clear_step",
    "step_context",
    "current_step",
    "banner",
    # Executors
    "get_plink",
    "get_plink2",
    "get_king",
    "run_command",
    "run_plink",
    "run_plink2",
    "run_king",
    "temp_plink_files",
    "CommandResult",
    "PLINK2_PSAM_COLS",
    # Config
    "BaseConfig",
    "PathConfig",
    "ThresholdConfig",
    "PipelineConfig",
]
