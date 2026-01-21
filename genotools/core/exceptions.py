"""Custom exception hierarchy for GenoTools.

This module provides a structured exception hierarchy that enables:
- Consistent error handling across the codebase
- Informative error messages with context
- Easy distinction between QC, ancestry, and external tool errors
"""

from typing import Optional


class GenoToolsError(Exception):
    """Base exception for all GenoTools errors.

    All custom exceptions in GenoTools inherit from this class,
    allowing callers to catch all GenoTools-specific errors with
    a single except clause if desired.
    """
    pass


class ValidationError(GenoToolsError):
    """Raised when input validation fails.

    Examples:
        - Invalid parameter values (e.g., threshold out of range)
        - Missing required files
        - Malformed input data
    """
    pass


class QCError(GenoToolsError):
    """QC-specific errors.

    Raised when a QC step fails in a way that cannot be recovered,
    such as all samples being removed or a required operation failing.
    """
    pass


class AncestryError(GenoToolsError):
    """Ancestry prediction errors.

    Raised when ancestry prediction fails, such as:
        - Model not fitted before prediction
        - Insufficient SNP overlap with reference panel
        - Invalid reference panel format
    """
    pass


class ExternalToolError(GenoToolsError):
    """External tool (PLINK, PLINK2, KING) execution failures.

    Provides detailed context about the failed command including:
        - Tool name
        - Full command that was executed
        - Return code
        - stderr output (truncated for readability)

    Attributes:
        tool: Name of the external tool (e.g., "plink2", "king")
        command: The full command that was executed
        returncode: The exit code from the process
        stderr: Standard error output from the process
        stdout: Standard output from the process (optional)
    """

    def __init__(
        self,
        tool: str,
        command: str,
        returncode: int,
        stderr: str,
        stdout: Optional[str] = None
    ):
        self.tool = tool
        self.command = command
        self.returncode = returncode
        self.stderr = stderr
        self.stdout = stdout

        # Truncate stderr for the message but keep full version in attribute
        stderr_preview = stderr[:500] if len(stderr) > 500 else stderr
        super().__init__(
            f"{tool} failed with exit code {returncode}:\n"
            f"Command: {command}\n"
            f"Error: {stderr_preview}"
        )


class FileFormatError(GenoToolsError):
    """Raised when file format is invalid or unsupported.

    Examples:
        - Unrecognized genotype file format
        - Corrupted pfile/bfile
        - Missing required columns in psam/fam files
    """
    pass


class DependencyError(GenoToolsError):
    """Raised when a required external dependency is unavailable.

    Examples:
        - PLINK/PLINK2 not found and cannot be downloaded
        - KING requested on non-Linux platform
    """
    pass
