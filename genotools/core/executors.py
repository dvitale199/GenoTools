"""External tool execution wrapper.

This module provides safe, consistent execution of external tools
(PLINK, PLINK2, KING) with:
- Lazy initialization (no imports trigger downloads)
- Proper error handling and reporting
- Platform-aware tool availability
- Context managers for temporary file cleanup

Usage:
    from genotools.core.executors import get_plink2, run_command, run_plink2

    # Get executable path (lazy-loaded)
    plink2 = get_plink2()

    # Run arbitrary command with error handling
    result = run_command(["plink2", "--version"], tool_name="plink2")

    # Run PLINK2 with automatic --make-pgen and psam-cols preservation
    result = run_plink2(
        input_path=Path("/data/input"),
        output_path=Path("/data/output"),
        extra_args=["--mind", "0.02"]
    )
"""

import platform
import shlex
import subprocess
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path
from typing import Generator, List, Optional, Sequence, Union

from .exceptions import DependencyError, ExternalToolError
from .logging import get_logger, raw_sink

logger = get_logger(__name__)


def _harvest_raw_log(cmd_list: Sequence[str], cmd_str: str) -> None:
    """Harvest an external tool's native ``.log`` into the run-scoped raw sink.

    Every PLINK/PLINK2 invocation writes ``{--out}.log`` and KING writes
    ``{--prefix}.log``. When a :class:`RunLog` is registered as the raw sink
    (during a pipeline run), read that file and hand it to the sink so it lands
    in the consolidated log and the per-step raw file. Best-effort: harvesting
    must never break execution, so all failures are swallowed. A no-op when no
    sink is registered (library/unit callers).
    """
    try:
        sink = raw_sink.get()
        if sink is None:
            return

        prefix: Optional[str] = None
        args = list(cmd_list)
        for flag in ("--out", "--prefix"):
            if flag in args:
                idx = args.index(flag)
                if idx + 1 < len(args):
                    prefix = args[idx + 1]
                break
        if prefix is None:
            return

        log_path = Path(f"{prefix}.log")
        if not log_path.is_file():
            return

        text = log_path.read_text(errors="replace")
        sink.append_raw(cmd_str, text)
    except Exception:  # pragma: no cover - harvesting is strictly best-effort
        pass

# Lazy-loaded executable paths
_plink: Optional[Path] = None
_plink2: Optional[Path] = None
_king: Optional[Path] = None

# CRITICAL: Always use this flag to preserve sample metadata
PLINK2_PSAM_COLS = "psam-cols=fid,parents,sex,pheno1,phenos"


@dataclass(frozen=True)
class CommandResult:
    """Result of running an external command.

    Attributes:
        returncode: Exit code from the process
        stdout: Standard output (decoded as UTF-8)
        stderr: Standard error (decoded as UTF-8)
        command: The command that was executed (as string)
        log_file: Path to log file if one was created
    """

    returncode: int
    stdout: str
    stderr: str
    command: str
    log_file: Optional[Path] = None


def get_plink() -> Path:
    """Get the path to PLINK 1.9 executable.

    Lazily initializes - only downloads/finds on first call.

    Returns:
        Path to the plink executable.

    Raises:
        DependencyError: If PLINK cannot be found or downloaded.
    """
    global _plink
    if _plink is None:
        try:
            from genotools.dependencies import check_plink
            plink_path = check_plink()  # type: ignore[no-untyped-call]
            if plink_path is None:
                raise DependencyError("PLINK 1.9 could not be found or downloaded")
            _plink = Path(plink_path)
        except Exception as e:
            raise DependencyError(f"Failed to initialize PLINK 1.9: {e}") from e

    return _plink


def get_plink2() -> Path:
    """Get the path to PLINK2 executable.

    Lazily initializes - only downloads/finds on first call.

    Returns:
        Path to the plink2 executable.

    Raises:
        DependencyError: If PLINK2 cannot be found or downloaded.
    """
    global _plink2
    if _plink2 is None:
        try:
            from genotools.dependencies import check_plink2
            plink2_path = check_plink2()  # type: ignore[no-untyped-call]
            if plink2_path is None:
                raise DependencyError("PLINK2 could not be found or downloaded")
            _plink2 = Path(plink2_path)
        except Exception as e:
            raise DependencyError(f"Failed to initialize PLINK2: {e}") from e

    return _plink2


def get_king() -> Optional[Path]:
    """Get the path to KING executable.

    KING is only available on Linux. Returns None on other platforms.

    Returns:
        Path to the king executable, or None if not available.

    Raises:
        DependencyError: If on Linux but KING cannot be found/downloaded.
    """
    global _king

    # KING only works on Linux
    if platform.system() != "Linux":
        logger.debug("KING is only available on Linux")
        return None

    if _king is None:
        try:
            from genotools.dependencies import check_king
            king_path = check_king()  # type: ignore[no-untyped-call]
            if king_path is None:
                logger.warning("KING could not be found or downloaded")
                return None
            _king = Path(king_path)
        except Exception as e:
            raise DependencyError(f"Failed to initialize KING: {e}") from e

    return _king


def run_command(
    command: Union[str, Sequence[str]],
    tool_name: str = "command",
    check: bool = True,
    capture_output: bool = True,
    timeout: Optional[float] = None,
) -> CommandResult:
    """Run an external command with proper error handling.

    Args:
        command: Command to run. Can be a string (will be split with shlex)
                or a sequence of arguments.
        tool_name: Name of the tool for error reporting.
        check: If True, raise ExternalToolError on non-zero exit code.
        capture_output: If True, capture stdout and stderr.
        timeout: Optional timeout in seconds.

    Returns:
        CommandResult with returncode, stdout, stderr, and command string.

    Raises:
        ExternalToolError: If check=True and command returns non-zero.

    Example:
        result = run_command(["plink2", "--version"], tool_name="plink2")
        print(result.stdout)
    """
    # Convert string command to list
    if isinstance(command, str):
        cmd_list = shlex.split(command)
    else:
        cmd_list = list(command)

    cmd_str = " ".join(cmd_list)
    logger.debug(f"Running {tool_name}: {cmd_str}")

    try:
        result = subprocess.run(
            cmd_list,
            capture_output=capture_output,
            text=True,
            timeout=timeout,
        )
    except subprocess.TimeoutExpired as e:
        raise ExternalToolError(
            tool=tool_name,
            command=cmd_str,
            returncode=-1,
            stderr=f"Command timed out after {timeout} seconds",
        ) from e
    except FileNotFoundError as e:
        raise ExternalToolError(
            tool=tool_name,
            command=cmd_str,
            returncode=-1,
            stderr=f"Executable not found: {cmd_list[0]}",
        ) from e

    cmd_result = CommandResult(
        returncode=result.returncode,
        stdout=result.stdout if capture_output else "",
        stderr=result.stderr if capture_output else "",
        command=cmd_str,
    )

    # Harvest the tool's native .log into the run-scoped raw sink (best-effort).
    # Done before the failure raise so a failing step's PLINK log is captured.
    _harvest_raw_log(cmd_list, cmd_str)

    if check and result.returncode != 0:
        raise ExternalToolError(
            tool=tool_name,
            command=cmd_str,
            returncode=result.returncode,
            stderr=result.stderr,
            stdout=result.stdout,
        )

    return cmd_result


def run_plink2(
    input_path: Path,
    output_path: Path,
    input_format: str = "pfile",
    extra_args: Optional[List[str]] = None,
    make_pgen: bool = True,
) -> CommandResult:
    """Run PLINK2 with standard options.

    This is a convenience wrapper that:
    - Automatically uses --pfile or --bfile based on input_format
    - Includes --make-pgen with psam-cols preservation by default
    - Handles common argument patterns

    Args:
        input_path: Input genotype file path (without extension).
        output_path: Output path (without extension).
        input_format: Either "pfile" or "bfile".
        extra_args: Additional PLINK2 arguments.
        make_pgen: If True, include --make-pgen with psam-cols.
                   Set to False for operations that produce other output.

    Returns:
        CommandResult from the PLINK2 execution.

    Raises:
        ExternalToolError: If PLINK2 fails.

    Example:
        result = run_plink2(
            input_path=Path("/data/input"),
            output_path=Path("/data/output"),
            extra_args=["--mind", "0.02", "--geno", "0.02"]
        )
    """
    plink2 = get_plink2()

    # Build command
    input_flag = "--pfile" if input_format == "pfile" else "--bfile"

    cmd: List[str] = [str(plink2), input_flag, str(input_path)]

    # Add extra arguments
    if extra_args:
        cmd.extend(extra_args)

    # Add make-pgen with psam-cols preservation
    if make_pgen:
        cmd.extend(["--make-pgen", PLINK2_PSAM_COLS])

    # Add output
    cmd.extend(["--out", str(output_path)])

    return run_command(cmd, tool_name="plink2")


def run_plink(
    input_path: Path,
    output_path: Path,
    input_format: str = "bfile",
    extra_args: Optional[List[str]] = None,
    make_bed: bool = True,
) -> CommandResult:
    """Run PLINK 1.9 with standard options.

    Args:
        input_path: Input genotype file path (without extension).
        output_path: Output path (without extension).
        input_format: Either "pfile" or "bfile".
        extra_args: Additional PLINK arguments.
        make_bed: If True, include --make-bed.

    Returns:
        CommandResult from the PLINK execution.

    Raises:
        ExternalToolError: If PLINK fails.
    """
    plink = get_plink()

    # Build command
    input_flag = "--pfile" if input_format == "pfile" else "--bfile"

    cmd: List[str] = [str(plink), input_flag, str(input_path)]

    # Add extra arguments
    if extra_args:
        cmd.extend(extra_args)

    # Add make-bed
    if make_bed:
        cmd.append("--make-bed")

    # Add output
    cmd.extend(["--out", str(output_path)])

    return run_command(cmd, tool_name="plink")


def run_king(
    input_path: Path,
    output_path: Path,
    extra_args: Optional[List[str]] = None,
) -> CommandResult:
    """Run KING relatedness analysis.

    Note: KING is only available on Linux.

    Args:
        input_path: Input genotype file path (without extension).
                   Must be in bfile format.
        output_path: Output path prefix.
        extra_args: Additional KING arguments.

    Returns:
        CommandResult from the KING execution.

    Raises:
        DependencyError: If KING is not available (non-Linux platform).
        ExternalToolError: If KING fails.
    """
    king = get_king()
    if king is None:
        raise DependencyError(
            "KING is only available on Linux. "
            "This operation cannot be performed on the current platform."
        )

    cmd: List[str] = [str(king), "-b", f"{input_path}.bed"]

    if extra_args:
        cmd.extend(extra_args)

    cmd.extend(["--prefix", str(output_path)])

    return run_command(cmd, tool_name="king")


@contextmanager
def temp_plink_files(*prefixes: Path) -> Generator[None, None, None]:
    """Context manager for cleaning up temporary PLINK files.

    Automatically removes all PLINK-related files for the given prefixes
    when the context exits, even if an exception occurred.

    Args:
        *prefixes: One or more path prefixes to clean up.

    Example:
        with temp_plink_files(Path("/tmp/temp1"), Path("/tmp/temp2")):
            run_plink2(input_path, Path("/tmp/temp1"), ...)
            run_plink2(Path("/tmp/temp1"), Path("/tmp/temp2"), ...)
        # All temp1.* and temp2.* files are automatically cleaned up
    """
    # All possible PLINK file extensions
    extensions = [
        # PLINK2 pfile
        ".pgen", ".pvar", ".psam",
        # PLINK1.9 bfile
        ".bed", ".bim", ".fam",
        # Common output files
        ".log", ".nosex", ".hh",
        # QC output files
        ".sexcheck", ".het", ".missing", ".smiss", ".vmiss",
        ".kin", ".kin0", ".grm", ".grm.id", ".grm.N.bin",
        # Prune files
        ".prune.in", ".prune.out",
        # Other
        ".snplist", ".id", ".eigenvec", ".eigenval",
    ]

    try:
        yield
    finally:
        for prefix in prefixes:
            for ext in extensions:
                filepath = prefix.parent / f"{prefix.name}{ext}"
                if filepath.exists():
                    try:
                        filepath.unlink()
                        logger.debug(f"Cleaned up: {filepath}")
                    except OSError as e:
                        logger.warning(f"Failed to clean up {filepath}: {e}")
