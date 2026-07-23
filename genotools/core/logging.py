"""Structured logging setup for GenoTools.

This module provides:
- Context-aware logging with current QC step tracking
- Configurable log levels and output destinations
- Consistent log formatting across the codebase

Usage:
    from genotools.core.logging import setup_logging, get_logger, current_step

    # Setup logging at application start
    setup_logging(level="INFO", log_file=Path("output/genotools.log"))

    # Get a logger for your module
    logger = get_logger(__name__)

    # Set current step context (automatically included in log messages)
    token = current_step.set("callrate_prune")
    logger.info("Starting callrate filtering")
    current_step.reset(token)
"""

import logging
import sys
from contextvars import ContextVar, Token
from pathlib import Path
from typing import List, Optional, Sequence, Tuple, Union


# Context variable for tracking the current QC step
# This allows log messages to automatically include which step is running
current_step: ContextVar[str] = ContextVar("current_step", default="")

# Full detailed format used for the consolidated on-disk log.
FILE_FORMAT = "%(asctime)s [%(levelname)s] %(step)s %(name)s: %(message)s"
DATE_FORMAT = "%Y-%m-%d %H:%M:%S"


class StepContextFilter(logging.Filter):
    """Logging filter that adds the current step to log records.

    This filter reads from the current_step context variable and adds
    it to each log record, allowing the formatter to include it.
    """

    def filter(self, record: logging.LogRecord) -> bool:
        """Add step context to the log record."""
        step = current_step.get()
        record.step = f"[{step}]" if step else ""
        return True


def setup_logging(
    level: Union[str, int] = "INFO",
    log_file: Optional[Path] = None,
    console: bool = True
) -> None:
    """Configure structured logging with optional file output.

    Args:
        level: Log level (DEBUG, INFO, WARNING, ERROR, CRITICAL) or int
        log_file: Optional path to write logs to file
        console: Whether to output to console (default True)

    Example:
        setup_logging(level="DEBUG", log_file=Path("run.log"))
    """
    # Convert string level to int if needed
    if isinstance(level, str):
        level = getattr(logging, level.upper(), logging.INFO)

    # Get the root genotools logger
    logger = logging.getLogger("genotools")
    logger.setLevel(level)

    # Clear any existing handlers, closing them first so repeated
    # setup_logging() calls in one process (e.g. successive pipeline runs)
    # don't leak file descriptors.
    for handler in logger.handlers[:]:
        handler.close()
    logger.handlers.clear()

    # Create formatter with step context placeholder
    formatter = logging.Formatter(
        "%(asctime)s [%(levelname)s] %(step)s %(name)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S"
    )

    # Add context filter to include step in all messages
    context_filter = StepContextFilter()

    if console:
        console_handler = logging.StreamHandler(sys.stderr)
        console_handler.setLevel(level)
        console_handler.setFormatter(formatter)
        console_handler.addFilter(context_filter)
        logger.addHandler(console_handler)

    if log_file is not None:
        # Ensure parent directory exists
        log_file.parent.mkdir(parents=True, exist_ok=True)

        file_handler = logging.FileHandler(log_file)
        file_handler.setLevel(level)
        file_handler.setFormatter(formatter)
        file_handler.addFilter(context_filter)
        logger.addHandler(file_handler)

    # Prevent propagation to root logger to avoid duplicate messages
    logger.propagate = False


def get_logger(name: str) -> logging.Logger:
    """Get a logger for a module.

    Args:
        name: Module name, typically __name__

    Returns:
        Configured logger instance

    Example:
        logger = get_logger(__name__)
        logger.info("Processing started")
    """
    # Ensure the logger is under the genotools namespace
    if not name.startswith("genotools"):
        name = f"genotools.{name}"
    return logging.getLogger(name)


def set_step(step_name: str) -> Token[str]:
    """Set the current step context.

    Args:
        step_name: Name of the current QC/processing step

    Returns:
        Token that can be used to reset the context

    Example:
        token = set_step("callrate_prune")
        try:
            # Do work - all logs will include [callrate_prune]
            logger.info("Filtering samples")
        finally:
            current_step.reset(token)
    """
    return current_step.set(step_name)


def clear_step(token: Token[str]) -> None:
    """Clear the step context using the token from set_step.

    Args:
        token: Token returned from set_step()
    """
    current_step.reset(token)


class step_context:
    """Context manager for setting the current step.

    Usage:
        with step_context("callrate_prune"):
            logger.info("This message includes [callrate_prune]")
        # Context automatically cleared after the block
    """

    def __init__(self, step_name: str):
        self.step_name = step_name
        self._token: Optional[Token[str]] = None

    def __enter__(self) -> "step_context":
        self._token = current_step.set(self.step_name)
        return self

    def __exit__(
        self,
        exc_type: Optional[type],
        exc_val: Optional[BaseException],
        exc_tb: Optional[object]
    ) -> None:
        if self._token is not None:
            current_step.reset(self._token)


def banner() -> str:
    """Return the GenoTools ASCII banner (for the consolidated log header)."""
    return """
     ██████╗ ███████╗███╗  ██╗ █████╗ ████████╗ █████╗  █████╗ ██╗      ██████╗
    ██╔════╝ ██╔════╝████╗ ██║██╔══██╗╚══██╔══╝██╔══██╗██╔══██╗██║     ██╔════╝
    ██║  ██╗ █████╗  ██╔██╗██║██║  ██║   ██║   ██║  ██║██║  ██║██║     ╚█████╗
    ██║  ╚██╗██╔══╝  ██║╚████║██║  ██║   ██║   ██║  ██║██║  ██║██║      ╚═══██╗
    ╚██████╔╝███████╗██║ ╚███║╚█████╔╝   ██║   ╚█████╔╝╚█████╔╝███████╗██████╔╝
    ╚═════╝ ╚══════╝╚═╝  ╚══╝ ╚════╝    ╚═╝    ╚════╝  ╚════╝ ╚══════╝╚═════╝
    """


# ---------------------------------------------------------------------------
# Consolidated run log (round 7: logging / pipeline-visibility redesign)
# ---------------------------------------------------------------------------

# A summary row: (step_label, outliers_removed, passed).
SummaryRow = Tuple[str, int, bool]

# Run-scoped raw-output sink. The executor reads this after every PLINK/KING
# invocation and, if set, hands the harvested native .log to the RunLog. Unset
# (None) => harvesting is a silent no-op, so library/unit callers are unaffected.
raw_sink: ContextVar[Optional["RunLog"]] = ContextVar("raw_sink", default=None)


class RunLog:
    """Single writer that owns the consolidated ``{out}_all_logs.log`` file.

    One object owns the file handle so structured records, section headers, and
    buffered raw PLINK blocks land in a deterministic order. Per step the layout
    is: a ``===== step =====`` header, the structured summary lines (written live
    via :meth:`write_record`), then the verbatim harvested PLINK output for that
    step (buffered via :meth:`append_raw`, flushed by :meth:`end_section`). The
    same buffered raw is also written to a per-step ``{out}_{step}.log`` file.
    """

    def __init__(self, path: Union[str, Path]) -> None:
        self.path = Path(path)
        self.path.parent.mkdir(parents=True, exist_ok=True)
        # Truncate: each run starts a fresh consolidated log. The runner's
        # re-run guard already prevents clobbering a prior run's output.
        self._fh: Optional[object] = open(self.path, "w", encoding="utf-8")
        self._formatter = logging.Formatter(FILE_FORMAT, datefmt=DATE_FORMAT)
        self._raw_buffer: List[str] = []

    # -- header / sections --------------------------------------------------

    def write_banner(self) -> None:
        self._write(banner() + "\n")

    def begin_section(self, title: str) -> None:
        """Open a section: flush any stray raw, write the header, reset buffer."""
        self._raw_buffer = []
        self._write(f"\n===== {title} =====\n")

    def write_record(self, record: logging.LogRecord) -> None:
        """Write a structured log record live (full detailed format)."""
        if not hasattr(record, "step"):
            record.step = ""  # type: ignore[attr-defined]
        self._write(self._formatter.format(record) + "\n")

    def append_raw(self, command: str, text: str) -> None:
        """Buffer a harvested raw PLINK block; flushed at :meth:`end_section`."""
        block = f"$ {command}\n{text.rstrip()}\n" if command else f"{text.rstrip()}\n"
        self._raw_buffer.append(block)

    def end_section(self, raw_log_path: Optional[Union[str, Path]] = None) -> None:
        """Flush buffered raw into the consolidated log and per-step file."""
        if self._raw_buffer:
            raw = "".join(self._raw_buffer)
            self._write(raw)
            if raw_log_path is not None:
                try:
                    Path(raw_log_path).write_text(raw, encoding="utf-8")
                except OSError:
                    pass
        self._raw_buffer = []

    # -- summary ------------------------------------------------------------

    def write_summary(self, rows: Sequence[SummaryRow]) -> None:
        """Append the end-of-run summary table."""
        lines = ["\n===== run summary =====\n"]
        for label, removed, passed in rows:
            status = "PASS" if passed else "FAIL"
            lines.append(f"  {label:<28} {removed:>10} removed   {status}\n")
        self._write("".join(lines))

    # -- lifecycle ----------------------------------------------------------

    def _write(self, text: str) -> None:
        if self._fh is not None:
            self._fh.write(text)
            self._fh.flush()

    def close(self) -> None:
        """Close the file handle. Idempotent."""
        if self._fh is not None:
            try:
                self._fh.close()
            finally:
                self._fh = None


class _RunLogHandler(logging.Handler):
    """Routes ``genotools.*`` structured records into a :class:`RunLog`."""

    def __init__(self, runlog: RunLog, level: int = logging.INFO) -> None:
        super().__init__(level)
        self._runlog = runlog

    def emit(self, record: logging.LogRecord) -> None:
        try:
            self._runlog.write_record(record)
        except Exception:  # pragma: no cover - never let logging crash a run
            self.handleError(record)


class _ConsoleFormatter(logging.Formatter):
    """Concise, step-tagged console format: ``[step] message`` (``! `` on WARN+)."""

    def format(self, record: logging.LogRecord) -> str:
        if not hasattr(record, "step"):
            record.step = ""  # type: ignore[attr-defined]
        prefix = "! " if record.levelno >= logging.WARNING else ""
        tag = f"{record.step} " if record.step else ""
        return f"{prefix}{tag}{record.getMessage()}"


def install_run_logging(
    out_path: Union[str, Path],
    level: Union[str, int] = "INFO",
    console: bool = True,
) -> RunLog:
    """Configure runtime logging around a :class:`RunLog` and return it.

    Builds the consolidated ``{out}_all_logs.log`` (banner written), attaches a
    :class:`_RunLogHandler` (+ an optional concise console handler) to the
    ``genotools`` logger, and registers the RunLog as the run-scoped raw sink so
    the executor harvests PLINK ``.log`` output into it.

    Args:
        out_path: Output prefix; the log is ``{out_path}_all_logs.log``.
        level: Log level for both handlers.
        console: Whether to also emit the curated console stream.

    Returns:
        The RunLog (the caller owns its lifecycle; call ``close()`` when done).
    """
    if isinstance(level, str):
        level = getattr(logging, level.upper(), logging.INFO)

    runlog = RunLog(f"{out_path}_all_logs.log")
    runlog.write_banner()

    logger = logging.getLogger("genotools")
    logger.setLevel(level)
    for handler in logger.handlers[:]:
        handler.close()
    logger.handlers.clear()

    context_filter = StepContextFilter()

    run_handler = _RunLogHandler(runlog, level)
    run_handler.addFilter(context_filter)
    logger.addHandler(run_handler)

    if console:
        console_handler = logging.StreamHandler(sys.stderr)
        console_handler.setLevel(level)
        console_handler.setFormatter(_ConsoleFormatter())
        console_handler.addFilter(context_filter)
        logger.addHandler(console_handler)

    logger.propagate = False
    raw_sink.set(runlog)
    return runlog
