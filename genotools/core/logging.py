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
from typing import Optional, Union


# Context variable for tracking the current QC step
# This allows log messages to automatically include which step is running
current_step: ContextVar[str] = ContextVar("current_step", default="")


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
