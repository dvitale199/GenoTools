"""Unit tests for the executor raw-log harvest (round 7)."""

from pathlib import Path

import pytest

from genotools.core.executors import _harvest_raw_log
from genotools.core.logging import raw_sink


class _FakeSink:
    """Captures append_raw calls in place of a RunLog."""

    def __init__(self) -> None:
        self.calls: list[tuple[str, str]] = []

    def append_raw(self, command: str, text: str) -> None:
        self.calls.append((command, text))


@pytest.fixture
def sink():
    s = _FakeSink()
    token = raw_sink.set(s)  # type: ignore[arg-type]
    try:
        yield s
    finally:
        raw_sink.reset(token)


def test_harvest_reads_out_log(tmp_path: Path, sink: _FakeSink):
    """--out <prefix> => reads <prefix>.log and hands text to the sink."""
    prefix = tmp_path / "step_output"
    Path(f"{prefix}.log").write_text("PLINK2 v2.0\n112 variants removed.\n")

    cmd = ["plink2", "--pfile", "in", "--geno", "0.02", "--out", str(prefix)]
    _harvest_raw_log(cmd, " ".join(cmd))

    assert len(sink.calls) == 1
    command, text = sink.calls[0]
    assert "112 variants removed." in text
    assert "--out" in command


def test_harvest_reads_prefix_log_for_king(tmp_path: Path, sink: _FakeSink):
    """KING uses --prefix; harvest reads <prefix>.log too."""
    prefix = tmp_path / "king_out"
    Path(f"{prefix}.log").write_text("KING relatedness output\n")

    cmd = ["king", "-b", "x.bed", "--prefix", str(prefix)]
    _harvest_raw_log(cmd, " ".join(cmd))

    assert len(sink.calls) == 1
    assert "KING relatedness output" in sink.calls[0][1]


def test_harvest_noop_when_no_log_file(tmp_path: Path, sink: _FakeSink):
    """No .log on disk => nothing handed to the sink, no error."""
    cmd = ["plink2", "--out", str(tmp_path / "missing")]
    _harvest_raw_log(cmd, " ".join(cmd))
    assert sink.calls == []


def test_harvest_noop_when_no_out_flag(tmp_path: Path, sink: _FakeSink):
    """No --out/--prefix (e.g. --version) => no-op."""
    _harvest_raw_log(["plink2", "--version"], "plink2 --version")
    assert sink.calls == []


def test_harvest_noop_without_sink(tmp_path: Path):
    """Unset sink => harvesting is a silent no-op (no raise)."""
    prefix = tmp_path / "s"
    Path(f"{prefix}.log").write_text("some log")
    # raw_sink defaults to None here (no fixture).
    _harvest_raw_log(["plink2", "--out", str(prefix)], "cmd")  # must not raise


def test_harvest_never_raises_on_bad_read(tmp_path: Path, sink: _FakeSink):
    """A directory at the .log path (unreadable) must not propagate."""
    prefix = tmp_path / "d"
    Path(f"{prefix}.log").mkdir()  # reading this as a file fails
    _harvest_raw_log(["plink2", "--out", str(prefix)], "cmd")  # must not raise
    assert sink.calls == []
