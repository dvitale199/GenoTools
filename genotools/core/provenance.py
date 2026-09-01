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

"""What software produced a result.

The previous round recorded *what settings* a run used. This one records *what
software*, on the two paths where the answer was previously unrecoverable:

- **A trained ancestry model** (:func:`package_versions`, :func:`version_drift`).
  ``metadata.json`` held hyperparameters but not the library versions behind the
  fit. A model fit under ``umap-learn==0.5.3`` and loaded under a newer umap
  unpickles *successfully* and embeds differently, so the run yields different
  ancestry calls with no error and no warning. Recording the versions at fit
  time makes that drift visible at load time.

- **The external executables** (:func:`describe_tool`, :func:`tool_lines`).
  ``dependencies.__check_package`` resolves plink/plink2/KING from
  ``$GENOTOOLS_DEP_DIR`` or ``~/.genotools/misc/executables`` and never consults
  ``PATH`` — good for reproducibility, but nothing recorded *which* build
  produced a result, and a dev box commonly has a different plink2 on ``PATH``
  that the pipeline quietly ignores.

Both halves are best-effort: provenance must never be the reason a run fails.
Every lookup here degrades to a string rather than raising.
"""

import os
import pathlib
import platform
import shutil
import subprocess
from dataclasses import dataclass
from importlib import metadata
from typing import Dict, List, Optional, Sequence, Tuple

from .. import __version__

__all__ = [
    "MODEL_PACKAGES",
    "ToolInfo",
    "describe_tool",
    "package_versions",
    "requirements_lines",
    "tool_lines",
    "version_drift",
]


# The libraries that shape a fitted ancestry model. umap-learn is the one that
# actually motivated this: it is the only pinned dependency in setup.py, and
# unpinning it is tracked separately. scipy is in the list because umap and
# scikit-learn both compute through it.
MODEL_PACKAGES: Tuple[str, ...] = (
    "umap-learn",
    "scikit-learn",
    "xgboost",
    "numpy",
    "pandas",
    "scipy",
)

# What a package resolves to when it is not installed at all. Recorded rather
# than omitted so that "was absent then, present now" is still visible drift.
NOT_INSTALLED = "not installed"

_UNKNOWN = "unknown"


def package_versions(
    packages: Sequence[str] = MODEL_PACKAGES,
) -> Dict[str, str]:
    """Return installed versions of ``packages``, plus genotools and python.

    Distribution names are used as given (``umap-learn``, not ``umap``), since
    that is what :mod:`importlib.metadata` and ``pip`` both speak.

    Args:
        packages: Distribution names to look up.

    Returns:
        Ordered mapping of name to version string. A package that is not
        installed maps to :data:`NOT_INSTALLED`; the mapping always carries
        ``genotools`` and ``python`` keys.
    """
    versions: Dict[str, str] = {}
    for name in packages:
        try:
            versions[name] = metadata.version(name)
        except Exception:
            # PackageNotFoundError, but also a broken dist-info: provenance is
            # never worth failing a fit over.
            versions[name] = NOT_INSTALLED
    versions["genotools"] = __version__
    versions["python"] = platform.python_version()
    return versions


def version_drift(
    recorded: Dict[str, str],
    current: Optional[Dict[str, str]] = None,
) -> List[Tuple[str, str, str]]:
    """Compare recorded versions against the current environment.

    Only packages present in *both* mappings are compared, so an explicitly
    passed ``current`` that omits a key treats it as uncomparable rather than
    as drift. The default ``current`` looks up exactly the keys ``recorded``
    holds, so a package that was recorded and has since been uninstalled does
    show up, as drift to :data:`NOT_INSTALLED` — that is a real change to the
    environment, not a change in what gets recorded.

    Args:
        recorded: Versions stored when the model was fitted.
        current: Versions to compare against. Defaults to this environment's,
            looked up for exactly the keys ``recorded`` holds (``genotools``
            and ``python`` included, from the interpreter rather than from
            distribution metadata).

    Returns:
        ``(package, recorded_version, current_version)`` for each package whose
        version differs, in the order ``recorded`` lists them.
    """
    if current is None:
        keys = [k for k in recorded if k not in ("genotools", "python")]
        current = package_versions(keys)
    return [
        (name, was, current[name])
        for name, was in recorded.items()
        if name in current and current[name] != was
    ]


def requirements_lines(recorded: Dict[str, str]) -> List[str]:
    """Render recorded versions as installable ``pip`` requirement specifiers.

    Detecting drift is not the same as being able to undo it. A warning that
    says "you are on umap-learn 0.5.12, this model was fitted on 0.5.3" leaves
    the reader to reconstruct an environment by hand; this turns the same
    record into something they can feed to ``pip install -r``.

    ``python`` is emitted as a comment rather than a requirement — pip cannot
    install an interpreter, and silently dropping it would hide a version
    difference that does affect unpickling.

    Args:
        recorded: A model's ``versions`` mapping.

    Returns:
        Lines suitable for a requirements file, interpreter first as a comment.
        Packages recorded as :data:`NOT_INSTALLED` are emitted as comments too,
        since "absent" is not a version pip can resolve.
    """
    lines: List[str] = []
    if "python" in recorded:
        lines.append(f"# python {recorded['python']}")
    for name, version in recorded.items():
        if name == "python":
            continue
        if version == NOT_INSTALLED:
            lines.append(f"# {name} was not installed when this model was fitted")
            continue
        lines.append(f"{name}=={version}")

    # umap-learn below 0.5.5 does `import pkg_resources` at import time, and
    # setuptools deleted pkg_resources in 82.0.0. Reinstalling the recorded
    # umap into a current environment therefore produces a model that cannot
    # be imported at all -- a reproducibility file that does not reproduce.
    # Pin the one thing that makes the old pin work again.
    if _older_than(recorded.get("umap-learn"), (0, 5, 5)):
        lines.append("setuptools<82  # umap-learn<0.5.5 needs pkg_resources")
    return lines


def _older_than(version: Optional[str], floor: Tuple[int, ...]) -> bool:
    """Whether ``version`` parses to something below ``floor``.

    Deliberately forgiving: an unparseable or absent version is not treated as
    old, because guessing wrong here would add a bogus pin to a requirements
    file rather than merely omit a helpful one.
    """
    if not version or version == NOT_INSTALLED:
        return False
    parts: List[int] = []
    for chunk in version.split(".")[: len(floor)]:
        digits = ""
        for char in chunk:
            if not char.isdigit():
                break
            digits += char
        if not digits:
            return False
        parts.append(int(digits))
    return tuple(parts) < floor


# ---------------------------------------------------------------------------
# External executables
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class ToolInfo:
    """Where an external tool came from and what it says it is.

    Attributes:
        name: Tool name as GenoTools refers to it ("plink2").
        path: Resolved executable, or None if not present in the executable
            folder. None means "not fetched yet", not "unavailable" —
            GenoTools downloads a pinned build on first use.
        version: First line of the tool's ``--version`` output, or a short
            explanation of why there isn't one.
        shadowed_by: A *different* executable of the same name found on
            ``PATH``. GenoTools does not use it; it is recorded because a dev
            box commonly has one and it explains a surprising result.
    """

    name: str
    path: Optional[pathlib.Path]
    version: str
    shadowed_by: Optional[pathlib.Path] = None


def executable_folder() -> pathlib.Path:
    """Return the folder GenoTools resolves its executables from.

    Mirrors ``dependencies.__get_executable_folder``. Kept here as a read-only
    probe so provenance can report a tool that has not been fetched without
    triggering the download that resolving it would.
    """
    override = os.environ.get("GENOTOOLS_DEP_DIR")
    if override:
        return pathlib.Path(os.path.abspath(override))
    return pathlib.Path.home() / ".genotools" / "misc" / "executables"


def tool_version(path: pathlib.Path) -> str:
    """Return the first line of ``path --version``, or why it isn't available."""
    try:
        result = subprocess.run(
            [str(path), "--version"],
            capture_output=True,
            text=True,
            timeout=30,
        )
    except (OSError, subprocess.SubprocessError):
        return _UNKNOWN
    # PLINK 1.9 and KING both answer on stdout; fall back to stderr for
    # anything that does not.
    output = result.stdout.strip() or result.stderr.strip()
    if not output:
        return _UNKNOWN
    return output.splitlines()[0].strip()


def describe_tool(name: str) -> ToolInfo:
    """Describe one external tool without fetching it.

    Args:
        name: Executable name as it lives in the executable folder
            ("plink", "plink2", "king").

    Returns:
        A :class:`ToolInfo`. ``path`` is None when the tool has not been
        downloaded yet, in which case GenoTools will fetch a pinned build the
        first time a step needs it.
    """
    binary = executable_folder() / name
    on_path = shutil.which(name)
    shadowed_by: Optional[pathlib.Path] = None
    if on_path:
        candidate = pathlib.Path(on_path)
        try:
            different = not binary.exists() or not candidate.samefile(binary)
        except OSError:  # pragma: no cover - racing filesystem
            different = True
        if different:
            shadowed_by = candidate

    if not binary.exists():
        return ToolInfo(
            name=name,
            path=None,
            version="not fetched yet; a pinned build is downloaded on first use",
            shadowed_by=shadowed_by,
        )
    return ToolInfo(
        name=name,
        path=binary,
        version=tool_version(binary),
        shadowed_by=shadowed_by,
    )


def tool_lines(names: Sequence[str] = ("plink2", "plink", "king")) -> List[str]:
    """Render one provenance line per tool, for the consolidated run log.

    A tool GenoTools ignores in favour of its own copy gets a second line
    naming the one on ``PATH``, because "I upgraded plink2 and nothing changed"
    is otherwise a genuinely hard thing to work out.
    """
    lines: List[str] = []
    for name in names:
        info = describe_tool(name)
        where = str(info.path) if info.path is not None else "not resolved"
        lines.append(f"{info.name}: {where} -- {info.version}")
        if info.shadowed_by is not None:
            lines.append(
                f"    note: a different {info.name} is on PATH "
                f"({info.shadowed_by}); GenoTools does not use it"
            )
    return lines
