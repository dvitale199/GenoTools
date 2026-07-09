#!/usr/bin/env bash
# ==============================================================================
# Build .venv-stable containing the PRE-REFACTOR genotools, for old-vs-new
# behavioral parity testing (tests/regression/test_parity.py).
#
# The refactor branch DELETED the old qc.py/pipeline.py, so the "old" baseline
# must be installed from a pre-refactor ref (default: origin/main).
#
# Usage:
#   bash tests/scripts/setup_stable_venv.sh [PRE_REFACTOR_REF]
#
# Examples:
#   bash tests/scripts/setup_stable_venv.sh                 # uses origin/main
#   bash tests/scripts/setup_stable_venv.sh v1.3.6          # a tagged release
# ==============================================================================
set -euo pipefail

REF="${1:-origin/main}"
# Interpreter used to create .venv-stable. Override to pin a version that the
# pre-refactor deps (e.g. umap-learn==0.5.3 / numba) still support, e.g.
#   PYTHON=python3.11 bash tests/scripts/setup_stable_venv.sh
PYTHON="${PYTHON:-python3}"
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
WORKTREE="$(mktemp -d)/genotools-stable-src"

echo ">> Building .venv-stable from ref: ${REF}"
echo ">> Repo root: ${ROOT}"

# Check the ref actually predates the refactor (has the old flat qc.py).
if git -C "${ROOT}" cat-file -e "${REF}:genotools/qc.py" 2>/dev/null; then
    echo ">> Confirmed ${REF} has legacy genotools/qc.py (pre-refactor)."
else
    echo "!! WARNING: ${REF} has no legacy genotools/qc.py -- it may not be a"
    echo "!! pre-refactor baseline. Parity testing needs the OLD implementation."
fi

# Materialize the old source in a throwaway worktree and install it.
git -C "${ROOT}" worktree add --detach "${WORKTREE}" "${REF}"
trap 'git -C "${ROOT}" worktree remove --force "${WORKTREE}" 2>/dev/null || true' EXIT

"${PYTHON}" -m venv "${ROOT}/.venv-stable"
"${ROOT}/.venv-stable/bin/pip" install --quiet --upgrade pip
"${ROOT}/.venv-stable/bin/pip" install "${WORKTREE}"
# Runtime deps the old code needs that its setup.py omits / that fresh venvs on
# modern Python (3.12+) no longer bundle: pkg_resources (via setuptools, needed
# by umap-learn==0.5.3) and psutil (imported by the old ancestry.py).
"${ROOT}/.venv-stable/bin/pip" install --quiet "setuptools<81" psutil

echo ">> Done. Old baseline installed at:"
echo "   ${ROOT}/.venv-stable/bin/genotools"
echo ">> Run parity tests with:"
echo "   pytest tests/regression/test_parity.py -v"
