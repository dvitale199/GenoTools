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

"""What AncestryModel.load accepts, and what it says when it refuses.

``--model`` is exactly where a GenoTools 1.x model gets pointed at 2.0, and
1.x pickled an ``sklearn.pipeline.Pipeline`` rather than an ``AncestryModel``.
The refusal has to name that, or the user is left reading a type mismatch.
"""

import pickle
from pathlib import Path

import pytest

from genotools.ancestry.model import AncestryModel
from genotools.core.exceptions import AncestryError


def _pickle(obj: object, path: Path) -> Path:
    with open(path, "wb") as f:
        pickle.dump(obj, f)
    return path


def test_legacy_1x_model_is_refused_with_a_next_step(tmp_path: Path) -> None:
    """A 1.x model is an sklearn Pipeline. Say so, and say what to do."""
    from sklearn.pipeline import Pipeline
    from sklearn.preprocessing import StandardScaler

    legacy = _pickle(
        Pipeline([("scale", StandardScaler())]), tmp_path / "pipeline_1x.pkl"
    )

    with pytest.raises(AncestryError) as excinfo:
        AncestryModel.load(legacy)

    message = str(excinfo.value)
    assert "1.x" in message, "the user cannot act on a bare type mismatch"
    assert "--ref-panel" in message, "name the way out, not just the problem"
    # Still says what it actually found, for anyone debugging a stranger file.
    assert "sklearn.pipeline" in message


def test_unrelated_pickle_is_refused_without_the_1x_hint(tmp_path: Path) -> None:
    """The hint must not fire on any old wrong file - that would send someone
    chasing a 1.x model they never had.
    """
    other = _pickle({"not": "a model"}, tmp_path / "random.pkl")

    with pytest.raises(AncestryError) as excinfo:
        AncestryModel.load(other)

    message = str(excinfo.value)
    assert "expected AncestryModel" in message
    assert "1.x" not in message


def test_missing_model_path_raises_file_not_found(tmp_path: Path) -> None:
    with pytest.raises(FileNotFoundError):
        AncestryModel.load(tmp_path / "nope.pkl")


def test_directory_without_pipeline_pkl_says_which_file_is_missing(
    tmp_path: Path,
) -> None:
    empty = tmp_path / "model_dir"
    empty.mkdir()

    with pytest.raises(FileNotFoundError, match="pipeline.pkl"):
        AncestryModel.load(empty)
