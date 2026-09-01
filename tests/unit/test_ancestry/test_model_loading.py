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


# ---------------------------------------------------------------------------
# Version provenance
#
# A model fitted under umap-learn 0.5.3 and loaded under a newer umap unpickles
# *successfully* and produces a different embedding: wrong ancestry calls, no
# error, no warning. Recording the versions at fit time is what makes that
# visible, so these tests assert on the warning, not on the stored dict.
# ---------------------------------------------------------------------------


def _model_with_versions(
    tmp_path: Path, versions: object, name: str = "model.pkl"
) -> Path:
    """Pickle a bare AncestryModel carrying ``versions`` (no fit required)."""
    model = AncestryModel()
    model.versions = versions  # type: ignore[assignment]
    return _pickle(model, tmp_path / name)


def _warnings(caplog: "pytest.LogCaptureFixture") -> str:
    return "\n".join(
        r.getMessage() for r in caplog.records if r.levelname == "WARNING"
    )


def test_drifted_versions_warn_and_name_the_package(
    tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    """The warning has to say which package moved and to what - "versions
    differ" leaves the user with nothing to act on."""
    path = _model_with_versions(tmp_path, {"umap-learn": "0.0.1-ancient"})

    with caplog.at_level("WARNING", logger="genotools"):
        AncestryModel.load(path)

    message = _warnings(caplog)
    assert "umap-learn" in message
    assert "0.0.1-ancient" in message, "must name the version the fit ran under"
    assert "ancestry calls" in message.lower(), "say what is actually at risk"


def test_matching_versions_load_quietly(
    tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    """A false drift warning on every load would train users to ignore it."""
    from genotools.core.provenance import package_versions

    path = _model_with_versions(tmp_path, package_versions())

    with caplog.at_level("WARNING", logger="genotools"):
        AncestryModel.load(path)

    assert _warnings(caplog) == ""


def test_a_model_without_versions_says_it_cannot_tell(
    tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    """Models predating version recording must not read as "no drift":
    "cannot tell" and "no drift" are different answers."""
    path = _model_with_versions(tmp_path, None)

    with caplog.at_level("WARNING", logger="genotools"):
        AncestryModel.load(path)

    message = _warnings(caplog)
    assert "provenance unknown" in message.lower()
    assert "retrain" in message.lower(), "name the way out"


def test_drift_never_blocks_the_load(
    tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    """Drift is often harmless, and hard-failing would strand every existing
    model the moment a floating dependency moved."""
    path = _model_with_versions(
        tmp_path, {"umap-learn": "0.0.1", "numpy": "0.0.1", "scipy": "0.0.1"}
    )

    with caplog.at_level("WARNING", logger="genotools"):
        model = AncestryModel.load(path)

    assert isinstance(model, AncestryModel)


def test_saved_metadata_carries_the_versions(tmp_path: Path) -> None:
    """metadata.json is the human-readable half; provenance has to reach it or
    nobody reads it without unpickling."""
    import json

    model = AncestryModel()
    model.versions = {"umap-learn": "0.5.3"}
    model._is_fitted = True

    out = model.save(tmp_path / "saved")
    metadata = json.loads((out / "metadata.json").read_text())

    assert metadata["versions"] == {"umap-learn": "0.5.3"}


def test_fit_is_what_records_the_versions() -> None:
    """Captured at fit, not at save: these are the versions that produced
    *this* fit, and they must travel with the pickle for the single-file
    format to carry provenance too.
    """
    import inspect

    from genotools.ancestry import model as model_module

    source = inspect.getsource(model_module.AncestryModel.fit)
    assert "package_versions()" in source, (
        "fit() no longer records versions; every model trained from here on "
        "would load with unknown provenance"
    )
