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

"""Tests for the downloadable-artifact catalogue.

The catalogue serves two mutually incompatible ancestry model formats, so the
thing worth pinning is not the resolver's mechanics but the invariant behind it:
whatever `genotools-download` hands a user by default has to be loadable by the
GenoTools that shipped it. That invariant was broken before 2.1.0 -- the default
was a 1.x model that 2.x rejects.
"""

import re

import pytest

from genotools.download_refs import (
    DEFAULT_MODEL,
    DEFAULT_REF,
    MODELS,
    REF_PANELS,
    resolve_name,
)


class TestResolveName:
    """Name resolution, including the two forms of "give me the default"."""

    def test_none_resolves_to_default(self) -> None:
        assert resolve_name(None, MODELS, DEFAULT_MODEL, "model") == DEFAULT_MODEL

    def test_literal_default_resolves_to_default(self) -> None:
        """The CLI help has always advertised `--model default`.

        It used to raise KeyError('default'), since the string was looked up in
        the checksum table like any other name.
        """
        assert resolve_name("default", MODELS, DEFAULT_MODEL, "model") == DEFAULT_MODEL

    def test_explicit_name_is_returned(self) -> None:
        assert resolve_name("nba_v1", MODELS, DEFAULT_MODEL, "model") == "nba_v1"

    def test_unknown_name_lists_what_is_available(self) -> None:
        with pytest.raises(SystemExit) as excinfo:
            resolve_name("nba_v9", MODELS, DEFAULT_MODEL, "model")
        message = str(excinfo.value)
        assert "nba_v9" in message
        for name in MODELS:
            assert name in message

    def test_unknown_name_names_the_kind(self) -> None:
        """The resolver serves panels and models; the error must say which."""
        with pytest.raises(SystemExit) as excinfo:
            resolve_name("nope", REF_PANELS, DEFAULT_REF, "reference panel")
        assert "reference panel" in str(excinfo.value)


class TestCatalogue:
    """Invariants the catalogue has to hold for a release to be coherent."""

    def test_default_model_is_a_2x_model(self) -> None:
        """The default must load in the code that ships it.

        This is the regression: through 2.0.x the default was `nba_v2`, a 1.x
        pickle that `AncestryModel.load` rejects by design, so the documented
        getting-started path ended in a load error.
        """
        _, model_format, _ = MODELS[DEFAULT_MODEL]
        assert model_format == "2.x"

    def test_defaults_are_in_their_catalogues(self) -> None:
        assert DEFAULT_MODEL in MODELS
        assert DEFAULT_REF in REF_PANELS

    @pytest.mark.parametrize("name", sorted(MODELS))
    def test_model_entries_are_well_formed(self, name: str) -> None:
        checksum, model_format, description = MODELS[name]
        assert re.fullmatch(r"[0-9a-f]{32}", checksum), "md5 expected"
        assert model_format in ("1.x", "2.x")
        assert description

    @pytest.mark.parametrize("name", sorted(REF_PANELS))
    def test_ref_panel_checksums_are_md5(self, name: str) -> None:
        assert re.fullmatch(r"[0-9a-f]{32}", REF_PANELS[name])

    def test_a_2x_model_is_offered_at_all(self) -> None:
        """Serving only 1.x models would leave 2.x users with no model."""
        assert any(fmt == "2.x" for _, fmt, _ in MODELS.values())
