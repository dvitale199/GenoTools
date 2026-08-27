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

"""Shared fixtures for the unit suite."""

import logging

import pytest


@pytest.fixture(autouse=True)
def restore_genotools_logger():
    """Undo any logging configuration a test leaves on the "genotools" logger.

    setup_logging() and install_run_logging() deliberately set
    ``propagate = False`` so library output does not escape to the root logger.
    A test that calls either one and then only clears handlers leaves
    propagation off for the rest of the session, and pytest's ``caplog``
    captures via a root handler - so every later test asserting on a warning
    silently sees nothing and fails for a reason that has nothing to do with
    the code under test.

    Snapshot and restore the logger's state around every unit test so the order
    tests run in cannot change their outcome.
    """
    logger = logging.getLogger("genotools")
    handlers = list(logger.handlers)
    propagate = logger.propagate
    level = logger.level
    try:
        yield
    finally:
        logger.handlers[:] = handlers
        logger.propagate = propagate
        logger.setLevel(level)
