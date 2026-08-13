#----------------------------------------------------------------------------#
# License
#     This file is part of cardiacFoam.
#
#     cardiacFoam is free software: you can redistribute it and/or modify it
#     under the terms of the GNU General Public License as published by the
#     Free Software Foundation, either version 3 of the License, or (at your
#     option) any later version.
#
#     cardiacFoam is distributed in the hope that it will be useful, but
#     WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#     General Public License for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with cardiacFoam.  If not, see <http://www.gnu.org/licenses/>.
#
# Module
#     test_ionic_catalog_live_verification
#
# Description
#     Runs listCellModelsVariables for every catalogued ionic model and
#     asserts the catalog matches. Skipped when OpenFOAM is not sourced.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

"""The live ratchet: does IONIC_MODEL_CATALOG match the built solver?

This is the check that keeps the catalog honest. An agent writes
``ionicConstantOverrides`` from the catalog, and an unknown name is a
**FatalError** at solver startup (``src/genericWriter/ionicModelIO.C:245-255``).
Constant naming follows no static rule -- ``AC_``-prefixed for most models,
unprefixed for TNNP and BuenoOrovio, and *mixed* inside TWorld -- so only the
solver can answer.

It skips when ``listCellModelsVariables`` is not on ``PATH``, which is the
normal state on a machine without OpenFOAM sourced. To run it::

    source /Volumes/OpenFOAM-v2412/etc/bashrc
    cd applications/scripts/driverFoam
    uv run pytest openfoam_driver/tests/plugins/cardiacfoam/\\
test_ionic_catalog_live_verification.py -q -s

**A skip is not a pass.** The predecessor of this module,
``scripts/generate_catalog.py``, was never wired to anything and rotted until
it no longer imported at all. A test that silently skips forever would repeat
that, so the skip reason names exactly what to run.
"""

from __future__ import annotations

import pytest

from openfoam_driver.plugins.cardiacfoam.ionic_catalog_verification import (
    find_listCellModelsVariables_binary,
    verify_ionic_catalog,
)
from openfoam_driver.plugins.cardiacfoam.ionic_model_catalog import (
    IONIC_MODEL_CATALOG,
)

requires_utility = pytest.mark.skipif(
    find_listCellModelsVariables_binary() is None,
    reason=(
        "listCellModelsVariables not on PATH -- source the OpenFOAM bashrc to "
        "run the live catalog check. This is a SKIP, not a pass."
    ),
)


@requires_utility
def test_every_catalogued_ionic_model_matches_the_built_solver():
    """The whole point: catalog names are what the solver actually exposes.

    On failure the message names the model and the direction, because the two
    directions have different consequences. A name the runtime has and the
    catalog lacks is invisible to an agent -- it never learns the override
    exists. A name the catalog has and the runtime lacks is worse in a
    different way: the agent writes it and the solver fatals at startup.
    """
    result = verify_ionic_catalog()

    assert result.utility_available is True
    assert set(result.results) == set(IONIC_MODEL_CATALOG)

    problems = {
        name: {
            "status": model.status,
            "runtime_has_catalog_lacks": {
                k: list(v) for k, v in model.missing_from_catalog.items()
            },
            "catalog_has_runtime_lacks": {
                k: list(v) for k, v in model.extra_in_catalog.items()
            },
            "reason": model.reason,
        }
        for name, model in sorted(result.results.items())
        if model.status != "match"
    }
    assert not problems, (
        "IONIC_MODEL_CATALOG disagrees with the built solver.\n"
        "Regenerate the affected entries from listCellModelsVariables rather "
        "than editing them by hand -- the naming convention differs per model "
        "and cannot be inferred.\n"
        f"{problems}"
    )
    assert result.all_match is True


@requires_utility
def test_every_model_can_actually_be_instantiated():
    """A model the driver cannot configure is a catalog problem in its own
    right, distinct from a name mismatch. Reported separately so a broken
    case-synthesis path is not mistaken for catalog drift."""
    result = verify_ionic_catalog()
    errored = {
        name: model.reason
        for name, model in sorted(result.results.items())
        if model.status == "error"
    }
    assert not errored, f"models the solver would not instantiate: {errored}"
