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
#     test_tissue_scoped_overrides
#
# Description
#     Tissue-scoped ionic constant overrides must be writable, not just global.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

"""Region-specific ionic overrides.

The solver applies `global` first, then the ONE tissue scope matching each
cell's tissue flag (`ionicModelIO.C:824-848`; flag mapping at `:173-180`).
Valid scopes are global / epicardialCells / mCells / endocardialCells /
myocyte (`:182-189`).

Only `global` used to be catalogued, so an agent asked for an epicardium-only
channel block could only write a whole-tissue override -- which is a different
experiment. `tutorials/PATHOS/BrugadaSyndrome/singleCell` uses exactly this
composition, so it was expressible by hand but not through the driver.
"""

from __future__ import annotations

import pytest

from openfoam_driver.specs.dict_builder import build_electro_properties

_SELECTORS = {
    "myocardiumSolver": "singleCellSolver",
    "ionicModel": "TNNP",
    "tissue": "epicardialCells",
}


@pytest.mark.parametrize(
    "scope", ["global", "epicardialCells", "mCells", "endocardialCells", "myocyte"]
)
def test_every_solver_recognised_scope_is_writable(scope):
    text = build_electro_properties(
        selectors=_SELECTORS,
        overrides={
            f"$ELECTRO_MODEL_COEFFS.ionicConstantOverrides.{scope}.scale.g_Kr": "0.5"
        },
    )
    assert scope in text
    assert "g_Kr" in text


def test_global_and_a_tissue_scope_compose_as_the_solver_applies_them():
    """The Brugada shape: a whole-tissue INa reduction plus epicardium-only
    Ito/ICaL changes. The solver applies global first, then the tissue scope."""
    text = build_electro_properties(
        selectors=_SELECTORS,
        overrides={
            "$ELECTRO_MODEL_COEFFS.ionicConstantOverrides.global.scale.g_Na": "0.2",
            "$ELECTRO_MODEL_COEFFS.ionicConstantOverrides.epicardialCells.scale.g_to": "6.0",
            "$ELECTRO_MODEL_COEFFS.ionicConstantOverrides.epicardialCells.scale.g_CaL": "0.5",
        },
    )
    assert "global" in text and "epicardialCells" in text
    for name in ("g_Na", "g_to", "g_CaL"):
        assert name in text


def test_set_is_writable_per_scope_too():
    text = build_electro_properties(
        selectors=_SELECTORS,
        overrides={
            "$ELECTRO_MODEL_COEFFS.ionicConstantOverrides.endocardialCells.set.g_Ks": "0.098"
        },
    )
    assert "endocardialCells" in text and "g_Ks" in text


def test_unprefixed_tnnp_constant_names_are_accepted():
    """TNNP's runtime names carry no AC_ prefix. The catalog must not force a
    convention the model does not use -- that would be a solver FatalError."""
    text = build_electro_properties(
        selectors=_SELECTORS,
        overrides={
            "$ELECTRO_MODEL_COEFFS.ionicConstantOverrides.global.scale.g_Kr": "0.5"
        },
    )
    assert "g_Kr" in text
    assert "AC_g_Kr" not in text
