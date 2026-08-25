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
#     test_dict_value_quoting
#
# Description
#     A value OpenFOAM cannot lex as a bare word must be emitted quoted.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

"""OpenFOAM cannot lex a bare token that starts with a digit but is not a
number: `dimension 3D;` raises a FatalIOError ("expected word, found label 3").
Tutorials write `dimension "3D";`. The emitter must do the same, and must NOT
quote anything else -- quoting a scalar, vector or dimension set would break
dictionaries that work today.
"""

from __future__ import annotations

import pytest

from openfoam_driver.core.specs.dict_builder import _openfoam_value_token


@pytest.mark.parametrize("value", ["3D", "1D", "2D", "3Dfoo"])
def test_a_word_starting_with_a_digit_is_quoted(value):
    assert _openfoam_value_token(value) == f'"{value}"'


@pytest.mark.parametrize(
    "value",
    [
        "0.5", "1e-5", "-3", "140000", "1.0E+06", ".5", "+2",   # numbers
        "epicardialCells", "TNNP", "godunov", "yes", "no",       # plain words
        "(1 0 0)", "(bath organ)",                               # lists/vectors
        "[ -1 -3 3 0 0 2 0 ] ( 0.11 0 0 )",                      # dimension set
        '"3D"',                                                  # already quoted
        "",                                                      # empty
    ],
)
def test_everything_else_is_left_untouched(value):
    assert _openfoam_value_token(value) == value


def test_dimension_reaches_the_dict_in_a_form_openfoam_can_parse():
    """The end-to-end regression: this exact value produced an unparseable
    dictionary, making $ELECTRO_MODEL_COEFFS.dimension unusable through the
    driver. Found by running listCellModelsVariables for real."""
    from openfoam_driver.plugins.cardiacfoam.dict_builder import build_electro_properties

    text = build_electro_properties(
        selectors={
            "myocardiumSolver": "singleCellSolver",
            "ionicModel": "monodomainFDAManufactured",
            "tissue": "myocyte",
        },
        overrides={"$ELECTRO_MODEL_COEFFS.dimension": "3D"},
    )
    assert 'dimension "3D";' in text
    assert "dimension 3D;" not in text
