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
#     test_sweep_derivation_catalog
#
# Description
#     Tests the fixed sweep derivation catalog.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

import pytest
from openfoam_driver.sweep_derivation_catalog import get_derivation, SWEEP_DERIVATION_CATALOG
from openfoam_driver.sweep_expansion import SweepValidationError


def test_case_id_template_joins_named_values():
    fn = get_derivation("case_id_template")
    result = fn({"ionicModel": "TNNP", "deltaT": 1e-6})
    assert result == {"caseId": "TNNP_1e-06"}


def test_case_id_template_is_registered():
    assert "case_id_template" in SWEEP_DERIVATION_CATALOG


def test_get_derivation_rejects_unknown_name():
    with pytest.raises(SweepValidationError, match="notARealDerivation"):
        get_derivation("notARealDerivation")


def test_case_id_template_rejects_path_unsafe_values():
    fn = get_derivation("case_id_template")
    with pytest.raises(SweepValidationError, match="path-safe|caseId"):
        fn({"ionicModel": "../TNNP"})


def test_get_derivation_rejects_non_string_name():
    with pytest.raises(SweepValidationError, match="derive"):
        get_derivation(["x"])
