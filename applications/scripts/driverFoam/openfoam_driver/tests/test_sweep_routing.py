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
#     test_sweep_routing
#
# Description
#     Tests routing of resolved sweep values to build_and_launch parameters.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

import pytest

from openfoam_driver.sweep_expansion import SweepValidationError
from openfoam_driver.sweep_routing import route_case_values


def test_selector_keys_route_to_electro_selectors():
    routed = route_case_values(base={}, resolved_axis_values={"ionicModel": "TNNP", "tissue": "epicardialCells"})
    assert routed["electro_selectors"] == {"ionicModel": "TNNP", "tissue": "epicardialCells"}
    assert routed["electro_overrides"] == {}


def test_type_routes_to_physics_selectors():
    routed = route_case_values(base={}, resolved_axis_values={"type": "electroModel"})
    assert routed["physics_selectors"] == {"type": "electroModel"}


def test_delta_t_and_end_time_route_to_dedicated_kwargs():
    routed = route_case_values(base={}, resolved_axis_values={"deltaT": 1e-6, "endTime": 0.5})
    assert routed["delta_t"] == 1e-6
    assert routed["end_time"] == 0.5


def test_other_keys_route_to_electro_overrides():
    routed = route_case_values(
        base={}, resolved_axis_values={"$ELECTRO_MODEL_COEFFS.singleCellStimulus.stim_amplitude": "80"},
    )
    assert routed["electro_overrides"] == {"$ELECTRO_MODEL_COEFFS.singleCellStimulus.stim_amplitude": "80"}


def test_base_selectors_are_preserved_and_extended():
    routed = route_case_values(
        base={"electro_selectors": {"myocardiumSolver": "singleCellSolver", "tissue": "epicardialCells"},
              "physics_selectors": {"type": "electroModel"}},
        resolved_axis_values={"ionicModel": "TNNP"},
    )
    assert routed["electro_selectors"] == {
        "myocardiumSolver": "singleCellSolver", "tissue": "epicardialCells", "ionicModel": "TNNP",
    }
    assert routed["physics_selectors"] == {"type": "electroModel"}


def test_derived_extra_keys_do_not_break_routing():
    # caseId is a bookkeeping value the expander adds; it must not be routed anywhere.
    routed = route_case_values(base={}, resolved_axis_values={"ionicModel": "TNNP", "caseId": "TNNP-label"})
    assert "caseId" not in routed["electro_selectors"]
    assert "caseId" not in routed["electro_overrides"]


def test_unsupported_control_dict_axis_is_rejected():
    # build_and_launch only exposes delta_t/end_time today. Other controlDict
    # entries must fail loudly instead of being misrouted as electro overrides.
    with pytest.raises(SweepValidationError, match="startTime|controlDict"):
        route_case_values(base={}, resolved_axis_values={"startTime": 0.0})
