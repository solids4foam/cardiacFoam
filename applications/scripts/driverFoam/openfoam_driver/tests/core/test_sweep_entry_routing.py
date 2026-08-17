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
#     test_sweep_entry_routing
#
# Description
#     Tests routing of resolved sweep values into a registered tutorial's
#     own make_spec kwargs (entry-based sweep mode).
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

from openfoam_driver.sweep_routing import route_entry_case_values


def test_axis_values_pass_through_unchanged_as_make_spec_kwargs():
    # Entry-based sweeps target an existing registered tutorial's own
    # make_spec(**kwargs); there is no fixed vocabulary to route into like
    # build_and_launch's electro_selectors/physics_selectors -- make_spec
    # validates its own kwargs, so the routing here is a pass-through.
    routed = route_entry_case_values(
        base={}, resolved_axis_values={"dx_values": [0.5], "dt_values": [0.01, 0.005]},
    )
    assert routed == {"dx_values": [0.5], "dt_values": [0.01, 0.005]}


def test_case_id_bookkeeping_key_is_stripped():
    # caseId is sweep_expansion's bookkeeping label for the on-disk case
    # directory name in generic (case_folder) mode; it is not a make_spec
    # kwarg for any registered tutorial and must not be forwarded.
    routed = route_entry_case_values(
        base={}, resolved_axis_values={"dx_values": [0.5], "caseId": "dx05"},
    )
    assert routed == {"dx_values": [0.5]}
    assert "caseId" not in routed


def test_nested_object_values_pass_through_unchanged():
    # niederer_2012.py's own make_spec accepts end_time_by_dx as a
    # {dx: end_time} mapping -- entry-mode routing must not flatten or
    # otherwise mutate nested object values.
    routed = route_entry_case_values(
        base={}, resolved_axis_values={"end_time_by_dx": {"0.5": 0.2}},
    )
    assert routed == {"end_time_by_dx": {"0.5": 0.2}}


def test_base_overrides_shared_across_every_case_are_included():
    # A real multi-case sweep (e.g. the Paper I Niederer sweep) has kwargs
    # that are fixed across every axis combination -- solvers/end_time_by_dx
    # aren't swept per case, they're supplied once in "base" and merged into
    # every case's routed overrides, same as generic mode's base electro/
    # physics selectors.
    routed = route_entry_case_values(
        base={"solvers": ["implicit"], "end_time_by_dx": {"0.5": 0.2, "0.1": 0.055}},
        resolved_axis_values={"dx_values": [0.5], "dt_values": [0.01]},
    )
    assert routed == {
        "solvers": ["implicit"],
        "end_time_by_dx": {"0.5": 0.2, "0.1": 0.055},
        "dx_values": [0.5],
        "dt_values": [0.01],
    }


def test_per_case_axis_values_override_base_on_conflict():
    routed = route_entry_case_values(
        base={"dx_values": [0.5]},
        resolved_axis_values={"dx_values": [0.1]},
    )
    assert routed["dx_values"] == [0.1]


def test_entry_key_in_base_is_not_forwarded():
    # "entry" selects which tutorial to target; it is sweep_runner's own
    # dispatch key, not a make_spec kwarg for that tutorial.
    routed = route_entry_case_values(
        base={"entry": "niederer2012", "solvers": ["implicit"]},
        resolved_axis_values={"dx_values": [0.5]},
    )
    assert "entry" not in routed
    assert routed["solvers"] == ["implicit"]


def test_archive_dir_name_key_in_base_is_not_forwarded():
    # archive_dir_name opts a sweep into output_collection.py's generic
    # snapshot/diff archiving (sweep_runner's own bookkeeping) -- it is not a
    # make_spec kwarg for any registered tutorial and must not be forwarded,
    # same as entry/caseId.
    routed = route_entry_case_values(
        base={"entry": "manufacturedBidomain", "archive_dir_name": "driverPostProcessingArchive_postProcessing"},
        resolved_axis_values={"number_cells": [10]},
    )
    assert "archive_dir_name" not in routed
    assert routed["number_cells"] == [10]
