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
#     test_generic_case
#
# Description
#     Tests generic_case.make_spec's workflow_dag: it must reflect what
#     run_case actually executes (solver_command/pre_solve_commands via
#     _run_direct, or the run-script/Allrun convention via _run_case),
#     since build_and_launch's non-dry path now runs this workflow_dag
#     through the strict executor instead of DriverEngine.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

from pathlib import Path

from openfoam_driver.specs.tutorials.generic_case import make_spec


def _make_spec(tmp_path, **overrides):
    kwargs = {
        "tutorials_root": tmp_path,
        "case_dir_name": "myCase",
        **overrides,
    }
    return make_spec(**kwargs)


def test_default_run_script_path_uses_allrun_step(tmp_path):
    # No solver_command: this is the registry's generic case-folder fallback
    # (an arbitrary discovered case directory with its own Allrun script),
    # unrelated to build_and_launch -- unchanged from before.
    spec = _make_spec(tmp_path)
    steps = spec.metadata["workflow_dag"]["steps"]
    assert steps == [{"id": "run", "command": "Allrun", "depends_on": []}]


def test_solver_command_only_builds_single_solve_step(tmp_path):
    spec = _make_spec(tmp_path, solver_command="cardiacFoam")
    steps = spec.metadata["workflow_dag"]["steps"]
    assert steps == [
        {"id": "solve", "command": "cardiacFoam", "args": [], "depends_on": []},
    ]


def test_solver_command_with_pre_solve_commands_builds_sequential_steps(tmp_path):
    spec = _make_spec(
        tmp_path,
        solver_command="cardiacFoam",
        pre_solve_commands=["vtkUnstructuredToFoam mesh.vtk", ["setTorsoOrganConductivityField"]],
    )
    steps = spec.metadata["workflow_dag"]["steps"]
    assert steps == [
        {"id": "pre_0", "command": "vtkUnstructuredToFoam", "args": ["mesh.vtk"], "depends_on": []},
        {"id": "pre_1", "command": "setTorsoOrganConductivityField", "args": [], "depends_on": ["pre_0"]},
        {"id": "solve", "command": "cardiacFoam", "args": [], "depends_on": ["pre_1"]},
    ]


def test_solver_command_as_list_is_split_into_command_and_args(tmp_path):
    spec = _make_spec(tmp_path, solver_command=["cardiacFoam", "-parallel"])
    steps = spec.metadata["workflow_dag"]["steps"]
    assert steps == [
        {"id": "solve", "command": "cardiacFoam", "args": ["-parallel"], "depends_on": []},
    ]
