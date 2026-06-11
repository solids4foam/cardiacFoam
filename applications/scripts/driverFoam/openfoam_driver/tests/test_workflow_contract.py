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
#     test_workflow_contract
#
# Description
#     Tests workflow contract logic and specification contracts.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

from __future__ import annotations

from openfoam_driver.core.runtime.models import DataArtifact
from openfoam_driver.core.runtime.workflow import normalize_workflow_dag


def test_normalize_workflow_dag_splits_command_string_into_argv() -> None:
    dag, diagnostics = normalize_workflow_dag({
        "steps": [
            {
                "id": "sample",
                "command": "postProcess -func Niedererpoints -latestTime",
                "depends_on": [],
            }
        ]
    })

    assert diagnostics == ()
    assert dag is not None
    step = dag["steps"][0]
    assert step["command"] == "postProcess"
    assert step["args"] == ["-func", "Niedererpoints", "-latestTime"]
    assert step["command_display"] == "postProcess -func Niedererpoints -latestTime"


def test_normalize_workflow_dag_attaches_unclaimed_artifacts_to_solve_step() -> None:
    expected = (
        DataArtifact(
            artifact_id="vm_series",
            path_pattern="{time}/Vm",
            format="openfoam_time_dirs",
        ),
    )
    dag, diagnostics = normalize_workflow_dag(
        {"steps": [{"id": "solve", "command": "cardiacFoam", "depends_on": []}]},
        expected_artifacts=expected,
    )

    assert diagnostics == ()
    assert dag is not None
    assert dag["steps"][0]["produces"] == ["vm_series"]


def test_normalize_workflow_dag_attaches_unclaimed_artifacts_to_allrun_step() -> None:
    expected = (
        DataArtifact(
            artifact_id="trace",
            path_pattern="postProcessing/*.txt",
            format="csv_sweep",
        ),
    )
    dag, diagnostics = normalize_workflow_dag(
        {"steps": [{"id": "run", "command": "Allrun", "depends_on": []}]},
        expected_artifacts=expected,
    )

    assert diagnostics == ()
    assert dag is not None
    assert dag["steps"][0]["produces"] == ["trace"]


def test_normalize_workflow_dag_reports_dependency_cycle() -> None:
    _dag, diagnostics = normalize_workflow_dag({
        "steps": [
            {"id": "a", "command": "Allrun", "depends_on": ["b"]},
            {"id": "b", "command": "Allrun", "depends_on": ["a"]},
        ]
    })

    assert any(
        diagnostic.code == "workflow_dependency_cycle"
        for diagnostic in diagnostics
    )


def test_normalize_workflow_dag_rejects_cwd_escape() -> None:
    _dag, diagnostics = normalize_workflow_dag({
        "steps": [
            {"id": "run", "command": "Allrun", "cwd": "../outside", "depends_on": []},
        ]
    })

    assert any(
        diagnostic.code == "workflow_cwd_not_case_relative"
        for diagnostic in diagnostics
    )
