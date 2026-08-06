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
#     test_workflow_state
#
# Description
#     Tests workflow state logic and specification contracts.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

from __future__ import annotations

from openfoam_driver.core.runtime.workflow_state import (
    initial_workflow_state,
    workflow_state_from_json,
)


def test_initial_workflow_state_marks_all_steps_pending() -> None:
    state = initial_workflow_state({
        "schema_version": "1",
        "step_status_values": ["pending", "running", "completed", "failed", "skipped"],
        "steps": [
            {
                "id": "mesh",
                "command": "blockMesh",
                "args": [],
                "cwd": ".",
                "depends_on": [],
                "produces": [],
                "consumes": [],
                "retry_policy": {"max_attempts": 1},
                "command_display": "blockMesh",
            },
            {
                "id": "solve",
                "command": "cardiacFoam",
                "args": [],
                "cwd": ".",
                "depends_on": ["mesh"],
                "produces": ["vm_series"],
                "consumes": [],
                "retry_policy": {"max_attempts": 1},
                "command_display": "cardiacFoam",
            },
        ],
    })

    assert state is not None
    payload = state.to_json()
    assert payload["status"] == "pending"
    assert payload["current_step_id"] == "mesh"
    assert payload["completed_steps"] == []
    assert payload["failed_step_id"] is None
    assert [step["status"] for step in payload["steps"]] == ["pending", "pending"]
    assert [step["attempt"] for step in payload["steps"]] == [0, 0]
    assert [step["step_id"] for step in payload["steps"]] == ["mesh", "solve"]
    assert workflow_state_from_json(payload).to_json() == payload


def test_initial_workflow_state_returns_none_without_steps() -> None:
    assert initial_workflow_state(None) is None
    assert initial_workflow_state({"steps": []}) is None
