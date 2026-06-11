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
#     test_workflow_runner
#
# Description
#     Tests workflow runner logic and specification contracts.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

from __future__ import annotations

import json
import sys
import tempfile
from pathlib import Path

import pytest

from openfoam_driver.core.runtime.workflow_runner import run_workflow_step
from openfoam_driver.core.runtime.workflow_state import initial_workflow_state


def _dag(command: str, args: list[str], *, produces: list[str] | None = None) -> dict:
    return {
        "schema_version": "1",
        "step_status_values": ["pending", "running", "completed", "failed", "skipped"],
        "steps": [
            {
                "id": "run",
                "command": command,
                "args": args,
                "cwd": ".",
                "depends_on": [],
                "produces": produces or [],
                "consumes": [],
                "retry_policy": {"max_attempts": 1},
                "command_display": " ".join([command, *args]),
            }
        ],
    }


def test_run_workflow_step_completes_and_records_logs_and_state_file() -> None:
    with tempfile.TemporaryDirectory() as temp_dir:
        root = Path(temp_dir)
        state_path = root / "workflow_state.json"
        code = (
            "import json, pathlib, sys; "
            "state=json.loads(pathlib.Path(sys.argv[1]).read_text()); "
            "print(state['status']); "
            "print(state['steps'][0]['status'])"
        )
        dag = _dag(sys.executable, ["-c", code, str(state_path)], produces=["result_csv"])
        state = initial_workflow_state(dag)
        assert state is not None

        result = run_workflow_step(
            dag,
            state,
            "run",
            case_root=root,
            log_dir=root / "logs",
            state_path=state_path,
        )

        payload = result.state.to_json()
        assert payload["status"] == "completed"
        assert payload["current_step_id"] is None
        assert payload["completed_steps"] == ["run"]
        assert payload["failed_step_id"] is None
        assert payload["steps"][0]["status"] == "completed"
        assert payload["steps"][0]["attempt"] == 1
        assert payload["steps"][0]["exit_code"] == 0
        assert payload["steps"][0]["produced_artifacts"] == ["result_csv"]
        assert Path(result.stdout_log).read_text().splitlines() == ["running", "running"]
        assert Path(result.stderr_log).read_text() == ""
        assert json.loads(state_path.read_text()) == payload


def test_run_workflow_step_marks_nonzero_exit_failed() -> None:
    with tempfile.TemporaryDirectory() as temp_dir:
        root = Path(temp_dir)
        code = "import sys; print('failure text', file=sys.stderr); sys.exit(7)"
        dag = _dag(sys.executable, ["-c", code])
        state = initial_workflow_state(dag)
        assert state is not None

        result = run_workflow_step(
            dag,
            state,
            "run",
            case_root=root,
            log_dir=root / "logs",
        )

        payload = result.state.to_json()
        assert payload["status"] == "failed"
        assert payload["current_step_id"] == "run"
        assert payload["failed_step_id"] == "run"
        assert payload["completed_steps"] == []
        assert payload["steps"][0]["status"] == "failed"
        assert payload["steps"][0]["exit_code"] == 7
        assert Path(result.stderr_log).read_text().strip() == "failure text"


def test_run_workflow_step_rejects_incomplete_dependencies() -> None:
    dag = {
        "schema_version": "1",
        "step_status_values": ["pending", "running", "completed", "failed", "skipped"],
        "steps": [
            {
                "id": "mesh",
                "command": sys.executable,
                "args": ["-c", "pass"],
                "cwd": ".",
                "depends_on": [],
                "produces": [],
                "consumes": [],
                "retry_policy": {"max_attempts": 1},
                "command_display": sys.executable,
            },
            {
                "id": "solve",
                "command": sys.executable,
                "args": ["-c", "pass"],
                "cwd": ".",
                "depends_on": ["mesh"],
                "produces": [],
                "consumes": [],
                "retry_policy": {"max_attempts": 1},
                "command_display": sys.executable,
            },
        ],
    }
    state = initial_workflow_state(dag)
    assert state is not None

    with tempfile.TemporaryDirectory() as temp_dir:
        with pytest.raises(ValueError, match="incomplete dependencies"):
            run_workflow_step(dag, state, "solve", case_root=Path(temp_dir))


def test_run_workflow_step_rejects_cwd_escape() -> None:
    dag = {
        "schema_version": "1",
        "step_status_values": ["pending", "running", "completed", "failed", "skipped"],
        "steps": [
            {
                "id": "run",
                "command": sys.executable,
                "args": ["-c", "pass"],
                "cwd": "..",
                "depends_on": [],
                "produces": [],
                "consumes": [],
                "retry_policy": {"max_attempts": 1},
                "command_display": sys.executable,
            },
        ],
    }
    state = initial_workflow_state(dag)
    assert state is not None

    with tempfile.TemporaryDirectory() as temp_dir:
        with pytest.raises(ValueError, match="escapes case root"):
            run_workflow_step(dag, state, "run", case_root=Path(temp_dir))
