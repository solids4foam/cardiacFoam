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
#     workflow_state
#
# Description
#     Tracks and mutates internal state of executing workflow steps.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import Any

from .workflow import STEP_STATUS_VALUES


WORKFLOW_STATE_STATUS_VALUES = STEP_STATUS_VALUES


@dataclass(frozen=True)
class WorkflowStepState:
    step_id: str
    status: str
    attempt: int
    command: str
    args: tuple[str, ...]
    cwd: str
    started_at: str | None = None
    finished_at: str | None = None
    exit_code: int | None = None
    stdout_log: str | None = None
    stderr_log: str | None = None
    produced_artifacts: tuple[str, ...] = ()
    diagnostics: tuple[dict[str, Any], ...] = ()

    def to_json(self) -> dict[str, Any]:
        payload = asdict(self)
        payload["args"] = list(self.args)
        payload["produced_artifacts"] = list(self.produced_artifacts)
        payload["diagnostics"] = list(self.diagnostics)
        return payload


def workflow_step_state_from_json(data: dict[str, Any]) -> WorkflowStepState:
    return WorkflowStepState(
        step_id=str(data["step_id"]),
        status=str(data["status"]),
        attempt=int(data["attempt"]),
        command=str(data["command"]),
        args=tuple(str(arg) for arg in data.get("args", ())),
        cwd=str(data["cwd"]),
        started_at=data.get("started_at"),
        finished_at=data.get("finished_at"),
        exit_code=data.get("exit_code"),
        stdout_log=data.get("stdout_log"),
        stderr_log=data.get("stderr_log"),
        produced_artifacts=tuple(str(item) for item in data.get("produced_artifacts", ())),
        diagnostics=tuple(dict(item) for item in data.get("diagnostics", ())),
    )


@dataclass(frozen=True)
class WorkflowRunState:
    status: str
    current_step_id: str | None
    completed_steps: tuple[str, ...]
    failed_step_id: str | None
    steps: tuple[WorkflowStepState, ...] = field(default_factory=tuple)

    def to_json(self) -> dict[str, Any]:
        return {
            "status": self.status,
            "current_step_id": self.current_step_id,
            "completed_steps": list(self.completed_steps),
            "failed_step_id": self.failed_step_id,
            "steps": [step.to_json() for step in self.steps],
        }


def workflow_state_from_json(data: dict[str, Any]) -> WorkflowRunState:
    return WorkflowRunState(
        status=str(data["status"]),
        current_step_id=data.get("current_step_id"),
        completed_steps=tuple(str(step_id) for step_id in data.get("completed_steps", ())),
        failed_step_id=data.get("failed_step_id"),
        steps=tuple(
            workflow_step_state_from_json(step)
            for step in data.get("steps", ())
        ),
    )


def replace_step_state(
    state: WorkflowRunState,
    updated_step: WorkflowStepState,
    *,
    status: str | None = None,
    current_step_id: str | None = None,
    completed_steps: tuple[str, ...] | None = None,
    failed_step_id: str | None = None,
) -> WorkflowRunState:
    return WorkflowRunState(
        status=status if status is not None else state.status,
        current_step_id=current_step_id,
        completed_steps=completed_steps if completed_steps is not None else state.completed_steps,
        failed_step_id=failed_step_id,
        steps=tuple(
            updated_step if step.step_id == updated_step.step_id else step
            for step in state.steps
        ),
    )


def initial_workflow_state(workflow_dag: dict[str, Any] | None) -> WorkflowRunState | None:
    """Build deterministic pending state from a normalized workflow DAG."""
    if workflow_dag is None:
        return None
    raw_steps = workflow_dag.get("steps")
    if not isinstance(raw_steps, list) or not raw_steps:
        return None

    step_states: list[WorkflowStepState] = []
    first_root_step_id: str | None = None
    for raw_step in raw_steps:
        if not isinstance(raw_step, dict):
            continue
        step_id = str(raw_step.get("id", ""))
        depends_on = raw_step.get("depends_on", [])
        if first_root_step_id is None and isinstance(depends_on, list) and not depends_on:
            first_root_step_id = step_id
        args = raw_step.get("args", [])
        step_states.append(WorkflowStepState(
            step_id=step_id,
            status="pending",
            attempt=0,
            command=str(raw_step.get("command", "")),
            args=tuple(str(arg) for arg in args) if isinstance(args, list) else (),
            cwd=str(raw_step.get("cwd", ".")),
        ))

    return WorkflowRunState(
        status="pending",
        current_step_id=first_root_step_id,
        completed_steps=(),
        failed_step_id=None,
        steps=tuple(step_states),
    )
