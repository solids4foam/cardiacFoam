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
#     workflow
#
# Description
#     Defines models and logic for parsing acyclic workflow dependencies.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

from __future__ import annotations

import shlex
from dataclasses import asdict, dataclass, field
from pathlib import PurePath
from typing import Any, Iterable

from .models import DataArtifact


STEP_STATUS_VALUES = ("pending", "running", "completed", "failed", "skipped")


@dataclass(frozen=True)
class WorkflowDiagnostic:
    level: str
    code: str
    message: str
    field: str = ""


@dataclass(frozen=True)
class WorkflowStep:
    id: str
    command: str
    args: tuple[str, ...] = ()
    cwd: str = "."
    depends_on: tuple[str, ...] = ()
    produces: tuple[str, ...] = ()
    consumes: tuple[str, ...] = ()
    timeout_s: int | None = None
    retry_policy: dict[str, Any] = field(default_factory=dict)
    command_display: str = ""

    def to_json(self) -> dict[str, Any]:
        payload = asdict(self)
        payload["args"] = list(self.args)
        payload["depends_on"] = list(self.depends_on)
        payload["produces"] = list(self.produces)
        payload["consumes"] = list(self.consumes)
        if payload["timeout_s"] is None:
            payload.pop("timeout_s")
        return payload


def _string_list(value: Any, *, field_name: str) -> tuple[tuple[str, ...], WorkflowDiagnostic | None]:
    if value is None:
        return (), None
    if not isinstance(value, list):
        return (), WorkflowDiagnostic(
            level="error",
            code="invalid_workflow_field",
            message=f"Workflow field {field_name!r} must be a list of strings.",
            field=field_name,
        )
    strings: list[str] = []
    for item in value:
        if not isinstance(item, str):
            return (), WorkflowDiagnostic(
                level="error",
                code="invalid_workflow_field",
                message=f"Workflow field {field_name!r} must be a list of strings.",
                field=field_name,
            )
        strings.append(item)
    return tuple(strings), None


def _command_parts(raw_command: Any) -> tuple[str, tuple[str, ...], str, WorkflowDiagnostic | None]:
    if isinstance(raw_command, str):
        try:
            parts = shlex.split(raw_command)
        except ValueError as exc:
            return "", (), raw_command, WorkflowDiagnostic(
                level="error",
                code="invalid_workflow_command",
                message=str(exc),
                field="command",
            )
        if not parts:
            return "", (), raw_command, WorkflowDiagnostic(
                level="error",
                code="workflow_step_without_command",
                message="Workflow step has an empty command.",
                field="command",
            )
        return parts[0], tuple(parts[1:]), raw_command, None
    if isinstance(raw_command, list) and raw_command and all(isinstance(item, str) for item in raw_command):
        return raw_command[0], tuple(raw_command[1:]), shlex.join(raw_command), None
    return "", (), "", WorkflowDiagnostic(
        level="error",
        code="invalid_workflow_command",
        message="Workflow command must be a string or non-empty argv list.",
        field="command",
    )


def _cwd_is_case_relative(cwd: str) -> bool:
    path = PurePath(cwd)
    return not path.is_absolute() and ".." not in path.parts


def normalize_workflow_dag(
    raw_dag: dict[str, Any] | None,
    *,
    expected_artifacts: Iterable[DataArtifact] = (),
    utility_produces: dict[str, tuple[str, ...]] | None = None,
) -> tuple[dict[str, Any] | None, tuple[WorkflowDiagnostic, ...]]:
    """Return an executable-shaped workflow DAG without executing it.

    Existing specs still author the compact form, for example
    ``{"command": "postProcess -func points"}``. Strict planning uses this
    normalizer to expose a stable argv-like contract for a future step runner.
    """
    if raw_dag is None:
        return None, (
            WorkflowDiagnostic(
                level="error",
                code="missing_workflow_dag",
                message="Strict planning requires a workflow_dag with steps.",
                field="workflow_dag",
            ),
        )
    raw_steps = raw_dag.get("steps")
    if not isinstance(raw_steps, list) or not raw_steps:
        return None, (
            WorkflowDiagnostic(
                level="error",
                code="missing_workflow_steps",
                message="Strict planning requires workflow_dag.steps to be a non-empty list.",
                field="workflow_dag.steps",
            ),
        )

    diagnostics: list[WorkflowDiagnostic] = []
    steps: list[WorkflowStep] = []
    seen_ids: set[str] = set()
    utility_produces = utility_produces or {}
    claimed_artifacts: set[str] = set()
    artifact_ids = tuple(artifact.artifact_id for artifact in expected_artifacts)

    for index, raw_step in enumerate(raw_steps):
        if not isinstance(raw_step, dict):
            diagnostics.append(WorkflowDiagnostic(
                level="error",
                code="invalid_workflow_step",
                message="Workflow step must be an object.",
                field=f"workflow_dag.steps[{index}]",
            ))
            continue

        step_id = str(raw_step.get("id", "")).strip()
        if not step_id:
            diagnostics.append(WorkflowDiagnostic(
                level="error",
                code="workflow_step_without_id",
                message="Workflow step id must be a non-empty string.",
                field=f"workflow_dag.steps[{index}].id",
            ))
            step_id = f"step-{index + 1}"
        if step_id in seen_ids:
            diagnostics.append(WorkflowDiagnostic(
                level="error",
                code="duplicate_workflow_step_id",
                message=f"Workflow step id {step_id!r} is duplicated.",
                field=f"workflow_dag.steps[{index}].id",
            ))
        seen_ids.add(step_id)

        command, parsed_args, command_display, command_error = _command_parts(raw_step.get("command"))
        if command_error is not None:
            diagnostics.append(command_error)

        explicit_args, args_error = _string_list(raw_step.get("args"), field_name="args")
        if args_error is not None:
            diagnostics.append(args_error)
        depends_on, depends_error = _string_list(raw_step.get("depends_on"), field_name="depends_on")
        if depends_error is not None:
            diagnostics.append(depends_error)
        produces, produces_error = _string_list(raw_step.get("produces"), field_name="produces")
        if produces_error is not None:
            diagnostics.append(produces_error)
        consumes, consumes_error = _string_list(raw_step.get("consumes"), field_name="consumes")
        if consumes_error is not None:
            diagnostics.append(consumes_error)

        if not produces and command in utility_produces:
            produces = utility_produces[command]
        if produces:
            claimed_artifacts.update(produces)

        retry_policy = raw_step.get("retry_policy", {})
        if not isinstance(retry_policy, dict):
            diagnostics.append(WorkflowDiagnostic(
                level="error",
                code="invalid_workflow_field",
                message="Workflow field 'retry_policy' must be an object.",
                field="retry_policy",
            ))
            retry_policy = {}
        else:
            retry_policy = dict(retry_policy)
            max_attempts = retry_policy.get("max_attempts")
            if max_attempts is not None and (
                isinstance(max_attempts, bool)
                or not isinstance(max_attempts, int)
                or max_attempts < 1
            ):
                diagnostics.append(WorkflowDiagnostic(
                    level="error",
                    code="invalid_workflow_field",
                    message="Workflow field 'retry_policy.max_attempts' must be an integer >= 1.",
                    field="retry_policy.max_attempts",
                ))
                del retry_policy["max_attempts"]
            backoff_seconds = retry_policy.get("backoff_seconds")
            if backoff_seconds is not None and (
                isinstance(backoff_seconds, bool)
                or not isinstance(backoff_seconds, (int, float))
                or backoff_seconds < 0
            ):
                diagnostics.append(WorkflowDiagnostic(
                    level="error",
                    code="invalid_workflow_field",
                    message="Workflow field 'retry_policy.backoff_seconds' must be a number >= 0.",
                    field="retry_policy.backoff_seconds",
                ))
                del retry_policy["backoff_seconds"]

        timeout_s = raw_step.get("timeout_s")
        if timeout_s is not None and (
            isinstance(timeout_s, bool)
            or not isinstance(timeout_s, int)
            or timeout_s <= 0
        ):
            diagnostics.append(WorkflowDiagnostic(
                level="error",
                code="invalid_workflow_field",
                message="Workflow field 'timeout_s' must be a positive integer.",
                field="timeout_s",
            ))
            timeout_s = None

        cwd = raw_step.get("cwd", ".")
        if not isinstance(cwd, str) or not cwd:
            diagnostics.append(WorkflowDiagnostic(
                level="error",
                code="invalid_workflow_field",
                message="Workflow field 'cwd' must be a non-empty string.",
                field="cwd",
            ))
            cwd = "."
        elif not _cwd_is_case_relative(cwd):
            diagnostics.append(WorkflowDiagnostic(
                level="error",
                code="workflow_cwd_not_case_relative",
                message="Workflow field 'cwd' must stay inside the case directory.",
                field="cwd",
            ))

        steps.append(WorkflowStep(
            id=step_id,
            command=command,
            args=parsed_args + explicit_args,
            cwd=cwd,
            depends_on=depends_on,
            produces=produces,
            consumes=consumes,
            timeout_s=timeout_s,
            retry_policy=retry_policy,
            command_display=command_display or command,
        ))

    known_ids = {step.id for step in steps}
    dependencies_by_id = {step.id: step.depends_on for step in steps}
    for step in steps:
        for dependency in step.depends_on:
            if dependency == step.id:
                diagnostics.append(WorkflowDiagnostic(
                    level="error",
                    code="workflow_step_self_dependency",
                    message=f"Workflow step {step.id!r} depends on itself.",
                    field=step.id,
                ))
            elif dependency not in known_ids:
                diagnostics.append(WorkflowDiagnostic(
                    level="error",
                    code="unknown_workflow_dependency",
                    message=f"Workflow step {step.id!r} depends on unknown step {dependency!r}.",
                    field=step.id,
                ))

    visited: set[str] = set()
    visiting: set[str] = set()
    cycle_reported = False

    def visit(step_id: str) -> None:
        nonlocal cycle_reported
        if step_id in visited:
            return
        if step_id in visiting:
            if not cycle_reported:
                diagnostics.append(WorkflowDiagnostic(
                    level="error",
                    code="workflow_dependency_cycle",
                    message="Workflow dependencies must form an acyclic graph.",
                    field=step_id,
                ))
                cycle_reported = True
            return
        visiting.add(step_id)
        for dependency in dependencies_by_id.get(step_id, ()):
            if dependency in dependencies_by_id:
                visit(dependency)
        visiting.remove(step_id)
        visited.add(step_id)

    for step in steps:
        visit(step.id)

    unclaimed_artifacts = tuple(artifact_id for artifact_id in artifact_ids if artifact_id not in claimed_artifacts)
    if unclaimed_artifacts:
        artifact_producer_steps = [
            step for step in steps if step.command in {"cardiacFoam", "Allrun"}
        ]
        if artifact_producer_steps:
            target_id = artifact_producer_steps[-1].id
            steps = [
                (
                    WorkflowStep(
                        id=step.id,
                        command=step.command,
                        args=step.args,
                        cwd=step.cwd,
                        depends_on=step.depends_on,
                        produces=tuple(dict.fromkeys((*step.produces, *unclaimed_artifacts))),
                        consumes=step.consumes,
                        timeout_s=step.timeout_s,
                        retry_policy=step.retry_policy,
                        command_display=step.command_display,
                    )
                    if step.id == target_id else step
                )
                for step in steps
            ]

    return {
        "schema_version": "1",
        "step_status_values": list(STEP_STATUS_VALUES),
        "steps": [step.to_json() for step in steps],
    }, tuple(diagnostics)
