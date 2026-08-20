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
#     Workflow DAG models, normalization, and the command allowlist
#     (validate_workflow_commands) shared by strict_planning and the
#     run-document adapter.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

from __future__ import annotations

import os
import shlex
import shutil
from dataclasses import asdict, dataclass, field
from pathlib import Path, PurePath
from typing import Any, Iterable

from .models import DataArtifact


STEP_STATUS_VALUES = ("pending", "running", "completed", "failed", "skipped")


# Single owner of the command allowlist, shared by strict_planning and the
# run-document adapter (see threat model in the RunDocument execution plan).
# Solver-NEUTRAL OpenFOAM and driver binaries, always allowed, resolved via
# PATH. Solver binaries (e.g. cardiacFoam) are plugin-owned and arrive through
# the CommandAuthorizationCapability; core must not name any solver here.
# Case-local scripts (Allrun-family) live in CASE_SCRIPT_COMMANDS instead.
CORE_NEUTRAL_COMMANDS = frozenset(
    {
        "blockMesh",
        "checkMesh",
        "decomposePar",
        "gmsh",
        "gmshToFoam",
        "mpirun",
        "postProcess",
        "reconstructPar",
        "setExprFields",
        "topoSet",
        "vtkUnstructuredToFoam",
    }
)

# Compatibility alias for pre-Phase-1 importers. Prefer CORE_NEUTRAL_COMMANDS
# plus the active context's plugin commands.
OPENFOAM_OR_DRIVER_COMMANDS = CORE_NEUTRAL_COMMANDS

# Bare command names that may resolve to a case-LOCAL executable. Every other
# bare name resolves via PATH only, so a case dir cannot shadow a trusted
# binary. Imported by workflow_runner for _resolve_command.
CASE_SCRIPT_COMMANDS = frozenset({"Allrun", "Allclean", "Allrun.pre", "Allrun.post"})


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
    driver_context: Any,
) -> tuple[dict[str, Any] | None, tuple[WorkflowDiagnostic, ...]]:
    """Return an executable-shaped workflow DAG without executing it.

    Existing specs still author the compact form, for example
    ``{"command": "postProcess -func points"}``. Strict planning uses this
    normalizer to expose a stable argv-like contract for a future step runner.

    ``driver_context`` is required, not defaulted: it selects which steps may
    be credited with unclaimed expected artifacts. A forgotten kwarg would
    silently change the returned document (artifacts attached to a different
    step, or to none) with no diagnostic, so the caller must pass it — ``None``
    is still accepted, and means "no plugin solver is authorized".
    """
    if raw_dag is None:
        return None, (
            WorkflowDiagnostic(
                level="error",
                code="missing_workflow_dag",
                message=(
                    "Strict planning requires a workflow_dag with steps. A case "
                    "folder needs an executable Allrun."
                ),
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
        # Only run-style steps may be credited with producing artifacts: the
        # case run script plus whatever solver binaries the context authorizes.
        # The rest of CASE_SCRIPT_COMMANDS is deliberately excluded -- Allclean
        # deletes output rather than producing it. Note solver_commands() only,
        # NOT the full authorized set: auxiliary_commands() are authorized to
        # run but are post-processing, and ``produces`` is enforced per-step
        # (workflow_runner's missing_artifacts check), so crediting one would
        # make a silent solver fail the wrong step.
        producer_commands = {"Allrun"}
        if driver_context is not None:
            producer_commands |= (
                driver_context.capabilities.command_authorization.solver_commands()
            )
        artifact_producer_steps = [
            step for step in steps if step.command in producer_commands
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


def workflow_output_artifacts(
    artifacts: Iterable[DataArtifact],
) -> tuple[DataArtifact, ...]:
    """Return artifacts whose existence is a responsibility of a workflow step.

    ``expectedArtifacts`` also records driverFOAM's own state and log files.
    Those files are created by the executor around a step transition, rather
    than by the solver command itself, so assigning them to a normalized step
    would make the solver incorrectly responsible for driver bookkeeping.
    """
    return tuple(artifact for artifact in artifacts if artifact.produced_by != "driverFOAM")


def _is_installed_openfoam_app(command: str) -> bool:
    """True when ``command`` resolves to an executable installed under
    ``$FOAM_APPBIN`` or ``$FOAM_USER_APPBIN`` — i.e. a real OpenFOAM
    application (core or user-compiled). When neither env var is set
    (OpenFOAM not sourced, e.g. the test suite) this returns ``False``, so
    the allowlist falls back to the always-on core set + utilities.
    """
    roots = []
    for var in ("FOAM_APPBIN", "FOAM_USER_APPBIN"):
        value = os.environ.get(var)
        if value:
            roots.append(Path(value).resolve())
    if not roots:
        return False
    resolved = shutil.which(command)
    if not resolved:
        return False
    resolved_path = Path(resolved).resolve()
    return any(resolved_path.is_relative_to(root) for root in roots)


def validate_workflow_commands(
    workflow_dag: dict[str, Any] | None,
    *,
    driver_context: Any = None,
) -> tuple[WorkflowDiagnostic, ...]:
    """Reject DAG steps whose command is not on the allowlist.

    The allowlist is the union of :data:`CORE_NEUTRAL_COMMANDS`, the commands
    ``driver_context`` authorizes through its
    ``CommandAuthorizationCapability``, :data:`CASE_SCRIPT_COMMANDS`, that
    context's utility manifests that declare ``produces``, and executables
    installed under ``$FOAM_APPBIN`` / ``$FOAM_USER_APPBIN`` (see
    :func:`_is_installed_openfoam_app`). Without a ``driver_context`` no plugin
    command and no utility is authorized, leaving only the core-neutral
    commands, case scripts, and installed OpenFOAM apps. An explicit path form
    (``command`` containing ``/``) is allowed only as ``./<name>`` where
    ``<name>`` is a case script — this keeps the gate in parity with
    ``_resolve_command`` (which lets ``./Allrun`` through) while still refusing
    arbitrary ``./script`` and absolute paths. This is the one owner of the
    command allowlist; both ``strict_plan`` and the run-document adapter call
    it so neither can drift. Runs on the *normalized* DAG, where ``command``
    is the bare executable (args already split out).
    """
    if driver_context is not None:
        authorization = driver_context.capabilities.command_authorization
        # Authorization is the union: both kinds of plugin command may run.
        # The split matters only to the artifact-producer heuristic above.
        plugin_commands = (
            authorization.solver_commands() | authorization.auxiliary_commands()
        )
        utilities = authorization.utility_manifests()
    else:
        plugin_commands = frozenset()
        utilities = {}

    diagnostics: list[WorkflowDiagnostic] = []
    for step in (workflow_dag or {}).get("steps", ()):
        command = step.get("command", "")
        step_id = str(step.get("id", ""))
        if not command:
            diagnostics.append(WorkflowDiagnostic(
                level="error",
                code="workflow_step_without_command",
                message=f"Workflow step {step_id or '<unknown>'!r} has no command.",
                field=step_id,
            ))
            continue
        if "/" in command:
            # Explicit path form. Only ``./<case-script>`` is permitted (gate
            # parity with _resolve_command); an arbitrary ``./script`` or any
            # absolute path is refused so it cannot bypass the allowlist.
            if command.startswith("./") and command[2:] in CASE_SCRIPT_COMMANDS:
                continue
            diagnostics.append(WorkflowDiagnostic(
                level="error",
                code="unknown_workflow_command",
                message=(
                    f"Workflow command {command!r} is an explicit path; only "
                    "./Allrun-family case scripts may be given as a path."
                ),
                field=step_id,
            ))
            continue
        if (
            command in CORE_NEUTRAL_COMMANDS
            or command in plugin_commands
            or command in CASE_SCRIPT_COMMANDS
        ):
            continue
        manifest = utilities.get(command)
        if manifest is not None:
            if manifest.produces:
                continue
            diagnostics.append(WorkflowDiagnostic(
                level="error",
                code="utility_without_produces",
                message=f"Utility {command!r} has no authoritative produces entries.",
                field=step_id,
            ))
            continue
        if _is_installed_openfoam_app(command):
            continue
        diagnostics.append(WorkflowDiagnostic(
            level="error",
            code="unknown_workflow_command",
            message=(
                f"Workflow command {command!r} is not a known OpenFOAM command, "
                "case script, registered utility, or installed OpenFOAM application."
            ),
            field=step_id,
        ))
    return tuple(diagnostics)
