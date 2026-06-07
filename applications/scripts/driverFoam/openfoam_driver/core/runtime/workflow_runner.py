from __future__ import annotations

import json
import os
import re
import subprocess
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Mapping

from .workflow_state import (
    WorkflowRunState,
    WorkflowStepState,
    replace_step_state,
)
from .models import DataArtifact


@dataclass(frozen=True)
class WorkflowStepRunResult:
    state: WorkflowRunState
    step_id: str
    exit_code: int | None
    stdout_log: str
    stderr_log: str


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _safe_step_id(step_id: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", step_id).strip("_") or "step"


def _atomic_write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(f"{path.name}.tmp")
    tmp.write_text(json.dumps(payload, indent=2) + "\n")
    os.replace(tmp, path)


def _step_by_id(workflow_dag: dict[str, Any], step_id: str) -> dict[str, Any]:
    for step in workflow_dag.get("steps", ()):
        if isinstance(step, dict) and step.get("id") == step_id:
            return step
    raise KeyError(f"Workflow step {step_id!r} does not exist")


def _step_state_by_id(state: WorkflowRunState, step_id: str) -> WorkflowStepState:
    for step_state in state.steps:
        if step_state.step_id == step_id:
            return step_state
    raise KeyError(f"Workflow state has no step {step_id!r}")


def _next_runnable_step_id(
    workflow_dag: dict[str, Any],
    state: WorkflowRunState,
) -> str | None:
    completed = set(state.completed_steps)
    states_by_id = {step.step_id: step for step in state.steps}
    for step in workflow_dag.get("steps", ()):
        if not isinstance(step, dict):
            continue
        step_id = str(step.get("id", ""))
        step_state = states_by_id.get(step_id)
        if step_state is None or step_state.status != "pending":
            continue
        depends_on = step.get("depends_on", [])
        if isinstance(depends_on, list) and all(str(dep) in completed for dep in depends_on):
            return step_id
    return None


def _dependencies_completed(
    step: dict[str, Any],
    state: WorkflowRunState,
) -> bool:
    completed = set(state.completed_steps)
    depends_on = step.get("depends_on", [])
    return isinstance(depends_on, list) and all(str(dep) in completed for dep in depends_on)


def _resolve_command(command: str, cwd: Path) -> str:
    if "/" in command:
        return command
    local_command = cwd / command
    if local_command.is_file() and os.access(local_command, os.X_OK):
        return str(local_command)
    return command


def _resolve_case_cwd(case_root: Path, cwd: str) -> Path:
    root = Path(case_root).resolve()
    resolved = (root / cwd).resolve()
    try:
        resolved.relative_to(root)
    except ValueError as exc:
        raise ValueError(f"Workflow cwd {cwd!r} escapes case root {root}") from exc
    return resolved


def run_workflow_step(
    workflow_dag: dict[str, Any],
    workflow_state: WorkflowRunState,
    step_id: str,
    *,
    case_root: Path,
    log_dir: Path | None = None,
    state_path: Path | None = None,
    env: Mapping[str, str] | None = None,
    expected_artifacts: tuple[DataArtifact, ...] = (),
) -> WorkflowStepRunResult:
    """Execute one normalized workflow step and return the updated state.

    This intentionally does not implement resume, retry loops, or multi-step
    orchestration. It only performs one subprocess transition and records logs.
    """
    step = _step_by_id(workflow_dag, step_id)
    previous_step_state = _step_state_by_id(workflow_state, step_id)
    if previous_step_state.status not in {"pending", "failed"}:
        raise ValueError(
            f"Workflow step {step_id!r} is {previous_step_state.status!r}; "
            "only pending or failed steps can be run by this low-level runner"
        )
    if not _dependencies_completed(step, workflow_state):
        raise ValueError(f"Workflow step {step_id!r} has incomplete dependencies")

    attempt = previous_step_state.attempt + 1
    resolved_log_dir = log_dir or (Path(case_root) / "postProcessing" / "workflow_logs")
    resolved_log_dir.mkdir(parents=True, exist_ok=True)
    safe_id = _safe_step_id(step_id)
    stdout_log = resolved_log_dir / f"{safe_id}.attempt{attempt}.stdout.log"
    stderr_log = resolved_log_dir / f"{safe_id}.attempt{attempt}.stderr.log"

    args = tuple(str(arg) for arg in step.get("args", ()))
    cwd = str(step.get("cwd", "."))
    command = str(step["command"])
    resolved_cwd = _resolve_case_cwd(Path(case_root), cwd)
    executable = _resolve_command(command, resolved_cwd)
    running_step = WorkflowStepState(
        step_id=step_id,
        status="running",
        attempt=attempt,
        command=command,
        args=args,
        cwd=cwd,
        started_at=_utc_now(),
        finished_at=None,
        exit_code=None,
        stdout_log=str(stdout_log),
        stderr_log=str(stderr_log),
        produced_artifacts=(),
        diagnostics=(),
    )
    running_state = replace_step_state(
        workflow_state,
        running_step,
        status="running",
        current_step_id=step_id,
        completed_steps=workflow_state.completed_steps,
        failed_step_id=None,
    )
    if state_path is not None:
        _atomic_write_json(Path(state_path), running_state.to_json())

    exit_code: int | None = None
    diagnostics: tuple[dict[str, Any], ...] = ()
    try:
        with stdout_log.open("w") as stdout_handle, stderr_log.open("w") as stderr_handle:
            completed = subprocess.run(
                (executable, *args),
                cwd=resolved_cwd,
                stdout=stdout_handle,
                stderr=stderr_handle,
                env=dict(env) if env is not None else None,
                text=True,
                timeout=step.get("timeout_s"),
                check=False,
            )
        exit_code = completed.returncode
    except subprocess.TimeoutExpired as exc:
        diagnostics = ({
            "level": "error",
            "code": "workflow_step_timeout",
            "message": f"Workflow step {step_id!r} timed out after {exc.timeout} seconds.",
            "field": step_id,
        },)
    except OSError as exc:
        diagnostics = ({
            "level": "error",
            "code": "workflow_step_exec_error",
            "message": str(exc),
            "field": step_id,
        },)

    status = "completed" if exit_code == 0 and not diagnostics else "failed"
    produced_artifacts = tuple(str(item) for item in step.get("produces", ())) if status == "completed" else ()

    if status == "completed" and produced_artifacts:
        import glob
        missing_artifacts = []
        for artifact_id in produced_artifacts:
            for artifact in expected_artifacts:
                if artifact.artifact_id == artifact_id:
                    pattern = str(case_root / artifact.path_pattern.format(case_id=case_root.name, time="*"))
                    if not glob.glob(pattern):
                        missing_artifacts.append(artifact_id)
        if missing_artifacts:
            status = "failed"
            produced_artifacts = ()
            diagnostics = (*diagnostics, {
                "level": "error",
                "code": "missing_artifacts",
                "message": f"Step {step_id!r} completed successfully but missing expected artifacts: {', '.join(missing_artifacts)}",
                "field": step_id,
            })

    final_step = WorkflowStepState(
        step_id=step_id,
        status=status,
        attempt=attempt,
        command=command,
        args=args,
        cwd=cwd,
        started_at=running_step.started_at,
        finished_at=_utc_now(),
        exit_code=exit_code,
        stdout_log=str(stdout_log),
        stderr_log=str(stderr_log),
        produced_artifacts=produced_artifacts,
        diagnostics=diagnostics,
    )

    completed_steps = workflow_state.completed_steps
    if status == "completed" and step_id not in completed_steps:
        completed_steps = (*completed_steps, step_id)

    provisional_state = replace_step_state(
        running_state,
        final_step,
        status="failed" if status == "failed" else "pending",
        current_step_id=step_id if status == "failed" else None,
        completed_steps=completed_steps,
        failed_step_id=step_id if status == "failed" else None,
    )
    if status == "completed":
        next_step_id = _next_runnable_step_id(workflow_dag, provisional_state)
        run_status = "pending" if next_step_id is not None else "completed"
        final_state = replace_step_state(
            provisional_state,
            final_step,
            status=run_status,
            current_step_id=next_step_id,
            completed_steps=completed_steps,
            failed_step_id=None,
        )
    else:
        final_state = provisional_state

    if state_path is not None:
        _atomic_write_json(Path(state_path), final_state.to_json())

    return WorkflowStepRunResult(
        state=final_state,
        step_id=step_id,
        exit_code=exit_code,
        stdout_log=str(stdout_log),
        stderr_log=str(stderr_log),
    )
