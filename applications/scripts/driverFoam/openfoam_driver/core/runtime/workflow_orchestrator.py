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
#     workflow_orchestrator
#
# Description
#     Drives a strict workflow to completion, retrying retryable step failures
#     up to a per-step bound with exponential backoff. Decisions derive from
#     WorkflowStepState.status, never exit_code. No config mutation.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

from __future__ import annotations

import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable

from .failure_classification import classify_failure
from .workflow_runner import _atomic_write_json, _step_by_id, _step_state_by_id, run_workflow_step
from .workflow_state import WorkflowRunState, WorkflowStepState, replace_step_state


@dataclass(frozen=True)
class WorkflowRunOutcome:
    state: WorkflowRunState
    steps: tuple[dict[str, Any], ...]


def backoff_delay(attempt: int, backoff_seconds: float, *, cap_seconds: float = 60.0) -> float:
    """Exponential backoff for the (1-based) attempt that just failed, capped."""
    return min(backoff_seconds * (2 ** (attempt - 1)), cap_seconds)


def _resolve_policy(step: dict[str, Any], default_max_attempts: int) -> tuple[int, float]:
    policy = step.get("retry_policy") or {}
    return (
        policy.get("max_attempts", default_max_attempts),
        policy.get("backoff_seconds", 0),
    )


def run_workflow(
    workflow_dag: dict[str, Any],
    workflow_state: WorkflowRunState,
    *,
    case_root: Path,
    output_dir: Path,
    expected_artifacts: tuple = (),
    default_max_attempts: int = 1,
    max_total_attempts: int | None = None,
    classification_overrides: dict[str, str] | None = None,
    runner: Callable[..., Any] = run_workflow_step,
    sleep: Callable[[float], None] = time.sleep,
    state_path: Path | None = None,
    env: dict[str, str] | None = None,
) -> WorkflowRunOutcome:
    """Run pending steps to completion, retrying retryable failures.

    The retry decision block (classify -> policy -> backoff -> re-run) is the
    extension point for a future remediation callback, which would drop in
    immediately before the re-run. This path performs no remediation.

    ``max_total_attempts`` caps the number of step executions across the whole
    run (retry-storm ceiling). ``None`` disables the ceiling, leaving only the
    per-step ``max_attempts`` bound in force.
    """
    resolved_state_path = state_path or (output_dir / "workflow_state.json")
    log_dir = output_dir / "workflow_logs"
    summaries: dict[str, dict[str, Any]] = {}
    total_attempts = 0

    while workflow_state.current_step_id is not None and workflow_state.status == "pending":
        step_id = workflow_state.current_step_id
        total_attempts += 1
        result = runner(
            workflow_dag,
            workflow_state,
            step_id,
            case_root=case_root,
            log_dir=log_dir,
            state_path=resolved_state_path,
            expected_artifacts=expected_artifacts,
            env=env,
        )
        workflow_state = result.state
        step_state = _step_state_by_id(workflow_state, step_id)
        summaries[step_id] = {
            "step": step_id,
            "status": "ok" if step_state.status == "completed" else "failed",
            "exit_code": result.exit_code,
            "stdout_log": result.stdout_log,
            "stderr_log": result.stderr_log,
            "attempts": step_state.attempt,
        }
        if step_state.status != "failed":
            continue

        classification = classify_failure(step_state, overrides=classification_overrides)
        max_attempts, backoff_seconds = _resolve_policy(
            _step_by_id(workflow_dag, step_id), default_max_attempts
        )
        budget_available = max_total_attempts is None or total_attempts < max_total_attempts
        if classification == "retryable" and step_state.attempt < max_attempts and budget_available:
            # Persist a resumable state so a crash during backoff resumes into a
            # retry rather than a refused "failed" state.
            resumable = replace_step_state(
                workflow_state,
                step_state,
                status="pending",
                current_step_id=step_id,
                completed_steps=workflow_state.completed_steps,
                failed_step_id=None,
            )
            _atomic_write_json(resolved_state_path, resumable.to_json())
            workflow_state = resumable
            sleep(backoff_delay(step_state.attempt, backoff_seconds))
            continue
        break  # terminal failure: persisted state remains "failed"

    return WorkflowRunOutcome(state=workflow_state, steps=tuple(summaries.values()))
