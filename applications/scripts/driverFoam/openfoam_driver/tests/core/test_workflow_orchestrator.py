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
#     test_workflow_orchestrator
#
# Description
#     Unit tests for run_workflow, backoff_delay, and resumable-state behaviour.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

import json
from dataclasses import dataclass
from pathlib import Path

from openfoam_driver.core.runtime.workflow_orchestrator import (
    backoff_delay,
    run_workflow,
)
from openfoam_driver.core.runtime.workflow_state import (
    WorkflowRunState,
    WorkflowStepState,
    replace_step_state,
)


def _initial_state(step_id="solve"):
    step = WorkflowStepState(
        step_id=step_id, status="pending", attempt=0,
        command="cardiacFoam", args=(), cwd=".",
    )
    return WorkflowRunState(
        status="pending", current_step_id=step_id,
        completed_steps=(), failed_step_id=None, steps=(step,),
    )


def _dag(step_id="solve", retry_policy=None):
    step = {"id": step_id, "command": "cardiacFoam", "args": [], "cwd": "."}
    step["retry_policy"] = retry_policy if retry_policy is not None else {}
    return {"steps": [step]}


@dataclass
class _FakeResult:
    state: WorkflowRunState
    step_id: str
    exit_code: int | None
    stdout_log: str = "out.log"
    stderr_log: str = "err.log"


def _make_runner(outcomes):
    """outcomes: list of (status, exit_code, codes) consumed per call, in order."""
    state = {"n": 0}

    def runner(workflow_dag, run_state, step_id, *, case_root, log_dir,
               state_path=None, expected_artifacts=(), env=None):
        status, exit_code, codes = outcomes[state["n"]]
        state["n"] += 1
        prev = next((s for s in run_state.steps if s.step_id == step_id), None)
        attempt = (prev.attempt if prev else 0) + 1
        diagnostics = tuple({"level": "error", "code": c, "message": "x"} for c in codes)
        new_step = WorkflowStepState(
            step_id=step_id, status=status, attempt=attempt,
            command="cardiacFoam", args=(), cwd=".",
            exit_code=exit_code, diagnostics=diagnostics,
        )
        if status == "completed":
            new_state = replace_step_state(
                run_state, new_step, status="completed",
                current_step_id=None,
                completed_steps=run_state.completed_steps + (step_id,),
                failed_step_id=None,
            )
        else:
            new_state = replace_step_state(
                run_state, new_step, status="failed",
                current_step_id=step_id,
                completed_steps=run_state.completed_steps,
                failed_step_id=step_id,
            )
        result = _FakeResult(state=new_state, step_id=step_id, exit_code=exit_code)
        if state_path is not None:
            Path(state_path).write_text(json.dumps(new_state.to_json()))
        return result

    return runner, state


def test_backoff_delay_exponential_with_cap():
    assert backoff_delay(1, 2) == 2
    assert backoff_delay(2, 2) == 4
    assert backoff_delay(3, 2) == 8
    assert backoff_delay(10, 2, cap_seconds=60.0) == 60.0
    assert backoff_delay(5, 0) == 0


def test_timeout_retries_then_completes(tmp_path):
    runner, _ = _make_runner([
        ("failed", 1, ["workflow_step_timeout"]),
        ("completed", 0, []),
    ])
    sleeps = []
    outcome = run_workflow(
        _dag(retry_policy={"max_attempts": 2, "backoff_seconds": 1}),
        _initial_state(),
        case_root=tmp_path, output_dir=tmp_path,
        runner=runner, sleep=sleeps.append,
    )
    assert outcome.state.status == "completed"
    summary = outcome.steps[0]
    assert summary["status"] == "ok"
    assert summary["attempts"] == 2
    assert sleeps == [1]  # backoff_delay(1, 1) == 1


def test_retryable_exhausts_attempts(tmp_path):
    runner, _ = _make_runner([
        ("failed", 1, ["workflow_step_timeout"]),
        ("failed", 1, ["workflow_step_timeout"]),
    ])
    outcome = run_workflow(
        _dag(retry_policy={"max_attempts": 2}),
        _initial_state(),
        case_root=tmp_path, output_dir=tmp_path,
        runner=runner, sleep=lambda s: None,
    )
    assert outcome.state.status == "failed"
    assert outcome.steps[0]["attempts"] == 2
    assert outcome.steps[0]["status"] == "failed"


def test_fatal_failure_does_not_retry(tmp_path):
    runner, calls = _make_runner([
        ("failed", 0, ["missing_artifacts"]),  # exit 0 but failed
    ])
    outcome = run_workflow(
        _dag(retry_policy={"max_attempts": 5}),
        _initial_state(),
        case_root=tmp_path, output_dir=tmp_path,
        runner=runner, sleep=lambda s: None,
    )
    assert calls["n"] == 1  # no retry
    assert outcome.state.status == "failed"
    assert outcome.steps[0]["status"] == "failed"  # NOT "ok" despite exit_code 0
    assert outcome.steps[0]["attempts"] == 1


def test_default_max_attempts_knob_enables_retry(tmp_path):
    runner, _ = _make_runner([
        ("failed", 1, ["workflow_step_timeout"]),
        ("completed", 0, []),
    ])
    outcome = run_workflow(
        _dag(retry_policy={}),  # empty policy -> falls back to default
        _initial_state(),
        case_root=tmp_path, output_dir=tmp_path,
        default_max_attempts=2,
        runner=runner, sleep=lambda s: None,
    )
    assert outcome.state.status == "completed"
    assert outcome.steps[0]["attempts"] == 2


def test_max_attempts_one_bails_on_first_failure(tmp_path):
    runner, calls = _make_runner([
        ("failed", 1, ["workflow_step_timeout"]),
    ])
    outcome = run_workflow(
        _dag(retry_policy={}),  # default_max_attempts default is 1
        _initial_state(),
        case_root=tmp_path, output_dir=tmp_path,
        runner=runner, sleep=lambda s: None,
    )
    assert calls["n"] == 1
    assert outcome.state.status == "failed"


def test_max_total_attempts_caps_retries_below_per_step_budget(tmp_path):
    # Per-step budget is generous (5), but the whole-run ceiling is 2: the run
    # must stop after 2 total attempts even though the step would keep retrying.
    runner, calls = _make_runner([
        ("failed", 1, ["workflow_step_timeout"]),
        ("failed", 1, ["workflow_step_timeout"]),
        ("failed", 1, ["workflow_step_timeout"]),
    ])
    outcome = run_workflow(
        _dag(retry_policy={"max_attempts": 5}),
        _initial_state(),
        case_root=tmp_path, output_dir=tmp_path,
        max_total_attempts=2,
        runner=runner, sleep=lambda s: None,
    )
    assert calls["n"] == 2  # stopped by the whole-run ceiling, not per-step
    assert outcome.state.status == "failed"
    assert outcome.steps[0]["attempts"] == 2


def test_max_total_attempts_none_preserves_per_step_behavior(tmp_path):
    # Default (None) must not change existing behavior: per-step max_attempts wins.
    runner, _ = _make_runner([
        ("failed", 1, ["workflow_step_timeout"]),
        ("completed", 0, []),
    ])
    outcome = run_workflow(
        _dag(retry_policy={"max_attempts": 2, "backoff_seconds": 1}),
        _initial_state(),
        case_root=tmp_path, output_dir=tmp_path,
        max_total_attempts=None,
        runner=runner, sleep=lambda s: None,
    )
    assert outcome.state.status == "completed"
    assert outcome.steps[0]["attempts"] == 2


def test_persisted_state_is_resumable_between_retries(tmp_path):
    state_path = tmp_path / "workflow_state.json"
    captured = {}

    base_runner, _ = _make_runner([
        ("failed", 1, ["workflow_step_timeout"]),
        ("completed", 0, []),
    ])

    def runner(*args, **kwargs):
        result = base_runner(*args, **kwargs)
        # Snapshot the on-disk state right after the FIRST (failed) call.
        if "first" not in captured and result.exit_code != 0:
            captured["first"] = json.loads(state_path.read_text())
        return result

    def sleep(_seconds):
        # During backoff the persisted state must be resumable (pending), not failed.
        captured["during_backoff"] = json.loads(state_path.read_text())

    outcome = run_workflow(
        _dag(retry_policy={"max_attempts": 2, "backoff_seconds": 1}),
        _initial_state(),
        case_root=tmp_path, output_dir=tmp_path,
        state_path=state_path,
        runner=runner, sleep=sleep,
    )
    assert outcome.state.status == "completed"
    assert captured["during_backoff"]["status"] == "pending"
    assert captured["during_backoff"]["current_step_id"] == "solve"
