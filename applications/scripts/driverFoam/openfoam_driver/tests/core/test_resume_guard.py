"""Interim guard: warn when a completed step predates its executable's mtime."""

from __future__ import annotations

import json
import os
from contextlib import redirect_stdout
from io import StringIO
from pathlib import Path

from openfoam_driver.cli import _execute_run
from openfoam_driver.core.runtime.resume_guard import stale_resume_warnings
from openfoam_driver.core.runtime.workflow_state import (
    WorkflowRunState,
    WorkflowStepState,
)


def _state(finished_at: str, command: str = "myTestSolver") -> WorkflowRunState:
    return WorkflowRunState(
        status="completed",
        current_step_id=None,
        completed_steps=("solve",),
        failed_step_id=None,
        steps=(
            WorkflowStepState(
                step_id="solve",
                status="completed",
                attempt=1,
                command=command,
                args=(),
                cwd=".",
                finished_at=finished_at,
            ),
        ),
    )


def _fake_executable(tmp_path: Path, name: str, mtime_epoch: float) -> None:
    """Put an executable on PATH with a controlled mtime."""
    binary = tmp_path / name
    binary.write_text("#!/bin/sh\nexit 0\n")
    binary.chmod(0o755)
    os.utime(binary, (mtime_epoch, mtime_epoch))
    os.environ["PATH"] = f"{tmp_path}{os.pathsep}" + os.environ["PATH"]


def test_warns_when_executable_is_newer_than_the_completed_step(
    tmp_path: Path, monkeypatch
) -> None:
    monkeypatch.setenv("PATH", os.environ["PATH"])
    # Step finished 2026-01-01; binary rebuilt 2026-06-01.
    _fake_executable(tmp_path, "myTestSolver", 1_780_000_000.0)
    warnings = stale_resume_warnings(
        _state("2026-01-01T00:00:00+00:00"), case_root=tmp_path
    )
    assert len(warnings) == 1
    assert warnings[0]["code"] == "possibly_stale_resume"
    assert warnings[0]["field"] == "solve"
    assert "--fresh" in warnings[0]["message"]


def test_silent_when_executable_is_older_than_the_completed_step(
    tmp_path: Path, monkeypatch
) -> None:
    monkeypatch.setenv("PATH", os.environ["PATH"])
    # Binary built 2020; step finished 2026.
    _fake_executable(tmp_path, "myOldSolver", 1_580_000_000.0)
    warnings = stale_resume_warnings(
        _state("2026-01-01T00:00:00+00:00", command="myOldSolver"),
        case_root=tmp_path,
    )
    assert warnings == ()


def test_silent_when_the_command_cannot_be_resolved(tmp_path: Path) -> None:
    warnings = stale_resume_warnings(
        _state("2026-01-01T00:00:00+00:00", command="definitelyNotOnPath"),
        case_root=tmp_path,
    )
    assert warnings == ()


def test_silent_for_a_step_that_never_completed(tmp_path: Path, monkeypatch) -> None:
    monkeypatch.setenv("PATH", os.environ["PATH"])
    _fake_executable(tmp_path, "myPendingSolver", 1_780_000_000.0)
    state = WorkflowRunState(
        status="pending",
        current_step_id="solve",
        completed_steps=(),
        failed_step_id=None,
        steps=(
            WorkflowStepState(
                step_id="solve",
                status="pending",
                attempt=0,
                command="myPendingSolver",
                args=(),
                cwd=".",
                finished_at=None,
            ),
        ),
    )
    assert stale_resume_warnings(state, case_root=tmp_path) == ()


def test_execute_run_payload_carries_or_omits_resume_warnings_per_cli_wiring(
    tmp_path: Path, monkeypatch
) -> None:
    """Integration test of the actual CLI wiring, not the pure function.

    Drives `cli._execute_run` directly (see cli.py:250-333) -- the failed-status
    early return at cli.py:283-292 and the success payload at cli.py:318-328 both
    thread `resume_warnings` through in the same way, so both are covered here via
    the (only reachable in this scenario) success-payload branch:

    1. A pre-existing workflow_state.json records a completed step whose
       executable was rebuilt after the step finished: the printed JSON payload
       must carry the `resume_warnings` key, populated by `stale_resume_warnings`.
    2. A pre-existing workflow_state.json with nothing to warn about: the printed
       payload must omit the `resume_warnings` key entirely -- not print it as an
       empty list -- so output is byte-identical to before this guard existed.
    """
    monkeypatch.setenv("PATH", os.environ["PATH"])
    # Step finished 2026-01-01; binary rebuilt 2026-06-01.
    _fake_executable(tmp_path, "resumeGuardCliTestSolver", 1_780_000_000.0)

    stale_output_dir = tmp_path / "postProcessingStale"
    stale_output_dir.mkdir()
    stale_state_path = stale_output_dir / "workflow_state.json"
    stale_state = _state(
        "2026-01-01T00:00:00+00:00", command="resumeGuardCliTestSolver"
    )
    stale_state_path.write_text(json.dumps(stale_state.to_json()))

    out = StringIO()
    with redirect_stdout(out):
        code = _execute_run(
            entry_label="resume-guard-cli-test-stale",
            workflow_dag={"steps": []},
            planned_state=None,
            case_root=tmp_path,
            output_dir=stale_output_dir,
            expected_artifacts=(),
            tail_lines=50,
        )

    payload = json.loads(out.getvalue())
    assert code == 0
    assert payload["status"] == "ok"
    assert "resume_warnings" in payload
    assert payload["resume_warnings"][0]["code"] == "possibly_stale_resume"
    assert payload["resume_warnings"][0]["field"] == "solve"
    assert "--fresh" in payload["resume_warnings"][0]["message"]

    clean_output_dir = tmp_path / "postProcessingClean"
    clean_output_dir.mkdir()
    clean_state_path = clean_output_dir / "workflow_state.json"
    clean_state = WorkflowRunState(
        status="completed",
        current_step_id=None,
        completed_steps=(),
        failed_step_id=None,
        steps=(),
    )
    clean_state_path.write_text(json.dumps(clean_state.to_json()))

    out2 = StringIO()
    with redirect_stdout(out2):
        code2 = _execute_run(
            entry_label="resume-guard-cli-test-clean",
            workflow_dag={"steps": []},
            planned_state=None,
            case_root=tmp_path,
            output_dir=clean_output_dir,
            expected_artifacts=(),
            tail_lines=50,
        )

    payload2 = json.loads(out2.getvalue())
    assert code2 == 0
    assert payload2["status"] == "ok"
    assert "resume_warnings" not in payload2
