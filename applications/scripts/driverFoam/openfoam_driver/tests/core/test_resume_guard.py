"""Interim guard: warn when a completed step predates its executable's mtime."""

from __future__ import annotations

import os
from pathlib import Path

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


def test_run_payload_omits_the_key_when_there_is_nothing_to_warn_about(
    tmp_path: Path,
) -> None:
    """An unaffected run must produce byte-identical output to before."""
    from openfoam_driver.core.runtime.resume_guard import stale_resume_warnings

    state = WorkflowRunState(
        status="completed",
        current_step_id=None,
        completed_steps=(),
        failed_step_id=None,
        steps=(),
    )
    assert stale_resume_warnings(state, case_root=tmp_path) == ()
