from __future__ import annotations

import dataclasses

from openfoam_driver.core.runtime.failure_context import (
    DEFAULT_TAIL_LINES,
    build_failure_context,
)
from openfoam_driver.core.runtime.workflow_state import WorkflowStepState


def _step_state(tmp_path, *, stdout="", stderr="", write_stdout=True, write_stderr=True):
    stdout_log = tmp_path / "step.attempt1.stdout.log"
    stderr_log = tmp_path / "step.attempt1.stderr.log"
    if write_stdout:
        stdout_log.write_bytes(stdout if isinstance(stdout, bytes) else stdout.encode())
    if write_stderr:
        stderr_log.write_bytes(stderr if isinstance(stderr, bytes) else stderr.encode())
    return WorkflowStepState(
        step_id="solve",
        status="failed",
        attempt=2,
        command="cardiacFoam",
        args=(),
        cwd=".",
        exit_code=1,
        stdout_log=str(stdout_log),
        stderr_log=str(stderr_log),
        diagnostics=({"level": "error", "code": "x", "message": "m"},),
    )


def test_tail_returns_last_lines(tmp_path):
    body = "\n".join(f"line{i}" for i in range(500)) + "\n"
    ctx = build_failure_context(_step_state(tmp_path, stdout=body), max_lines=10)
    tail_lines = ctx["stdout_tail"].splitlines()
    assert tail_lines == [f"line{i}" for i in range(490, 500)]
    assert ctx["stdout_truncated"] is True


def test_short_file_not_truncated(tmp_path):
    ctx = build_failure_context(_step_state(tmp_path, stdout="a\nb\nc\n"), max_lines=10)
    assert ctx["stdout_tail"] == "a\nb\nc"
    assert ctx["stdout_truncated"] is False


def test_truncates_on_max_bytes(tmp_path):
    body = "x" * 5000 + "\n"
    ctx = build_failure_context(_step_state(tmp_path, stdout=body), max_lines=1000, max_bytes=1000)
    assert len(ctx["stdout_tail"]) <= 1000
    assert ctx["stdout_truncated"] is True


def test_missing_path_yields_empty_tail(tmp_path):
    state = _step_state(tmp_path, stderr="boom\n", write_stdout=False)
    state = dataclasses.replace(state, stdout_log=str(tmp_path / "nope.log"))
    ctx = build_failure_context(state)
    assert ctx["stdout_tail"] == ""
    assert ctx["stdout_truncated"] is False
    assert ctx["stderr_tail"] == "boom"


def test_empty_file_yields_empty_tail(tmp_path):
    ctx = build_failure_context(_step_state(tmp_path, stdout=""))
    assert ctx["stdout_tail"] == ""
    assert ctx["stdout_truncated"] is False


def test_max_lines_override(tmp_path):
    body = "\n".join(f"l{i}" for i in range(50)) + "\n"
    ctx = build_failure_context(_step_state(tmp_path, stdout=body), max_lines=3)
    assert ctx["stdout_tail"].splitlines() == ["l47", "l48", "l49"]


def test_non_utf8_does_not_raise(tmp_path):
    ctx = build_failure_context(_step_state(tmp_path, stdout=b"ok\n\xff\xfe bad\n"))
    assert "ok" in ctx["stdout_tail"]


def test_echoes_identifying_fields(tmp_path):
    ctx = build_failure_context(_step_state(tmp_path, stdout="hi\n"))
    assert ctx["step_id"] == "solve"
    assert ctx["attempt"] == 2
    assert ctx["exit_code"] == 1
    assert ctx["diagnostics"] == [{"level": "error", "code": "x", "message": "m"}]
    assert ctx["stdout_log"].endswith("step.attempt1.stdout.log")


def test_default_max_lines_applied_when_omitted(tmp_path):
    total = DEFAULT_TAIL_LINES + 25
    body = "\n".join(f"d{i}" for i in range(total)) + "\n"
    ctx = build_failure_context(_step_state(tmp_path, stdout=body))
    tail_lines = ctx["stdout_tail"].splitlines()
    assert len(tail_lines) == DEFAULT_TAIL_LINES
    assert tail_lines == [f"d{i}" for i in range(25, total)]
    assert ctx["stdout_truncated"] is True


def test_none_path_yields_empty_tail(tmp_path):
    state = dataclasses.replace(_step_state(tmp_path, stderr="boom\n"), stdout_log=None)
    ctx = build_failure_context(state)
    assert ctx["stdout_tail"] == ""
    assert ctx["stdout_truncated"] is False
    assert ctx["stderr_tail"] == "boom"


def test_midfile_seek_drops_partial_leading_line(tmp_path):
    body = "AAAAAAAAAA\nBBBBBBBBBB\nCCCCCCCCCC\nDDDDDDDDDD\n"  # 44 bytes
    ctx = build_failure_context(
        _step_state(tmp_path, stdout=body), max_lines=1000, max_bytes=25
    )
    # The seek lands mid-"B"; the partial "BB" fragment must be dropped.
    assert ctx["stdout_tail"] == "CCCCCCCCCC\nDDDDDDDDDD"
    assert ctx["stdout_truncated"] is True


def test_attach_failure_context_includes_candidate_remediations(tmp_path):
    from openfoam_driver import cli
    from openfoam_driver.core.runtime.workflow_state import WorkflowStepState, WorkflowRunState

    step = WorkflowStepState(
        step_id="solve", status="failed", attempt=1, command="cardiacFoam",
        args=(), cwd=str(tmp_path), exit_code=1,
        stdout_log=None, stderr_log=None,
        diagnostics=({"level": "error", "code": "missing_artifacts",
                      "message": "x", "field": "solve"},),
    )
    state = WorkflowRunState(status="failed", current_step_id="solve",
                             completed_steps=(), failed_step_id="solve", steps=(step,))
    payload: dict = {}
    cli._attach_failure_context(payload, state, "solve", tail_lines=50)

    assert "failure_context" in payload
    cr = payload["failure_context"]["candidate_remediations"]
    assert cr and cr[0]["diagnostic_code"] == "missing_artifacts"
    assert cr[0]["source"] == "static"
