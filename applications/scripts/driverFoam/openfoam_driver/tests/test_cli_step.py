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
#     test_cli_step
#
# Description
#     Tests cli step logic and specification contracts.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

from __future__ import annotations

import json
import os
import tempfile
from contextlib import redirect_stdout
from io import StringIO
from pathlib import Path

from openfoam_driver.cli import main


REPO_ROOT = Path(__file__).resolve().parents[5]
SINGLE_CELL_ROOT = REPO_ROOT / "tutorials" / "electrophysiologyProtocols" / "singleCell"


def _write_case(root: Path, *, allrun: str, steps: list[dict]) -> Path:
    case_root = root / "cliStepCase"
    (case_root / "constant").mkdir(parents=True)
    (case_root / "system").mkdir()
    (case_root / "constant" / "electroProperties").write_text(
        (SINGLE_CELL_ROOT / "constant" / "electroProperties").read_text()
    )
    (case_root / "constant" / "physicsProperties").write_text(
        (SINGLE_CELL_ROOT / "constant" / "physicsProperties").read_text()
    )
    for name in ("controlDict", "fvSchemes", "fvSolution"):
        (case_root / "system" / name).write_text("\n")
    allrun_path = case_root / "Allrun"
    allrun_path.write_text(allrun)
    os.chmod(allrun_path, 0o755)
    (case_root / "workflow_contract.json").write_text(json.dumps({
        "tutorial_family": "cli-step-test",
        "steps": steps,
    }))
    return case_root


def test_cli_step_executes_single_step_and_writes_state_and_logs() -> None:
    with tempfile.TemporaryDirectory() as temp_dir:
        tutorials_root = Path(temp_dir)
        case_root = _write_case(
            tutorials_root,
            allrun="#!/bin/sh\nmkdir -p postProcessing 0.001\ntouch postProcessing/cliStepCase_1.txt 0.001/Vm 0.001/AV_Ta\nprintf 'step ok\\n'\n",
            steps=[{"id": "run", "command": "Allrun", "depends_on": []}],
        )

        out = StringIO()
        with redirect_stdout(out):
            code = main([
                "step",
                "--strict",
                "--entry",
                "cliStepCase",
                "--step",
                "run",
                "--tutorials-root",
                str(tutorials_root),
            ])

        payload = json.loads(out.getvalue())
        state_path = case_root / "postProcessing" / "workflow_state.json"
        assert code == 0
        assert payload["status"] == "ok"
        assert payload["exit_code"] == 0
        assert Path(payload["stdout_log"]).read_text() == "step ok\n"
        assert Path(payload["stderr_log"]).read_text() == ""
        assert json.loads(state_path.read_text()) == payload["workflow_state"]
        assert payload["workflow_state"]["status"] == "completed"
        assert payload["workflow_state"]["steps"][0]["status"] == "completed"


def test_cli_step_returns_nonzero_for_failing_step() -> None:
    with tempfile.TemporaryDirectory() as temp_dir:
        tutorials_root = Path(temp_dir)
        _write_case(
            tutorials_root,
            allrun="#!/bin/sh\nmkdir -p postProcessing 0.001\ntouch postProcessing/cliStepCase_1.txt 0.001/Vm 0.001/AV_Ta\nprintf 'step failed\\n' >&2\nexit 7\n",
            steps=[{"id": "run", "command": "Allrun", "depends_on": []}],
        )

        out = StringIO()
        with redirect_stdout(out):
            code = main([
                "step",
                "--strict",
                "--entry",
                "cliStepCase",
                "--step",
                "run",
                "--tutorials-root",
                str(tutorials_root),
            ])

        payload = json.loads(out.getvalue())
        assert code == 1
        assert payload["status"] == "failed"
        assert payload["exit_code"] == 7
        assert Path(payload["stderr_log"]).read_text() == "step failed\n"
        assert payload["workflow_state"]["status"] == "failed"
        assert payload["workflow_state"]["failed_step_id"] == "run"


def test_cli_step_refuses_dependency_incomplete_step() -> None:
    with tempfile.TemporaryDirectory() as temp_dir:
        tutorials_root = Path(temp_dir)
        _write_case(
            tutorials_root,
            allrun="#!/bin/sh\nmkdir -p postProcessing 0.001\ntouch postProcessing/cliStepCase_1.txt 0.001/Vm 0.001/AV_Ta\nexit 0\n",
            steps=[
                {"id": "mesh", "command": "Allrun", "depends_on": []},
                {"id": "solve", "command": "Allrun", "depends_on": ["mesh"]},
            ],
        )

        out = StringIO()
        with redirect_stdout(out):
            code = main([
                "step",
                "--strict",
                "--entry",
                "cliStepCase",
                "--step",
                "solve",
                "--tutorials-root",
                str(tutorials_root),
            ])

        payload = json.loads(out.getvalue())
        assert code == 1
        assert payload["status"] == "failed"
        assert "incomplete dependencies" in payload["error"]


def test_cli_step_continues_from_existing_workflow_state() -> None:
    with tempfile.TemporaryDirectory() as temp_dir:
        tutorials_root = Path(temp_dir)
        case_root = _write_case(
            tutorials_root,
            allrun="#!/bin/sh\nmkdir -p postProcessing 0.001\ntouch postProcessing/cliStepCase_1.txt 0.001/Vm 0.001/AV_Ta\nprintf 'ran %s\\n' \"$1\"\n",
            steps=[
                {"id": "mesh", "command": "Allrun", "args": ["mesh"], "depends_on": []},
                {"id": "solve", "command": "Allrun", "args": ["solve"], "depends_on": ["mesh"]},
            ],
        )

        first_out = StringIO()
        with redirect_stdout(first_out):
            first_code = main([
                "step",
                "--strict",
                "--entry",
                "cliStepCase",
                "--step",
                "mesh",
                "--tutorials-root",
                str(tutorials_root),
            ])
        assert first_code == 0

        second_out = StringIO()
        with redirect_stdout(second_out):
            second_code = main([
                "step",
                "--strict",
                "--entry",
                "cliStepCase",
                "--step",
                "solve",
                "--tutorials-root",
                str(tutorials_root),
            ])

        payload = json.loads(second_out.getvalue())
        state_path = case_root / "postProcessing" / "workflow_state.json"
        assert second_code == 0
        assert payload["workflow_state"]["status"] == "completed"
        assert payload["workflow_state"]["completed_steps"] == ["mesh", "solve"]
        assert [step["status"] for step in payload["workflow_state"]["steps"]] == [
            "completed",
            "completed",
        ]
        assert json.loads(state_path.read_text()) == payload["workflow_state"]


def test_cli_run_executes_all_runnable_steps_in_order() -> None:
    with tempfile.TemporaryDirectory() as temp_dir:
        tutorials_root = Path(temp_dir)
        case_root = _write_case(
            tutorials_root,
            allrun="#!/bin/sh\nmkdir -p postProcessing 0.001\ntouch postProcessing/cliStepCase_1.txt 0.001/Vm 0.001/AV_Ta\nprintf 'ran %s\\n' \"$1\"\n",
            steps=[
                {"id": "mesh", "command": "Allrun", "args": ["mesh"], "depends_on": []},
                {"id": "solve", "command": "Allrun", "args": ["solve"], "depends_on": ["mesh"]},
            ],
        )

        out = StringIO()
        with redirect_stdout(out):
            code = main([
                "run",
                "--strict",
                "--entry",
                "cliStepCase",
                "--tutorials-root",
                str(tutorials_root),
            ])

        payload = json.loads(out.getvalue())
        state_path = case_root / "postProcessing" / "workflow_state.json"
        assert code == 0
        assert payload["status"] == "ok"
        assert [step["step"] for step in payload["steps"]] == ["mesh", "solve"]
        assert payload["workflow_state"]["status"] == "completed"
        assert payload["workflow_state"]["completed_steps"] == ["mesh", "solve"]
        assert json.loads(state_path.read_text()) == payload["workflow_state"]


def test_cli_run_stops_on_failed_step() -> None:
    with tempfile.TemporaryDirectory() as temp_dir:
        tutorials_root = Path(temp_dir)
        _write_case(
            tutorials_root,
            allrun=(
                "#!/bin/sh\n"
                    "mkdir -p postProcessing 0.001\n"
                    "touch postProcessing/cliStepCase_1.txt 0.001/Vm 0.001/AV_Ta\n"
                "if [ \"$1\" = mesh ]; then printf 'mesh\\n'; exit 0; fi\n"
                "printf 'solve failed\\n' >&2\n"
                "exit 9\n"
            ),
            steps=[
                {"id": "mesh", "command": "Allrun", "args": ["mesh"], "depends_on": []},
                {"id": "solve", "command": "Allrun", "args": ["solve"], "depends_on": ["mesh"]},
            ],
        )

        out = StringIO()
        with redirect_stdout(out):
            code = main([
                "run",
                "--strict",
                "--entry",
                "cliStepCase",
                "--tutorials-root",
                str(tutorials_root),
            ])

        payload = json.loads(out.getvalue())
        assert code == 1
        assert payload["status"] == "failed"
        assert [step["step"] for step in payload["steps"]] == ["mesh", "solve"]
        assert payload["steps"][-1]["exit_code"] == 9
        assert payload["workflow_state"]["status"] == "failed"
        assert payload["workflow_state"]["failed_step_id"] == "solve"


def test_cli_run_does_not_retry_failed_saved_state() -> None:
    with tempfile.TemporaryDirectory() as temp_dir:
        tutorials_root = Path(temp_dir)
        case_root = _write_case(
            tutorials_root,
            allrun="#!/bin/sh\nmkdir -p postProcessing 0.001\ntouch postProcessing/cliStepCase_1.txt 0.001/Vm 0.001/AV_Ta\nexit 4\n",
            steps=[{"id": "run", "command": "Allrun", "depends_on": []}],
        )

        first_out = StringIO()
        with redirect_stdout(first_out):
            first_code = main([
                "run",
                "--strict",
                "--entry",
                "cliStepCase",
                "--tutorials-root",
                str(tutorials_root),
            ])
        assert first_code == 1

        (case_root / "Allrun").write_text("#!/bin/sh\nexit 0\n")
        os.chmod(case_root / "Allrun", 0o755)
        second_out = StringIO()
        with redirect_stdout(second_out):
            second_code = main([
                "run",
                "--strict",
                "--entry",
                "cliStepCase",
                "--tutorials-root",
                str(tutorials_root),
            ])

        payload = json.loads(second_out.getvalue())
        assert second_code == 1
        assert payload["status"] == "failed"
        assert "failed; use action=step" in payload["error"]
        assert payload["workflow_state"]["steps"][0]["attempt"] == 1


def _failed_exit0_runner(
    workflow_dag,
    workflow_state,
    step_id,
    *,
    case_root,
    log_dir,
    state_path,
    expected_artifacts=(),
    env=None,
):
    """Stub for run_workflow_step: marks the step failed with exit_code == 0.

    Simulates the missing_artifacts case (command 'succeeded' but produced
    nothing). Writes real stdout/stderr log files and persists state.
    """
    from openfoam_driver.core.runtime.workflow_runner import (
        WorkflowStepRunResult,
        _step_state_by_id,
    )
    from openfoam_driver.core.runtime.workflow_state import (
        WorkflowStepState,
        replace_step_state,
    )

    log_path = Path(log_dir)
    log_path.mkdir(parents=True, exist_ok=True)
    previous = _step_state_by_id(workflow_state, step_id)
    attempt = previous.attempt + 1
    stdout_log = log_path / f"{step_id}.attempt{attempt}.stdout.log"
    stderr_log = log_path / f"{step_id}.attempt{attempt}.stderr.log"
    stdout_log.write_text("starting solve\n")
    stderr_log.write_text("FOAM FATAL ERROR: missing expected artifacts\n")

    failed_step = WorkflowStepState(
        step_id=step_id,
        status="failed",
        attempt=attempt,
        command=previous.command,
        args=previous.args,
        cwd=previous.cwd,
        exit_code=0,
        stdout_log=str(stdout_log),
        stderr_log=str(stderr_log),
        diagnostics=(
            {"level": "error", "code": "missing_artifacts", "message": "missing"},
        ),
    )
    state = replace_step_state(
        workflow_state,
        failed_step,
        status="failed",
        current_step_id=step_id,
        completed_steps=workflow_state.completed_steps,
        failed_step_id=step_id,
    )
    if state_path is not None:
        Path(state_path).write_text(json.dumps(state.to_json()))
    return WorkflowStepRunResult(
        state=state,
        step_id=step_id,
        exit_code=0,
        stdout_log=str(stdout_log),
        stderr_log=str(stderr_log),
    )


def test_cli_step_attaches_failure_context_on_failure(monkeypatch) -> None:
    with tempfile.TemporaryDirectory() as temp_dir:
        tutorials_root = Path(temp_dir)
        _write_case(
            tutorials_root,
            allrun="#!/bin/sh\nmkdir -p postProcessing 0.001\ntouch postProcessing/cliStepCase_1.txt 0.001/Vm 0.001/AV_Ta\nexit 0\n",
            steps=[{"id": "run", "command": "Allrun", "depends_on": []}],
        )
        monkeypatch.setattr("openfoam_driver.cli.run_workflow_step", _failed_exit0_runner)

        out = StringIO()
        with redirect_stdout(out):
            code = main([
                "step", "--strict",
                "--entry", "cliStepCase",
                "--step", "run",
                "--tutorials-root", str(tutorials_root),
            ])

        payload = json.loads(out.getvalue())
        assert code == 1
        assert payload["status"] == "failed"
        ctx = payload["failure_context"]
        assert ctx["step_id"] == "run"
        assert ctx["exit_code"] == 0
        assert "FOAM FATAL ERROR" in ctx["stderr_tail"]


def test_cli_run_attaches_failure_context_for_failed_step() -> None:
    with tempfile.TemporaryDirectory() as temp_dir:
        tutorials_root = Path(temp_dir)
        _write_case(
            tutorials_root,
            allrun=(
                "#!/bin/sh\n"
                "mkdir -p postProcessing 0.001\n"
                "touch postProcessing/cliStepCase_1.txt 0.001/Vm 0.001/AV_Ta\n"
                "if [ \"$1\" = mesh ]; then printf 'mesh\\n'; exit 0; fi\n"
                "printf 'solve blew up\\n' >&2\n"
                "exit 9\n"
            ),
            steps=[
                {"id": "mesh", "command": "Allrun", "args": ["mesh"], "depends_on": []},
                {"id": "solve", "command": "Allrun", "args": ["solve"], "depends_on": ["mesh"]},
            ],
        )

        out = StringIO()
        with redirect_stdout(out):
            code = main([
                "run", "--strict",
                "--entry", "cliStepCase",
                "--tutorials-root", str(tutorials_root),
                "--tail-lines", "50",
            ])

        payload = json.loads(out.getvalue())
        assert code == 1
        assert payload["status"] == "failed"
        ctx = payload["failure_context"]
        assert ctx["step_id"] == "solve"
        assert "solve blew up" in ctx["stderr_tail"]


def test_cli_step_reports_failed_when_status_failed_with_exit_code_zero(monkeypatch) -> None:
    with tempfile.TemporaryDirectory() as temp_dir:
        tutorials_root = Path(temp_dir)
        _write_case(
            tutorials_root,
            allrun="#!/bin/sh\nmkdir -p postProcessing 0.001\ntouch postProcessing/cliStepCase_1.txt 0.001/Vm 0.001/AV_Ta\nexit 0\n",
            steps=[{"id": "run", "command": "Allrun", "depends_on": []}],
        )
        monkeypatch.setattr("openfoam_driver.cli.run_workflow_step", _failed_exit0_runner)

        out = StringIO()
        with redirect_stdout(out):
            code = main([
                "step", "--strict",
                "--entry", "cliStepCase",
                "--step", "run",
                "--tutorials-root", str(tutorials_root),
            ])

        payload = json.loads(out.getvalue())
        assert code == 1
        assert payload["status"] == "failed"
        assert payload["exit_code"] == 0
        assert payload["workflow_state"]["status"] == "failed"


def test_cli_step_apply_invalid_override_does_not_rerun() -> None:
    with tempfile.TemporaryDirectory() as temp_dir:
        tutorials_root = Path(temp_dir)
        _write_case(
            tutorials_root,
            allrun="#!/bin/sh\nmkdir -p postProcessing 0.001\ntouch postProcessing/cliStepCase_1.txt 0.001/Vm 0.001/AV_Ta\nexit 0\n",
            steps=[{"id": "run", "command": "Allrun", "depends_on": []}],
        )
        bad = tutorials_root / "ov.json"
        bad.write_text('[{"driver_path": "notAKey", "value": "1"}]')

        out = StringIO()
        with redirect_stdout(out):
            code = main([
                "step", "--strict", "--entry", "cliStepCase", "--step", "run",
                "--tutorials-root", str(tutorials_root),
                "--apply", str(bad),
            ])

        payload = json.loads(out.getvalue())
        assert code == 1
        assert payload["status"] == "failed"
        assert "notAKey" in payload["error"]
        # rejected before any rerun: no audit line written.
        assert not (tutorials_root / "cliStepCase" / "postProcessing"
                    / "remediation_history.jsonl").exists()


def test_cli_step_apply_valid_override_mutates_reruns_and_audits() -> None:
    with tempfile.TemporaryDirectory() as temp_dir:
        tutorials_root = Path(temp_dir)
        case_root = _write_case(
            tutorials_root,
            allrun="#!/bin/sh\nmkdir -p postProcessing 0.001\ntouch postProcessing/cliStepCase_1.txt 0.001/Vm 0.001/AV_Ta\nexit 0\n",
            steps=[{"id": "run", "command": "Allrun", "depends_on": []}],
        )
        (case_root / "system" / "controlDict").write_text("deltaT    0.001;\nendTime    1;\n")
        good = tutorials_root / "ov.json"
        good.write_text('[{"driver_path": "deltaT", "value": "0.0005"}]')

        out = StringIO()
        with redirect_stdout(out):
            code = main([
                "step", "--strict", "--entry", "cliStepCase", "--step", "run",
                "--tutorials-root", str(tutorials_root),
                "--apply", str(good),
            ])

        payload = json.loads(out.getvalue())
        assert code == 0
        assert payload["status"] == "ok"
        assert "0.0005" in (case_root / "system" / "controlDict").read_text()
        audit = case_root / "postProcessing" / "remediation_history.jsonl"
        rec = json.loads(audit.read_text().splitlines()[0])
        assert rec["applied_overrides"][0]["driver_path"] == "deltaT"
        assert rec["resulting_status"] == "ok"


def test_cli_step_apply_audits_rerun_error(monkeypatch) -> None:
    # A valid override is applied, then the rerun itself raises. The applied mutation must
    # still be recorded so it is never silently lost. monkeypatch the runner on the cli
    # module (where it is bound: `from .core.runtime.workflow_runner import run_workflow_step`).
    def _boom(*args, **kwargs):
        raise RuntimeError("rerun blew up")

    monkeypatch.setattr("openfoam_driver.cli.run_workflow_step", _boom)

    with tempfile.TemporaryDirectory() as temp_dir:
        tutorials_root = Path(temp_dir)
        case_root = _write_case(
            tutorials_root,
            allrun="#!/bin/sh\nexit 0\n",
            steps=[{"id": "run", "command": "Allrun", "depends_on": []}],
        )
        (case_root / "system" / "controlDict").write_text("deltaT    0.001;\nendTime    1;\n")
        # A real --apply is always a *rerun*, so the output dir already exists from the
        # prior attempt; the monkeypatched runner never creates it, so pre-create it here.
        (case_root / "postProcessing").mkdir()
        good = tutorials_root / "ov.json"
        good.write_text('[{"driver_path": "deltaT", "value": "0.0005"}]')

        out = StringIO()
        with redirect_stdout(out):
            code = main([
                "step", "--strict", "--entry", "cliStepCase", "--step", "run",
                "--tutorials-root", str(tutorials_root),
                "--apply", str(good),
            ])

        assert code == 1
        # overrides were applied before the rerun raised -> file mutated and audited.
        assert "0.0005" in (case_root / "system" / "controlDict").read_text()
        audit = case_root / "postProcessing" / "remediation_history.jsonl"
        rec = json.loads(audit.read_text().splitlines()[0])
        assert rec["resulting_status"] == "rerun_error"
        assert rec["applied_overrides"][0]["driver_path"] == "deltaT"
