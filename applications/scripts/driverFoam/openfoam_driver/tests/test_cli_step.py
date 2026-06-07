from __future__ import annotations

import json
import os
import tempfile
from contextlib import redirect_stdout
from io import StringIO
from pathlib import Path

from openfoam_driver.cli import main


REPO_ROOT = Path(__file__).resolve().parents[5]
SINGLE_CELL_ROOT = REPO_ROOT / "tutorials" / "singleCellprotocols" / "singleCell"


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
            allrun="#!/bin/sh\nmkdir -p postProcessing 0.001\ntouch postProcessing/cliStepCase_1.txt 0.001/Vm\nprintf 'step ok\\n'\n",
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
            allrun="#!/bin/sh\nmkdir -p postProcessing 0.001\ntouch postProcessing/cliStepCase_1.txt 0.001/Vm\nprintf 'step failed\\n' >&2\nexit 7\n",
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
            allrun="#!/bin/sh\nmkdir -p postProcessing 0.001\ntouch postProcessing/cliStepCase_1.txt 0.001/Vm\nexit 0\n",
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
            allrun="#!/bin/sh\nmkdir -p postProcessing 0.001\ntouch postProcessing/cliStepCase_1.txt 0.001/Vm\nprintf 'ran %s\\n' \"$1\"\n",
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
            allrun="#!/bin/sh\nmkdir -p postProcessing 0.001\ntouch postProcessing/cliStepCase_1.txt 0.001/Vm\nprintf 'ran %s\\n' \"$1\"\n",
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
                    "touch postProcessing/cliStepCase_1.txt 0.001/Vm\n"
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
            allrun="#!/bin/sh\nmkdir -p postProcessing 0.001\ntouch postProcessing/cliStepCase_1.txt 0.001/Vm\nexit 4\n",
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
