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
#     test_workflow_runner
#
# Description
#     Tests workflow runner logic and specification contracts.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

from __future__ import annotations

import json
import sys
import tempfile
from pathlib import Path

import pytest

from openfoam_driver.core.runtime.models import DataArtifact
from openfoam_driver.core.runtime.workflow_runner import run_workflow_step
from openfoam_driver.core.runtime.workflow_state import initial_workflow_state


def _dag(command: str, args: list[str], *, produces: list[str] | None = None) -> dict:
    return {
        "schema_version": "1",
        "step_status_values": ["pending", "running", "completed", "failed", "skipped"],
        "steps": [
            {
                "id": "run",
                "command": command,
                "args": args,
                "cwd": ".",
                "depends_on": [],
                "produces": produces or [],
                "consumes": [],
                "retry_policy": {"max_attempts": 1},
                "command_display": " ".join([command, *args]),
            }
        ],
    }


def test_run_workflow_step_completes_and_records_logs_and_state_file() -> None:
    with tempfile.TemporaryDirectory() as temp_dir:
        root = Path(temp_dir)
        state_path = root / "workflow_state.json"
        code = (
            "import json, pathlib, sys; "
            "state=json.loads(pathlib.Path(sys.argv[1]).read_text()); "
            "print(state['status']); "
            "print(state['steps'][0]['status'])"
        )
        dag = _dag(sys.executable, ["-c", code, str(state_path)], produces=["result_csv"])
        state = initial_workflow_state(dag)
        assert state is not None

        result = run_workflow_step(
            dag,
            state,
            "run",
            case_root=root,
            log_dir=root / "logs",
            state_path=state_path,
        )

        payload = result.state.to_json()
        assert payload["status"] == "completed"
        assert payload["current_step_id"] is None
        assert payload["completed_steps"] == ["run"]
        assert payload["failed_step_id"] is None
        assert payload["steps"][0]["status"] == "completed"
        assert payload["steps"][0]["attempt"] == 1
        assert payload["steps"][0]["exit_code"] == 0
        assert payload["steps"][0]["produced_artifacts"] == ["result_csv"]
        assert Path(result.stdout_log).read_text().splitlines() == ["running", "running"]
        assert Path(result.stderr_log).read_text() == ""
        assert json.loads(state_path.read_text()) == payload


def test_run_workflow_step_marks_nonzero_exit_failed() -> None:
    with tempfile.TemporaryDirectory() as temp_dir:
        root = Path(temp_dir)
        code = "import sys; print('failure text', file=sys.stderr); sys.exit(7)"
        dag = _dag(sys.executable, ["-c", code])
        state = initial_workflow_state(dag)
        assert state is not None

        result = run_workflow_step(
            dag,
            state,
            "run",
            case_root=root,
            log_dir=root / "logs",
        )

        payload = result.state.to_json()
        assert payload["status"] == "failed"
        assert payload["current_step_id"] == "run"
        assert payload["failed_step_id"] == "run"
        assert payload["completed_steps"] == []
        assert payload["steps"][0]["status"] == "failed"
        assert payload["steps"][0]["exit_code"] == 7
        assert Path(result.stderr_log).read_text().strip() == "failure text"


def test_run_workflow_step_allows_missing_optional_artifacts() -> None:
    with tempfile.TemporaryDirectory() as temp_dir:
        root = Path(temp_dir)
        dag = _dag(sys.executable, ["-c", "print('ok')"], produces=["optional_vm"])
        state = initial_workflow_state(dag)
        assert state is not None

        result = run_workflow_step(
            dag,
            state,
            "run",
            case_root=root,
            log_dir=root / "logs",
            expected_artifacts=(
                DataArtifact(
                    artifact_id="optional_vm",
                    path_pattern="{time}/Vm",
                    format="openfoam_time_dirs",
                    optional=True,
                ),
            ),
        )

        payload = result.state.to_json()
        assert payload["status"] == "completed"
        assert payload["steps"][0]["status"] == "completed"


def test_run_workflow_step_accepts_decomposed_time_artifact() -> None:
    # A parallel, not-yet-reconstructed run writes processor0/<time>/<field>.
    # The required time_indexed artifact must be considered produced.
    with tempfile.TemporaryDirectory() as temp_dir:
        root = Path(temp_dir)
        script = (
            "import pathlib; "
            "d=pathlib.Path('processor0/0.001'); d.mkdir(parents=True); "
            "(d/'Vm').write_text('x')"
        )
        dag = _dag(sys.executable, ["-c", script], produces=["vm_field"])
        state = initial_workflow_state(dag)
        assert state is not None

        result = run_workflow_step(
            dag, state, "run",
            case_root=root, log_dir=root / "logs",
            expected_artifacts=(
                DataArtifact(
                    artifact_id="vm_field",
                    path_pattern="{time}/Vm",
                    format="openfoam_time_dirs",
                    time_indexed=True,
                ),
            ),
        )

        payload = result.state.to_json()
        assert payload["status"] == "completed"
        assert payload["steps"][0]["status"] == "completed"
        assert payload["steps"][0]["produced_artifacts"] == ["vm_field"]


def test_run_workflow_step_accepts_reconstructed_time_artifact() -> None:
    # Serial / reconstructed location <time>/<field> still satisfies the gate.
    with tempfile.TemporaryDirectory() as temp_dir:
        root = Path(temp_dir)
        script = (
            "import pathlib; "
            "d=pathlib.Path('0.001'); d.mkdir(parents=True); "
            "(d/'Vm').write_text('x')"
        )
        dag = _dag(sys.executable, ["-c", script], produces=["vm_field"])
        state = initial_workflow_state(dag)
        assert state is not None

        result = run_workflow_step(
            dag, state, "run",
            case_root=root, log_dir=root / "logs",
            expected_artifacts=(
                DataArtifact(
                    artifact_id="vm_field",
                    path_pattern="{time}/Vm",
                    format="openfoam_time_dirs",
                    time_indexed=True,
                ),
            ),
        )
        assert result.state.to_json()["status"] == "completed"


def test_run_workflow_step_missing_time_artifact_still_fails() -> None:
    # Neither reconstructed nor decomposed location exists -> missing_artifacts.
    with tempfile.TemporaryDirectory() as temp_dir:
        root = Path(temp_dir)
        dag = _dag(sys.executable, ["-c", "print('ok')"], produces=["vm_field"])
        state = initial_workflow_state(dag)
        assert state is not None

        result = run_workflow_step(
            dag, state, "run",
            case_root=root, log_dir=root / "logs",
            expected_artifacts=(
                DataArtifact(
                    artifact_id="vm_field",
                    path_pattern="{time}/Vm",
                    format="openfoam_time_dirs",
                    time_indexed=True,
                ),
            ),
        )
        payload = result.state.to_json()
        assert payload["status"] == "failed"
        codes = {d["code"] for d in payload["steps"][0]["diagnostics"]}
        assert "missing_artifacts" in codes


def test_run_workflow_step_rejects_incomplete_dependencies() -> None:
    dag = {
        "schema_version": "1",
        "step_status_values": ["pending", "running", "completed", "failed", "skipped"],
        "steps": [
            {
                "id": "mesh",
                "command": sys.executable,
                "args": ["-c", "pass"],
                "cwd": ".",
                "depends_on": [],
                "produces": [],
                "consumes": [],
                "retry_policy": {"max_attempts": 1},
                "command_display": sys.executable,
            },
            {
                "id": "solve",
                "command": sys.executable,
                "args": ["-c", "pass"],
                "cwd": ".",
                "depends_on": ["mesh"],
                "produces": [],
                "consumes": [],
                "retry_policy": {"max_attempts": 1},
                "command_display": sys.executable,
            },
        ],
    }
    state = initial_workflow_state(dag)
    assert state is not None

    with tempfile.TemporaryDirectory() as temp_dir:
        with pytest.raises(ValueError, match="incomplete dependencies"):
            run_workflow_step(dag, state, "solve", case_root=Path(temp_dir))


def test_run_workflow_step_rejects_cwd_escape() -> None:
    dag = {
        "schema_version": "1",
        "step_status_values": ["pending", "running", "completed", "failed", "skipped"],
        "steps": [
            {
                "id": "run",
                "command": sys.executable,
                "args": ["-c", "pass"],
                "cwd": "..",
                "depends_on": [],
                "produces": [],
                "consumes": [],
                "retry_policy": {"max_attempts": 1},
                "command_display": sys.executable,
            },
        ],
    }
    state = initial_workflow_state(dag)
    assert state is not None

    with tempfile.TemporaryDirectory() as temp_dir:
        with pytest.raises(ValueError, match="escapes case root"):
            run_workflow_step(dag, state, "run", case_root=Path(temp_dir))


def test_case_script_step_preserves_dyld_vars_through_shell_hop() -> None:
    # On macOS, /bin/sh is SIP-protected: the OS silently strips inherited
    # DYLD_* env vars before a shebang-interpreted script's own body runs,
    # even though `env=` correctly carried them into the subprocess call.
    # A case-local Allrun-family script is exactly such a shebang script, so
    # invoking it directly with env=execution_env used to lose
    # DYLD_LIBRARY_PATH silently, crashing cardiacFoam with "Library not
    # loaded" deep inside the script. This test proves the value the real
    # script sees matches what was passed in `env`, on whatever platform CI
    # runs on -- the macOS-specific failure mode this guards against can
    # only be observed by actually running on macOS (verified manually).
    with tempfile.TemporaryDirectory() as temp_dir:
        root = Path(temp_dir)
        allrun = root / "Allrun"
        allrun.write_text("#!/bin/sh\necho \"$DYLD_LIBRARY_PATH\"\n")
        allrun.chmod(0o755)

        dag = _dag("Allrun", [])
        state = initial_workflow_state(dag)
        assert state is not None

        marker = "/marker/path/for/regression/test"
        result = run_workflow_step(
            dag,
            state,
            "run",
            case_root=root,
            log_dir=root / "logs",
            env={"PATH": __import__("os").environ.get("PATH", ""), "DYLD_LIBRARY_PATH": marker},
        )

        assert result.state.to_json()["steps"][0]["exit_code"] == 0
        assert Path(result.stdout_log).read_text().strip() == marker


def test_case_script_invocation_embeds_dyld_vars_literally_in_argv() -> None:
    # Mechanism-level check independent of macOS SIP actually being active:
    # for a CASE_SCRIPT_COMMANDS-family step, the argv passed to
    # subprocess.run must carry DYLD_* values as literal text (surviving even
    # if the OS strips them from the *inherited* environment of the shell
    # that's about to exec them), not rely solely on `env=`.
    import subprocess as subprocess_module

    captured: dict[str, object] = {}
    real_run = subprocess_module.run

    def fake_run(argv, **kwargs):
        captured["argv"] = argv
        return real_run((sys.executable, "-c", "pass"), **{k: v for k, v in kwargs.items() if k != "timeout"})

    with tempfile.TemporaryDirectory() as temp_dir:
        root = Path(temp_dir)
        allrun = root / "Allrun"
        allrun.write_text("#!/bin/sh\ntrue\n")
        allrun.chmod(0o755)

        dag = _dag("Allrun", [])
        state = initial_workflow_state(dag)
        assert state is not None

        import openfoam_driver.core.runtime.workflow_runner as workflow_runner_module
        original = workflow_runner_module.subprocess.run
        workflow_runner_module.subprocess.run = fake_run
        try:
            run_workflow_step(
                dag, state, "run", case_root=root, log_dir=root / "logs",
                env={"PATH": __import__("os").environ.get("PATH", ""), "DYLD_LIBRARY_PATH": "/marker/xyz"},
            )
        finally:
            workflow_runner_module.subprocess.run = original

        argv = captured["argv"]
        assert any("/marker/xyz" in str(part) for part in argv), argv
