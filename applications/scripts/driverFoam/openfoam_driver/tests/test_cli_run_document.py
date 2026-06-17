"""CLI integration tests for the --run-document execution path (Task 5).

Proves the planning→execution loop: `plan --strict --entry` produces a
RunDocument; `run --run-document <file>` executes that document. Uses a local
Allrun script, so no cardiacFoam binary is required (SKIP_ENV_DIAGNOSTICS is
set suite-wide by tests/conftest.py).
"""
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
    case_root = root / "runDocCase"
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
        "tutorial_family": "run-doc-test",
        "steps": steps,
    }))
    return case_root


def _plan_to_file(tutorials_root: Path, doc_path: Path) -> dict:
    """Run `plan --strict` and write its run_document to doc_path."""
    out = StringIO()
    with redirect_stdout(out):
        code = main([
            "plan", "--strict", "--entry", "runDocCase",
            "--tutorials-root", str(tutorials_root),
        ])
    report = json.loads(out.getvalue())
    assert code == 0, report
    run_document = report["run_document"]
    assert run_document is not None
    doc_path.write_text(json.dumps(run_document))
    return run_document


def test_plan_then_run_document_round_trip_executes() -> None:
    with tempfile.TemporaryDirectory() as temp_dir:
        tutorials_root = Path(temp_dir)
        _write_case(
            tutorials_root,
            allrun="#!/bin/sh\nmkdir -p postProcessing 0.001\ntouch postProcessing/runDocCase_1.txt 0.001/Vm\nprintf 'ran\\n'\n",
            steps=[{"id": "run", "command": "Allrun", "depends_on": []}],
        )
        doc_path = tutorials_root / "run.json"
        _plan_to_file(tutorials_root, doc_path)

        out = StringIO()
        with redirect_stdout(out):
            code = main(["run", "--run-document", str(doc_path)])

        payload = json.loads(out.getvalue())
        assert code == 0, payload
        assert payload["status"] == "ok"
        assert payload["workflow_state"]["status"] == "completed"
        assert Path(payload["workflow_state_path"]).exists()
        assert payload["steps"]
        assert payload["steps"][0]["status"] == "ok"


def test_step_via_run_document_executes_named_step() -> None:
    with tempfile.TemporaryDirectory() as temp_dir:
        tutorials_root = Path(temp_dir)
        _write_case(
            tutorials_root,
            allrun="#!/bin/sh\nmkdir -p postProcessing 0.001\ntouch postProcessing/runDocCase_1.txt 0.001/Vm\nprintf 'ran\\n'\n",
            steps=[{"id": "run", "command": "Allrun", "depends_on": []}],
        )
        doc_path = tutorials_root / "run.json"
        _plan_to_file(tutorials_root, doc_path)

        out = StringIO()
        with redirect_stdout(out):
            code = main(["step", "--run-document", str(doc_path), "--step", "run"])

        payload = json.loads(out.getvalue())
        assert code == 0, payload
        assert payload["status"] == "ok"
        assert payload["workflow_state"]["steps"][0]["status"] == "completed"


def test_run_document_with_bad_dag_surfaces_diagnostics() -> None:
    with tempfile.TemporaryDirectory() as temp_dir:
        doc_path = Path(temp_dir) / "run.json"
        doc_path.write_text(json.dumps({
            "version": "2",
            "id": "d",
            "name": "bad",
            "status": "planned",
            "config": {"anatomy": {}, "physics": {}, "stimulus": {}, "solver": {}},
            "launch": {},
            "workflowDag": None,
        }))

        out = StringIO()
        with redirect_stdout(out):
            code = main(["run", "--run-document", str(doc_path)])

        payload = json.loads(out.getvalue())
        assert code == 1
        assert payload["status"] == "failed"
        codes = {d["code"] for d in payload["diagnostics"]}
        assert "missing_workflow_dag" in codes
        assert "missing_case_root" in codes


def test_run_document_rejects_unknown_command() -> None:
    with tempfile.TemporaryDirectory() as temp_dir:
        doc_path = Path(temp_dir) / "run.json"
        doc_path.write_text(json.dumps({
            "version": "2",
            "id": "d",
            "name": "danger",
            "status": "planned",
            "config": {"anatomy": {}, "physics": {}, "stimulus": {}, "solver": {}},
            "launch": {"caseRoot": temp_dir, "outputDir": temp_dir},
            "workflowDag": {
                "schema_version": "1",
                "step_status_values": [
                    "pending", "running", "completed", "failed", "skipped",
                ],
                "steps": [
                    {"id": "s", "command": "rm", "args": ["-rf", "/"], "cwd": ".",
                     "depends_on": [], "produces": [], "consumes": [],
                     "retry_policy": {}, "command_display": "rm -rf /"},
                ],
            },
        }))

        out = StringIO()
        with redirect_stdout(out):
            code = main(["run", "--run-document", str(doc_path)])

        payload = json.loads(out.getvalue())
        assert code == 1
        codes = {d["code"] for d in payload["diagnostics"]}
        assert "unknown_workflow_command" in codes


def test_run_document_and_entry_are_mutually_exclusive() -> None:
    with tempfile.TemporaryDirectory() as temp_dir:
        doc_path = Path(temp_dir) / "run.json"
        doc_path.write_text("{}")
        try:
            main(["run", "--run-document", str(doc_path), "--entry", "singleCell"])
        except SystemExit as exc:
            assert exc.code == 2  # argparse parser.error
            return
        raise AssertionError("expected SystemExit from mutually-exclusive args")


def test_run_document_only_valid_for_run_and_step() -> None:
    with tempfile.TemporaryDirectory() as temp_dir:
        doc_path = Path(temp_dir) / "run.json"
        doc_path.write_text("{}")
        try:
            main(["describe", "--run-document", str(doc_path)])
        except SystemExit as exc:
            assert exc.code == 2
            return
        raise AssertionError("expected SystemExit for --run-document with describe")
