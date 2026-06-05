from __future__ import annotations

import json
import tempfile
from contextlib import redirect_stdout
from io import StringIO
from pathlib import Path

from openfoam_driver.cli import main
from openfoam_driver.scripts._dict_keys_scanner import (
    compute_dict_key_drift,
    strict_dict_key_report,
)
from openfoam_driver.strict_planning import strict_plan


REPO_ROOT = Path(__file__).resolve().parents[5]


def test_strict_plan_succeeds_for_single_cell() -> None:
    report = strict_plan("singleCell")
    payload = report.to_json()

    assert payload["status"] == "ok"
    assert payload["resolved_entry"]["entry_kind"] == "registered_tutorial"
    assert payload["expected_artifacts"]
    assert payload["run_document"]["version"] == "2"
    assert payload["run_document"]["validation"]["status"] == "ok"
    assert payload["workflow_diagnostics"] == []
    assert payload["workflow_dag"]["schema_version"] == "1"
    assert payload["workflow_dag"]["step_status_values"] == [
        "pending",
        "running",
        "completed",
        "failed",
        "skipped",
    ]
    solve_step = payload["workflow_dag"]["steps"][0]
    assert solve_step["command"] == "cardiacFoam"
    assert solve_step["args"] == []
    assert solve_step["cwd"] == "."
    assert solve_step["retry_policy"] == {"max_attempts": 1}
    assert {artifact["artifact_id"] for artifact in payload["expected_artifacts"]} <= set(
        solve_step["produces"]
    )
    assert payload["run_document"]["workflowDag"] == payload["workflow_dag"]
    assert payload["workflow_state"]["status"] == "pending"
    assert payload["workflow_state"]["current_step_id"] == "solve"
    assert payload["workflow_state"]["completed_steps"] == []
    assert payload["workflow_state"]["failed_step_id"] is None
    assert payload["workflow_state"]["steps"][0]["step_id"] == "solve"
    assert payload["workflow_state"]["steps"][0]["status"] == "pending"
    assert payload["workflow_state"]["steps"][0]["attempt"] == 0
    assert payload["workflow_state"]["steps"][0]["command"] == "cardiacFoam"
    assert payload["workflow_state"]["steps"][0]["args"] == []
    assert payload["workflow_state"]["steps"][0]["cwd"] == "."
    assert payload["workflow_state"]["steps"][0]["exit_code"] is None
    assert payload["workflow_state"]["steps"][0]["stdout_log"] is None
    assert payload["workflow_state"]["steps"][0]["stderr_log"] is None
    assert payload["run_document"]["workflowState"] == payload["workflow_state"]


def test_strict_plan_succeeds_for_manufactured_tutorial() -> None:
    report = strict_plan("manufacturedFDA")
    payload = report.to_json()

    assert payload["status"] == "ok"
    assert payload["workflow_dag"]["steps"]
    assert payload["workflow_state"]["current_step_id"] == "mesh"
    assert [step["status"] for step in payload["workflow_state"]["steps"]] == [
        "pending",
        "pending",
    ]
    assert {
        step["step_id"] for step in payload["workflow_state"]["steps"]
    } == {
        step["id"] for step in payload["workflow_dag"]["steps"]
    }
    assert any(
        artifact["artifact_id"] == "verification_error_summary"
        for artifact in payload["expected_artifacts"]
    )


def test_cli_plan_strict_prints_json_and_returns_zero() -> None:
    out = StringIO()
    with redirect_stdout(out):
        code = main(["plan", "--strict", "--entry", "singleCell"])

    payload = json.loads(out.getvalue())
    assert code == 0
    assert payload["status"] == "ok"
    assert payload["launch"]["command"]


def test_strict_plan_fails_on_unknown_workflow_command() -> None:
    with tempfile.TemporaryDirectory() as temp_dir:
        tutorials_root = Path(temp_dir)
        case_root = tutorials_root / "badCase"
        (case_root / "constant").mkdir(parents=True)
        (case_root / "system").mkdir()
        (case_root / "constant" / "physicsProperties").write_text("type electroModel;\n")
        (case_root / "constant" / "electroProperties").write_text(
            "myocardiumSolver singleCellSolver;\n"
            "singleCellSolverCoeffs\n"
            "{\n"
            "    ionicModel AlievPanfilov;\n"
            "    tissue myocyte;\n"
            "    solutionAlgorithm explicit;\n"
            "}\n"
        )
        for name in ("controlDict", "fvSchemes", "fvSolution"):
            (case_root / "system" / name).write_text("\n")
        (case_root / "workflow_contract.json").write_text(json.dumps({
            "tutorial_family": "strict-test",
            "steps": [
                {"id": "unknown", "command": "notARealUtility", "depends_on": []}
            ],
        }))

        report = strict_plan(
            "badCase",
            overrides={"tutorials_root": str(tutorials_root)},
        )

    payload = report.to_json()
    assert payload["status"] == "failed"
    assert any(
        item["code"] == "unknown_workflow_command"
        for item in payload["artifact_diagnostics"]
    )


def test_strict_plan_fails_on_unknown_workflow_dependency() -> None:
    with tempfile.TemporaryDirectory() as temp_dir:
        tutorials_root = Path(temp_dir)
        case_root = tutorials_root / "badDependency"
        (case_root / "constant").mkdir(parents=True)
        (case_root / "system").mkdir()
        (case_root / "constant" / "physicsProperties").write_text("type electroModel;\n")
        (case_root / "constant" / "electroProperties").write_text(
            "myocardiumSolver singleCellSolver;\n"
            "singleCellSolverCoeffs\n"
            "{\n"
            "    ionicModel AlievPanfilov;\n"
            "    tissue myocyte;\n"
            "    solutionAlgorithm explicit;\n"
            "}\n"
        )
        for name in ("controlDict", "fvSchemes", "fvSolution"):
            (case_root / "system" / name).write_text("\n")
        (case_root / "workflow_contract.json").write_text(json.dumps({
            "tutorial_family": "strict-test",
            "steps": [
                {"id": "solve", "command": "cardiacFoam", "depends_on": ["mesh"]}
            ],
        }))

        report = strict_plan(
            "badDependency",
            overrides={"tutorials_root": str(tutorials_root)},
        )

    payload = report.to_json()
    assert payload["status"] == "failed"
    assert payload["run_document"]["status"] == "failed"
    assert payload["run_document"]["validation"]["status"] == "failed"
    assert any(
        item["code"] == "unknown_workflow_dependency"
        for item in payload["workflow_diagnostics"]
    )


def test_strict_plan_fails_when_artifact_prediction_is_empty() -> None:
    with tempfile.TemporaryDirectory() as temp_dir:
        tutorials_root = Path(temp_dir)
        case_root = tutorials_root / "missingArtifacts"
        (case_root / "constant").mkdir(parents=True)
        (case_root / "system").mkdir()
        (case_root / "constant" / "physicsProperties").write_text("type electroModel;\n")
        (case_root / "constant" / "electroProperties").write_text(
            "myocardiumSolver futureSolver;\n"
            "futureSolverCoeffs\n"
            "{\n"
            "    ionicModel AlievPanfilov;\n"
            "}\n"
        )
        for name in ("controlDict", "fvSchemes", "fvSolution"):
            (case_root / "system" / name).write_text("\n")

        report = strict_plan(
            "missingArtifacts",
            overrides={"tutorials_root": str(tutorials_root)},
        )

    payload = report.to_json()
    assert payload["status"] == "failed"
    assert any(
        item["code"] == "empty_artifact_prediction"
        for item in payload["artifact_diagnostics"]
    )


def test_strict_dict_key_scanner_allowlist_is_current() -> None:
    report = strict_dict_key_report(REPO_ROOT / "src")
    assert report.status == "ok"
    assert report.to_json()["unused_allowlist"] == []


def test_strict_dict_key_scanner_fails_on_unallowlisted_key() -> None:
    with tempfile.TemporaryDirectory() as temp_dir:
        src_root = Path(temp_dir) / "src"
        src_root.mkdir()
        (src_root / "reader.C").write_text(
            'void read(const Foam::dictionary& dict) { dict.lookup("unlistedStrictKey"); }\n'
        )
        drift = compute_dict_key_drift(src_root)
        allowlist_path = Path(temp_dir) / "allowlist.json"
        allowlist_path.write_text(json.dumps({
            "absent_keys": sorted(drift["absent_keys"] - {"unlistedStrictKey"}),
            "stale_paths": sorted(drift["stale_paths"]),
            "unmatched_subdicts": sorted(drift["unmatched_subdicts"]),
        }))

        report = strict_dict_key_report(src_root, allowlist_path=allowlist_path)

    payload = report.to_json()
    assert payload["status"] == "failed"
    assert payload["absent_keys"] == ["unlistedStrictKey"]
