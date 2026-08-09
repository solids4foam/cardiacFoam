from __future__ import annotations

import os
import stat
import json
from pathlib import Path

from openfoam_driver.cli import main
from openfoam_driver.core.plugin_interface import driver_context
from openfoam_driver.core.runtime.run_document_exec import build_execution_inputs
from openfoam_driver.core.runtime.run_model import RunDocument
from openfoam_driver.tests.plugins.minimal_plugin import MinimalOpenFOAMPlugin


def test_generic_plugin_executes_plain_allrun_case(tmp_path: Path) -> None:
    case_root = tmp_path / "plainOpenFoamCase"
    case_root.mkdir()
    allrun = case_root / "Allrun"
    allrun.write_text("#!/bin/sh\nprintf complete > generic-proof.txt\n")
    allrun.chmod(allrun.stat().st_mode | stat.S_IXUSR)

    exit_code = main([
        "run",
        "--strict",
        "--plugin", "none",
        "--entry", "plainOpenFoamCase",
        "--tutorials-root", str(tmp_path),
    ])

    assert exit_code == 0
    assert (case_root / "generic-proof.txt").read_text() == "complete"
    assert (case_root / "postProcessing" / "workflow_state.json").is_file()


def test_trusted_minimal_plugin_executes_plain_allrun_case(
    tmp_path: Path,
    capsys,
) -> None:
    case_root = tmp_path / "plainOpenFoamCase"
    case_root.mkdir()
    allrun = case_root / "Allrun"
    allrun.write_text("#!/bin/sh\nprintf minimal > minimal-proof.txt\n")
    allrun.chmod(allrun.stat().st_mode | stat.S_IXUSR)

    exit_code = main([
        "run",
        "--strict",
        "--plugin",
        "openfoam_driver.tests.plugins.minimal_plugin:MinimalOpenFOAMPlugin",
        "--entry", "plainOpenFoamCase",
        "--tutorials-root", str(tmp_path),
    ])

    assert exit_code == 0
    assert (case_root / "minimal-proof.txt").read_text() == "minimal"
    payload = json.loads(capsys.readouterr().out)
    assert payload["status"] == "ok"


def test_run_document_validation_uses_the_selected_plugin(tmp_path: Path) -> None:
    case_root = tmp_path / "plainOpenFoamCase"
    case_root.mkdir()
    allrun = case_root / "Allrun"
    allrun.write_text("#!/bin/sh\nexit 0\n")
    allrun.chmod(allrun.stat().st_mode | stat.S_IXUSR)

    context = driver_context(MinimalOpenFOAMPlugin(), source="test")
    run_doc = RunDocument(
        id="minimal-plugin-document",
        name="plainOpenFoamCase",
        status="planned",
        plugin=context.identity.to_json(),
        config={"anatomy": {}, "physics": {}, "stimulus": {}, "solver": {}},
        workflowDag={
            "steps": [{"id": "run", "command": "Allrun", "depends_on": []}],
        },
        launch={
            "caseRoot": str(case_root),
            "outputDir": str(case_root / "postProcessing"),
        },
    )

    inputs, diagnostics = build_execution_inputs(
        run_doc, driver_context=context,
    )

    assert inputs is not None, diagnostics
    assert not [item for item in diagnostics if item["code"] == "run_validation"]


def test_run_document_rejects_a_mismatched_supplied_plugin(tmp_path: Path) -> None:
    case_root = tmp_path / "plainOpenFoamCase"
    case_root.mkdir()
    allrun = case_root / "Allrun"
    allrun.write_text("#!/bin/sh\nexit 0\n")
    allrun.chmod(allrun.stat().st_mode | stat.S_IXUSR)

    context = driver_context(MinimalOpenFOAMPlugin(), source="test")
    planned_plugin = context.identity.to_json() | {"capability_digest": "sha256:wrong"}
    run_doc = RunDocument(
        id="mismatched-plugin-document",
        name="plainOpenFoamCase",
        status="planned",
        plugin=planned_plugin,
        config={"anatomy": {}, "physics": {}, "stimulus": {}, "solver": {}},
        workflowDag={"steps": [{"id": "run", "command": "Allrun", "depends_on": []}]},
        launch={
            "caseRoot": str(case_root),
            "outputDir": str(case_root / "postProcessing"),
        },
    )

    inputs, diagnostics = build_execution_inputs(
        run_doc, driver_context=context,
    )

    assert inputs is None
    assert [item["code"] for item in diagnostics] == ["plugin_identity_mismatch"]
