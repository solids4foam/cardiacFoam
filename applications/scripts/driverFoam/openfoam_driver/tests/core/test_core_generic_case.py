from __future__ import annotations

from pathlib import Path

from openfoam_driver.core.plugin_interface import generic_openfoam_context
from openfoam_driver.strict_planning import strict_plan


def test_plain_allrun_case_plans_without_cardiac_dictionaries(tmp_path: Path) -> None:
    case_root = tmp_path / "plainOpenFoamCase"
    case_root.mkdir()
    (case_root / "Allrun").write_text("#!/bin/sh\nexit 0\n")

    report = strict_plan(
        "plainOpenFoamCase",
        overrides={"tutorials_root": str(tmp_path)},
    )

    assert report.status == "ok"
    assert report.resolved_entry["entry_kind"] == "case_folder"
    assert report.validation_diagnostics == ()
    assert report.readiness_score["blocked_stages"] == []
    assert {artifact.artifact_id for artifact in report.expected_artifacts} == {
        "core.workflow_state",
        "core.workflow_logs",
    }


def test_plain_allrun_case_works_with_the_no_domain_context(tmp_path: Path) -> None:
    case_root = tmp_path / "plainOpenFoamCase"
    case_root.mkdir()
    (case_root / "Allrun").write_text("#!/bin/sh\nexit 0\n")

    report = strict_plan(
        "plainOpenFoamCase",
        overrides={"tutorials_root": str(tmp_path)},
        driver_context=generic_openfoam_context(),
    )

    assert report.status == "ok"
    assert report.plugin["id"] == "org.driverfoam.generic-openfoam"
    assert report.run_document is not None
    assert report.run_document.plugin == report.plugin
