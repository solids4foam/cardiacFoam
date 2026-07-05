"""Tests for the capability manifest — the machine-readable surface of what the
driver will accept (allowed commands + samplable field names)."""

from openfoam_driver.capability_manifest import build_capability_manifest
from openfoam_driver.core.runtime.workflow import (
    OPENFOAM_OR_DRIVER_COMMANDS,
    CASE_SCRIPT_COMMANDS,
    validate_workflow_commands,
)


def test_core_commands_match_enforcer():
    manifest = build_capability_manifest()
    assert set(manifest["allowed_commands"]["core"]) == set(OPENFOAM_OR_DRIVER_COMMANDS)
    assert set(manifest["allowed_commands"]["case_scripts"]) == set(CASE_SCRIPT_COMMANDS)


def test_manifest_utilities_are_accepted_by_validator():
    manifest = build_capability_manifest()
    for cmd in manifest["allowed_commands"]["utilities"]:
        dag = {"steps": [{"id": "s", "command": cmd}]}
        errors = [d for d in validate_workflow_commands(dag) if d.level == "error"]
        assert errors == [], f"utility {cmd!r} in manifest but rejected by validator: {errors}"


def test_samplable_fields_for_tnnp_single_cell():
    manifest = build_capability_manifest(
        resolved_solver="singleCellSolver", resolved_ionic_model="TNNP"
    )
    electro = manifest["samplable_fields"]["electro"]
    assert "membrane_V" in electro
    assert "Vm" in electro
    assert "bananas" not in electro
    # single-cell has no mechanics region
    assert manifest["samplable_fields"]["solid"] == []


def test_samplable_fields_multi_region_tags_solid():
    manifest = build_capability_manifest(
        resolved_solver="monodomainSolver",
        resolved_ionic_model="TNNP",
        resolved_active_tension="LandNiederer",
    )
    solid = manifest["samplable_fields"]["solid"]
    assert "Ta" in solid
    assert "lambda" in solid
    # active-tension state variables are included
    assert "XW" in solid


def test_unknown_model_is_not_an_error():
    # An unresolved / unknown model just yields the fixed solver fields, no crash.
    manifest = build_capability_manifest(resolved_ionic_model="NotARealModel")
    assert "Vm" in manifest["samplable_fields"]["electro"]


def test_resolve_case_models_missing_file_is_none():
    from openfoam_driver.capability_manifest import resolve_case_models

    assert resolve_case_models("/nonexistent/case") == (None, None, None)


def test_describe_entry_includes_capability_manifest():
    from openfoam_driver.introspection import describe_entry

    payload = describe_entry("singleCell")
    manifest = payload["capability_manifest"]
    assert "cardiacFoam" in manifest["allowed_commands"]["core"]
    assert "electro" in manifest["samplable_fields"]


def test_strict_plan_carries_capability_manifest(monkeypatch):
    monkeypatch.setenv("SKIP_ENV_DIAGNOSTICS", "1")
    from openfoam_driver.strict_planning import strict_plan

    report = strict_plan("singleCell").to_json()
    assert "cardiacFoam" in report["capability_manifest"]["allowed_commands"]["core"]
