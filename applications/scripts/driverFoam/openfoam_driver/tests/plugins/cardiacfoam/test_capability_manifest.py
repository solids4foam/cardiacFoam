"""Tests for the capability manifest — the machine-readable surface of what the
driver will accept (allowed commands + samplable field names)."""

from openfoam_driver.capability_manifest import build_capability_manifest

from openfoam_driver.plugins.cardiacfoam.ionic_model_catalog import IONIC_MODEL_CATALOG
from openfoam_driver.plugins.cardiacfoam.active_tension_catalog import ACTIVE_TENSION_MODEL_CATALOG
from openfoam_driver.plugins.cardiacfoam.command_authorization import (
    CARDIAC_AUXILIARY_COMMANDS,
    CARDIAC_SOLVER_COMMANDS,
    utility_manifests,
)
import functools

# The manifest advertises the accept-surface, which is the union of both kinds
# of authorized plugin command -- matching CardiacFoamPlugin.get_capabilities.
CARDIAC_AUTHORIZED_COMMANDS = CARDIAC_SOLVER_COMMANDS | CARDIAC_AUXILIARY_COMMANDS

build_capability_manifest = functools.partial(
    build_capability_manifest,
    ionic_model_catalog=IONIC_MODEL_CATALOG,
    active_tension_model_catalog=ACTIVE_TENSION_MODEL_CATALOG,
    plugin_commands=CARDIAC_AUTHORIZED_COMMANDS,
    utility_manifests=dict(utility_manifests()),
)

from openfoam_driver.core.plugin_interface import default_driver_context
from openfoam_driver.core.runtime.workflow import (
    CORE_NEUTRAL_COMMANDS,
    CASE_SCRIPT_COMMANDS,
    validate_workflow_commands,
)


def test_core_commands_match_enforcer():
    manifest = build_capability_manifest()
    # The enforcer accepts the core-neutral set plus whatever the active
    # plugin authorizes; the manifest must advertise exactly that union.
    assert set(manifest["allowed_commands"]["core"]) == (
        set(CORE_NEUTRAL_COMMANDS) | set(CARDIAC_AUTHORIZED_COMMANDS)
    )
    assert set(manifest["allowed_commands"]["case_scripts"]) == set(CASE_SCRIPT_COMMANDS)


def test_manifest_utilities_are_accepted_by_validator():
    manifest = build_capability_manifest()
    context = default_driver_context()
    for cmd in manifest["allowed_commands"]["utilities"]:
        dag = {"steps": [{"id": "s", "command": cmd}]}
        errors = [
            d for d in validate_workflow_commands(dag, driver_context=context)
            if d.level == "error"
        ]
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


def test_species_labels_are_not_samplable_fields():
    manifest = build_capability_manifest(
        resolved_solver="monodomainSolver", resolved_ionic_model="TNNP"
    )
    electro = manifest["samplable_fields"]["electro"]
    assert "human" not in electro
    assert "pig" not in electro
    assert "generic" not in electro


def test_plain_spatial_ep_has_no_solid_region():
    for solver in ("monodomainSolver", "bidomainSolver", "eikonalSolver"):
        manifest = build_capability_manifest(resolved_solver=solver)
        assert manifest["samplable_fields"]["solid"] == []


def test_single_cell_active_tension_does_not_imply_solid_region():
    manifest = build_capability_manifest(
        resolved_solver="singleCellSolver", resolved_active_tension="LandNiederer"
    )
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
