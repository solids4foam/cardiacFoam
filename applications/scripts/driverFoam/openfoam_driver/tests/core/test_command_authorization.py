"""Command authorization is plugin-owned, not baked into core."""

from __future__ import annotations

from openfoam_driver.core.plugin_interface import (
    generic_openfoam_context,
    default_driver_context,
)
from openfoam_driver.core.runtime.workflow import (
    CORE_NEUTRAL_COMMANDS,
    validate_workflow_commands,
)


def _dag(command: str) -> dict:
    return {"steps": [{"id": "s", "command": command, "depends_on": []}]}


def test_core_neutral_commands_contain_no_solver_names() -> None:
    assert "cardiacFoam" not in CORE_NEUTRAL_COMMANDS
    assert "bathBidomainInterfaceMetrics" not in CORE_NEUTRAL_COMMANDS
    # Solver-neutral OpenFOAM tooling stays in core.
    assert "blockMesh" in CORE_NEUTRAL_COMMANDS
    assert "decomposePar" in CORE_NEUTRAL_COMMANDS
    assert "mpirun" in CORE_NEUTRAL_COMMANDS


def test_cardiac_plugin_authorizes_its_unmanifested_utility() -> None:
    """bathBidomainInterfaceMetrics has no utility.manifest.toml, so it cannot
    come through utility_manifests(); the plugin must authorize it directly or
    the manufacturedFDABathBidomain workflow stops validating."""
    context = default_driver_context()
    errors = [
        d for d in validate_workflow_commands(
            _dag("bathBidomainInterfaceMetrics"), driver_context=context
        )
        if d.level == "error"
    ]
    assert errors == []


def test_cardiac_plugin_authorizes_its_own_solver() -> None:
    context = default_driver_context()
    errors = [
        d for d in validate_workflow_commands(_dag("cardiacFoam"), driver_context=context)
        if d.level == "error"
    ]
    assert errors == []


def test_generic_plugin_does_not_authorize_the_cardiac_solver() -> None:
    context = generic_openfoam_context()
    codes = {
        d.code for d in validate_workflow_commands(
            _dag("cardiacFoam"), driver_context=context
        )
    }
    assert "unknown_workflow_command" in codes


def test_no_context_accepts_only_core_neutral_commands() -> None:
    assert validate_workflow_commands(_dag("blockMesh")) == ()
    codes = {d.code for d in validate_workflow_commands(_dag("cardiacFoam"))}
    assert "unknown_workflow_command" in codes


def test_case_scripts_remain_core_owned() -> None:
    context = generic_openfoam_context()
    assert validate_workflow_commands(_dag("Allrun"), driver_context=context) == ()
    assert validate_workflow_commands(_dag("./Allrun"), driver_context=context) == ()


def test_cardiac_utilities_come_from_the_plugin() -> None:
    context = default_driver_context()
    manifests = context.capabilities.command_authorization.utility_manifests()
    assert "listCellModelsVariables" in manifests
    generic = generic_openfoam_context()
    assert generic.capabilities.command_authorization.utility_manifests() == {}
