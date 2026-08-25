"""Command authorization is plugin-owned, not baked into core."""

from __future__ import annotations

import pytest

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
    the manufacturedBathBidomain workflow stops validating."""
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


def _without_installed_openfoam_apps(monkeypatch) -> None:
    """Isolate the PLUGIN authorization rule from the $FOAM_APPBIN rule.

    The allowlist deliberately accepts any executable installed under
    $FOAM_APPBIN / $FOAM_USER_APPBIN, so that a user's own compiled utility
    runs without being catalogued. With OpenFOAM sourced, cardiacFoam IS such
    an executable -- so it is accepted by that rule regardless of which plugin
    is active.

    These two tests originally asserted outright rejection and passed only
    because OpenFOAM happened not to be sourced in the authoring environment:
    they encoded an environment accident as a security property. Suppressing
    the installed-app rule pins what they actually mean -- that the generic
    plugin does not authorize cardiacFoam *as a plugin command*.
    """
    monkeypatch.setattr(
        "openfoam_driver.core.runtime.workflow._is_installed_openfoam_app",
        lambda command: False,
    )


def test_generic_plugin_does_not_authorize_the_cardiac_solver(monkeypatch) -> None:
    _without_installed_openfoam_apps(monkeypatch)
    context = generic_openfoam_context()
    codes = {
        d.code for d in validate_workflow_commands(
            _dag("cardiacFoam"), driver_context=context
        )
    }
    assert "unknown_workflow_command" in codes


def test_no_context_accepts_only_core_neutral_commands(monkeypatch) -> None:
    _without_installed_openfoam_apps(monkeypatch)
    assert validate_workflow_commands(_dag("blockMesh")) == ()
    codes = {d.code for d in validate_workflow_commands(_dag("cardiacFoam"))}
    assert "unknown_workflow_command" in codes


def test_an_installed_openfoam_app_is_authorized_whatever_the_plugin(monkeypatch) -> None:
    """The other half, pinned deliberately rather than left to the environment.

    This is documented behaviour, not a leak: an executable present under
    $FOAM_APPBIN / $FOAM_USER_APPBIN is accepted so a user's own compiled
    utility can run. Stating it here means the boundary is described by tests
    in both directions instead of only the one the environment happened to
    exercise.
    """
    monkeypatch.setattr(
        "openfoam_driver.core.runtime.workflow._is_installed_openfoam_app",
        lambda command: command == "someInstalledApp",
    )
    context = generic_openfoam_context()
    errors = [
        d for d in validate_workflow_commands(
            _dag("someInstalledApp"), driver_context=context
        )
        if d.level == "error"
    ]
    assert errors == []


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


def test_solver_and_auxiliary_commands_are_distinct() -> None:
    """Both are authorized, but only solver_commands() may be credited with a
    run's artifacts (see normalize_workflow_dag's producer heuristic)."""
    auth = default_driver_context().capabilities.command_authorization
    assert auth.solver_commands() == frozenset({"cardiacFoam"})
    assert auth.auxiliary_commands() == frozenset(
        {"bathBidomainInterfaceMetrics", "gradientReconstructionOrder"}
    )
    assert not (auth.solver_commands() & auth.auxiliary_commands())


def test_generic_plugin_authorizes_neither_kind_of_command() -> None:
    auth = generic_openfoam_context().capabilities.command_authorization
    assert auth.solver_commands() == frozenset()
    assert auth.auxiliary_commands() == frozenset()


def test_utility_manifests_are_not_a_shared_mutable_dict() -> None:
    """The cache hands the same object to every caller, so no consumer may be
    able to corrupt the authorization input of all the others."""
    from openfoam_driver.plugins.cardiacfoam.command_authorization import (
        utility_manifests,
    )

    cached = utility_manifests()
    with pytest.raises(TypeError):
        cached["injected"] = object()  # type: ignore[index]

    plugin = default_driver_context().plugin
    handed_out = plugin.get_utility_manifests()
    handed_out["injected"] = object()
    assert "injected" not in plugin.get_utility_manifests()


def test_plugin_utilities_root_matches_the_utility_catalog_root() -> None:
    """Derived from one constant, not recomputed -- a drift would silently
    degrade to no authorized utilities at all."""
    from openfoam_driver.core.utility_catalog import UTILITIES_ROOT
    from openfoam_driver.plugins.cardiacfoam.command_authorization import utility_roots

    assert utility_roots() == (UTILITIES_ROOT,)
