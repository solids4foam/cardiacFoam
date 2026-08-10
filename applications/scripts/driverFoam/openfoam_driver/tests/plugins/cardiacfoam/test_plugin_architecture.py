from __future__ import annotations

from openfoam_driver.core.generic_plugin import GenericOpenFOAMPlugin
from openfoam_driver.core.plugin_interface import (
    SolverPlugin,
    driver_context,
    validate_plugin,
)
from openfoam_driver.plugins.cardiacfoam_plugin import CardiacFoamPlugin
from openfoam_driver.tests.plugins.minimal_plugin import MinimalOpenFOAMPlugin


def test_cardiacfoam_plugin_satisfies_runtime_contract() -> None:
    plugin = validate_plugin(CardiacFoamPlugin())
    assert isinstance(plugin, SolverPlugin)

    ctx = driver_context(plugin, source="test")

    assert ctx.identity.id == "org.cardiacfoam"
    assert plugin.plugin_name == "cardiacFoam"
    assert plugin.get_dict_entries()
    assert "registered_tutorials" in plugin.get_tutorial_catalog()
    assert "spec_factories" in plugin.get_tutorial_catalog()
    assert {"deltaT", "endTime"} <= {
        entry.driver_path
        for entry in plugin.get_dictionary_catalog().documents["controlDict"]
    }


def test_generic_openfoam_plugin_satisfies_runtime_contract() -> None:
    plugin = validate_plugin(GenericOpenFOAMPlugin())
    ctx = driver_context(plugin, source="test")

    assert ctx.identity.id == "org.driverfoam.generic-openfoam"
    assert plugin.get_dict_entries() == ()
    assert plugin.get_tutorial_catalog() == {"registered_tutorials": (), "spec_factories": {}}


def test_minimal_plugin_proves_non_cardiac_solver_contract() -> None:
    plugin = validate_plugin(MinimalOpenFOAMPlugin())
    ctx = driver_context(plugin, source="test")

    assert ctx.identity.id == "org.driverfoam.test-minimal"
    assert plugin.get_dict_entries() == ()
    assert plugin.get_capabilities() == {}
