"""A plugin owns its dictionary phase vocabulary, not just its assignments."""
from __future__ import annotations

from openfoam_driver.core.contracts.dictionary import DictEntry
from openfoam_driver.core.contracts.dictionary_catalog import DictionaryCatalog
from openfoam_driver.core.plugin_interface import driver_context
from openfoam_driver.plugins.cardiacfoam_plugin import CardiacFoamPlugin
from openfoam_driver.tests.plugins.minimal_plugin import MinimalOpenFOAMPlugin


class _PhasedPlugin(MinimalOpenFOAMPlugin):
    """A plugin with a non-cardiac phase vocabulary."""

    _ENTRIES = (
        DictEntry(
            driver_path="solidProperties.rho",
            description="Density",
            required=True,
            phases=frozenset({"material"}),
        ),
    )

    @property
    def plugin_id(self) -> str:
        return "org.test.phased"

    def get_phases(self) -> tuple[str, ...]:
        return ("mesh", "material", "loading", "solver")

    def get_dict_entries(self):
        return self._ENTRIES

    def get_dictionary_catalog(self):
        return DictionaryCatalog({"solidProperties": self._ENTRIES})

    def get_dict_groups(self):
        return {"material": self._ENTRIES}


def test_cardiac_phases_are_unchanged() -> None:
    ctx = driver_context(CardiacFoamPlugin(), source="test")
    assert ctx.capabilities.dictionaries.phases() == (
        "anatomy", "physics", "stimulus", "solver",
    )


def test_a_plugin_may_declare_its_own_phase_vocabulary() -> None:
    ctx = driver_context(_PhasedPlugin(), source="test")
    assert ctx.capabilities.dictionaries.phases() == (
        "mesh", "material", "loading", "solver",
    )


def test_a_plugin_without_the_hook_gets_the_phases_its_entries_declare() -> None:
    ctx = driver_context(MinimalOpenFOAMPlugin(), source="test")
    # No entries, no hook -> no phases, and specifically not cardiac's four.
    assert ctx.capabilities.dictionaries.phases() == ()
