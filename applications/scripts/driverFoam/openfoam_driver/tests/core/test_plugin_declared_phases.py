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


def test_entries_in_a_plugin_declared_phase_are_still_required_checked() -> None:
    """The silent defect: a non-cardiac phase made primary_phase() return None,
    and both the required-field and enum checks did `continue`."""
    from openfoam_driver.core.specs.validation import primary_phase

    entry = _PhasedPlugin._ENTRIES[0]
    assert primary_phase(entry, ("mesh", "material", "loading", "solver")) == "material"


def test_an_entry_outside_the_declared_order_is_reported_not_skipped() -> None:
    from openfoam_driver.core.specs.validation import primary_phase

    entry = _PhasedPlugin._ENTRIES[0]
    assert primary_phase(entry, ("anatomy", "physics")) is None


class _MisdeclaredPhasePlugin(_PhasedPlugin):
    """Declares a phase order that does not cover its own entry's phase.

    This is the exact shape that used to vanish: `primary_phase()` returned
    None and the required-field check did `continue`, so a REQUIRED entry was
    never checked and nothing was reported.
    """

    @property
    def plugin_id(self) -> str:
        return "org.test.misdeclared"

    def get_phases(self) -> tuple[str, ...]:
        return ("mesh", "loading")  # note: no "material"


def test_an_unvalidatable_entry_is_reported_rather_than_silently_skipped() -> None:
    from openfoam_driver.core.runtime.run_model import RunDocument
    from openfoam_driver.core.specs.validation import validate_run

    ctx = driver_context(_MisdeclaredPhasePlugin(), source="test")
    run = RunDocument(
        id="t", name="t", status="draft", config={"mesh": {}, "loading": {}},
    )

    errors = validate_run(run, driver_context=ctx)
    messages = [e.message for e in errors]

    assert any("solidProperties.rho" in m for m in messages), (
        "a required entry whose phase is outside the declared order was not "
        f"reported at all -- this is the silent skip. Got: {messages}"
    )
    assert any("cannot be validated" in m for m in messages), messages
