"""Coverage tests for ``DictEntry.phases`` classification.

Enforces that every catalogued entry declares at least one workflow phase, and
that every declared value is a valid ``Phase`` literal.
"""

from __future__ import annotations

import typing

from openfoam_driver.dict_entries import (
    ELECTRO_PROPERTY_ENTRY_GROUPS,
    PHYSICS_PROPERTY_ENTRIES,
    Phase,
)

VALID_PHASES = set(typing.get_args(Phase))


def _all_entries():
    yield from PHYSICS_PROPERTY_ENTRIES
    for group in ELECTRO_PROPERTY_ENTRY_GROUPS.values():
        yield from group


def test_every_dict_entry_has_at_least_one_phase():
    unclassified = [e for e in _all_entries() if not e.phases]
    assert not unclassified, (
        f"{len(unclassified)} entries have no phases: "
        + ", ".join(e.driver_path for e in unclassified[:10])
    )


def test_every_phase_value_is_a_valid_literal():
    invalid = []
    for e in _all_entries():
        bad = [p for p in e.phases if p not in VALID_PHASES]
        if bad:
            invalid.append((e.driver_path, bad))
    assert not invalid, f"invalid phases: {invalid}"
