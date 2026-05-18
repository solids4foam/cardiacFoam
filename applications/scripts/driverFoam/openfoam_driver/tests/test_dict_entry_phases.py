"""Tests for the ``Phase`` literal and multi-phase ``DictEntry.phases`` field.

Covers the DictEntry dataclass contract that lets a single entry belong to more
than one workflow phase.
"""

from __future__ import annotations

import typing

from openfoam_driver.dict_entries import DictEntry, Phase

VALID_PHASES = {"anatomy", "physics", "stimulus", "solver"}


def test_dict_entry_has_phases_field_accepting_a_frozenset():
    entry = DictEntry(
        driver_path="foo",
        description="x",
        source_refs=("bar",),
        phases=frozenset({"physics"}),
    )
    assert entry.phases == frozenset({"physics"})


def test_dict_entry_phases_supports_multi_phase_ownership():
    entry = DictEntry(
        driver_path="nRegions",
        description="number of regions",
        source_refs=("bar",),
        phases=frozenset({"anatomy", "solver"}),
    )
    assert entry.phases == frozenset({"anatomy", "solver"})


def test_dict_entry_phases_default_is_empty_frozenset():
    entry = DictEntry(driver_path="foo", description="x", source_refs=("bar",))
    assert entry.phases == frozenset()


def test_phase_literal_values():
    args = typing.get_args(Phase)
    assert set(args) == VALID_PHASES
