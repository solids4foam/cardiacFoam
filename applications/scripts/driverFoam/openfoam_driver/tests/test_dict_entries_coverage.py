#----------------------------------------------------------------------------#
# License
#     This file is part of cardiacFoam.
#
#     cardiacFoam is free software: you can redistribute it and/or modify it
#     under the terms of the GNU General Public License as published by the
#     Free Software Foundation, either version 3 of the License, or (at your
#     option) any later version.
#
#     cardiacFoam is distributed in the hope that it will be useful, but
#     WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#     General Public License for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with cardiacFoam.  If not, see <http://www.gnu.org/licenses/>.
#
# Module
#     test_dict_entries_coverage
#
# Description
#     Tests dict entries coverage logic and specification contracts.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

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
