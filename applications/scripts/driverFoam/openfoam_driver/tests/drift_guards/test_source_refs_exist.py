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
#     test_source_refs_exist
#
# Description
#     Drift guard: every source_refs path in every DictEntry must resolve to
#     a real file in the monorepo src/ tree. Fails when a C++ source file is
#     renamed or deleted without updating the catalogue.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

"""Drift guard: catalogue ``source_refs`` must resolve to real files on disk.

Every ``DictEntry`` in the driverFOAM catalogue carries a ``source_refs``
tuple of relative paths that document *where in the C++ source tree* the
dictionary key is actually read. This test verifies that every path in every
``source_refs`` tuple points to a file that actually exists in the monorepo.

Why this matters
----------------
When a C++ file is renamed, refactored, or deleted, the catalogue entry that
references it silently becomes stale. The old path stays in ``source_refs``
pointing at nothing, and driverFOAM may open it to extract validation rules —
crashing at runtime instead of at CI time.

This is precisely what happened with the ``eikonalSolver`` refactor: the old
``myocardiumModels/eikonalSolver/eikonalSolver.C`` was deleted and a new path
(``eikonalMyocardiumDomain.C``) took over, but the catalogue reference to the
old file was not removed.  This test would have caught that the moment the
deletion PR landed.

What is checked
---------------
* ``PHYSICS_PROPERTY_ENTRIES`` (physicsModel / controlDict section)
* ``CONTROL_DICT_ENTRIES`` (controlDict section)
* All groups returned by ``get_electro_property_entry_groups()`` (the main
  ``electroProperties`` catalogue)

What is NOT checked
-------------------
* ``source_refs`` that are empty tuples — no path, nothing to validate.
* Paths inside ``notes`` or ``description`` fields (free text, not a contract).

Failure message
---------------
The test lists every missing file together with the ``driver_path`` of the
entry that references it, so the fix is unambiguous: either update the path
in the catalogue to the new file location, or remove the stale reference.
"""

from __future__ import annotations

import pytest
from pathlib import Path

from openfoam_driver.dict_entries import (
    get_electro_property_entry_groups,
    PHYSICS_PROPERTY_ENTRIES,
)
from openfoam_driver.plugins.cardiacfoam.common_dict_entries import (
    CONTROL_DICT_ENTRIES,
)
from openfoam_driver.tests.conftest import monorepo_root

# ---------------------------------------------------------------------------
# Module-level skip (inherited from drift_guards/conftest.py pytestmark, but
# also expressed explicitly here so the skip reason is self-contained when
# this file is run in isolation).
# ---------------------------------------------------------------------------

pytestmark = pytest.mark.skipif(
    monorepo_root is None,
    reason=(
        "Requires the full cardiacFoam monorepo tree (src/). "
        "Clone the full repository to enable drift guard tests."
    ),
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _collect_broken_refs(repo_root: Path) -> list[tuple[str, str]]:
    """Return ``[(missing_path, driver_path), ...]`` for every stale ref."""
    broken: list[tuple[str, str]] = []

    def _check_group(entries):
        for entry in entries:
            for ref in entry.source_refs:
                if not (repo_root / ref).is_file():
                    broken.append((ref, entry.driver_path))

    _check_group(PHYSICS_PROPERTY_ENTRIES)
    _check_group(CONTROL_DICT_ENTRIES)
    for group_entries in get_electro_property_entry_groups().values():
        _check_group(group_entries)

    return broken


def _collect_duplicate_refs_within_entry() -> list[tuple[str, str, int]]:
    """Return ``[(driver_path, ref, count), ...]`` for in-tuple duplicates."""
    from collections import Counter

    dups: list[tuple[str, str, int]] = []

    def _check_group(entries):
        for entry in entries:
            c = Counter(entry.source_refs)
            for ref, n in c.items():
                if n > 1:
                    dups.append((entry.driver_path, ref, n))

    _check_group(PHYSICS_PROPERTY_ENTRIES)
    _check_group(CONTROL_DICT_ENTRIES)
    for group_entries in get_electro_property_entry_groups().values():
        _check_group(group_entries)

    return dups


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------


def test_all_source_refs_resolve_to_real_files():
    """Every path in every DictEntry.source_refs must exist on disk.

    Failure means a C++ file was renamed or deleted without updating the
    catalogue. Fix: update the stale path to the new file location, or remove
    the reference if the file no longer exists.

    The failure message lists every missing path and the catalogue entry that
    references it.
    """
    assert monorepo_root is not None  # guard; skip should have fired first
    broken = _collect_broken_refs(monorepo_root)

    if not broken:
        return  # all good

    lines = [
        "Catalogue source_refs point to files that do not exist on disk.",
        "Update the path in the catalogue to the new location, or remove it.",
        "",
    ]
    for ref, driver_path in sorted(broken):
        lines.append(f"  MISSING  {ref}")
        lines.append(f"    referenced by driver_path: {driver_path!r}")

    pytest.fail("\n".join(lines))


def test_no_duplicate_source_refs_within_a_single_entry():
    """No DictEntry should list the same file twice in its source_refs tuple.

    Duplicate refs are harmless but indicate a copy-paste error and make the
    catalogue harder to audit. Remove the duplicate occurrence.
    """
    dups = _collect_duplicate_refs_within_entry()

    if not dups:
        return  # all good

    lines = [
        "DictEntry.source_refs tuples contain duplicate paths:",
        "",
    ]
    for driver_path, ref, count in sorted(dups):
        lines.append(f"  driver_path: {driver_path!r}")
        lines.append(f"    duplicated {count}x: {ref}")

    pytest.fail("\n".join(lines))
