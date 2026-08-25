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
#     test_tutorial_keys_are_catalog_addressable
#
# Description
#     Every key a tracked tutorial actually sets must be settable through the
#     driver, at the scope the tutorial writes it.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

"""Catalog coverage measured against real tutorial dicts.

The C++ key scanner (``_dict_keys_scanner``) compares *bare leaf names*, so a
key catalogued at one scope masks its absence at another: ``solver`` counted
as present because a top-level entry existed, while
``…purkinjeGraphModelCoeffs.solver`` was unsettable. Making that scanner
path-aware would need interprocedural analysis -- measured on this tree, only
4 of 340 C++ reads have a scope resolvable from a local
``member_(parent.subDict("X"))`` binding; the rest read from a dict handed in
by a caller.

This test closes the same gap from the data side instead. It walks the
tracked tutorial dicts, reconstructs the full scope path of every entry, and
asserts the driver can address it. That catches a key catalogued at the wrong
scope, a key added to a tutorial with no catalog entry, and a catalog entry
deleted while tutorials still set it -- with no C++ parsing at all.

Keys that are genuinely unaddressable are waived below *with a reason*, so
the set stays an explicit, reviewed list rather than silent drift.
"""

from __future__ import annotations

import re
import subprocess
from pathlib import Path

import pytest

from openfoam_driver.core.specs.apply_overrides import OverrideError, validate_overrides
from openfoam_driver.tests.conftest import monorepo_root, skip_without_monorepo

pytestmark = skip_without_monorepo

#: ``<solver>Coeffs`` blocks all map onto the $ELECTRO_MODEL_COEFFS token.
_COEFFS_BLOCK = re.compile(
    r"^(monodomain|bidomain|singleCell|eikonal|monodomain1D)Solver(1D)?Coeffs$"
)
_ENTRY = re.compile(r"([A-Za-z_]\w*)\s+[^{};]*;\s*$")
_BLOCK_HEADER = re.compile(r"([A-Za-z_][\w.|\"]*)\s*$")

#: Unaddressable on purpose. Each entry needs a reason; delete the entry when
#: the reason stops holding, rather than widening the waiver.
WAIVED: dict[str, str] = {
    # Empty, and worth keeping that way. Every key a tracked tutorial sets is
    # a key something reads. Dead entries were deleted from the dicts rather
    # than waived here -- a waiver hides drift, deleting the key removes it.
}


#: Sub-dict names that are user-chosen instance names, not fixed keys.
_DYNAMIC_PARENTS = ("conductionNetworkDomains", "domainCouplings", "ecgDomains")


#: ``/* ... */`` spans. tutorials/template/constant/electroProperties keeps a
#: large commented-out example of an ecgDomains block; parsing it as live
#: content invents entries that no case actually sets.
_BLOCK_COMMENT = re.compile(r"/\*.*?\*/", re.DOTALL)


def _scoped_keys(path: Path) -> set[str]:
    """Full driver_path of every ``key value;`` entry in one dict file."""
    found: set[str] = set()
    stack: list[str | None] = []
    pending: str | None = None
    text = _BLOCK_COMMENT.sub("", path.read_text(errors="ignore"))
    for line in text.splitlines():
        code = line.split("//", 1)[0]
        stripped = code.strip()

        entry = _ENTRY.match(stripped)
        if entry and stack:
            segments = [s for s in stack if s]
            if segments and _COEFFS_BLOCK.match(segments[0]):
                segments = ["$ELECTRO_MODEL_COEFFS"] + segments[1:]
                # user-named instances become the catalog's <name> wildcard
                segments = [
                    "<name>" if i and segments[i - 1] in _DYNAMIC_PARENTS else s
                    for i, s in enumerate(segments)
                ]
                found.add(".".join(segments + [entry.group(1)]))

        header = _BLOCK_HEADER.match(stripped)
        if header:
            pending = header.group(1)
        for char in code:
            if char == "{":
                stack.append(pending)
                pending = None
            elif char == "}":
                if stack:
                    stack.pop()
    return found


def _tracked_electro_dicts() -> list[Path]:
    listing = subprocess.run(
        ["git", "ls-files", "tutorials"],
        cwd=monorepo_root, capture_output=True, text=True, check=True,
    ).stdout.split()
    return [
        monorepo_root / f for f in listing if "electroProperties" in Path(f).name
    ]


def test_every_tutorial_dict_key_is_catalog_addressable() -> None:
    keys: set[str] = set()
    for dict_file in _tracked_electro_dicts():
        keys |= _scoped_keys(dict_file)

    assert keys, "found no tutorial dict entries -- the parser is broken"

    unaddressable = []
    for driver_path in sorted(keys):
        concrete = driver_path.replace("<name>", "X")
        try:
            validate_overrides([{"driver_path": concrete, "value": "1"}])
        except OverrideError as exc:
            if "not catalog-addressable" in str(exc):
                unaddressable.append(driver_path)
        except Exception:  # enum/type rejection still means it is addressable
            pass

    unexpected = [p for p in unaddressable if p not in WAIVED]
    stale_waivers = [p for p in WAIVED if p not in unaddressable]

    assert not unexpected, (
        "tutorial keys the driver cannot set:\n  "
        + "\n  ".join(unexpected)
        + "\nAdd a catalog entry at that scope, or waive it with a reason."
    )
    assert not stale_waivers, (
        "these waivers are no longer needed -- delete them:\n  "
        + "\n  ".join(stale_waivers)
    )
