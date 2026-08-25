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
#     case_dict_keys
#
# Description
#     Warns about case-dictionary keys absent from the plugin catalogue.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

"""Warn about case-dictionary keys the active plugin's catalogue does not know.

This is the opposite direction to ``scripts/_dict_keys_scanner.py``. That one
asks *what keys does the C++ accept?* -- information that lives only in C++
source, because a defaulted read nobody sets appears in no dictionary at all.
This one asks *what keys did the author actually write?* -- information that
lives only in the case file, because driverFOAM's builder only ever writes
keys it already knows and the solver only ever asks "is key X present?".

Nothing else closes this gap. OpenFOAM ignores unrecognised keys by design, so
a misspelled optional key is silently dropped: misspelling
``activeTensionModel`` in singleCell yields exit 0, a clean solver log, and an
output set quietly missing the active-tension trace. driverFOAM's artifact
reconciliation cannot catch it either -- it reads the same misspelled file, so
it never expects the missing artifact and the two agree on the wrong thing.

**Warn, never error.** cardiacFoam does not own every key that may legitimately
appear in these dictionaries: OpenFOAM's own machinery reads some of them, and
the catalogue deliberately documents only keys subject to programmatic
override. An unmatched key is a question for a human, not a failed plan.
"""

from __future__ import annotations

import os
import re
from collections.abc import Iterable, Iterator, Mapping, Sequence
from pathlib import Path

from openfoam_driver.core.planning_types import StrictDiagnostic, diagnostic

_WILDCARD = re.compile(r"<[^>]+>")

SKIP_ENV_VAR = "SKIP_CASE_DICT_KEY_DIAGNOSTICS"

# OpenFOAM's runtime-selection convention: a model selected by <key> reads its
# settings from a sibling <modelName>Coeffs sub-dictionary. Catalogue paths
# address the inside of that dictionary through a stripped $SCOPE_TOKEN.
# prefix, so the literal name never enters the known set. It belongs to
# OpenFOAM, not to any plugin's catalogue, and warning about it would fire on
# every case. A misspelling ("...Coefs") still fails this test and is reported.
_RTS_COEFFS_SUFFIX = "Coeffs"


def _scope_relative(trail: tuple[str, ...]) -> tuple[str, ...]:
    """Drop a leading runtime-selection scope dict from a trail.

    Catalogue paths carry a ``$SCOPE_TOKEN.`` prefix that parsing strips, so
    they are expressed *relative to the inside* of the ``<model>Coeffs``
    dictionary. A trail walked from the file root is absolute. Align them, or
    every single key inside the coeffs dict reports as unknown.
    """
    if trail and trail[0].endswith(_RTS_COEFFS_SUFFIX):
        return trail[1:]
    return trail


def _prefixes(catalogued_paths: Iterable[str]) -> set[tuple[str, ...]]:
    """Every catalogue path and every prefix of one, as segment tuples.

    Prefixes matter because a container (``ecgDomains``) is legitimate even
    though no catalogue path ends there.
    """
    out: set[tuple[str, ...]] = set()
    for path in catalogued_paths:
        segments = tuple(path.split("."))
        for i in range(1, len(segments) + 1):
            out.add(segments[:i])
    return out


def _matches(trail: tuple[str, ...], known: set[tuple[str, ...]]) -> bool:
    """Does ``trail`` match a catalogue prefix, ``<placeholder>`` matching any?"""
    for candidate in known:
        if len(candidate) != len(trail):
            continue
        if all(
            _WILDCARD.fullmatch(c) or c == t
            for c, t in zip(candidate, trail)
        ):
            return True
    return False


def case_dict_key_diagnostics(
    case_root: str | Path,
    *,
    catalogued_paths: Iterable[str],
    dict_relpaths: Sequence[str],
) -> tuple[StrictDiagnostic, ...]:
    """Warn (never error) about keys in ``dict_relpaths`` absent from the catalogue.

    ``catalogued_paths`` are scope-stripped catalogue ``driver_path`` strings.
    Matching is by *position*, not by bare name: a key is known when its full
    trail matches a catalogue path or a prefix of one, with ``<placeholder>``
    segments matching any single name. That is what distinguishes an author's
    own instance label (``ecgDomains { ECG { ... } }``, where the catalogue
    says ``ecgDomains.<name>...``) from a genuine misspelling -- a flat set of
    names cannot tell them apart, and would warn on every case that names an
    ECG domain.

    Degrades to silence on any parse or IO failure: a parser limitation must
    never surface as a spurious key warning. Honors ``SKIP_ENV_VAR``.
    """
    if os.environ.get(SKIP_ENV_VAR):
        return ()

    known = _prefixes(catalogued_paths)
    root = Path(case_root)
    diagnostics: list[StrictDiagnostic] = []

    for relpath in dict_relpaths:
        path = root / relpath
        if not path.is_file():
            continue
        try:
            from foamlib import FoamFile

            parsed = FoamFile(path)
            unmatched = _unmatched(parsed, known)
        except Exception:
            continue
        for trail in unmatched:
            where = ".".join(trail)
            diagnostics.append(
                diagnostic(
                    "warning",
                    "uncatalogued_case_dict_key",
                    (
                        f"{relpath}: key {where!r} is not in the active "
                        "plugin's dictionary catalogue. OpenFOAM ignores "
                        "unrecognised keys, so if this is a misspelling the "
                        "solver will silently fall back to its default."
                    ),
                    source=relpath,
                    field=trail[-1],
                )
            )
    return tuple(diagnostics)


def _unmatched(
    node: Mapping,
    known: set[tuple[str, ...]],
    trail: tuple[str, ...] = (),
) -> list[tuple[str, ...]]:
    """Outermost unmatched keys, depth-first.

    An unmatched container is reported once and NOT descended into: if
    ``singleCellSolverCoefs`` is misspelled then every key beneath it is
    unreachable too, and reporting all of them buries the one that matters.
    """
    found: list[tuple[str, ...]] = []
    for key in node:
        try:
            value = node[key]
        except Exception:
            continue
        full = trail + (str(key),)
        if not _matches(_scope_relative(full), known) and not str(key).endswith(
            _RTS_COEFFS_SUFFIX
        ):
            found.append(full)
            continue
        if hasattr(value, "keys"):
            found.extend(_unmatched(value, known, full))
    return found
