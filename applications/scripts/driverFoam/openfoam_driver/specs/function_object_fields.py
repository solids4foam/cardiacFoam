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
#     function_object_fields
#
# Description
#     Plan-time validation of the field names sampled by controlDict function
#     objects. Function objects themselves are OpenFOAM's (not driverFOAM's),
#     but the *fields* they can sample are this solver's — a name the solver
#     does not expose is silently dropped at solve time. This module turns that
#     silent drop into a non-blocking warning during strict planning.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

from __future__ import annotations

import os
import re
from pathlib import Path
from typing import Mapping, Sequence

from ..planning_types import StrictDiagnostic, diagnostic as _diagnostic

# controlDict locations to scan: top-level and the electro sub-region used by
# multi-region electromechanical cases.
_CONTROLDICT_RELPATHS = ("system/controlDict", "system/electro/controlDict")

_REGION_RE = re.compile(r"\bregion\s+(\w+)\s*;")
_FIELDS_RE = re.compile(r"\bfields\s*\(([^)]*)\)", re.DOTALL)


def _strip_comments(text: str) -> str:
    text = re.sub(r"/\*.*?\*/", " ", text, flags=re.DOTALL)
    text = re.sub(r"//[^\n]*", " ", text)
    return text


def _balanced_block(text: str, open_index: int) -> tuple[str, int] | None:
    """Return (inner_text, index_after_close) for the brace opening at
    ``open_index`` (which must point at ``{``). ``None`` if unbalanced."""
    depth = 0
    for i in range(open_index, len(text)):
        char = text[i]
        if char == "{":
            depth += 1
        elif char == "}":
            depth -= 1
            if depth == 0:
                return text[open_index + 1 : i], i + 1
    return None


def _functions_block(text: str) -> str | None:
    match = re.search(r"\bfunctions\b\s*", text)
    if match is None:
        return None
    brace = text.find("{", match.end())
    if brace == -1:
        return None
    result = _balanced_block(text, brace)
    return result[0] if result else None


def _iter_subdicts(block: str):
    """Yield the inner text of each top-level ``name { ... }`` sub-dictionary in
    a functions block. ``#includeFunc`` lines carry no braces and are skipped."""
    i = 0
    length = len(block)
    while i < length:
        brace = block.find("{", i)
        if brace == -1:
            return
        result = _balanced_block(block, brace)
        if result is None:
            return
        inner, after = result
        yield inner
        i = after


def _sampled_fields(subdict: str) -> list[str]:
    fields: list[str] = []
    for match in _FIELDS_RE.finditer(subdict):
        fields.extend(match.group(1).split())
    return fields


def _region_of(subdict: str) -> str:
    match = _REGION_RE.search(subdict)
    return match.group(1) if match else "electro"


def _diagnostics_for_text(
    text: str, samplable: Mapping[str, set[str]], source: str
) -> list[StrictDiagnostic]:
    block = _functions_block(_strip_comments(text))
    if block is None:
        return []
    diagnostics: list[StrictDiagnostic] = []
    for subdict in _iter_subdicts(block):
        region = _region_of(subdict)
        if region not in samplable:
            # A region this plugin's samplable-fields map doesn't name
            # (e.g. a not-yet-cataloged bath/torso domain) -- skip rather
            # than guess. Forcing it into "electro" would fabricate a
            # warning for a region electro never claimed to cover; the
            # same "never surface a spurious warning from a limitation"
            # rule this module already applies to parse/IO failures.
            continue
        allowed = samplable[region]
        for field_name in _sampled_fields(subdict):
            if field_name not in allowed:
                diagnostics.append(
                    _diagnostic(
                        "warning",
                        "unknown_sampled_field",
                        (
                            f"Function object samples field {field_name!r} "
                            f"(region {region!r}) which the resolved model does "
                            "not expose; the solver will silently drop it."
                        ),
                        source=source,
                        field=field_name,
                    )
                )
    return diagnostics


def function_object_field_diagnostics(
    case_root: str | Path,
    *,
    samplable: Mapping[str, Sequence[str] | set[str]],
) -> tuple[StrictDiagnostic, ...]:
    """Warn (never error) about controlDict function objects sampling fields
    absent from ``samplable`` (an open ``{region_name: {field, ...}}`` map,
    typically from :func:`capability_manifest.build_capability_manifest` --
    the built-in cardiac plugin currently declares ``"electro"`` and
    ``"solid"``, but core imposes no fixed key set). A function object whose
    ``region`` isn't a key in ``samplable`` at all is skipped rather than
    checked against a guessed bucket -- see :func:`_diagnostics_for_text`.

    Degrades to silence on any parse or IO failure — a parser limitation must
    never surface as a spurious field warning. Honors
    ``SKIP_FUNCTION_OBJECT_DIAGNOSTICS`` to bypass the check entirely.
    """
    if os.environ.get("SKIP_FUNCTION_OBJECT_DIAGNOSTICS"):
        return ()
    normalized = {region: set(fields) for region, fields in samplable.items()}
    root = Path(case_root)
    diagnostics: list[StrictDiagnostic] = []
    for relpath in _CONTROLDICT_RELPATHS:
        path = root / relpath
        try:
            if not path.is_file():
                continue
            text = path.read_text()
        except OSError:
            continue
        try:
            diagnostics.extend(_diagnostics_for_text(text, normalized, str(path)))
        except Exception:
            # A parse failure must not fabricate a warning.
            continue
    return tuple(diagnostics)
