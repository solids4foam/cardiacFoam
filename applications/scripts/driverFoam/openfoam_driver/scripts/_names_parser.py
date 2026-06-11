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
#     _names_parser
#
# Description
#     Parses string identifiers and names from definitions.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

"""Parser for ionic-model `*_Names.H` headers.

Shared by the regenerate-ionic-catalog script and the catalogue contract test
so both sides interpret the headers identically.

Returns the names that appear in each of the three enums (states, algebraic,
constants), with the `NUM_*` sentinel stripped. Identifiers are returned in
declaration order and with any `= <value>` initialiser removed.

The variant enum names observed in this repo are matched by case-insensitive
substring:
    states     -> any enum whose name contains "STATE"
    algebraic  -> any enum whose name contains "ALGEBRAIC"
    constants  -> any enum whose name contains "CONSTANT"
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class ParsedNames:
    states: tuple[str, ...]
    algebraic: tuple[str, ...]
    constants: tuple[str, ...]


_CATEGORY_KEYWORDS = {
    "states": "STATE",
    "algebraic": "ALGEBRAIC",
    "constants": "CONSTANT",
}

_BLOCK_COMMENT = re.compile(r"/\*.*?\*/", re.DOTALL)
_LINE_COMMENT = re.compile(r"//[^\n]*")
_ENUM_HEADER = re.compile(r"enum\s+([A-Za-z_][A-Za-z0-9_]*)\s*\{", re.MULTILINE)


def _strip_comments(text: str) -> str:
    text = _BLOCK_COMMENT.sub("", text)
    text = _LINE_COMMENT.sub("", text)
    return text


def _is_sentinel(identifier: str) -> bool:
    """Return True for `NUM_STATES`, `BATH_BIDOMAIN_NUM_STATES`, etc."""
    upper = identifier.upper()
    return upper.startswith("NUM_") or "_NUM_" in upper


def _parse_enum_body(body: str) -> tuple[str, ...]:
    body = _strip_comments(body)
    items: list[str] = []
    for raw in body.split(","):
        token = raw.strip()
        if not token:
            continue
        # Drop `= <value>` initialiser if present.
        if "=" in token:
            token = token.split("=", 1)[0].strip()
        if not token:
            continue
        if _is_sentinel(token):
            continue
        items.append(token)
    return tuple(items)


def _extract_enum_body(text: str, name_keyword: str) -> tuple[str, ...] | None:
    """Find the first enum whose name contains `name_keyword` (case-insensitive)
    and return its parsed identifier list. Returns None if no such enum exists.
    """
    text_nocomment = _strip_comments(text)
    keyword_upper = name_keyword.upper()
    for match in _ENUM_HEADER.finditer(text_nocomment):
        enum_name = match.group(1)
        if keyword_upper not in enum_name.upper():
            continue
        # Locate the matching closing brace from the position after `{`.
        start = match.end()
        depth = 1
        i = start
        while i < len(text_nocomment) and depth > 0:
            ch = text_nocomment[i]
            if ch == "{":
                depth += 1
            elif ch == "}":
                depth -= 1
            i += 1
        if depth != 0:
            return None
        body = text_nocomment[start : i - 1]
        return _parse_enum_body(body)
    return None


def parse_names_header(path: Path) -> ParsedNames:
    """Parse a `<Model>_<year>Names.H` file into (states, algebraic, constants).

    Missing categories are returned as empty tuples (some models have no
    algebraic or no constants list).
    """
    text = path.read_text(encoding="utf-8", errors="replace")
    states = _extract_enum_body(text, _CATEGORY_KEYWORDS["states"]) or ()
    algebraic = _extract_enum_body(text, _CATEGORY_KEYWORDS["algebraic"]) or ()
    constants = _extract_enum_body(text, _CATEGORY_KEYWORDS["constants"]) or ()
    return ParsedNames(states=states, algebraic=algebraic, constants=constants)


def find_names_header(model_dir: Path) -> Path | None:
    """Locate the `<Model>_<year>Names.H` file under a model directory.

    Returns the first match; returns None if the model has no Names header.
    Skips `lnInclude/` (symlink mirror) which would yield duplicates.
    """
    if not model_dir.is_dir():
        return None
    candidates = [
        p
        for p in model_dir.glob("*Names.H")
        if "lnInclude" not in p.parts
    ]
    return candidates[0] if candidates else None


# Models whose catalogue entries are intentionally user-facing semantic labels
# (e.g. monodomain manufactured maps a 4-state enum down to a single "phi"
# field, plus "Iion_cm"). These are NOT regenerated from headers and are NOT
# compared against headers by the contract test. The C++ enum members exist
# for internal indexing only.
EXCLUDED_FROM_HEADER_SYNC: frozenset[str] = frozenset({
    "monodomainFDAManufactured",
    "bidomainFDAManufactured",
    "bathBidomainFDAManufactured",
    # compactBatched GPU models live under <Model>Batched/ directories and
    # have no own *Names.H — they reuse the parent CPU model's header.
    "AlievPanfilovcompactBatched",
    "BuenoOroviocompactBatched",
    "CourtemanchecompactBatched",
    "FabbricompactBatched",
    "GaurcompactBatched",
    "GrandicompactBatched",
    "PerisYaguecompactBatched",
    "StewartcompactBatched",
    "TNNPcompactBatched",
    "ToRORd_dynClcompactBatched",
    "TrovatocompactBatched",
    "TWorldcompactBatched",
})
