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
#     _rtst_scanner
#
# Description
#     Scans runtime selection tables for available implementations.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

"""Scanner for OpenFOAM runtime-selection-table registrations.

For each `addToRunTimeSelectionTable(base, derived, ctor)` call in
`src/**/*.C` (skipping `lnInclude/` and `Make/`), this resolves the registered
type name used as the dictionary value (which may differ from the C++ class
name via `OverrideTypeName(...)`) by reading the derived class's header.

Used by the RTST contract test to verify that catalogue `enum_values` match
the names actually selectable at runtime.
"""

from __future__ import annotations

import re
from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path

from openfoam_driver.dict_entries import (
    get_electro_property_entry_groups,
    PHYSICS_PROPERTY_ENTRIES,
    DictEntry,
)


# `addToRunTimeSelectionTable(base, derived, ctor);` — argument layout is
# whitespace-tolerant including newlines.
_RTST_RE = re.compile(
    r"addToRunTimeSelectionTable\s*\(\s*"
    r"(?P<base>[A-Za-z_]\w*)\s*,\s*"
    r"(?P<derived>[A-Za-z_]\w*)\s*,\s*"
    r"(?P<ctor>[A-Za-z_]\w*)\s*\)",
    re.DOTALL,
)

_OVERRIDE_TYPENAME_RE = re.compile(
    r'OverrideTypeName\s*\(\s*"([^"]+)"\s*\)'
)

# `class Foo`, `class Foo final`, `class Foo : public Bar`, etc.
# Skip forward declarations (`class Foo;` with no `{`).
_CLASS_DECL_RE = re.compile(
    r"\bclass\s+([A-Za-z_]\w*)\b(?![^;{]*;)"
)


@dataclass(frozen=True)
class RtstReg:
    base: str
    registered_name: str
    derived_class: str
    source_file: Path


def _iter_src_files(src_root: Path, suffix: str) -> Iterable[Path]:
    for path in src_root.rglob(f"*{suffix}"):
        parts = path.parts
        if "lnInclude" in parts or "Make" in parts:
            continue
        yield path


def _index_override_type_names(src_root: Path) -> dict[str, str]:
    """Map *C++ class name* -> registered name (via `OverrideTypeName`).

    Each header is scanned for `class Foo { ... OverrideTypeName("X") ... }`
    patterns; each `OverrideTypeName` is paired with the nearest preceding
    `class <Name>` declaration in the same file. This handles the general
    case where a model's C++ class name could differ from its registered
    type string, and works uniformly for the common case (`BuenoOrovio`,
    `ORd`, `PerisYague`, …) where the class name and the override name are
    the same.

    Classes without an `OverrideTypeName` simply do not appear in the
    index; the caller falls back to the class name.
    """
    index: dict[str, str] = {}
    for header in _iter_src_files(src_root, ".H"):
        text = header.read_text(encoding="utf-8", errors="replace")
        # Strip comments first so commented-out code never matches.
        text = re.sub(r"/\*.*?\*/", "", text, flags=re.DOTALL)
        text = re.sub(r"//[^\n]*", "", text)
        class_positions: list[tuple[int, str]] = [
            (m.start(), m.group(1)) for m in _CLASS_DECL_RE.finditer(text)
        ]
        if not class_positions:
            continue
        for override in _OVERRIDE_TYPENAME_RE.finditer(text):
            pos = override.start()
            # Nearest preceding class declaration.
            enclosing = None
            for cls_pos, cls_name in class_positions:
                if cls_pos < pos:
                    enclosing = cls_name
                else:
                    break
            if enclosing is not None:
                index[enclosing] = override.group(1)
    return index


def scan_rtst_registrations(
    src_root: Path,
) -> dict[str, dict[str, RtstReg]]:
    """Return `{base_class: {registered_name: RtstReg}}` for every RTST
    registration found under `src_root`.
    """
    override_map = _index_override_type_names(src_root)
    out: dict[str, dict[str, RtstReg]] = {}
    for source in _iter_src_files(src_root, ".C"):
        text = source.read_text(encoding="utf-8", errors="replace")
        # Strip block comments so we don't pick up commented-out
        # registrations.
        text = re.sub(r"/\*.*?\*/", "", text, flags=re.DOTALL)
        text = re.sub(r"//[^\n]*", "", text)
        for match in _RTST_RE.finditer(text):
            base = match.group("base")
            derived = match.group("derived")
            registered = override_map.get(derived, derived)
            out.setdefault(base, {})[registered] = RtstReg(
                base=base,
                registered_name=registered,
                derived_class=derived,
                source_file=source,
            )
    return out


# ---------------------------------------------------------------------------
# Catalogue side


@dataclass(frozen=True)
class CatalogueEnum:
    driver_path: str
    enum_values: tuple[str, ...]
    source_refs: tuple[str, ...]


def iter_catalogue_enums() -> Iterable[CatalogueEnum]:
    """Yield every `DictEntry` with `value_kind == "enum"` and a non-empty
    `enum_values` tuple, across both the physics and electro catalogues.
    """
    def from_entry(e: DictEntry) -> CatalogueEnum | None:
        if e.value_kind != "enum" or not e.enum_values:
            return None
        return CatalogueEnum(
            driver_path=e.driver_path,
            enum_values=e.enum_values,
            source_refs=e.source_refs,
        )

    for e in PHYSICS_PROPERTY_ENTRIES:
        c = from_entry(e)
        if c is not None:
            yield c
    for group in get_electro_property_entry_groups().values():
        for e in group:
            c = from_entry(e)
            if c is not None:
                yield c
