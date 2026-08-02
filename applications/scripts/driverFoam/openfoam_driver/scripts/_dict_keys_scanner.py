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
#     _dict_keys_scanner
#
# Description
#     Scans dictionary definitions for schema validation.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

"""Scanner for OpenFOAM dictionary-read call sites in C++ source.

For each `.C` / `.H` file under a given `src_root` (skipping `lnInclude/`,
`Make/`, and `*_Names.H` files), this module detects patterns of the form:

    <receiver>.lookup("key")
    <receiver>.lookupOrDefault<T>("key", default)
    <receiver>.get<T>("key")
    <receiver>.getOrDefault<T>("key", default)
    <receiver>.found("key")
    <receiver>.readEntry("key", out)
    <receiver>.subDict("name")
    <receiver>.subOrEmptyDict("name")
    <receiver>.optionalSubDict("name")
    readScalar(<receiver>.lookup("key"))
    readLabel(<receiver>.lookup("key"))
    readBool(<receiver>.lookup("key"))

Returns a flat list of `DictRead` records.  Sub-dict opens are flagged with
``kind="subdict"``; all other patterns use ``kind="key"``.

Comments are stripped before scanning so commented-out code is never matched.

Catalogue-side helpers (`CataloguePath`, `iter_catalogue_paths`) parse every
`driver_path` in `PHYSICS_PROPERTY_ENTRIES` and `get_electro_property_entry_groups()`
into a structured form for comparison against the scanner output.

Accuracy is ~80%; false positives/negatives are expected.  The output is for
human review only.
"""

from __future__ import annotations

import re
import json
from collections import defaultdict
from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path

from openfoam_driver.dict_entries import (
    get_electro_property_entry_groups,
    PHYSICS_PROPERTY_ENTRIES,
)


# ---------------------------------------------------------------------------
# Comment stripping  (identical pattern to _rtst_scanner.py)

_BLOCK_COMMENT = re.compile(r"/\*.*?\*/", re.DOTALL)
_LINE_COMMENT = re.compile(r"//[^\n]*")


def _strip_comments(text: str) -> str:
    text = _BLOCK_COMMENT.sub("", text)
    text = _LINE_COMMENT.sub("", text)
    return text


# ---------------------------------------------------------------------------
# Patterns for dictionary key reads
#
# Strategy: one combined regex with named groups.  The receiver identifier is
# captured but not used for path reconstruction (see module docstring).
#
# The string literal is always a double-quoted token without embedded quotes.
# We allow arbitrary whitespace (including newlines) between the method name,
# the opening paren, and the first argument.

_STRING_LIT = r'"([^"]+)"'

# Methods that read a *key* from a dictionary.
_KEY_METHOD = (
    r"(?:"
    r"lookupOrDefault(?:\s*<[^>]+>)?"      # lookupOrDefault<T>( or lookupOrDefault(
    r"|lookup"                               # lookup(
    r"|getOrDefault(?:\s*<[^>]+>)?"         # getOrDefault<T>(
    r"|get(?:\s*<[^>]+>)"                   # get<T>(   (require type param to avoid false-positives on e.g. get())
    r"|found"                               # found(
    r"|readEntry"                           # readEntry(
    r")"
)

# Methods that open a *sub-dictionary*.
_SUBDICT_METHOD = r"(?:subDict|subOrEmptyDict|optionalSubDict)"

# readScalar/readLabel/readBool wrapping a .lookup("key")
_WRAP_FUNC = r"(?:readScalar|readLabel|readBool)"

# Rather than one combined regex (hard to maintain), use three focused ones.
# Each captures the string literal as the *last* group in the pattern.

_KEY_RE = re.compile(
    r"[A-Za-z_]\w*(?:\s*\.\s*[A-Za-z_]\w*)*"
    r"\s*\.\s*" + _KEY_METHOD + r"\s*\(\s*" + _STRING_LIT,
    re.DOTALL,
)

_WRAP_RE = re.compile(
    r"(?:" + _WRAP_FUNC + r")\s*\(\s*"
    r"[A-Za-z_]\w*(?:\s*\.\s*[A-Za-z_]\w*)*"
    r"\s*\.\s*lookup\s*\(\s*" + _STRING_LIT,
    re.DOTALL,
)

_SUB_RE = re.compile(
    r"[A-Za-z_]\w*(?:\s*\.\s*[A-Za-z_]\w*)*"
    r"\s*\.\s*(?P<meth>" + _SUBDICT_METHOD + r")\s*\(\s*" + _STRING_LIT,
    re.DOTALL,
)


# ---------------------------------------------------------------------------
# Public dataclass

@dataclass(frozen=True)
class DictRead:
    kind: str       # "key" | "subdict"
    name: str       # the string literal
    source_file: Path
    line: int       # 1-based


# ---------------------------------------------------------------------------
# Scanner implementation

def _iter_src_files(src_root: Path) -> Iterable[Path]:
    for path in src_root.rglob("*"):
        if not path.is_file():
            continue
        if path.suffix not in {".C", ".H"}:
            continue
        parts = path.parts
        if "lnInclude" in parts or "Make" in parts:
            continue
        # Skip *_Names.H headers — they contain enum identifiers, not dict keys.
        if path.name.endswith("_Names.H") or path.name.endswith("Names.H"):
            continue
        yield path


def _line_of(text: str, pos: int) -> int:
    """Return 1-based line number for character position *pos* in *text*."""
    return text.count("\n", 0, pos) + 1


def scan_dict_reads(src_root: Path) -> list[DictRead]:
    """Return all dictionary-read sites found under *src_root*.

    Files in ``lnInclude/``, ``Make/``, and ``*Names.H`` are skipped.
    Comments are stripped before scanning.
    """
    results: list[DictRead] = []

    for source in _iter_src_files(src_root):
        raw = source.read_text(encoding="utf-8", errors="replace")
        text = _strip_comments(raw)

        # Key reads
        for m in _KEY_RE.finditer(text):
            key = m.group(m.lastindex)   # last capture group = the string literal
            results.append(
                DictRead(
                    kind="key",
                    name=key,
                    source_file=source,
                    line=_line_of(text, m.start()),
                )
            )

        # Wrapped reads: readScalar/readLabel/readBool(recv.lookup("key"))
        for m in _WRAP_RE.finditer(text):
            key = m.group(m.lastindex)
            results.append(
                DictRead(
                    kind="key",
                    name=key,
                    source_file=source,
                    line=_line_of(text, m.start()),
                )
            )

        # Sub-dict opens
        for m in _SUB_RE.finditer(text):
            key = m.group(m.lastindex)
            results.append(
                DictRead(
                    kind="subdict",
                    name=key,
                    source_file=source,
                    line=_line_of(text, m.start()),
                )
            )

    return results


# ---------------------------------------------------------------------------
# Catalogue-side helper

@dataclass(frozen=True)
class CataloguePath:
    driver_path: str            # original value from DictEntry
    normalised: str             # driver_path with $ELECTRO_MODEL_COEFFS. stripped
    leaf: str                   # last dot-segment
    parents: tuple[str, ...]    # all segments before the leaf
    has_wildcard: bool          # True if any segment matches <...>

    # Whether the entry is flagged as dynamic_path=True in the catalogue.
    dynamic_path: bool


@dataclass(frozen=True)
class DictKeyStrictReport:
    """Allowlist-backed catalogue drift report used by strict planning."""

    status: str
    absent_keys: tuple[str, ...]
    stale_paths: tuple[str, ...]
    unmatched_subdicts: tuple[str, ...]
    unused_allowlist: tuple[str, ...]

    def to_json(self) -> dict[str, object]:
        return {
            "status": self.status,
            "absent_keys": list(self.absent_keys),
            "stale_paths": list(self.stale_paths),
            "unmatched_subdicts": list(self.unmatched_subdicts),
            "unused_allowlist": list(self.unused_allowlist),
        }


_WILDCARD_RE = re.compile(r"<[^>]+>")
_PREFIX = "$ELECTRO_MODEL_COEFFS."


def _parse_path(driver_path: str, is_dynamic: bool) -> CataloguePath:
    # Strip the common prefix.
    if driver_path.startswith(_PREFIX):
        normalised = driver_path[len(_PREFIX):]
    else:
        normalised = driver_path

    segments = normalised.split(".")
    leaf = segments[-1]
    parents = tuple(segments[:-1])
    has_wildcard = any(_WILDCARD_RE.search(s) for s in segments)

    return CataloguePath(
        driver_path=driver_path,
        normalised=normalised,
        leaf=leaf,
        parents=parents,
        has_wildcard=has_wildcard,
        dynamic_path=is_dynamic,
    )


def iter_catalogue_paths() -> Iterable[CataloguePath]:
    """Yield a `CataloguePath` for every entry in the two catalogues."""
    for entry in PHYSICS_PROPERTY_ENTRIES:
        yield _parse_path(entry.driver_path, entry.dynamic_path)
    for group in get_electro_property_entry_groups().values():
        for entry in group:
            yield _parse_path(entry.driver_path, entry.dynamic_path)


IGNORED_FOAMFILE_KEYS: frozenset[str] = frozenset(
    {
        "version",
        "format",
        "class",
        "object",
        "location",
        "dimensions",
        "internalField",
        "boundaryField",
        "FoamFile",
        "note",
        "arch",
        "root",
        "case",
        "time",
        "path",
    }
)


def _default_allowlist_path() -> Path:
    return Path(__file__).with_name("dict_key_allowlist.json")


def load_dict_key_allowlist(path: Path | None = None) -> dict[str, set[str]]:
    """Load the reviewed strict-scanner allowlist.

    The file is intentionally JSON so the strict scanner can be used from both
    tests and the CLI without importing project-specific test fixtures.
    """
    allowlist_path = path or _default_allowlist_path()
    payload = json.loads(allowlist_path.read_text())
    return {
        "absent_keys": set(payload.get("absent_keys", [])),
        "stale_paths": set(payload.get("stale_paths", [])),
        "unmatched_subdicts": set(payload.get("unmatched_subdicts", [])),
    }


def compute_dict_key_drift(src_root: Path) -> dict[str, set[str]]:
    """Compute approximate C++ dictionary-reader drift against dict_entries."""
    reads = scan_dict_reads(src_root)
    cat_paths = list(iter_catalogue_paths())

    key_reads: dict[str, list[DictRead]] = defaultdict(list)
    subdict_reads: dict[str, list[DictRead]] = defaultdict(list)
    for read in reads:
        if read.kind == "key":
            key_reads[read.name].append(read)
        else:
            subdict_reads[read.name].append(read)

    code_keys_set: set[str] = set(key_reads.keys())
    cat_leaves: set[str] = set()
    cat_parent_segs: set[str] = set()
    for path in cat_paths:
        if not (path.has_wildcard and path.dynamic_path):
            cat_leaves.add(path.leaf)
        for seg in path.parents:
            if not _WILDCARD_RE.fullmatch(seg):
                cat_parent_segs.add(seg)

    absent_keys = {
        key
        for key in code_keys_set
        if key not in cat_leaves and key not in IGNORED_FOAMFILE_KEYS
    }
    stale_paths = {
        path.driver_path
        for path in cat_paths
        if not (path.has_wildcard and path.dynamic_path)
        and path.leaf not in code_keys_set
    }
    unmatched_subdicts = {
        name
        for name in set(subdict_reads) | cat_parent_segs
        if (name in subdict_reads) != (name in cat_parent_segs)
    }

    return {
        "absent_keys": absent_keys,
        "stale_paths": stale_paths,
        "unmatched_subdicts": unmatched_subdicts,
    }


def strict_dict_key_report(
    src_root: Path,
    *,
    allowlist_path: Path | None = None,
) -> DictKeyStrictReport:
    """Return the allowlist-backed strict scanner result."""
    drift = compute_dict_key_drift(src_root)
    allowlist = load_dict_key_allowlist(allowlist_path)

    unexpected: dict[str, set[str]] = {}
    unused: set[str] = set()
    for key in ("absent_keys", "stale_paths", "unmatched_subdicts"):
        unexpected[key] = drift[key] - allowlist[key]
        unused.update(f"{key}:{item}" for item in sorted(allowlist[key] - drift[key]))

    status = "ok" if not any(unexpected.values()) and not unused else "failed"
    return DictKeyStrictReport(
        status=status,
        absent_keys=tuple(sorted(unexpected["absent_keys"])),
        stale_paths=tuple(sorted(unexpected["stale_paths"])),
        unmatched_subdicts=tuple(sorted(unexpected["unmatched_subdicts"])),
        unused_allowlist=tuple(sorted(unused)),
    )
