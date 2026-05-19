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
`driver_path` in `PHYSICS_PROPERTY_ENTRIES` and `ELECTRO_PROPERTY_ENTRY_GROUPS`
into a structured form for comparison against the scanner output.

Accuracy is ~80%; false positives/negatives are expected.  The output is for
human review only.
"""

from __future__ import annotations

import re
from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path

from openfoam_driver.dict_entries import (
    ELECTRO_PROPERTY_ENTRY_GROUPS,
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
    for group in ELECTRO_PROPERTY_ENTRY_GROUPS.values():
        for entry in group:
            yield _parse_path(entry.driver_path, entry.dynamic_path)
