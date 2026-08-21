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
#     mutators
#
# Description
#     Provides programmatic mutators for modifying OpenFOAM dictionary entries.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

from __future__ import annotations

import re
from pathlib import Path
from typing import Any

from . import foam_backend


def _format_value(value: Any) -> str:
    if isinstance(value, bool):
        return "yes" if value else "no"

    text = str(value)

    # Override values arrive verbatim from sweep.json and the CLI and are
    # written straight into a case dictionary, so a value carrying a `;` can
    # append a second entry, and a `#`-directive becomes code OpenFOAM will
    # compile and run. Neither is a legitimate scalar override; refuse both
    # rather than trusting the caller. See SECURITY.md.
    if ";" in text or "\n" in text:
        raise ValueError(
            f"override value {text!r} contains a statement separator; "
            "a value may not introduce additional dictionary entries"
        )
    if "#" in text:
        raise ValueError(
            f"override value {text!r} contains an OpenFOAM directive; "
            "directives are not permitted in override values"
        )

    return text


def _strip_inline_comment(line: str) -> str:
    return line.split("//", 1)[0]


def _normalize_scope(scope: str | list[str] | tuple[str, ...] | None) -> list[str]:
    if scope is None:
        return []
    if isinstance(scope, str):
        normalized = scope.strip()
        if not normalized:
            raise ValueError("scope cannot be an empty string")
        return [normalized]
    normalized = [str(item).strip() for item in scope]
    if not normalized or any(not item for item in normalized):
        raise ValueError("scope must contain one or more non-empty names")
    return normalized


def _explode_inline_blocks_with_spans(
    lines: list[str],
) -> list[tuple[str, int, int, int]]:
    """Rewrite ``a { b 1; }`` as one brace or entry per virtual line.

    The scope machinery reasons in whole lines, so a block written inline --
    legal OpenFOAM, and present in tracked tutorial dicts -- collapses to a
    degenerate line range and resolves to nothing. Splitting at braces and
    semicolons lets the existing line-based logic handle it unchanged.

    Each result is ``(text, line_index, start_col, end_col)``, so the writing
    path can splice a replacement back into the original line instead of
    reformatting the file the way foamDictionary does. Comments are dropped
    from the virtual text but survive in the untouched remainder of the line.
    """
    exploded: list[tuple[str, int, int, int]] = []
    for index, line in enumerate(lines):
        code = _strip_inline_comment(line)
        if ("{" not in code and "}" not in code) or code.strip() in ("{", "}"):
            exploded.append((line, index, 0, len(line)))
            continue

        buffer = ""
        start = 0
        for position, char in enumerate(code):
            if char in "{}":
                if buffer.strip():
                    exploded.append((buffer.strip() + "\n", index, start, position))
                exploded.append((char + "\n", index, position, position + 1))
                buffer = ""
                start = position + 1
            elif char == ";":
                buffer += char
                exploded.append((buffer.strip() + "\n", index, start, position + 1))
                buffer = ""
                start = position + 1
            else:
                buffer += char
        if buffer.strip():
            exploded.append((buffer.strip() + "\n", index, start, len(code)))
    return exploded


def _explode_inline_blocks(lines: list[str]) -> list[str]:
    """The virtual-line texts of :func:`_explode_inline_blocks_with_spans`."""
    return [text for text, _, _, _ in _explode_inline_blocks_with_spans(lines)]


def _iter_direct_child_lines(lines: list[str], start: int, end: int):
    """Yield the indices in ``[start, end)`` that sit at that span's own
    level, skipping lines owned by a nested sub-dictionary.

    A scope names one dictionary, not it and all its descendants -- so both
    the key scans and the block-header scan below must ignore nested content.
    Braces inside comments don't count; tracked dicts do contain ``// }``.
    """
    depth = 0
    for idx in range(start, end):
        # Yield before this line's braces, so a sub-dictionary's header and
        # its closing brace both count as part of the nested block.
        if depth == 0:
            yield idx
        for ch in _strip_inline_comment(lines[idx]):
            if ch == "{":
                depth += 1
            elif ch == "}":
                depth -= 1


def _quoted_pattern_headers(
    lines: list[str], start: int, end: int
) -> list[tuple[str, str]]:
    """Return ``(regex_source, on_disk_name)`` for quoted block headers.

    OpenFOAM lets a sub-dictionary be keyed by a quoted regex, e.g.
    ``"Vm|VmFinal|u|uFinal"``, which matches the field ``Vm``. Both this
    module's line scanner and foamlib match block names literally, so such a
    block was previously unreachable by member name.
    """
    headers: list[tuple[str, str]] = []
    for index in _iter_direct_child_lines(lines, start, end):
        candidate = _strip_inline_comment(lines[index]).strip()
        if not candidate.startswith('"'):
            continue
        closing = candidate.find('"', 1)
        if closing <= 0:
            continue
        headers.append((candidate[1:closing], candidate[: closing + 1]))
    return headers


def _resolve_pattern_scope(
    lines: list[str], dict_name: str, *, start: int, end: int
) -> str | None:
    """Map a member name onto the quoted-regex block header that matches it.

    OpenFOAM's precedence: an exact literal key wins; otherwise the
    *last-declared* matching pattern wins. Returns the on-disk header text
    (quotes included) so the caller can match it literally, or ``None``.
    """
    for regex_source, on_disk in reversed(
        _quoted_pattern_headers(lines, start, end)
    ):
        try:
            if re.fullmatch(regex_source, dict_name):
                return on_disk
        except re.error:
            continue
    return None


def _find_dict_block_bounds(
    lines: list[str],
    dict_name: str,
    *,
    start: int,
    end: int,
) -> tuple[int, int]:
    # A trailing \b fails to match a scope name ending in a non-word
    # character (e.g. a quoted regex-style block name like
    # "phiE|phiEFinal|phiI|phiIFinal" -- both the closing quote and whatever
    # follows it, whitespace or newline, are non-word, so there is no word
    # boundary there at all). Require whitespace, an opening brace, or
    # end-of-line instead, which also still correctly rejects a longer name
    # that merely has this one as a prefix (e.g. "singleCellSolverCoeffs"
    # must not match a line starting "singleCellSolverCoeffsExtra").
    header_pattern = re.compile(rf"^\s*{re.escape(dict_name)}(?=\s|\{{|$)")

    for i in _iter_direct_child_lines(lines, start, end):
        candidate = _strip_inline_comment(lines[i])
        if not header_pattern.match(candidate):
            continue

        stripped = candidate.strip()
        if stripped.endswith(";") and "{" not in stripped:
            # This candidate is shaped like a scalar entry (name value;)
            # sharing dict_name, not a block header. Scanning forward from
            # here for the first '{' would silently walk into an unrelated
            # sibling block and treat its contents as this scope's -- keep
            # looking for a genuine block header with this name instead.
            continue

        # OpenFOAM dicts commonly appear as:
        #   someDict
        #   {
        # or
        #   someDict {
        open_line = i
        while open_line < end and "{" not in _strip_inline_comment(lines[open_line]):
            open_line += 1

        if open_line >= end:
            raise KeyError(f"Scope '{dict_name}' has no opening brace")

        depth = 0
        saw_open = False
        close_line: int | None = None
        for j in range(open_line, end):
            text = _strip_inline_comment(lines[j])
            for ch in text:
                if ch == "{":
                    depth += 1
                    saw_open = True
                elif ch == "}" and saw_open:
                    depth -= 1
                    if depth == 0:
                        close_line = j
                        break
            if close_line is not None:
                break

        if close_line is None:
            raise KeyError(f"Scope '{dict_name}' has unbalanced braces")

        return open_line + 1, close_line

    if not dict_name.startswith('"'):
        resolved = _resolve_pattern_scope(lines, dict_name, start=start, end=end)
        if resolved is not None:
            return _find_dict_block_bounds(lines, resolved, start=start, end=end)

    raise KeyError(f"Scope '{dict_name}' not found")


def _resolve_search_region(
    lines: list[str],
    scope: str | list[str] | tuple[str, ...] | None,
) -> tuple[int, int]:
    scope_path = _normalize_scope(scope)
    if not scope_path:
        return 0, len(lines)

    start, end = 0, len(lines)
    for dict_name in scope_path:
        start, end = _find_dict_block_bounds(lines, dict_name, start=start, end=end)
    return start, end


def read_foam_entry(
    file_path: Path,
    key: str,
    *,
    scope: str | list[str] | tuple[str, ...] | None = None,
) -> str | None:
    """Read the value of a key from an OpenFOAM dictionary-like text file.

    Reuses the ``_resolve_search_region`` / ``_find_dict_block_bounds``
    infrastructure from :func:`update_foam_entry`. Returns the raw value
    string — trailing semicolon and inline comments stripped — or ``None``
    if the key or its scope block is absent.

    Deliberately does NOT shell out to foamDictionary, unlike its writing
    siblings. foamDictionary re-serialises what it reads (``0.0`` -> ``0``,
    ``5.5e-3`` -> ``0.0055``), and those values feed the dict builders, so
    preferring it made generated dicts and their provenance digests depend on
    whether OpenFOAM happened to be sourced. It also *evaluates* the file:
    a ``#calc``/``#codeStream`` entry is compiled and executed to produce the
    value, which is not an acceptable side effect of reading a case whose
    override values are written in verbatim. Returning the literal source
    text is both deterministic and inert.
    """
    if not file_path.exists():
        return None

    key_pattern = re.compile(rf"^\s*{re.escape(key)}\b")
    lines = _explode_inline_blocks(file_path.read_text().splitlines(keepends=True))
    try:
        search_start, search_end = _resolve_search_region(lines, scope)
    except KeyError:
        return None

    for idx in _iter_direct_child_lines(lines, search_start, search_end):
        line = lines[idx]
        stripped = _strip_inline_comment(line).strip()
        if stripped.startswith("//"):
            continue
        if not key_pattern.match(line):
            continue
        value_part = stripped[len(key):].strip().rstrip(";").strip()
        return value_part if value_part else None

    return None


def update_control_dict(
    control_dict_path: Path,
    *,
    delta_t: float | str | None = None,
    end_time: float | str | None = None,
    start_time: float | str | None = None,
    write_interval: float | str | None = None,
    write_control: str | None = None,
    write_format: str | None = None,
    purge_write: int | str | None = None,
) -> None:
    """Patch entries in an existing OpenFOAM ``controlDict``.

    Each parameter is optional — pass only the values you want to change.
    Raises ``FileNotFoundError`` (via :func:`update_foam_entry`) if the file
    does not exist.
    """
    patches: dict[str, float | str | int] = {
        "deltaT": delta_t,
        "endTime": end_time,
        "startTime": start_time,
        "writeInterval": write_interval,
        "writeControl": write_control,
        "writeFormat": write_format,
        "purgeWrite": purge_write,
    }

    for key, value in patches.items():
        if value is not None:
            update_foam_entry(control_dict_path, key, value)


def update_foam_entry(
    file_path: Path,
    key: str,
    value: Any,
    *,
    scope: str | list[str] | tuple[str, ...] | None = None,
    add_if_missing: bool = False,
) -> None:
    """
    Update a key in an OpenFOAM dictionary-like text file.

    Matches the first non-comment line starting with `key` and rewrites it as:
        <indent><key>    <value>;

    If `scope` is provided, the update is restricted to the named dictionary
    block (or nested path of blocks).

    add_if_missing=True appends `<key>    <value>;` as a new line at the end
    of the scoped block instead of raising when the key isn't already
    present -- e.g. an optional fvSolution PIMPLE control (nNonOrthogonal-
    Correctors) that a case's committed dict may or may not already declare.
    foamDictionary's own `-set` already does this implicitly (auto-creates a
    missing key); this mirrors that for the pure-Python fallback path so
    behaviour doesn't depend on whether OpenFOAM happens to be sourced.
    """
    if not file_path.exists():
        raise FileNotFoundError(f"Dictionary file not found: {file_path}")

    key_pattern = re.compile(rf"^\s*{re.escape(key)}\b")
    lines = file_path.read_text().splitlines(keepends=True)
    virtual = _explode_inline_blocks_with_spans(lines)
    try:
        search_start, search_end = _resolve_search_region(
            [t for t, _, _, _ in virtual], scope
        )
    except KeyError:
        # The line scanner couldn't even resolve the scope -- e.g. a brace
        # inside a quoted value defeats its brace counting. Let a real parser
        # have a go before giving up.
        return foam_backend.update_entry(
            file_path, key, value, scope=scope, add_if_missing=add_if_missing
        )
    direct = _iter_direct_child_lines([t for t, _, _, _ in virtual], search_start, search_end)

    target: tuple[int, int, int] | None = None
    for idx in direct:
        text, line_index, start, end = virtual[idx]
        if text.strip().startswith("//") or not key_pattern.match(text):
            continue
        target = (line_index, start, end)
        break

    replaced = target is not None
    if replaced:
        line_index, start, end = target
        line = lines[line_index]
        if start == 0 and end >= len(line.rstrip("\n")):
            # The entry owns the whole line: keep the original indentation.
            indent = line[: len(line) - len(line.lstrip())]
            lines[line_index] = f"{indent}{key}    {_format_value(value)};\n"
        else:
            # Inline block: splice in place so the rest of the line -- sibling
            # entries, closing braces, any trailing comment -- is preserved.
            lines[line_index] = (
                line[:start] + f"{key}    {_format_value(value)};" + line[end:]
            )
        file_path.write_text("".join(lines))

    if not replaced:
        if add_if_missing:
            if scope is None:
                raise ValueError("add_if_missing requires a scope")
            if search_end >= len(virtual):
                raise KeyError(f"Scope '{scope}' not found in {file_path}")
            insert_before_index = virtual[search_end][1]
            indent = "    "
            for idx in _iter_direct_child_lines(
                [t for t, _, _, _ in virtual], search_start, search_end
            ):
                text, line_index, _start, _end = virtual[idx]
                if text.strip() and not text.strip().startswith("//"):
                    sibling_line = lines[line_index]
                    indent = sibling_line[: len(sibling_line) - len(sibling_line.lstrip())]
                    break
            lines.insert(insert_before_index, f"{indent}{key}    {_format_value(value)};\n")
            file_path.write_text("".join(lines))
            return
        return foam_backend.update_entry(
            file_path, key, value, scope=scope, add_if_missing=add_if_missing
        )


def remove_foam_dict(
    file_path: Path,
    dict_name: str,
    *,
    scope: str | list[str] | tuple[str, ...] | None = None,
    missing_ok: bool = False,
) -> None:
    """Remove a dictionary block from an OpenFOAM dictionary-like text file."""
    if not file_path.exists():
        raise FileNotFoundError(f"Dictionary file not found: {file_path}")

    lines = file_path.read_text().splitlines(keepends=True)
    try:
        search_start, search_end = _resolve_search_region(lines, scope)
    except KeyError:
        return foam_backend.remove_dict(
            file_path, dict_name, scope=scope, missing_ok=missing_ok
        )
    # A trailing \b fails to match a scope name ending in a non-word
    # character (e.g. a quoted regex-style block name like
    # "phiE|phiEFinal|phiI|phiIFinal" -- both the closing quote and whatever
    # follows it, whitespace or newline, are non-word, so there is no word
    # boundary there at all). Require whitespace, an opening brace, or
    # end-of-line instead, which also still correctly rejects a longer name
    # that merely has this one as a prefix (e.g. "singleCellSolverCoeffs"
    # must not match a line starting "singleCellSolverCoeffsExtra").
    header_pattern = re.compile(rf"^\s*{re.escape(dict_name)}(?=\s|\{{|$)")

    remove_start: int | None = None
    remove_end: int | None = None

    i = search_start
    while i < search_end:
        candidate = _strip_inline_comment(lines[i])
        if not header_pattern.match(candidate):
            i += 1
            continue

        open_line = i
        while open_line < search_end and "{" not in _strip_inline_comment(lines[open_line]):
            open_line += 1

        if open_line >= search_end:
            return foam_backend.remove_dict(
                file_path, dict_name, scope=scope, missing_ok=missing_ok
            )

        depth = 0
        saw_open = False
        for j in range(open_line, search_end):
            text = _strip_inline_comment(lines[j])
            for ch in text:
                if ch == "{":
                    depth += 1
                    saw_open = True
                elif ch == "}" and saw_open:
                    depth -= 1
                    if depth == 0:
                        remove_start = i
                        remove_end = j + 1
                        break
            if remove_end is not None:
                break
        break

    if remove_start is None or remove_end is None:
        if not dict_name.startswith('"'):
            resolved = _resolve_pattern_scope(
                lines, dict_name, start=search_start, end=search_end
            )
            if resolved is not None:
                return remove_foam_dict(
                    file_path, resolved, scope=scope, missing_ok=missing_ok
                )
        return foam_backend.remove_dict(
            file_path, dict_name, scope=scope, missing_ok=missing_ok
        )

    del lines[remove_start:remove_end]
    file_path.write_text("".join(lines))


def remove_foam_entry(
    file_path: Path,
    entry_name: str,
    *,
    scope: str | list[str] | tuple[str, ...] | None = None,
    missing_ok: bool = False,
) -> None:
    """Remove a scalar entry (``name value;``) from an OpenFOAM dictionary file.

    The scalar counterpart of :func:`remove_foam_dict`. That one deletes a
    ``name { ... }`` block and raises if the matched name has no opening
    brace; this one deletes the entry's line (or lines, for a value that
    wraps before its terminating ``;``) and raises if the matched name turns
    out to open a block instead.

    Needed because a caller who wants a scalar gone has no way to say so with
    the block remover: ``missing_ok`` covers only *absence*, so a key that is
    present but the wrong shape raises regardless. Deleting
    ``surfaceCurrentPatches.xMin`` when switching a bath-bidomain case between
    boundary variants is exactly that case.
    """
    if not file_path.exists():
        raise FileNotFoundError(f"Dictionary file not found: {file_path}")

    lines = file_path.read_text().splitlines(keepends=True)
    try:
        search_start, search_end = _resolve_search_region(lines, scope)
    except KeyError:
        if missing_ok:
            return
        raise

    header_pattern = re.compile(rf"^\s*{re.escape(entry_name)}(?=\s|;|$)")

    for i in range(search_start, search_end):
        candidate = _strip_inline_comment(lines[i])
        if not header_pattern.match(candidate):
            continue

        # Distinguish a scalar from a block: a block's brace opens either on
        # the header line or before the first ``;``.
        end = i
        while end < search_end and ";" not in _strip_inline_comment(lines[end]):
            if "{" in _strip_inline_comment(lines[end]):
                raise KeyError(
                    f"'{entry_name}' is a dictionary, not a scalar entry; "
                    "use remove_foam_dict"
                )
            end += 1
        if end >= search_end:
            raise KeyError(f"Entry '{entry_name}' has no terminating ';'")
        if "{" in _strip_inline_comment(lines[i]):
            raise KeyError(
                f"'{entry_name}' is a dictionary, not a scalar entry; "
                "use remove_foam_dict"
            )

        del lines[i:end + 1]
        file_path.write_text("".join(lines))
        return

    if missing_ok:
        return
    if scope is None:
        raise KeyError(f"Entry '{entry_name}' not found in {file_path}")
    raise KeyError(f"Entry '{entry_name}' not found in scope '{scope}' in {file_path}")


def read_foam_dict_block(
    file_path: Path,
    dict_name: str,
    *,
    scope: str | list[str] | tuple[str, ...] | None = None,
) -> str | None:
    """Read the raw text of a named sub-dictionary block -- header line,
    braces, and body verbatim -- from an OpenFOAM dictionary-like text
    file, or ``None`` if the file, its scope, or the block itself is
    absent.

    Block-location logic is the read-only twin of :func:`remove_foam_dict`'s
    fallback scan (same header pattern, same brace-depth bookkeeping),
    deliberately operating on raw, unexploded lines rather than
    :func:`read_foam_entry`'s exploded ones: the goal here is to hand back
    exactly what is on disk so it can be replayed verbatim via
    :func:`ensure_foam_dict`'s ``block_text`` parameter elsewhere (e.g. a
    dict regenerator carrying a ``conductionNetworkDomains`` block forward
    unmodified across a top-level solver switch it has no way to
    reconstruct from the catalog alone), not to reformat or reason about
    individual leaf entries.

    Deliberately text-based, not foamDictionary-based, for the same reason
    as :func:`read_foam_entry`: foamDictionary re-serialises and evaluates
    the file it reads, which is not an acceptable side effect of a read.
    """
    if not file_path.exists():
        return None

    lines = file_path.read_text().splitlines(keepends=True)
    try:
        search_start, search_end = _resolve_search_region(lines, scope)
    except KeyError:
        return None

    header_pattern = re.compile(rf"^\s*{re.escape(dict_name)}(?=\s|\{{|$)")

    i = search_start
    while i < search_end:
        candidate = _strip_inline_comment(lines[i])
        if not header_pattern.match(candidate):
            i += 1
            continue

        open_line = i
        while open_line < search_end and "{" not in _strip_inline_comment(lines[open_line]):
            open_line += 1
        if open_line >= search_end:
            return None

        depth = 0
        saw_open = False
        for j in range(open_line, search_end):
            text = _strip_inline_comment(lines[j])
            for ch in text:
                if ch == "{":
                    depth += 1
                    saw_open = True
                elif ch == "}" and saw_open:
                    depth -= 1
                    if depth == 0:
                        return "".join(lines[i:j + 1])
        return None

    return None


# No foamlib fallback here, deliberately. block_text is inserted verbatim to
# preserve comments and formatting for the dynamic-container carry-forward
# in plugins/cardiacfoam/dict_builder.py; routing it through a parser would
# re-serialise exactly what this path exists to keep intact.
def ensure_foam_dict(
    file_path: Path,
    dict_name: str,
    block_text: str,
    *,
    scope: str | list[str] | tuple[str, ...] | None = None,
) -> bool:
    """Insert a dictionary block if it is missing from the selected scope."""
    if not file_path.exists():
        raise FileNotFoundError(f"Dictionary file not found: {file_path}")

    lines = file_path.read_text().splitlines(keepends=True)
    search_start, search_end = _resolve_search_region(lines, scope)
    # A trailing \b fails to match a scope name ending in a non-word
    # character (e.g. a quoted regex-style block name like
    # "phiE|phiEFinal|phiI|phiIFinal" -- both the closing quote and whatever
    # follows it, whitespace or newline, are non-word, so there is no word
    # boundary there at all). Require whitespace, an opening brace, or
    # end-of-line instead, which also still correctly rejects a longer name
    # that merely has this one as a prefix (e.g. "singleCellSolverCoeffs"
    # must not match a line starting "singleCellSolverCoeffsExtra").
    header_pattern = re.compile(rf"^\s*{re.escape(dict_name)}(?=\s|\{{|$)")

    for idx in range(search_start, search_end):
        candidate = _strip_inline_comment(lines[idx])
        if header_pattern.match(candidate):
            return False

    if not dict_name.startswith('"'):
        # dict_name may already be covered by a quoted-regex block header
        # (e.g. "Vm|VmFinal" already matches "Vm") even though no literal
        # header matched above. Without this check, ensure_foam_dict would
        # insert a duplicate literal block that OpenFOAM's own
        # literal-beats-pattern precedence then shadows the existing one
        # with -- a surprising side effect of a false "not found".
        if (
            _resolve_pattern_scope(lines, dict_name, start=search_start, end=search_end)
            is not None
        ):
            return False

    block_lines = block_text.splitlines(keepends=True)
    if not block_lines:
        raise ValueError("block_text cannot be empty")
    if not block_lines[-1].endswith("\n"):
        block_lines[-1] = f"{block_lines[-1]}\n"

    lines[search_end:search_end] = block_lines
    file_path.write_text("".join(lines))

    return True
