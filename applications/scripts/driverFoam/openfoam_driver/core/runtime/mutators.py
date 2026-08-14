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
import shutil
import subprocess
from pathlib import Path
from typing import Any


def _format_value(value: Any) -> str:
    if isinstance(value, bool):
        return "yes" if value else "no"
    return str(value)


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

    # Check if foamDictionary is available in the current environment
    has_foam_dict = shutil.which("foamDictionary") is not None

    for key, value in patches.items():
        if value is not None:
            success = False
            if has_foam_dict:
                try:
                    update_foam_entry_via_foamDictionary(control_dict_path, key, value)
                    success = True
                except Exception:
                    # Fallback to python string parsing on failure
                    pass

            if not success:
                update_foam_entry(control_dict_path, key, value)


_FOAM_ENTRY_LINE = re.compile(r"^[A-Za-z_][\w.]*\s+\S.*;\s*$")


def _count_foam_entries(text: str) -> int:
    """Count lines that look like a top-level ``key value;`` entry.

    A coarse, syntax-unaware count used only to sanity-check that
    foamDictionary didn't silently discard most of a file's content while
    still exiting 0 (see update_foam_entry_via_foamDictionary).
    """
    count = 0
    for line in text.splitlines():
        stripped = _strip_inline_comment(line).strip()
        if stripped and _FOAM_ENTRY_LINE.match(stripped):
            count += 1
    return count


def update_foam_entry_via_foamDictionary(
    file_path: Path,
    key: str,
    value: Any,
    *,
    scope: str | list[str] | tuple[str, ...] | None = None,
) -> None:
    """Update a key using OpenFOAM's foamDictionary utility."""
    if not file_path.exists():
        raise FileNotFoundError(f"Dictionary file not found: {file_path}")

    original_text = file_path.read_text()

    scope_path = _normalize_scope(scope)
    entry_path = "/".join(scope_path + [key]) if scope_path else key

    cmd = [
        "foamDictionary",
        str(file_path),
        "-entry",
        entry_path,
        "-set",
        _format_value(value),
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(
            f"foamDictionary failed to update '{entry_path}' in {file_path}:\n"
            f"STDOUT: {result.stdout}\nSTDERR: {result.stderr}"
        )

    # foamDictionary can exit 0 while having silently rewritten the file as a
    # near-empty dict (e.g. a malformed header comment makes it fail to parse
    # the existing content, then `-set` auto-creates the missing key). Byte
    # count alone doesn't catch this: foamDictionary always re-serializes its
    # full banner, which can pad a gutted file back up to a similar size. Count
    # recognizable `key value;` entries instead, and revert + raise if most of
    # them vanished, instead of leaving a gutted file behind.
    new_text = file_path.read_text()
    original_entries = _count_foam_entries(original_text)
    new_entries = _count_foam_entries(new_text)
    if original_entries >= 2 and new_entries < original_entries * 0.5:
        file_path.write_text(original_text)
        raise RuntimeError(
            f"foamDictionary reported success but the result has only "
            f"{new_entries} recognizable entries versus {original_entries} "
            f"before, for '{entry_path}' in {file_path}; the input likely "
            "failed to parse (e.g. a malformed header comment). Reverted the "
            "file to avoid silent data loss."
        )


def update_foam_entry(
    file_path: Path,
    key: str,
    value: Any,
    *,
    scope: str | list[str] | tuple[str, ...] | None = None,
) -> None:
    """
    Update a key in an OpenFOAM dictionary-like text file.

    Matches the first non-comment line starting with `key` and rewrites it as:
        <indent><key>    <value>;

    If `scope` is provided, the update is restricted to the named dictionary
    block (or nested path of blocks).
    """
    if not file_path.exists():
        raise FileNotFoundError(f"Dictionary file not found: {file_path}")

    has_foam_dict = shutil.which("foamDictionary") is not None
    if has_foam_dict:
        try:
            update_foam_entry_via_foamDictionary(file_path, key, value, scope=scope)
            return
        except Exception:
            pass

    key_pattern = re.compile(rf"^\s*{re.escape(key)}\b")
    lines = file_path.read_text().splitlines(keepends=True)
    virtual = _explode_inline_blocks_with_spans(lines)
    search_start, search_end = _resolve_search_region([t for t, _, _, _ in virtual], scope)
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
        if scope is None:
            raise KeyError(f"Key '{key}' not found in {file_path}")
        raise KeyError(f"Key '{key}' not found in scope '{scope}' in {file_path}")


def remove_foam_dict_via_foamDictionary(
    file_path: Path,
    dict_name: str,
    *,
    scope: str | list[str] | tuple[str, ...] | None = None,
    missing_ok: bool = False,
) -> None:
    """Remove a dictionary block using OpenFOAM's foamDictionary utility."""
    if not file_path.exists():
        raise FileNotFoundError(f"Dictionary file not found: {file_path}")

    scope_path = _normalize_scope(scope)
    entry_path = "/".join(scope_path + [dict_name]) if scope_path else dict_name

    cmd = [
        "foamDictionary",
        str(file_path),
        "-entry",
        entry_path,
        "-remove",
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        if missing_ok:
            return
        raise KeyError(f"Dictionary '{entry_path}' not found by foamDictionary in {file_path}")


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

    has_foam_dict = shutil.which("foamDictionary") is not None
    if has_foam_dict:
        try:
            remove_foam_dict_via_foamDictionary(file_path, dict_name, scope=scope, missing_ok=missing_ok)
            return
        except Exception:
            pass

    lines = file_path.read_text().splitlines(keepends=True)
    try:
        search_start, search_end = _resolve_search_region(lines, scope)
    except KeyError:
        if missing_ok:
            return
        raise
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
            raise KeyError(f"Dictionary '{dict_name}' has no opening brace")

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
        if missing_ok:
            return
        if scope is None:
            raise KeyError(f"Dictionary '{dict_name}' not found in {file_path}")
        raise KeyError(f"Dictionary '{dict_name}' not found in scope '{scope}' in {file_path}")

    del lines[remove_start:remove_end]
    file_path.write_text("".join(lines))


def ensure_foam_dict_via_foamDictionary(
    file_path: Path,
    dict_name: str,
    block_text: str,
    *,
    scope: str | list[str] | tuple[str, ...] | None = None,
) -> bool:
    """Insert a dictionary block using OpenFOAM's foamDictionary utility."""
    if not file_path.exists():
        raise FileNotFoundError(f"Dictionary file not found: {file_path}")

    scope_path = _normalize_scope(scope)
    entry_path = "/".join(scope_path + [dict_name]) if scope_path else dict_name

    # Check if the entry already exists
    cmd_check = ["foamDictionary", str(file_path), "-entry", entry_path]
    res_check = subprocess.run(cmd_check, capture_output=True)
    if res_check.returncode == 0:
        return False

    # Extract the internal brace content since -add only takes the value
    start_idx = block_text.find("{")
    end_idx = block_text.rfind("}")
    if start_idx == -1 or end_idx == -1 or start_idx > end_idx:
        raise ValueError("block_text must contain { and }")

    inner_val = block_text[start_idx : end_idx + 1]

    cmd_add = [
        "foamDictionary",
        str(file_path),
        "-entry",
        entry_path,
        "-add",
        inner_val,
    ]
    res_add = subprocess.run(cmd_add, capture_output=True, text=True)
    if res_add.returncode != 0:
        raise RuntimeError(f"foamDictionary failed to add block to '{entry_path}' in {file_path}:\n{res_add.stderr}")
    return True


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

    has_foam_dict = shutil.which("foamDictionary") is not None
    if has_foam_dict:
        try:
            return ensure_foam_dict_via_foamDictionary(file_path, dict_name, block_text, scope=scope)
        except Exception:
            pass

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

    block_lines = block_text.splitlines(keepends=True)
    if not block_lines:
        raise ValueError("block_text cannot be empty")
    if not block_lines[-1].endswith("\n"):
        block_lines[-1] = f"{block_lines[-1]}\n"

    lines[search_end:search_end] = block_lines
    file_path.write_text("".join(lines))

    return True
