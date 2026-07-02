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


def _find_dict_block_bounds(
    lines: list[str],
    dict_name: str,
    *,
    start: int,
    end: int,
) -> tuple[int, int]:
    header_pattern = re.compile(rf"^\s*{re.escape(dict_name)}\b")

    i = start
    while i < end:
        candidate = _strip_inline_comment(lines[i])
        if not header_pattern.match(candidate):
            i += 1
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


def validate_foam_entries(file_path: Path, entries: list[dict]) -> list[str]:
    """
    Validate that a list of override entries have corresponding keys in the file.
    Returns a list of error messages for any missing keys.
    """
    if not file_path.exists():
        return [f"Dictionary file not found: {file_path}"]

    lines = file_path.read_text().splitlines(keepends=True)
    errors = []

    for entry in entries:
        key = entry["key"]
        scope = entry.get("scope")
        key_pattern = re.compile(rf"^\s*{re.escape(key)}\b")
        try:
            search_start, search_end = _resolve_search_region(lines, scope)
        except KeyError as e:
            errors.append(str(e).strip("'"))
            continue

        found = False
        for idx in range(search_start, search_end):
            line = lines[idx].strip()
            if line.startswith("//"):
                continue
            if key_pattern.match(lines[idx]):
                found = True
                break
        
        if not found:
            if scope is None:
                errors.append(f"Key '{key}' not found in {file_path}")
            else:
                errors.append(f"Key '{key}' not found in scope '{scope}' in {file_path}")

    return errors


def read_foam_entry_via_foamDictionary(
    file_path: Path,
    key: str,
    *,
    scope: str | list[str] | tuple[str, ...] | None = None,
) -> str | None:
    """Read a key using OpenFOAM's foamDictionary utility."""
    if not file_path.exists():
        return None

    scope_path = _normalize_scope(scope)
    entry_path = "/".join(scope_path + [key]) if scope_path else key

    cmd = [
        "foamDictionary",
        str(file_path),
        "-entry",
        entry_path,
        "-value",
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise KeyError(f"Key not found by foamDictionary: {entry_path}")
    return result.stdout.strip()


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
    """
    if not file_path.exists():
        return None

    has_foam_dict = shutil.which("foamDictionary") is not None
    if has_foam_dict:
        try:
            return read_foam_entry_via_foamDictionary(file_path, key, scope=scope)
        except Exception:
            pass

    key_pattern = re.compile(rf"^\s*{re.escape(key)}\b")
    lines = file_path.read_text().splitlines(keepends=True)
    try:
        search_start, search_end = _resolve_search_region(lines, scope)
    except KeyError:
        return None

    for idx in range(search_start, search_end):
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

    key_pattern = re.compile(rf"^\s*{re.escape(key)}\b")
    lines = file_path.read_text().splitlines(keepends=True)
    search_start, search_end = _resolve_search_region(lines, scope)

    replaced = False
    with file_path.open("w") as handle:
        for idx, line in enumerate(lines):
            stripped = line.strip()

            if stripped.startswith("//"):
                handle.write(line)
                continue

            within_scope = search_start <= idx < search_end
            if within_scope and (not replaced) and key_pattern.match(line):
                indent = line[: len(line) - len(line.lstrip())]
                handle.write(f"{indent}{key}    {_format_value(value)};\n")
                replaced = True
            else:
                handle.write(line)

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
    search_start, search_end = _resolve_search_region(lines, scope)
    header_pattern = re.compile(rf"^\s*{re.escape(dict_name)}\b")

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
    header_pattern = re.compile(rf"^\s*{re.escape(dict_name)}\b")

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
