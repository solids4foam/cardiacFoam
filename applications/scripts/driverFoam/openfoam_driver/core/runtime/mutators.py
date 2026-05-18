from __future__ import annotations

import re
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
