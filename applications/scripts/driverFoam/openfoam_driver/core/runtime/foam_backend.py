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
#     foam_backend
#
# Description
#     Complex-syntax fallback for OpenFOAM dictionary mutation, backed by
#     foamlib.
#
#     This replaces the former `foamDictionary` subprocess tier. That tier made
#     written bytes depend on whether OpenFOAM happened to be sourced, and it
#     *executed* `#calc`/`#codeStream` to produce values. foamlib parses in
#     process, never evaluates directives, and is available unconditionally.
#
#     It is a fallback, not the primary path. `mutators.py`'s line-based tier
#     runs first because it writes values verbatim; foamlib is consulted only
#     when the line scanner cannot locate the target (for example a brace
#     inside a quoted string, which defeats brace counting).
#
#     Reads never come here. `read_foam_entry` must return the verbatim on-disk
#     token, and foamlib returns typed values (`5e-6` -> `5e-06`, `#calc` -> a
#     tuple) with no public API for source text.
#----------------------------------------------------------------------------#

from __future__ import annotations

import difflib
import re
import warnings
from pathlib import Path
from typing import Any

from foamlib import FoamFile

ScopeArg = str | list[str] | tuple[str, ...] | None

_TRUE_TOKENS = {"yes", "true", "on"}
_FALSE_TOKENS = {"no", "false", "off"}


def _normalize_scope(scope: ScopeArg) -> tuple[str, ...]:
    if scope is None:
        return ()
    if isinstance(scope, str):
        stripped = scope.strip()
        if not stripped:
            raise ValueError("scope cannot be an empty string")
        return (stripped,)
    parts = tuple(str(item).strip() for item in scope)
    if not parts or any(not part for part in parts):
        raise ValueError("scope must contain one or more non-empty names")
    return parts


def coerce_value(value: Any) -> Any:
    """Convert a driverFOAM override value into a type foamlib will store.

    Override values arrive from ``sweep.json`` and the CLI as strings, but
    foamlib is type-strict on write: it refuses a ``str`` that would be read
    back as something else (``"1e-6"``, ``"0.0"``, ``"uniform 0"`` all raise
    ``ValueError``). Anything that is not clearly numeric or boolean is left as
    a string and allowed to fail loudly in foamlib -- that refusal is a feature,
    because it is what rejects an injected ``"1e-6;  rogue  1"``.
    """
    if not isinstance(value, str):
        return value

    token = value.strip()
    lowered = token.lower()
    if lowered in _TRUE_TOKENS:
        return True
    if lowered in _FALSE_TOKENS:
        return False

    try:
        return int(token)
    except ValueError:
        pass
    try:
        return float(token)
    except ValueError:
        pass

    return token


_ENTRY_SEPARATOR_RE = re.compile(r'(?<![\w"])([A-Za-z_]\w*)( )([^\s{};"][^{};]*;)')


def _reformat_separators(line: str) -> str:
    """Widen ``key value;`` to ``key    value;`` (tier 1's four-space form).

    Matches ``identifier<one space>value;`` anywhere in the line, not just at
    the start, so it also fixes an inline entry like ``Vm { tolerance 1e-12; }``
    without touching the ``Vm {`` that precedes it. The negative lookbehind
    excludes a key immediately preceded by a quote, so it does not fire inside
    a quoted string like ``note "a value with { a brace";``.
    """

    def repl(match: re.Match[str]) -> str:
        key, _space, rest = match.groups()
        return f"{key}    {rest}"

    return _ENTRY_SEPARATOR_RE.sub(repl, line)


def normalize_output(before: str, after: str) -> str:
    """Reshape foamlib's write into the byte form the line-based tier emits.

    foamlib writes ``k 250;`` (one space) and, *measured directly against
    1.7.5*, inserts a blank line both before the edited entry and at
    end-of-file. Tier 1 writes ``k    250;`` (four spaces) and inserts
    nothing.

    A whole-file heuristic scan cannot fix this safely: reformatting every
    line matching ``key value;`` also rewrites lines foamlib never touched --
    e.g. an existing ``note  "a value with { a brace";`` two-space entry
    would be corrupted to four spaces even though foamlib left it
    byte-identical. Instead, diff ``before`` against ``after`` and act only
    on the lines the write actually changed: reformat their separator width,
    and drop any *purely inserted* blank line. Unchanged lines pass through
    verbatim, so untouched entries -- however they happen to be spaced --
    are never rewritten.
    """
    before_lines = before.splitlines(keepends=True)
    after_lines = after.splitlines(keepends=True)
    matcher = difflib.SequenceMatcher(a=before_lines, b=after_lines, autojunk=False)

    out: list[str] = []
    for tag, _i1, _i2, j1, j2 in matcher.get_opcodes():
        if tag == "equal":
            out.extend(after_lines[j1:j2])
            continue
        for line in after_lines[j1:j2]:
            if line.strip() == "":
                continue
            out.append(_reformat_separators(line))
    return "".join(out)


def _require_file(file_path: Path) -> None:
    if not file_path.exists():
        raise FileNotFoundError(f"Dictionary file not found: {file_path}")


def update_entry(
    file_path: Path,
    key: str,
    value: Any,
    *,
    scope: ScopeArg = None,
    add_if_missing: bool = False,
) -> None:
    """Set ``key`` within ``scope``, failing closed when the key is absent.

    ``add_if_missing`` with no ``scope`` is rejected here for the same
    reason tier 1 rejects it at ``mutators.py:434``: without a scope there
    is no well-defined insertion point. Mirroring the guard keeps the two
    tiers agreeing on a case that tier 1's own ``raise ValueError`` already
    intercepts *before* any fallback to this module would ever run --
    leaving this adapter more permissive here would document a capability
    that the real, composed ``update_foam_entry`` never actually exposes.
    """
    _require_file(file_path)
    if add_if_missing and scope is None:
        raise ValueError("add_if_missing requires a scope")
    path = tuple(_normalize_scope(scope)) + (key,)

    before = file_path.read_text()
    foam_file = FoamFile(file_path)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        try:
            if not add_if_missing:
                try:
                    foam_file[path]
                except KeyError as exc:
                    raise KeyError(
                        f"Key '{key}' not found in scope '{scope}' in {file_path}"
                        if scope is not None
                        else f"Key '{key}' not found in {file_path}"
                    ) from exc
            foam_file[path] = coerce_value(value)
        except KeyError:
            raise
        except (TypeError, ValueError) as exc:
            # Catches foamlib's own ValueError (type-inconsistent value) and
            # FoamFileDecodeError (a ValueError subclass raised by a parse
            # failure on either the pre-check read or the write). Neither
            # foamlib exception type is allowed to escape unmapped -- see
            # the spec's Exception contract.
            file_path.write_text(before)
            raise ValueError(f"cannot write value {value!r} to {key!r}: {exc}") from exc

    file_path.write_text(normalize_output(before, file_path.read_text()))


def remove_dict(
    file_path: Path,
    dict_name: str,
    *,
    scope: ScopeArg = None,
    missing_ok: bool = False,
) -> None:
    """Delete a sub-dictionary block.

    foamlib's ``del`` leaves the emptied block as a whitespace-only line
    between its braces (measured: ``solvers\\n{\\n    Vm {...}\\n}\\n`` becomes
    ``solvers\\n{\\n    \\n}\\n``), not a clean removal. Route through
    ``normalize_output`` so the stray line is dropped along with the same
    blank-line class ``update_entry`` already has to handle.
    """
    _require_file(file_path)
    path = tuple(_normalize_scope(scope)) + (dict_name,)

    before = file_path.read_text()
    foam_file = FoamFile(file_path)
    try:
        del foam_file[path]
    except KeyError:
        if missing_ok:
            return
        raise KeyError(f"Dictionary '{dict_name}' not found in {file_path}") from None

    file_path.write_text(normalize_output(before, file_path.read_text()))
