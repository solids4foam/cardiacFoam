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
#     apply_overrides
#
# Description
#     Mechanically apply an agent-chosen override set to a case's dicts for
#     `step --strict --apply`. Validates each override for *applyability*
#     (catalog-addressable AND writable by the router) before any write, then
#     routes controlDict leaves to update_control_dict and $ELECTRO_MODEL_COEFFS
#     keys through the existing solver-coeffs resolver. Lives in the specs layer
#     because routing needs that resolver (detect_myocardium_solver_name +
#     _entry_scope_and_key). The driver never decides *what* to change — the
#     agent authors the override set; this only applies a validated one.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

from __future__ import annotations

import re
import shutil
from pathlib import Path, PurePath
from typing import Any, Iterable

from .common import detect_myocardium_solver_name
from .dict_builder import _entry_scope_and_key
from ..core.runtime.mutators import update_foam_entry, update_foam_entry_via_foamDictionary
from ..dict_entries import ELECTRO_PROPERTY_ENTRY_GROUPS

_PREFIX = "$ELECTRO_MODEL_COEFFS."


def _is_safe_system_path(path_str: str) -> bool:
    """Validate that the path is strictly inside system/ and has no traversal segments."""
    if not path_str.startswith("system/"):
        return False
    path = PurePath(path_str)
    return not path.is_absolute() and ".." not in path.parts


class OverrideError(ValueError):
    """An override is malformed, non-applyable, out-of-enum, or failed to apply."""


def _electro_by_path() -> dict[str, Any]:
    """Map every electro entry's full driver_path ($ELECTRO_MODEL_COEFFS.<...>) -> entry."""
    out: dict[str, Any] = {}
    for group in ELECTRO_PROPERTY_ENTRY_GROUPS.values():
        for entry in group:
            out[entry.driver_path] = entry
    return out


def _match_dynamic_entry(dp: str, all_entries: Iterable[Any]) -> Any | None:
    """Return the dynamic catalog entry whose template matches concrete *dp*."""
    for entry in all_entries:
        if not getattr(entry, "dynamic_path", False):
            continue

        template = entry.driver_path
        pattern_parts: list[str] = []
        previous_end = 0
        for placeholder in re.finditer(r"<[^.<>]+>", template):
            pattern_parts.append(re.escape(template[previous_end:placeholder.start()]))
            pattern_parts.append(r"[^.]+")
            previous_end = placeholder.end()
        pattern_parts.append(re.escape(template[previous_end:]))

        if re.fullmatch("".join(pattern_parts), dp):
            return entry
    return None


def validate_overrides(overrides: Any) -> None:
    """Reject anything not safely applyable, *before* any write. Raises OverrideError."""
    if not isinstance(overrides, list):
        raise OverrideError(
            "overrides payload must be a JSON list of {driver_path, value} objects"
        )
    electro = _electro_by_path()
    for ov in overrides:
        if not isinstance(ov, dict) or "driver_path" not in ov or "value" not in ov:
            raise OverrideError(
                f"each override must be an object with 'driver_path' and 'value' (got {ov!r})"
            )
        dp = ov["driver_path"]
        if ":" in dp:
            file_path, _, entry_path = dp.partition(":")
            if not _is_safe_system_path(file_path):
                raise OverrideError(f"override file path {file_path!r} is not a safe system/ path")
            if not entry_path:
                raise OverrideError(f"override driver_path {dp!r} is missing an entry path after ':'")
            continue
        elif not dp.startswith("$"):
            # Backward compatibility: flat strings are treated as controlDict entries.
            continue

        if "<" in dp or ">" in dp:
            raise OverrideError(
                f"override driver_path {dp!r} contains a placeholder; substitute the "
                f"concrete name"
            )

        entry = electro.get(dp)
        if entry is None:
            entry = _match_dynamic_entry(dp, electro.values())
            if entry is None:
                raise OverrideError(
                    f"override driver_path {dp!r} is not catalog-addressable / applyable"
                )
        enum_values = getattr(entry, "enum_values", None)
        if enum_values and ov["value"] not in enum_values:
            raise OverrideError(
                f"override {dp!r} value {ov['value']!r} not in enum {tuple(enum_values)}"
            )


def apply_overrides(overrides: list[dict[str, Any]], *, case_root: Path) -> None:
    """Apply validated overrides to the case dicts.

    Raises OverrideError on any mutator failure (caught at the CLI boundary). Not
    transactional: a mid-list failure can leave earlier overrides applied.
    """
    electro_path = case_root / "constant" / "electroProperties"
    coeffs_scope: str | None = None
    for ov in overrides:
        dp, value = ov["driver_path"], ov["value"]
        try:
            if ":" in dp:
                file_path, _, entry_path = dp.partition(":")
                update_foam_entry_via_foamDictionary(case_root / file_path, entry_path, value)
            elif not dp.startswith("$"):
                if shutil.which("foamDictionary"):
                    update_foam_entry_via_foamDictionary(case_root / "system" / "controlDict", dp, value)
                else:
                    update_foam_entry(case_root / "system" / "controlDict", dp, value)
            else:
                if coeffs_scope is None:
                    coeffs_scope = f"{detect_myocardium_solver_name(electro_path)}Coeffs"
                scope_path, key = _entry_scope_and_key(dp, coeffs_scope)
                update_foam_entry(electro_path, key, value, scope=scope_path)
        except (OSError, KeyError, ValueError, RuntimeError) as exc:
            raise OverrideError(f"failed to apply override {dp!r}: {exc}") from exc
