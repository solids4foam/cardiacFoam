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

_PREFIX = "$ELECTRO_MODEL_COEFFS."


def _is_safe_system_path(path_str: str) -> bool:
    """Validate that the path is strictly inside system/ and has no traversal segments."""
    if not path_str.startswith("system/"):
        return False
    path = PurePath(path_str)
    return not path.is_absolute() and ".." not in path.parts


class OverrideError(ValueError):
    """An override is malformed, non-applyable, out-of-enum, or failed to apply."""


def _catalog_entries(driver_context=None) -> tuple[set[str], tuple[Any, ...]]:
    from openfoam_driver.core.compatibility import resolve_public_driver_context

    driver_context = resolve_public_driver_context(driver_context)
    catalog = driver_context.capabilities.dictionaries.catalog()
    return (
        {entry.driver_path for entry in catalog.entries_for("controlDict")},
        catalog.entries_for("electroProperties"),
    )


def _electro_by_path(entries: Iterable[Any]) -> dict[str, Any]:
    """Map every electro entry's full driver_path ($ELECTRO_MODEL_COEFFS.<...>) -> entry."""
    out: dict[str, Any] = {}
    for entry in entries:
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


def validate_overrides(overrides: Any, *, driver_context=None) -> None:
    """Reject anything not safely applyable, *before* any write. Raises OverrideError."""
    if not isinstance(overrides, list):
        raise OverrideError(
            "overrides payload must be a JSON list of {driver_path, value} objects"
        )
    control_dict_keys, electro_entries = _catalog_entries(driver_context)
    electro = _electro_by_path(electro_entries)
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
            # Still must be a real controlDict key -- otherwise this silently passes
            # validation and, at apply time, either raises a raw KeyError (no
            # foamDictionary) or silently writes a brand-new bogus key into
            # controlDict (foamDictionary auto-creates missing keys on `-set`).
            if dp not in control_dict_keys:
                known = ", ".join(sorted(control_dict_keys))
                raise OverrideError(
                    f"override driver_path {dp!r} is not a known controlDict entry. "
                    f"Known controlDict entries: {known}"
                )
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
                # foamDictionary spells this scope as "solvers/V/tolerance";
                # update_foam_entry takes it apart. Going through it rather
                # than straight to foamDictionary keeps this route usable
                # without a sourced OpenFOAM, like every other override path.
                *scope_path, key = entry_path.split("/")
                update_foam_entry(
                    case_root / file_path, key, value, scope=scope_path or None
                )
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
