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

from pathlib import Path
from typing import Any

from .common import detect_myocardium_solver_name
from .dict_builder import _entry_scope_and_key
from ..core.runtime.mutators import update_control_dict, update_foam_entry
from ..dict_entries import ELECTRO_PROPERTY_ENTRY_GROUPS

_PREFIX = "$ELECTRO_MODEL_COEFFS."

# controlDict leaf -> update_control_dict kwarg. Only these leaves are *applyable*
# (startFrom / stopAt are catalog-addressable but update_control_dict has no kwarg for
# them, so they are rejected at validation time).
_CONTROL_DICT_KWARG: dict[str, str] = {
    "deltaT": "delta_t",
    "endTime": "end_time",
    "startTime": "start_time",
    "writeInterval": "write_interval",
    "writeControl": "write_control",
    "writeFormat": "write_format",
    "purgeWrite": "purge_write",
}


class OverrideError(ValueError):
    """An override is malformed, non-applyable, out-of-enum, or failed to apply."""


def _electro_by_path() -> dict[str, Any]:
    """Map every electro entry's full driver_path ($ELECTRO_MODEL_COEFFS.<...>) -> entry."""
    out: dict[str, Any] = {}
    for group in ELECTRO_PROPERTY_ENTRY_GROUPS.values():
        for entry in group:
            out[entry.driver_path] = entry
    return out


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
        if dp in _CONTROL_DICT_KWARG:
            continue
        entry = electro.get(dp)
        if entry is None:
            raise OverrideError(
                f"override driver_path {dp!r} is not catalog-addressable / applyable"
            )
        if getattr(entry, "dynamic_path", False):
            raise OverrideError(
                f"override driver_path {dp!r} is a dynamic_path entry and is not applyable"
            )
        if "<" in dp:
            raise OverrideError(
                f"override driver_path {dp!r} contains a placeholder; substitute the "
                f"concrete name"
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
            kwarg = _CONTROL_DICT_KWARG.get(dp)
            if kwarg is not None:
                update_control_dict(case_root / "system" / "controlDict", **{kwarg: value})
            else:
                if coeffs_scope is None:
                    coeffs_scope = f"{detect_myocardium_solver_name(electro_path)}Coeffs"
                scope_path, key = _entry_scope_and_key(dp, coeffs_scope)
                update_foam_entry(electro_path, key, value, scope=scope_path)
        except (OSError, KeyError, ValueError) as exc:
            raise OverrideError(f"failed to apply override {dp!r}: {exc}") from exc
