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
#     sweep_routing
#
# Description
#     Routes resolved sweep values to build_and_launch's structured parameters.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

from __future__ import annotations

from typing import Any

from .dict_entries import CONTROL_DICT_ENTRIES, PHYSICS_PROPERTY_ENTRIES
from .sweep_expansion import SweepValidationError
from .specs.dict_builder import SELECTOR_KEYS

_CONTROL_DICT_KEYS: frozenset[str] = frozenset(entry.driver_path for entry in CONTROL_DICT_ENTRIES)
_SUPPORTED_CONTROL_DICT_KEYS: frozenset[str] = frozenset({"deltaT", "endTime"})
_UNSUPPORTED_CONTROL_DICT_KEYS: frozenset[str] = _CONTROL_DICT_KEYS - _SUPPORTED_CONTROL_DICT_KEYS
_PHYSICS_SELECTOR_KEYS: frozenset[str] = frozenset(entry.driver_path for entry in PHYSICS_PROPERTY_ENTRIES)

# Bookkeeping keys the expander may add that are never real dict values.
_NON_ROUTABLE_KEYS: frozenset[str] = frozenset({"caseId"})


def route_case_values(
    *, base: dict[str, Any], resolved_axis_values: dict[str, Any],
) -> dict[str, Any]:
    """Classify a resolved case's values into build_and_launch's parameters.

    `base` may carry pre-set "electro_selectors"/"physics_selectors"/
    "electro_overrides"/"physics_overrides" dicts (fixed across the whole
    sweep); resolved_axis_values are merged on top per the routing rules.
    """
    electro_selectors: dict[str, Any] = dict(base.get("electro_selectors", {}))
    physics_selectors: dict[str, Any] = dict(base.get("physics_selectors", {}))
    electro_overrides: dict[str, Any] = dict(base.get("electro_overrides", {}))
    physics_overrides: dict[str, Any] = dict(base.get("physics_overrides", {}))
    delta_t: Any = base.get("delta_t")
    end_time: Any = base.get("end_time")

    for key, value in resolved_axis_values.items():
        if key in _NON_ROUTABLE_KEYS:
            continue
        if key in SELECTOR_KEYS:
            electro_selectors[key] = value
        elif key in _PHYSICS_SELECTOR_KEYS:
            physics_selectors[key] = value
        elif key == "deltaT":
            delta_t = value
        elif key == "endTime":
            end_time = value
        elif key in _UNSUPPORTED_CONTROL_DICT_KEYS:
            known = ", ".join(sorted(_SUPPORTED_CONTROL_DICT_KEYS))
            raise SweepValidationError(
                f"controlDict sweep axis '{key}' is not supported by this plan; "
                f"supported controlDict axes are: {known}"
            )
        else:
            electro_overrides[key] = value

    return {
        "electro_selectors": electro_selectors,
        "physics_selectors": physics_selectors,
        "electro_overrides": electro_overrides,
        "physics_overrides": physics_overrides,
        "delta_t": delta_t,
        "end_time": end_time,
    }
