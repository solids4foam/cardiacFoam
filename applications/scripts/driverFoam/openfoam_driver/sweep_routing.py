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

_NON_ROUTABLE_KEYS: frozenset[str] = frozenset({"caseId"})


def _route_case_values_legacy(
    *, base: dict[str, Any], resolved_axis_values: dict[str, Any], driver_context=None,
) -> dict[str, Any]:
    """Legacy public route, implemented by the cardiacFoam plugin."""
    from .plugins.cardiacfoam.sweep import route_case_values

    return route_case_values(
        base=base,
        resolved_axis_values=resolved_axis_values,
        driver_context=driver_context,
    )


def route_case_values(
    *, base: dict[str, Any], resolved_axis_values: dict[str, Any], driver_context=None,
) -> dict[str, Any]:
    """Route through the selected plugin while preserving the public API."""

    from .core.compatibility import resolve_public_driver_context
    from .core.plugin_capabilities import SweepRoutingRequest

    driver_context = resolve_public_driver_context(driver_context)
    return driver_context.capabilities.sweep_materializer.route(
        SweepRoutingRequest(
            base=base,
            resolved_axis_values=resolved_axis_values,
        ),
        driver_context=driver_context,
    )


_ENTRY_NON_ROUTABLE_KEYS: frozenset[str] = _NON_ROUTABLE_KEYS | frozenset(
    {"entry", "archive_dir_name"}
)


def route_entry_case_values(
    *, base: dict[str, Any], resolved_axis_values: dict[str, Any],
) -> dict[str, Any]:
    """Route a resolved entry-based sweep case's values to make_spec kwargs.

    Entry-based sweeps target an existing registered tutorial's own
    make_spec(**kwargs) (e.g. niederer_2012.py's dx_values/dt_values/
    end_time_by_dx), which already validates its own keyword arguments --
    unlike route_case_values's build_and_launch target, there is no fixed
    vocabulary to classify values into here. `base` carries kwargs fixed
    across every case in the sweep (e.g. solvers, end_time_by_dx); per-case
    resolved_axis_values are merged on top, winning on conflict. Every value
    passes straight through except sweep_expansion's bookkeeping keys and
    "entry" itself (sweep_runner's own dispatch key, not a make_spec kwarg).
    """
    merged = {**base, **resolved_axis_values}
    return {
        key: value
        for key, value in merged.items()
        if key not in _ENTRY_NON_ROUTABLE_KEYS
    }
