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
#     plugin_discovery
#
# Description
#     Discovery of installed solver plugins through the 'driverfoam.plugins'
#     entry-point group. This is discovery, not sandboxing: loading a plugin
#     executes its Python code, exactly as the trusted module:Class form does.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

from __future__ import annotations

from importlib.metadata import entry_points
from typing import Any

ENTRY_POINT_GROUP = "driverfoam.plugins"


def _entry_points() -> tuple[Any, ...]:
    """Indirection seam so tests can inject entry points without installing."""
    return tuple(entry_points(group=ENTRY_POINT_GROUP))


def ambiguous_plugin_names() -> dict[str, tuple[str, ...]]:
    """Entry-point names claimed by more than one installed distribution.

    Two distributions exporting the same name is a packaging conflict, not
    something to resolve by dictionary insertion order -- which distribution
    won would depend on installation order and be invisible in the plan.
    """
    seen: dict[str, list[str]] = {}
    for entry_point in _entry_points():
        dist = getattr(entry_point, "dist", None)
        origin = f"{dist.name}={dist.version}" if dist is not None else "<unknown>"
        seen.setdefault(entry_point.name, []).append(origin)
    return {
        name: tuple(sorted(origins))
        for name, origins in seen.items()
        if len(origins) > 1
    }


def discover_plugins() -> dict[str, Any]:
    """Return unambiguously installed plugin entry points keyed by name.

    Never raises: a broken third-party distribution must not make the CLI
    unusable for everyone else. A name claimed by several distributions is
    omitted here and reported by :func:`ambiguous_plugin_names`, so it fails
    loudly at load time rather than silently resolving to whichever
    distribution happened to be enumerated last.
    """
    ambiguous = set(ambiguous_plugin_names())
    return {
        entry_point.name: entry_point
        for entry_point in _entry_points()
        if entry_point.name not in ambiguous
    }


def load_discovered_plugin(name: str):
    """Load and validate a discovered plugin by entry-point name.

    Loading executes the plugin's Python code. Identity provenance records the
    installing distribution's name and version so a plan states which package
    supplied the semantics it was built against.
    """
    from .plugin_interface import driver_context

    ambiguous = ambiguous_plugin_names().get(name)
    if ambiguous is not None:
        raise KeyError(
            f"driverFOAM plugin name {name!r} is claimed by more than one "
            f"installed distribution ({', '.join(ambiguous)}); uninstall one "
            "or select it with the module:Class form"
        )
    entry_point = discover_plugins().get(name)
    if entry_point is None:
        raise KeyError(
            f"No installed driverFOAM plugin named {name!r} in entry-point "
            f"group {ENTRY_POINT_GROUP!r}"
        )
    dist = getattr(entry_point, "dist", None)
    source = (
        f"entry-point:{dist.name}={dist.version}"
        if dist is not None
        else f"entry-point:{name}"
    )
    plugin_class = entry_point.load()
    return driver_context(plugin_class(), source=source)
