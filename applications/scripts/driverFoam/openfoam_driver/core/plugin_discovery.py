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


def discover_plugins() -> dict[str, Any]:
    """Return installed plugin entry points keyed by name.

    Never raises: a broken third-party distribution must not make the CLI
    unusable for everyone else.
    """
    return {entry_point.name: entry_point for entry_point in _entry_points()}


def load_discovered_plugin(name: str):
    """Load and validate a discovered plugin by entry-point name.

    Loading executes the plugin's Python code. Identity provenance records the
    installing distribution's name and version so a plan states which package
    supplied the semantics it was built against.
    """
    from .plugin_interface import driver_context

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
