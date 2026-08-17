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
#     dict_entries
#
# Description
#     Defines schema contracts for OpenFOAM dictionary parameters.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

from __future__ import annotations

from typing import TYPE_CHECKING
from .core.contracts.dictionary import DictEntry, Phase, build_group
if TYPE_CHECKING:
    from openfoam_driver.core.plugin_interface import DriverContext

# Ionic models that implement transmural tissue heterogeneity
# (configureIonicHeterogeneity, endo/M/epi blend and/or namedRegions) on
# CPU and/or GPU.

def get_heterogeneity_models(
    driver_context: "DriverContext | None" = None,
) -> tuple[str, ...]:
    from openfoam_driver.core.compatibility import resolve_public_driver_context

    driver_context = resolve_public_driver_context(driver_context)
    return driver_context.capabilities.manifest.manifest().get("heterogeneity_models", ())

def get_electro_property_entry_groups(
    driver_context: "DriverContext | None" = None,
) -> dict[str, tuple[DictEntry, ...]]:
    from openfoam_driver.core.compatibility import resolve_public_driver_context

    driver_context = resolve_public_driver_context(driver_context)
    return driver_context.capabilities.dictionaries.groups()

def all_documented_driver_paths(
    driver_context: "DriverContext | None" = None,
) -> tuple[str, ...]:
    from openfoam_driver.core.compatibility import resolve_public_driver_context

    driver_context = resolve_public_driver_context(driver_context)
    paths = [
        entry.driver_path
        for entry in driver_context.capabilities.dictionaries.entries()
    ]
    return tuple(dict.fromkeys(paths))


def __getattr__(name: str):
    """PEP 562 lazy resolution for deprecated cardiac-catalog re-exports.

    ``CONTROL_DICT_ENTRIES``/``PHYSICS_PROPERTY_ENTRIES`` used to be imported
    at module scope from the cardiacFoam plugin, which meant every consumer
    of this (solver-neutral) module transitively imported cardiac plugin
    internals just by importing ``dict_entries``. Resolving them lazily here
    keeps ``from openfoam_driver.dict_entries import CONTROL_DICT_ENTRIES``
    (and ``PHYSICS_PROPERTY_ENTRIES``) working unchanged for existing
    external callers while removing the module-scope import.
    """
    if name in ("CONTROL_DICT_ENTRIES", "PHYSICS_PROPERTY_ENTRIES"):
        from .plugins.cardiacfoam.common_dict_entries import (
            CONTROL_DICT_ENTRIES,
            PHYSICS_PROPERTY_ENTRIES,
        )

        return {
            "CONTROL_DICT_ENTRIES": CONTROL_DICT_ENTRIES,
            "PHYSICS_PROPERTY_ENTRIES": PHYSICS_PROPERTY_ENTRIES,
        }[name]
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
