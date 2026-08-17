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
#     named_catalogs
#
# Description
#     cardiacFoam's own named catalogs (ionic models, active-tension models),
#     exposed through the plugin's ``get_named_catalogs()`` hook and namespaced
#     generically under introspection's ``plugin_catalogs`` key -- core no
#     longer hardcodes the field names ``ionic_model_catalog``/
#     ``active_tension_catalog``. Kept as a plain function, not a method, so
#     both the plugin's real hook and the v1-compatibility fallback in
#     ``core/compatibility.py`` share one authored shape, the same split
#     ``override_schema.py`` uses for ``config_schema``/``dict_entry_catalog``.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

from __future__ import annotations

from typing import Any


def named_catalogs(capabilities: dict[str, Any]) -> dict[str, Any]:
    """Return the plugin's ionic-model and active-tension catalogs.

    ``capabilities`` is the plugin's own ``get_capabilities()`` manifest,
    which already carries ``ionic_models``/``active_tension_models``/
    ``solver_compatibility_rules``. Values are returned unserialized --
    core owns serialization.
    """
    return {
        "ionic_model_catalog": {
            "schema_version": "1.0",
            "ionic_models": dict(capabilities.get("ionic_models", {})),
            "solver_compatibility": list(
                capabilities.get("solver_compatibility_rules", ())
            ),
        },
        "active_tension_catalog": {
            "schema_version": "1.0",
            "active_tension_models": dict(
                capabilities.get("active_tension_models", {})
            ),
        },
    }
