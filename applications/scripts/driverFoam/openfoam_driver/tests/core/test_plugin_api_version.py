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
#     test_plugin_api_version
#
# Description
#     The plugin API version is an enforced contract, not a free-form string.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

from __future__ import annotations

import pytest

from openfoam_driver.core.generic_plugin import GenericOpenFOAMPlugin
from openfoam_driver.core.plugin_interface import (
    SUPPORTED_PLUGIN_API_VERSIONS,
    default_driver_context,
    driver_context,
    generic_openfoam_context,
)


def test_supported_version_is_two() -> None:
    assert SUPPORTED_PLUGIN_API_VERSIONS == frozenset({"2"})


def test_builtin_plugins_are_v2() -> None:
    assert default_driver_context().identity.api_version == "2"
    assert generic_openfoam_context().identity.api_version == "2"


def test_unsupported_version_is_rejected_before_any_catalog_runs() -> None:
    class FuturePlugin(GenericOpenFOAMPlugin):
        @property
        def plugin_api_version(self) -> str:
            return "99"

        def get_dict_entries(self):
            raise AssertionError("must be rejected before catalogs are read")

        def get_profile(self):
            raise AssertionError("must be rejected before the profile is read")

    with pytest.raises(TypeError, match="99"):
        driver_context(FuturePlugin(), source="test")


def test_declaring_the_contract_without_implementing_it_is_rejected() -> None:
    """A version string is not a contract unless the shape is checked. Without
    this, a plugin claiming to speak the contract while missing a required
    member would fail only much later, deep inside whichever core module
    first called the missing method."""

    class HalfMigratedPlugin(GenericOpenFOAMPlugin):
        # Drops one required member.
        get_artifact_value_reader = None

    with pytest.raises(TypeError, match="does not implement the plugin contract"):
        driver_context(HalfMigratedPlugin(), source="test")


def test_the_shape_check_names_what_is_missing() -> None:
    class MissingTwo(GenericOpenFOAMPlugin):
        get_solve_step_commands = None
        get_utility_roots = None

    with pytest.raises(TypeError) as excinfo:
        driver_context(MissingTwo(), source="test")
    message = str(excinfo.value)
    assert "get_solve_step_commands" in message
    assert "get_utility_roots" in message


def test_both_builtin_plugins_satisfy_the_full_protocol() -> None:
    """cardiac AND generic exercise every required capability. The generic
    plugin previously declared v2 while implementing 8 of 12, riding the
    adapter's degrade-to-empty fallback."""
    from openfoam_driver.core.plugin_interface import SolverPlugin
    from openfoam_driver.plugins.cardiacfoam_plugin import CardiacFoamPlugin

    for plugin in (CardiacFoamPlugin(), GenericOpenFOAMPlugin()):
        assert isinstance(plugin, SolverPlugin), type(plugin).__name__
