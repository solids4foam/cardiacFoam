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

from pathlib import Path

import pytest

from openfoam_driver.core.generic_plugin import GenericOpenFOAMPlugin
from openfoam_driver.core.plugin_interface import (
    SUPPORTED_PLUGIN_API_VERSIONS,
    default_driver_context,
    driver_context,
    generic_openfoam_context,
)
from openfoam_driver.core.plugin_profile import load_plugin_profile

_LEGACY_PROFILE = Path(__file__).parents[1] / "fixtures" / "legacy_v1_plugin.yaml"


class _LegacyV1Plugin(GenericOpenFOAMPlugin):
    """A third-party v1 plugin: no v2 members, own identity."""

    @property
    def plugin_id(self) -> str:
        return "org.example.legacy"

    @property
    def plugin_api_version(self) -> str:
        return "1"

    def get_profile(self):
        return load_plugin_profile(_LEGACY_PROFILE)


def test_supported_versions_are_one_and_two() -> None:
    assert SUPPORTED_PLUGIN_API_VERSIONS == frozenset({"1", "2"})


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


def test_a_v1_plugin_still_loads_through_compatibility() -> None:
    context = driver_context(_LegacyV1Plugin(), source="test")
    assert context.identity.api_version == "1"
    # A third-party v1 plugin declares no solver commands, and core must not
    # invent a cardiac-shaped default for it.
    assert context.capabilities.command_authorization.solver_commands() == frozenset()
    assert context.capabilities.override_schema.config_schema("x", {}) == {}
