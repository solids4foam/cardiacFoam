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
#     test_plugin_discovery
#
# Description
#     Installed plugins are discoverable; module:Class remains the dev path.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

from __future__ import annotations

import pytest

from openfoam_driver.core import plugin_discovery
from openfoam_driver.core.plugin_interface import load_plugin_context


class _FakeEntryPoint:
    name = "fakeplugin"
    value = "openfoam_driver.core.generic_plugin:GenericOpenFOAMPlugin"
    dist = type("D", (), {"name": "fake-dist", "version": "9.9"})()

    def load(self):
        from openfoam_driver.core.generic_plugin import GenericOpenFOAMPlugin

        return GenericOpenFOAMPlugin


def test_a_colon_still_means_a_trusted_local_import() -> None:
    context = load_plugin_context(
        "openfoam_driver.core.generic_plugin:GenericOpenFOAMPlugin"
    )
    assert context.identity.id == "org.driverfoam.generic-openfoam"
    assert context.identity.source.startswith("trusted-import:")


def test_an_unknown_discovered_id_raises_keyerror() -> None:
    with pytest.raises(KeyError, match="definitelyNotInstalled"):
        load_plugin_context("definitelyNotInstalled")


def test_discovery_reads_the_driverfoam_plugins_group(monkeypatch) -> None:
    monkeypatch.setattr(
        plugin_discovery, "_entry_points", lambda: (_FakeEntryPoint(),)
    )
    assert "fakeplugin" in plugin_discovery.discover_plugins()
    context = plugin_discovery.load_discovered_plugin("fakeplugin")
    # Identity provenance records the installing distribution, so a plan says
    # which package supplied the semantics it was built against.
    assert context.identity.source == "entry-point:fake-dist=9.9"


def test_a_discovered_id_wins_only_when_there_is_no_colon(monkeypatch) -> None:
    monkeypatch.setattr(
        plugin_discovery, "_entry_points", lambda: (_FakeEntryPoint(),)
    )
    # A colon always means the trusted import form, never discovery.
    context = load_plugin_context(
        "openfoam_driver.core.generic_plugin:GenericOpenFOAMPlugin"
    )
    assert context.identity.source.startswith("trusted-import:")


def test_discovery_is_empty_by_default_and_does_not_raise() -> None:
    # No third-party plugin is installed in this repository's environment.
    assert isinstance(plugin_discovery.discover_plugins(), dict)


class _RivalEntryPoint(_FakeEntryPoint):
    """A second distribution claiming the same entry-point name."""

    dist = type("D", (), {"name": "rival-dist", "version": "0.1"})()


def test_a_name_claimed_by_two_distributions_is_reported_not_resolved(
    monkeypatch,
) -> None:
    """Insertion order must not silently decide which plugin wins -- that
    would depend on installation order and be invisible in the plan."""
    monkeypatch.setattr(
        plugin_discovery,
        "_entry_points",
        lambda: (_FakeEntryPoint(), _RivalEntryPoint()),
    )
    ambiguous = plugin_discovery.ambiguous_plugin_names()
    assert ambiguous == {"fakeplugin": ("fake-dist=9.9", "rival-dist=0.1")}
    # The ambiguous name is withheld from discovery rather than resolved.
    assert "fakeplugin" not in plugin_discovery.discover_plugins()


def test_loading_an_ambiguous_name_fails_loudly(monkeypatch) -> None:
    monkeypatch.setattr(
        plugin_discovery,
        "_entry_points",
        lambda: (_FakeEntryPoint(), _RivalEntryPoint()),
    )
    with pytest.raises(KeyError, match="claimed by more than one"):
        plugin_discovery.load_discovered_plugin("fakeplugin")
