from __future__ import annotations

from pathlib import Path

import pytest

from openfoam_driver.core.generic_plugin import GenericOpenFOAMPlugin
from openfoam_driver.core.plugin_profile import PluginProfile, load_plugin_profile
from openfoam_driver.plugins.cardiacfoam_plugin import CardiacFoamPlugin


def test_cardiac_profile_declares_case_files_and_cxx_provenance() -> None:
    profile = CardiacFoamPlugin().get_profile()

    assert profile.plugin_id == "org.cardiacfoam"
    assert {rule.path for rule in profile.case_files} >= {
        "system/controlDict",
        "constant/physicsProperties",
        "constant/electroProperties",
    }
    assert profile.cxx_mapping is not None
    assert profile.cxx_mapping.allowlist_path.is_file()
    assert all(path.is_dir() for path in profile.cxx_mapping.source_roots)
    assert profile.digest.startswith("sha256:")


def test_cardiac_catalog_partitions_entries_by_document() -> None:
    catalog = CardiacFoamPlugin().get_dictionary_catalog()

    assert {"electroProperties", "physicsProperties", "controlDict"} <= set(catalog.documents)
    assert {entry.driver_path for entry in catalog.entries_for("physicsProperties")} == {"type"}
    assert {entry.driver_path for entry in catalog.entries_for("controlDict")} >= {"deltaT", "endTime"}


def test_generic_profile_declares_no_solver_specific_files() -> None:
    profile = GenericOpenFOAMPlugin().get_profile()

    assert profile.plugin_id == "org.driverfoam.generic-openfoam"
    assert profile.case_files == ()
    assert profile.cxx_mapping is None


def test_profile_rejects_case_path_escape(tmp_path: Path) -> None:
    path = tmp_path / "plugin.yaml"
    path.write_text(
        """schema_version: 1
plugin: {id: example.bad, api_version: '1'}
case_profile:
  dictionaries:
    - path: ../outside
      kind: openfoam_dictionary
      role: plugin.configuration
      required: always
"""
    )

    with pytest.raises(ValueError, match="escapes the case"):
        load_plugin_profile(path)


def test_profile_digest_is_stable_after_payload_mutation() -> None:
    payload = {
        "schema_version": 1,
        "plugin": {"id": "example.profile", "api_version": "1"},
        "case_profile": {"dictionaries": []},
    }
    profile = PluginProfile(
        path=Path("example-plugin.yaml"),
        plugin_id="example.profile",
        api_version="1",
        case_files=(),
        cxx_mapping=None,
        payload=payload,
    )

    digest = profile.digest
    payload["plugin"]["id"] = "example.mutated"

    assert profile.digest == digest
