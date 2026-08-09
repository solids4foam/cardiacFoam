from __future__ import annotations

import pytest
from pathlib import Path

from openfoam_driver.core.plugin_interface import driver_context, validate_plugin
from openfoam_driver.core.plugin_profile import PluginProfile
from openfoam_driver.core.runtime.registry import list_tutorials
from openfoam_driver.core.contracts.dictionary import DictEntry
from openfoam_driver.core.contracts.dictionary_catalog import DictionaryCatalog


class _Plugin:
    def __init__(self, plugin_id: str, tutorial_name: str) -> None:
        self._plugin_id = plugin_id
        self._tutorial_name = tutorial_name

    @property
    def plugin_name(self) -> str:
        return self._plugin_id

    @property
    def plugin_id(self) -> str:
        return self._plugin_id

    @property
    def plugin_version(self) -> str:
        return "1.0.0"

    @property
    def plugin_api_version(self) -> str:
        return "1"

    def get_profile(self):
        return PluginProfile(
            path=Path("test-plugin.yaml"),
            plugin_id=self._plugin_id,
            api_version="1",
            case_files=(),
            cxx_mapping=None,
            payload={
                "schema_version": 1,
                "plugin": {"id": self._plugin_id, "api_version": "1"},
                "case_profile": {"dictionaries": []},
            },
        )

    def get_dict_entries(self):
        return ()

    def get_dictionary_catalog(self):
        return DictionaryCatalog({})

    def get_dict_groups(self):
        return {}

    def get_capabilities(self):
        return {}

    def get_tutorial_catalog(self):
        return {"registered_tutorials": (self._tutorial_name,), "spec_factories": {}}

    def get_tutorial_displays(self):
        return ()

    def validate_configuration(self, spec):
        return ()

    def validate_run_semantics(self, context):
        return ()

    def predict_data_artifacts(self, case_root, spec):
        return ()


def test_contexts_do_not_share_plugin_selection() -> None:
    alpha = driver_context(_Plugin("example.alpha", "alpha"), source="test")
    beta = driver_context(_Plugin("example.beta", "beta"), source="test")

    assert list_tutorials(alpha) == ["alpha"]
    assert list_tutorials(beta) == ["beta"]
    assert alpha.identity.to_json()["id"] == "example.alpha"
    assert beta.identity.to_json()["id"] == "example.beta"


def test_plugin_contract_rejects_missing_members() -> None:
    with pytest.raises(TypeError, match="missing required members"):
        validate_plugin(object())


def test_plugin_contract_rejects_an_invalid_stable_id() -> None:
    with pytest.raises(TypeError, match="plugin_id must use lowercase"):
        validate_plugin(_Plugin("Example Plugin", "example"))


def test_driver_context_rejects_duplicate_catalog_paths() -> None:
    plugin = _Plugin("example.duplicates", "duplicates")
    plugin.get_dict_entries = lambda: (
        DictEntry(driver_path="shared", description="first"),
        DictEntry(driver_path="shared", description="second"),
    )

    with pytest.raises(TypeError, match="duplicate paths: shared"):
        driver_context(plugin, source="test")
