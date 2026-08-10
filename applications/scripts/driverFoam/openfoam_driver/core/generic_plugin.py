"""Built-in no-domain plugin for generic OpenFOAM case-folder orchestration."""

from __future__ import annotations

from functools import lru_cache
from pathlib import Path

from .plugin_profile import load_plugin_profile
from .contracts.dictionary_catalog import DictionaryCatalog


class GenericOpenFOAMPlugin:
    """Provides no solver semantics beyond the core OpenFOAM execution model."""

    @property
    def plugin_name(self) -> str:
        return "generic OpenFOAM"

    @property
    def plugin_id(self) -> str:
        return "org.driverfoam.generic-openfoam"

    @property
    def plugin_version(self) -> str:
        return "1"

    @property
    def plugin_api_version(self) -> str:
        return "1"

    @staticmethod
    @lru_cache(maxsize=1)
    def get_profile():
        return load_plugin_profile(Path(__file__).with_name("generic-plugin.yaml"))

    def get_dict_entries(self):
        return ()

    def get_dictionary_catalog(self):
        return DictionaryCatalog({})

    def get_dict_groups(self):
        return {}

    def get_capabilities(self):
        return {"samplable_fields": {"electro": [], "solid": []}}

    def get_tutorial_catalog(self):
        return {"registered_tutorials": (), "spec_factories": {}}

    def get_tutorial_displays(self):
        return ()

    def validate_configuration(self, spec):
        return ()

    def validate_run_semantics(self, context):
        return ()

    def predict_data_artifacts(self, case_root, spec):
        return ()

    def build_run_document_config(self, spec):
        del spec
        return {
            "anatomy": {},
            "physics": {},
            "stimulus": {},
            "solver": {},
        }, ()
