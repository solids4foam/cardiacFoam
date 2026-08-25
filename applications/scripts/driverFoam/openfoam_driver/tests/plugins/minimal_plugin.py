"""A genuine no-domain plugin used by the cross-plugin contract tests."""

from __future__ import annotations

from pathlib import Path

from openfoam_driver.core.plugin_profile import PluginProfile
from openfoam_driver.core.contracts.dictionary_catalog import DictionaryCatalog


class MinimalOpenFOAMPlugin:
    """Implements only the required plugin contract; adds no solver meaning."""

    @property
    def plugin_name(self) -> str:
        return "minimal OpenFOAM test plugin"

    @property
    def plugin_id(self) -> str:
        return "org.driverfoam.test-minimal"

    @property
    def plugin_version(self) -> str:
        return "1.0.0"

    @property
    def plugin_api_version(self) -> str:
        return "2"

    def get_profile(self) -> PluginProfile:
        return PluginProfile(
            path=Path(__file__),
            plugin_id=self.plugin_id,
            api_version=self.plugin_api_version,
            case_files=(),
            cxx_mapping=None,
            payload={
                "schema_version": 1,
                "plugin": {
                    "id": self.plugin_id,
                    "api_version": self.plugin_api_version,
                },
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
        return {"registered_tutorials": (), "spec_factories": {}}

    def get_tutorial_displays(self):
        return ()

    def validate_configuration(self, spec):
        return ()

    def validate_run_semantics(self, context):
        return ()

    def predict_data_artifacts(self, case_root, spec):
        return ()

    def get_solver_commands(self) -> frozenset[str]:
        return frozenset()

    def get_auxiliary_commands(self) -> frozenset[str]:
        return frozenset()

    def get_utility_manifests(self) -> dict:
        return {}

    def get_utility_roots(self) -> tuple[Path, ...]:
        return ()

    def resolve_case_models(self, case_root):
        del case_root
        return {}

    def get_samplable_fields(self, resolved):
        del resolved
        return {}

    def get_override_schema(self, tutorial_name, make_spec_info):
        del tutorial_name, make_spec_info
        return {}

    def get_dict_entry_catalog(self):
        return {}

    def get_solve_step_commands(self) -> frozenset:
        return frozenset()

    def get_telemetry_source_globs(self, command: str) -> tuple:
        del command
        return ()

    def get_extra_provenance_paths(self, case_root) -> tuple:
        del case_root
        return ()

    def get_artifact_value_reader(self, artifact_format: str):
        del artifact_format
        return None

    def get_run_document_config_schema(self) -> dict:
        """No solver semantics means no constraint on the config shape."""
        return {"type": "object", "additionalProperties": True}
