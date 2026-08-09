"""A genuine no-domain plugin used by the cross-plugin contract tests."""

from __future__ import annotations

from pathlib import Path

from openfoam_driver.core.plugin_profile import PluginProfile


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
        return "1"

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
