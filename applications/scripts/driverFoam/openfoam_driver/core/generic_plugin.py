"""Built-in no-domain plugin for generic OpenFOAM case-folder orchestration.

This module is the **canonical scaffold for new solver plugins**. To add support
for a new solver to driverFOAM:

  1. Copy this file to your own package as ``my_solver_plugin.py``.
  2. Fill in the identity properties (``plugin_name``, ``plugin_id``, etc.).
  3. Author a ``plugin.yaml`` (see ``generic-plugin.yaml`` for the template).
  4. Register an entry-point in your ``pyproject.toml`` under the
     ``[project.entry-points."driverfoam.plugins"]`` group.
  5. Follow the full step-by-step guide in:
     ``.agents/skills/driverfoam-plugin-builder/SKILL.md``

See also ``core/plugin_interface.py`` for the full ``SolverPlugin`` /
``SolverPluginOptionalHooks`` contracts.
"""

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
        return "2"

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
        """Advertise a real, if empty, accept-surface.

        Built through the same assembler the cardiac plugin uses, so a generic
        plan carries an ``allowed_commands`` block naming the solver-neutral
        OpenFOAM commands it may actually run. Returning a hand-written stub
        here previously omitted that block entirely and hardcoded cardiac
        region names (``electro``/``solid``) for a plugin that has neither.
        """
        from openfoam_driver.core.capability_manifest import build_capability_manifest

        return build_capability_manifest(
            plugin_commands=self.get_solver_commands() | self.get_auxiliary_commands(),
            utility_manifests=self.get_utility_manifests(),
            samplable_fields=self.get_samplable_fields({}),
        )

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
        """No solver semantics means no solver binary to authorize."""
        return frozenset()

    def get_auxiliary_commands(self) -> frozenset[str]:
        """No solver semantics means no plugin-specific helpers either."""
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

    def get_solve_step_commands(self) -> frozenset:
        """No solver semantics means no step is a solve step."""
        return frozenset()

    def get_telemetry_source_globs(self, command: str) -> tuple:
        """No declared solver log locations."""
        del command
        return ()

    def get_extra_provenance_paths(self, case_root) -> tuple:
        """No plugin-owned inputs beyond the case itself."""
        del case_root
        return ()

    def get_artifact_value_reader(self, artifact_format: str):
        """No solver-specific artifact formats to read."""
        del artifact_format
        return None

    def get_dict_entry_catalog(self):
        """No plugin documents, so no document names."""
        return {}

    def get_named_catalogs(self):
        """No solver semantics means no named catalogs to expose."""
        return {}

    def get_override_scopes(self):
        """No solver semantics means no $TOKEN. override scopes to declare."""
        return ()

    def build_run_document_config(self, spec):
        del spec
        # Return an empty config — the generic plugin imposes no key structure.
        # (Contrast with the cardiac plugin which returns anatomy/physics/stimulus/solver.)
        return {}, ()


    def get_run_document_config_schema(self) -> dict:
        """No solver semantics means no constraint on the config shape."""
        return {"type": "object", "additionalProperties": True}
