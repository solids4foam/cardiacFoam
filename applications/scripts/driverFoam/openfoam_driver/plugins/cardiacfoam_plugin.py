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
#     cardiacfoam_plugin
#
# Description
#     Implements the SolverPlugin interface for the cardiacFoam solver.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

from __future__ import annotations

from typing import TYPE_CHECKING
from functools import lru_cache
from pathlib import Path
from openfoam_driver.core.plugin_interface import SolverPlugin, CapabilityManifest

# Note: now imported from the local plugin catalog instead of dict_entries
from openfoam_driver.plugins.cardiacfoam.dict_entries_catalog import ELECTRO_PROPERTY_ENTRY_GROUPS, HETEROGENEITY_MODELS
from openfoam_driver.plugins.cardiacfoam.common_dict_entries import (
    CONTROL_DICT_ENTRIES,
    PHYSICS_PROPERTY_ENTRIES,
)
from openfoam_driver.core.contracts.dictionary_catalog import DictionaryCatalog
from openfoam_driver.plugins.cardiacfoam.active_tension_catalog import ACTIVE_TENSION_MODEL_CATALOG
from openfoam_driver.plugins.cardiacfoam.ionic_model_catalog import IONIC_MODEL_CATALOG
from openfoam_driver.core.capability_manifest import build_capability_manifest
from openfoam_driver.plugins.cardiacfoam.solver_coupling import SOLVER_COMPATIBILITY_RULES
from openfoam_driver.core.runtime.registry import list_tutorials
from openfoam_driver.core.planning_types import StrictDiagnostic, diagnostic

if TYPE_CHECKING:
    from openfoam_driver.core.runtime.models import TutorialSpec, CaseConfig, DataArtifact
    from openfoam_driver.core.tutorials_display import TutorialDisplay
    from pathlib import Path


class CardiacFoamPlugin:
    """
    Plugin implementation for cardiacFoam.
    Provides domain-specific dictionaries, tutorials, and capabilities 
    to the generic driverFOAM engine.
    """
    
    @property
    def plugin_name(self) -> str:
        return "cardiacFoam"

    @property
    def plugin_id(self) -> str:
        return "org.cardiacfoam"

    @property
    def plugin_version(self) -> str:
        return "0.1.0"

    @property
    def plugin_api_version(self) -> str:
        return "2"

    @staticmethod
    @lru_cache(maxsize=1)
    def get_profile():
        from openfoam_driver.core.plugin_profile import load_plugin_profile

        return load_plugin_profile(Path(__file__).parent / "cardiacfoam" / "plugin.yaml")

    def configure_execution_environment(self, env: dict[str, str]):
        """Apply the plugin's declared backend and build-manifest contract."""
        from openfoam_driver.plugins.cardiacfoam.runtime_profile import (
            configure_runtime_environment,
        )

        return configure_runtime_environment(env)

    def get_openfoam_bashrc(self, env: dict[str, str]) -> str | None:
        from openfoam_driver.plugins.cardiacfoam.runtime_profile import (
            configured_openfoam_bashrc,
        )

        return configured_openfoam_bashrc(env)
        
    def get_dict_groups(self) -> dict[str, tuple[DictEntry, ...]]:
        """
        Return the dictionary entries organized by logical group.
        """
        return ELECTRO_PROPERTY_ENTRY_GROUPS

    def get_dict_entries(self) -> tuple[DictEntry, ...]:
        """
        Aggregate and return all dictionary entries specific to cardiacFoam.
        """
        entries: list[DictEntry] = list(PHYSICS_PROPERTY_ENTRIES)
        for group in self.get_dict_groups().values():
            entries.extend(group)
        return tuple(entries)

    @staticmethod
    @lru_cache(maxsize=1)
    def get_dictionary_catalog() -> DictionaryCatalog:
        electro_entries: list[DictEntry] = []
        for group in ELECTRO_PROPERTY_ENTRY_GROUPS.values():
            electro_entries.extend(group)
        return DictionaryCatalog({
            "electroProperties": tuple(electro_entries),
            "physicsProperties": PHYSICS_PROPERTY_ENTRIES,
            "controlDict": CONTROL_DICT_ENTRIES,
        })

    def get_capabilities(self) -> CapabilityManifest:
        """
        Return the cardiacFoam capabilities (models, solvers, etc.).
        """
        # No case_root is available at this call site, so this resolves to
        # the fixed solver fields only (matches historical behaviour: no
        # resolved model, no ionic/active-tension-specific field names).
        resolved: dict = {}
        manifest = build_capability_manifest(
            # The manifest advertises the accept-surface, so it lists both
            # kinds of authorized plugin command -- the solver/auxiliary split
            # only governs who may be credited with a run's artifacts.
            plugin_commands=self.get_solver_commands() | self.get_auxiliary_commands(),
            utility_manifests=self.get_utility_manifests(),
            samplable_fields=self.get_samplable_fields(resolved),
        )
        manifest["heterogeneity_models"] = HETEROGENEITY_MODELS
        manifest["ionic_models"] = IONIC_MODEL_CATALOG
        manifest["active_tension_models"] = ACTIVE_TENSION_MODEL_CATALOG
        manifest["solver_compatibility_rules"] = SOLVER_COMPATIBILITY_RULES
        return manifest

    def resolve_case_models(self, case_root: Path) -> dict:
        """Best-effort ``{"solver", "ionic_model", "active_tension"}`` from a
        case's ``constant/electroProperties``. Never raises."""
        from openfoam_driver.plugins.cardiacfoam.case_introspection import (
            resolve_case_models,
        )

        return resolve_case_models(case_root)

    def get_samplable_fields(self, resolved: dict) -> dict:
        """Field names the resolved cardiac model exposes, by region."""
        from openfoam_driver.plugins.cardiacfoam.case_introspection import (
            samplable_fields,
        )

        return samplable_fields(resolved)

    def get_override_schema(self, tutorial_name: str, make_spec_info: dict) -> dict:
        """cardiacFoam's authored --config schema, including a worked example."""
        from openfoam_driver.plugins.cardiacfoam.override_schema import config_schema

        return config_schema(tutorial_name, make_spec_info)

    def get_solve_step_commands(self) -> frozenset:
        """Commands that actually run the solver (Phase 4 telemetry)."""
        from openfoam_driver.plugins.cardiacfoam.runtime_evidence import (
            solve_step_commands,
        )

        return solve_step_commands()

    def get_telemetry_source_globs(self, command: str) -> tuple:
        """Where this command's solver log lands beyond captured stdout."""
        from openfoam_driver.plugins.cardiacfoam.runtime_evidence import (
            telemetry_source_globs,
        )

        return telemetry_source_globs(command)

    def get_extra_provenance_paths(self, case_root) -> tuple:
        """Extra inputs Phase 2 must digest beyond system/ and constant/."""
        from openfoam_driver.plugins.cardiacfoam.runtime_evidence import (
            extra_provenance_paths,
        )

        return extra_provenance_paths(case_root)

    def get_artifact_value_reader(self, artifact_format: str):
        """Reader for a cardiac artifact format, or None (Phase 5)."""
        from openfoam_driver.plugins.cardiacfoam.runtime_evidence import (
            artifact_value_reader,
        )

        return artifact_value_reader(artifact_format)

    def get_required_inputs(self, case_root, resolved_case, selected_start_time) -> tuple:
        """Model-dependent required inputs (CaseProvenanceCapability). See
        ``case_provenance.py`` for why this defers to the safe default."""
        from openfoam_driver.plugins.cardiacfoam.case_provenance import (
            required_inputs,
        )

        return required_inputs(case_root, resolved_case, selected_start_time)

    def get_generated_output_globs(self, case_root, resolved_case, selected_start_time) -> tuple:
        """Fixed mesh-diagnostic outputs nothing in src/ or applications/ reads."""
        from openfoam_driver.plugins.cardiacfoam.case_provenance import (
            generated_output_globs,
        )

        return generated_output_globs(case_root, resolved_case, selected_start_time)

    def get_dict_entry_catalog(self) -> dict:
        """Dictionary entries arranged by cardiacFoam's own document names."""
        from openfoam_driver.plugins.cardiacfoam.override_schema import (
            dict_entry_catalog,
        )

        return dict_entry_catalog(
            self.get_dictionary_catalog(), self.get_dict_groups(),
        )

    def get_named_catalogs(self) -> dict:
        """This plugin's own catalogs -- ionic models and active-tension
        models -- namespaced under introspection's generic
        ``plugin_catalogs`` key instead of core-hardcoded field names."""
        from openfoam_driver.plugins.cardiacfoam.named_catalogs import (
            named_catalogs,
        )

        return named_catalogs(self.get_capabilities())

    def get_override_scopes(self) -> tuple:
        """This plugin's one `step --strict --apply` override scope:
        $ELECTRO_MODEL_COEFFS -> constant/electroProperties."""
        from openfoam_driver.plugins.cardiacfoam.overrides import (
            electro_model_coeffs_scope,
        )

        return (electro_model_coeffs_scope(),)

    def get_regeneration_scopes(self) -> tuple:
        """This plugin's one `step --strict --apply` regeneration scope:
        myocardiumSolver -> constant/electroProperties. Switching
        myocardiumSolver renames the active <solver>Coeffs sub-block and
        changes which sibling keys are legal, so it needs a full rebuild
        rather than the key-patch $ELECTRO_MODEL_COEFFS route above."""
        from openfoam_driver.plugins.cardiacfoam.overrides import (
            electro_properties_regeneration_scope,
        )

        return (electro_properties_regeneration_scope(),)

    def get_tutorial_catalog(self) -> dict:
        from openfoam_driver.plugins.cardiacfoam.tutorials.registry import SPEC_FACTORIES, REGISTERED_TUTORIALS
        from openfoam_driver.core.runtime.generic_case import make_spec as make_generic_case_spec
        return {
            "spec_factories": SPEC_FACTORIES,
            "registered_tutorials": REGISTERED_TUTORIALS,
            "make_generic_case_spec": make_generic_case_spec
        }

    def get_tutorial_displays(self) -> tuple[TutorialDisplay, ...]:
        from openfoam_driver.plugins.cardiacfoam.tutorials.display import TUTORIALS
        return TUTORIALS

    def validate_configuration(self, spec: TutorialSpec) -> tuple[StrictDiagnostic, ...]:
        from pathlib import Path
        from openfoam_driver.plugins.cardiacfoam.detection import (
            detect_myocardium_solver_name,
            detect_ionic_model_name,
        )
        from openfoam_driver.core.planning_types import diagnostic as _diagnostic

        diagnostics = []
        case_root = Path(spec.case_root)
        electro_path = case_root / "constant" / "electroProperties"

        if electro_path.exists():
            try:
                solver = detect_myocardium_solver_name(electro_path)
                if solver not in {"singleCellSolver", "monodomainSolver", "bidomainSolver", "eikonalSolver"}:
                    diagnostics.append(_diagnostic(
                        "error",
                        "unknown_solver",
                        f"No strict artifact handler is registered for myocardiumSolver {solver!r}.",
                        source=str(electro_path),
                        field="myocardiumSolver",
                    ))
            except KeyError as exc:
                diagnostics.append(_diagnostic("error", "missing_solver", str(exc), source=str(electro_path)))

            try:
                ionic_model = detect_ionic_model_name(electro_path)
            except KeyError:
                ionic_model = None
            
            capabilities = self.get_capabilities()
            if ionic_model is not None and ionic_model not in capabilities.get("ionic_models", {}):
                diagnostics.append(_diagnostic(
                    "error",
                    "unknown_ionic_model",
                    f"Ionic model {ionic_model!r} is not supported by the active plugin..",
                    source=str(electro_path),
                    field="ionicModel",
                ))

            from openfoam_driver.plugins.cardiacfoam.validation import (
                _evaluate_pvj_resistance_requirement,
            )
            diagnostics.extend(
                _evaluate_pvj_resistance_requirement(case_root, electro_path)
            )

        return tuple(diagnostics)

    def validate_run_semantics(self, context):
        """Apply cardiacFoam's cross-field rules after core validation."""
        from openfoam_driver.plugins.cardiacfoam.validation import (
            _evaluate_block_references,
            _evaluate_dynamic_required_fields,
            _evaluate_heterogeneity,
            _evaluate_solver_coupling,
            _evaluate_tissue_compatibility,
        )

        return tuple(
            _evaluate_solver_coupling(context)
            + _evaluate_block_references(context)
            + _evaluate_dynamic_required_fields(context)
            + _evaluate_heterogeneity(context)
            + _evaluate_tissue_compatibility(context)
        )

    def predict_data_artifacts(self, case_root: Path, spec: TutorialSpec) -> tuple[DataArtifact, ...]:
        from openfoam_driver.plugins.cardiacfoam.artifacts_predictor import predict_cardiac_artifacts
        return predict_cardiac_artifacts(case_root, spec)

    def build_run_document_config(self, spec):
        from openfoam_driver.plugins.cardiacfoam.run_document_config import build_config

        return build_config(spec)

    def get_run_document_config_schema(self) -> dict:
        """cardiacFoam's RunDocument.config JSON Schema (the phase vocabulary)."""
        from openfoam_driver.plugins.cardiacfoam.config_schema import (
            get_run_document_config_schema,
        )

        return get_run_document_config_schema()

    def has_case_marker(self, case_root: Path) -> bool:
        """Return the historical cardiac case-folder discovery evidence."""
        from openfoam_driver.plugins.cardiacfoam.case_compatibility import has_case_marker

        return has_case_marker(case_root)

    def is_case_runnable_without_workflow(self, case_root: Path) -> bool:
        """Preserve legacy runnability for uncontracted cardiac cases."""
        from openfoam_driver.plugins.cardiacfoam.case_compatibility import (
            is_runnable_without_workflow,
        )

        return is_runnable_without_workflow(case_root)

    def route_sweep_case_values(
        self,
        *,
        base,
        resolved_axis_values,
        driver_context,
    ):
        from openfoam_driver.plugins.cardiacfoam.sweep import route_case_values

        return route_case_values(
            base=base,
            resolved_axis_values=resolved_axis_values,
            driver_context=driver_context,
        )

    def materialize_sweep_case(self, *, case_dir: Path, routed) -> None:
        from openfoam_driver.plugins.cardiacfoam.sweep import materialize_case

        materialize_case(case_dir=case_dir, routed=routed)

    def get_solver_commands(self) -> frozenset[str]:
        """This plugin's artifact-producing solver commands."""
        from openfoam_driver.plugins.cardiacfoam.command_authorization import (
            solver_commands,
        )

        return solver_commands()

    def get_auxiliary_commands(self) -> frozenset[str]:
        """Authorized plugin commands that do not produce the run's artifacts."""
        from openfoam_driver.plugins.cardiacfoam.command_authorization import (
            auxiliary_commands,
        )

        return auxiliary_commands()

    def get_utility_manifests(self) -> dict:
        """This plugin's ``utility.manifest.toml`` sidecars, by command name."""
        from openfoam_driver.plugins.cardiacfoam.command_authorization import (
            utility_manifests,
        )

        # utility_manifests() is cached and returns a read-only view; copy so a
        # caller mutating what it gets back cannot reach the shared cache.
        return dict(utility_manifests())

    def get_utility_roots(self) -> tuple[Path, ...]:
        """Roots searched for this plugin's utility manifests."""
        from openfoam_driver.plugins.cardiacfoam.command_authorization import (
            utility_roots,
        )

        return utility_roots()

    def is_nondimensional_case(self, spec) -> bool:
        from openfoam_driver.plugins.cardiacfoam.planning_policy import (
            is_nondimensional_case,
        )

        return is_nondimensional_case(spec)

    def get_mesh_geometry_diagnostics(self, case_root: Path) -> tuple:
        """Scale checks for point sets core's region discovery cannot see.

        cardiacFoam cases may carry a Purkinje conduction tree in
        ``constant/purkinjeGraph*`` -- a Foam dictionary with its own point
        list, not a mesh region.
        """
        from openfoam_driver.plugins.cardiacfoam.mesh_geometry import (
            purkinje_graph_diagnostics,
        )

        return purkinje_graph_diagnostics(case_root)
