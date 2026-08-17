"""Internal, focused capability seams for solver plugins.

The public :class:`SolverPlugin` protocol remains the compatibility contract for
Plan 1.  Core code consumes this bundle instead of reaching through
``DriverContext.plugin`` directly.  The adapters deliberately preserve the
legacy method calls, return values, call order, and exception behaviour.

Optional case-compatibility and sweep hooks let a plugin take ownership of
solver-specific behaviour without adding new required members to the public
protocol.  Plugins that do not provide those hooks retain the historical
driverFOAM fallbacks.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Protocol, TYPE_CHECKING

if TYPE_CHECKING:
    from .plugin_interface import SolverPlugin
    from .plugin_profile import CaseFileRule
    from .runtime.models import DataArtifact, TutorialSpec
    from ..planning_types import StrictDiagnostic
    from ..report_catalog import ReportDefinition


@dataclass(frozen=True)
class ConfigurationValidationRequest:
    spec: "TutorialSpec"


@dataclass(frozen=True)
class RunSemanticValidationRequest:
    context: dict[str, Any]


@dataclass(frozen=True)
class ArtifactPredictionRequest:
    case_root: Path
    spec: "TutorialSpec"


@dataclass(frozen=True)
class RunDocumentConfigurationRequest:
    spec: "TutorialSpec"


@dataclass(frozen=True)
class CaseCompatibilityRequest:
    case_root: Path


@dataclass(frozen=True)
class SweepRoutingRequest:
    base: dict[str, Any]
    resolved_axis_values: dict[str, Any]


@dataclass(frozen=True)
class SweepMaterializationRequest:
    case_dir: Path
    routed: dict[str, Any]


@dataclass(frozen=True)
class ResolvedInput:
    """One field-level input a solver plugin's case model resolves to an
    actual on-disk path -- or fails to.

    Globs are insufficient here: field *names* are dictionary-configurable
    (``fieldName``, ``ionicHeterogeneity/field``, ``bathConductivityField``)
    and field *locations* resolve by a backward ``Time::findInstance``
    search with a ``constant/`` fallback, so a field's canonical path is not
    knowable from its name alone. ``consumer`` records which model/domain
    resolved it, for diagnostics -- never for severity.
    """

    name: str
    path: Path | None
    required: bool
    consumer: str


@dataclass(frozen=True)
class RuntimeDependency:
    """One thing the workflow's *executable* consumes at run time, outside
    the case tree: the solver binary itself, a library it links or loads,
    or a case-local shared object built from sources inside the case.

    Replaces the earlier ``extra_provenance_paths() -> tuple[Path, ...]``
    stub, which could not express "this was required and I could not find
    it" -- a tuple of paths can only omit, and omission reads as "nothing to
    check". ``path is None`` on a ``required=True`` dependency must surface
    as ``unavailable`` rather than silently vanishing from the list.

    Most cases run through an ``Allrun`` script, so a workflow step's
    command fingerprints the script, never the solver binary the script
    invokes -- the exact gap that let a rebuilt solver replay a resumed
    run's previous numbers as fresh. Declaring dependencies this way, apart
    from however a step happens to be launched, is the fix.
    """

    name: str
    path: Path | None
    required: bool


class TutorialCatalogCapability(Protocol):
    def catalog(self) -> dict[str, Any]: ...
    def displays(self) -> tuple[Any, ...]: ...


class DictionaryCatalogCapability(Protocol):
    def entries(self) -> tuple[Any, ...]: ...
    def catalog(self) -> Any: ...
    def groups(self) -> dict[str, tuple[Any, ...]]: ...


class CapabilityManifestCapability(Protocol):
    def manifest(self) -> Any: ...


class ConfigurationValidatorCapability(Protocol):
    def validate(
        self, request: ConfigurationValidationRequest,
    ) -> tuple["StrictDiagnostic", ...]: ...


class RunSemanticValidatorCapability(Protocol):
    def validate(self, request: RunSemanticValidationRequest) -> tuple[Any, ...]: ...


class ArtifactPredictorCapability(Protocol):
    def predict(self, request: ArtifactPredictionRequest) -> tuple["DataArtifact", ...]: ...


class RunDocumentConfigurationCapability(Protocol):
    def build(
        self, request: RunDocumentConfigurationRequest,
    ) -> tuple[dict[str, dict[str, Any]], tuple["StrictDiagnostic", ...]]: ...
    def schema(self) -> dict[str, Any]: ...


class CxxMappingCapability(Protocol):
    def profile(self) -> Any: ...


class MeshDiagnosticPolicyCapability(Protocol):
    def is_nondimensional(self, spec: "TutorialSpec") -> bool: ...


class CaseCompatibilityCapability(Protocol):
    def has_case_marker(self, request: CaseCompatibilityRequest) -> bool: ...
    def is_runnable_without_workflow(self, request: CaseCompatibilityRequest) -> bool: ...


class SweepMaterializerCapability(Protocol):
    def route(self, request: SweepRoutingRequest, *, driver_context: Any) -> dict[str, Any]: ...
    def materialize(self, request: SweepMaterializationRequest) -> None: ...


class CommandAuthorizationCapability(Protocol):
    """What the active plugin authorizes a workflow step to invoke.

    ``solver_commands`` and ``auxiliary_commands`` are both authorized, but
    only ``solver_commands`` names binaries that produce a run's artifacts;
    core's artifact-producer heuristic must consult that one alone.
    """

    def solver_commands(self) -> frozenset[str]: ...
    def auxiliary_commands(self) -> frozenset[str]: ...
    def utility_manifests(self) -> dict[str, Any]: ...
    def utility_roots(self) -> tuple[Path, ...]: ...


class CaseIntrospectionCapability(Protocol):
    """Solver-specific case-model resolution and the fields it exposes.

    ``resolve_case_models`` is a best-effort, never-raising read of a case's
    on-disk configuration; ``samplable_fields`` names the fields the resolved
    model exposes for sampling by function objects, split by region. A
    plugin with no solver semantics (the generic plugin) resolves nothing and
    exposes no fields.
    """

    def resolve_case_models(self, case_root: Path) -> dict[str, Any]: ...
    def samplable_fields(self, resolved: dict[str, Any]) -> dict[str, tuple[str, ...]]: ...


class CaseFileContractCapability(Protocol):
    """Which case files the active plugin's profile declares, and how strictly.

    Sourced directly from ``PluginProfile.case_files``: ``required_files``
    lists every rule whose ``required`` is ``"always"``; ``conditional_files``
    lists the rest. ``required_rules`` returns the same required rules with
    their ``role`` intact.

    **Roles are namespaced and the prefix is load-bearing.** ``openfoam.*``
    marks a file the OpenFOAM runtime itself requires (``openfoam.control_dict``,
    ``openfoam.discretisation``, ...); ``plugin.*`` marks one the solver plugin
    requires (``plugin.configuration``). Consumers split on that prefix -- a
    rule written as ``control_dict`` rather than ``openfoam.control_dict`` will
    be silently classified as plugin-owned. The profile loader does not yet
    validate the namespace. ``get_profile()`` is a required v1 plugin member and
    ``case_files`` is already part of ``PluginProfile``, so every plugin
    already carries this data -- no compatibility fallback is needed.

    ``describe_config_resolution`` is different: it is a human-readable
    sentence, not derived from ``case_files`` data, so it *does* need a
    compatibility fallback (``legacy_describe_config_resolution``) for v1
    plugins and plugins that never authored one -- cardiac-shaped only for
    the built-in cardiac plugin, plugin-neutral for everyone else.
    """

    def required_files(self) -> tuple[str, ...]: ...
    def conditional_files(self) -> tuple[str, ...]: ...
    def required_rules(self) -> tuple["CaseFileRule", ...]: ...
    def describe_config_resolution(self) -> str: ...


class OverrideSchemaCapability(Protocol):
    """The plugin's authored configuration vocabulary.

    ``config_schema`` is the machine-readable description of the ``--config``
    JSON an agent writes, including a worked example for the named tutorial.
    ``dict_entry_catalog`` returns the plugin's dictionary entries arranged by
    its own document names, **unserialized** -- core owns serialization, the
    plugin owns the vocabulary and the document shape.
    """

    def config_schema(
        self, tutorial_name: str, make_spec_info: dict[str, Any],
    ) -> dict[str, Any]: ...
    def dict_entry_catalog(self) -> dict[str, Any]: ...


class RuntimeEvidenceCapability(Protocol):
    """Where the plugin's runtime evidence lives.

    Phase 4 (telemetry) consumes ``solve_step_commands`` and
    ``telemetry_source_globs``; Phase 5 (observables) will consume
    ``artifact_value_reader``; both remain declaration-only for now. Phase 2
    (provenance) now consumes ``extra_provenance_paths`` for real.

    Every member degrades to empty for a plugin that declares nothing, which
    is the honest answer rather than a solver-shaped guess -- so this
    capability needs no compatibility fallback.
    """

    def solve_step_commands(self) -> frozenset[str]: ...
    def telemetry_source_globs(self, command: str) -> tuple[str, ...]: ...
    def extra_provenance_paths(self, case_root: Path) -> tuple[RuntimeDependency, ...]: ...
    def artifact_value_reader(self, artifact_format: str) -> Any | None: ...


class CaseProvenanceCapability(Protocol):
    """Solver-declared case classification for the provenance snapshot.

    ``required_inputs`` returns already-*resolved* paths, not patterns --
    field names are dictionary-configurable and field locations resolve by
    a backward ``Time::findInstance`` search with a ``constant/`` fallback,
    so a field's canonical path is not knowable from its name alone.
    ``generated_output_globs`` may stay globs: generated diagnostic outputs
    have fixed names.

    Both take the resolved case dictionaries (not just the model name) and
    the selected start time, because gating is by dictionary *value*: e.g.
    ``conductivitySource field`` vs ``uniform`` flips a mandatory read on
    and off, and an absent key silently defaults to ``uniform``.

    Routed through the capability adapter exactly like every other plugin
    capability -- deliberately **not** a mandatory ``SolverPluginV2``
    member, so existing v2 third-party plugins keep loading. The adapter's
    fallback returns empty for both, which under the resolution precedence
    (a DAG step's ``consumes``, then a plugin's ``required_inputs``, then
    ``generated_output_globs``, then: unknown files are ``required_input``)
    means "everything unknown is a required input" -- the safe default for
    a plugin that declares nothing.
    """

    def required_inputs(
        self,
        case_root: Path,
        resolved_case: dict[str, Any],
        selected_start_time: str,
    ) -> tuple[ResolvedInput, ...]: ...

    def generated_output_globs(
        self,
        case_root: Path,
        resolved_case: dict[str, Any],
        selected_start_time: str,
    ) -> tuple[str, ...]: ...


class ReportCatalogCapability(Protocol):
    """Post-run report definitions the active plugin wants offered.

    ``report_catalog`` (:mod:`openfoam_driver.report_catalog`) owns the
    solver-neutral machinery -- ``ReportDefinition``, the ``applicable_when``
    predicate evaluator, the JSON record shape -- but the *catalog itself*
    (which reports exist, e.g. "Vm field" or "activation map") is
    solver-specific data. Not a mandatory ``SolverPluginV2`` member, so
    existing v2 third-party plugins keep loading; the fallback
    (``legacy_report_catalog``) is cardiac-shaped only for the built-in
    cardiac plugin and empty for everyone else -- the honest answer for a
    plugin that declares no reports, matching the pattern already used by
    :class:`CaseProvenanceCapability`.
    """

    def reports(self) -> tuple["ReportDefinition", ...]: ...


class NamedCatalogsCapability(Protocol):
    """The plugin's own named catalogs, namespaced generically.

    ``catalogs`` returns a mapping from plugin-chosen catalog name to
    plugin-chosen catalog content (e.g. the cardiac plugin's
    ``ionic_model_catalog``/``active_tension_catalog``) -- core imposes no
    key set, it only namespaces the whole mapping under
    ``describe_entry``'s ``plugin_catalogs`` key and serializes it. Not a
    mandatory ``SolverPluginV2`` member, so existing v2 third-party plugins
    keep loading; the fallback (``legacy_named_catalogs``) is cardiac-shaped
    only for the built-in cardiac plugin and empty for everyone else,
    matching the pattern already used by :class:`ReportCatalogCapability`.
    """

    def catalogs(self) -> dict[str, Any]: ...


@dataclass(frozen=True)
class _TutorialCatalogAdapter:
    plugin: "SolverPlugin"

    def catalog(self) -> dict[str, Any]:
        return self.plugin.get_tutorial_catalog()

    def displays(self) -> tuple[Any, ...]:
        return self.plugin.get_tutorial_displays()


@dataclass(frozen=True)
class _DictionaryCatalogAdapter:
    plugin: "SolverPlugin"

    def entries(self) -> tuple[Any, ...]:
        return self.plugin.get_dict_entries()

    def catalog(self) -> Any:
        return self.plugin.get_dictionary_catalog()

    def groups(self) -> dict[str, tuple[Any, ...]]:
        return self.plugin.get_dict_groups()


@dataclass(frozen=True)
class _CapabilityManifestAdapter:
    plugin: "SolverPlugin"

    def manifest(self) -> Any:
        return self.plugin.get_capabilities()


@dataclass(frozen=True)
class _ConfigurationValidatorAdapter:
    plugin: "SolverPlugin"

    def validate(
        self, request: ConfigurationValidationRequest,
    ) -> tuple["StrictDiagnostic", ...]:
        return self.plugin.validate_configuration(request.spec)


@dataclass(frozen=True)
class _RunSemanticValidatorAdapter:
    plugin: "SolverPlugin"

    def validate(self, request: RunSemanticValidationRequest) -> tuple[Any, ...]:
        return self.plugin.validate_run_semantics(request.context)


@dataclass(frozen=True)
class _ArtifactPredictorAdapter:
    plugin: "SolverPlugin"

    def predict(self, request: ArtifactPredictionRequest) -> tuple["DataArtifact", ...]:
        return self.plugin.predict_data_artifacts(request.case_root, request.spec)


@dataclass(frozen=True)
class _RunDocumentConfigurationAdapter:
    plugin: "SolverPlugin"

    def build(
        self, request: RunDocumentConfigurationRequest,
    ) -> tuple[dict[str, dict[str, Any]], tuple["StrictDiagnostic", ...]]:
        hook = getattr(self.plugin, "build_run_document_config", None)
        if callable(hook):
            return hook(request.spec)
        # Existing plugins were interpreted through the cardiac-shaped v2
        # adapter.  Preserve that fallback until Plan 2 changes the document.
        from .compatibility import legacy_run_document_config

        return legacy_run_document_config(request.spec)

    def schema(self) -> dict[str, Any]:
        hook = getattr(self.plugin, "get_run_document_config_schema", None)
        if callable(hook):
            return hook()
        from .compatibility import legacy_run_document_config_schema

        return legacy_run_document_config_schema(self.plugin)


@dataclass(frozen=True)
class _CxxMappingAdapter:
    plugin: "SolverPlugin"

    def profile(self) -> Any:
        return self.plugin.get_profile()


@dataclass(frozen=True)
class _MeshDiagnosticPolicyAdapter:
    plugin: "SolverPlugin"

    def is_nondimensional(self, spec: "TutorialSpec") -> bool:
        hook = getattr(self.plugin, "is_nondimensional_case", None)
        if callable(hook):
            return bool(hook(spec))
        from .compatibility import legacy_nondimensional_case

        return legacy_nondimensional_case(spec)


@dataclass(frozen=True)
class _CaseCompatibilityAdapter:
    plugin: "SolverPlugin"

    def has_case_marker(self, request: CaseCompatibilityRequest) -> bool:
        hook = getattr(self.plugin, "has_case_marker", None)
        if callable(hook):
            return bool(hook(request.case_root))
        from .compatibility import legacy_case_marker

        return legacy_case_marker(request.case_root)

    def is_runnable_without_workflow(self, request: CaseCompatibilityRequest) -> bool:
        hook = getattr(self.plugin, "is_case_runnable_without_workflow", None)
        if callable(hook):
            return bool(hook(request.case_root))
        from .compatibility import legacy_case_runnable_without_workflow

        return legacy_case_runnable_without_workflow(request.case_root)


@dataclass(frozen=True)
class _SweepMaterializerAdapter:
    plugin: "SolverPlugin"

    def route(self, request: SweepRoutingRequest, *, driver_context: Any) -> dict[str, Any]:
        hook = getattr(self.plugin, "route_sweep_case_values", None)
        if callable(hook):
            return hook(
                base=request.base,
                resolved_axis_values=request.resolved_axis_values,
                driver_context=driver_context,
            )
        # Compatibility bridge for existing third-party-style plugins.  Plan 1
        # preserves the historical cardiac-shaped generic sweep fallback.
        from .compatibility import legacy_route_sweep_case

        return legacy_route_sweep_case(
            base=request.base,
            resolved_axis_values=request.resolved_axis_values,
            driver_context=driver_context,
        )

    def materialize(self, request: SweepMaterializationRequest) -> None:
        hook = getattr(self.plugin, "materialize_sweep_case", None)
        if callable(hook):
            hook(case_dir=request.case_dir, routed=request.routed)
            return
        from .compatibility import legacy_materialize_sweep_case

        legacy_materialize_sweep_case(case_dir=request.case_dir, routed=request.routed)


@dataclass(frozen=True)
class _CommandAuthorizationAdapter:
    plugin: "SolverPlugin"

    def solver_commands(self) -> frozenset[str]:
        hook = getattr(self.plugin, "get_solver_commands", None)
        if callable(hook):
            return frozenset(hook())
        from .compatibility import legacy_solver_commands

        return legacy_solver_commands(self.plugin)

    def auxiliary_commands(self) -> frozenset[str]:
        hook = getattr(self.plugin, "get_auxiliary_commands", None)
        if callable(hook):
            return frozenset(hook())
        from .compatibility import legacy_auxiliary_commands

        return legacy_auxiliary_commands(self.plugin)

    def utility_manifests(self) -> dict[str, Any]:
        hook = getattr(self.plugin, "get_utility_manifests", None)
        if callable(hook):
            return dict(hook())
        from .compatibility import legacy_utility_manifests

        return legacy_utility_manifests(self.plugin)

    def utility_roots(self) -> tuple[Path, ...]:
        hook = getattr(self.plugin, "get_utility_roots", None)
        if callable(hook):
            return tuple(hook())
        from .compatibility import legacy_utility_roots

        return legacy_utility_roots(self.plugin)


@dataclass(frozen=True)
class _CaseIntrospectionAdapter:
    plugin: "SolverPlugin"

    def resolve_case_models(self, case_root: Path) -> dict[str, Any]:
        hook = getattr(self.plugin, "resolve_case_models", None)
        if callable(hook):
            return dict(hook(case_root))
        from .compatibility import legacy_resolve_case_models

        return legacy_resolve_case_models(self.plugin, case_root)

    def samplable_fields(self, resolved: dict[str, Any]) -> dict[str, tuple[str, ...]]:
        hook = getattr(self.plugin, "get_samplable_fields", None)
        if callable(hook):
            return {k: tuple(v) for k, v in hook(resolved).items()}
        from .compatibility import legacy_samplable_fields

        return legacy_samplable_fields(self.plugin, resolved)


@dataclass(frozen=True)
class _CaseFileContractAdapter:
    plugin: "SolverPlugin"

    def _rules(self) -> tuple["CaseFileRule", ...]:
        return tuple(self.plugin.get_profile().case_files)

    def required_rules(self) -> tuple["CaseFileRule", ...]:
        """Required rules with their ``role`` intact, so a consumer need not
        re-derive plugin semantics from a path prefix."""
        return tuple(rule for rule in self._rules() if rule.required == "always")

    def required_files(self) -> tuple[str, ...]:
        return tuple(rule.path for rule in self.required_rules())

    def conditional_files(self) -> tuple[str, ...]:
        return tuple(rule.path for rule in self._rules() if rule.required != "always")

    def describe_config_resolution(self) -> str:
        hook = getattr(self.plugin, "get_config_resolution_description", None)
        if callable(hook):
            return str(hook())
        from .compatibility import legacy_describe_config_resolution

        return legacy_describe_config_resolution(self.plugin)


@dataclass(frozen=True)
class _OverrideSchemaAdapter:
    plugin: "SolverPlugin"

    def config_schema(
        self, tutorial_name: str, make_spec_info: dict[str, Any],
    ) -> dict[str, Any]:
        hook = getattr(self.plugin, "get_override_schema", None)
        if callable(hook):
            return dict(hook(tutorial_name, make_spec_info))
        from .compatibility import legacy_override_schema

        return legacy_override_schema(self.plugin, tutorial_name, make_spec_info)

    def dict_entry_catalog(self) -> dict[str, Any]:
        hook = getattr(self.plugin, "get_dict_entry_catalog", None)
        if callable(hook):
            return dict(hook())
        from .compatibility import legacy_dict_entry_catalog

        return legacy_dict_entry_catalog(self.plugin)


@dataclass(frozen=True)
class _RuntimeEvidenceAdapter:
    plugin: "SolverPlugin"

    def solve_step_commands(self) -> frozenset[str]:
        hook = getattr(self.plugin, "get_solve_step_commands", None)
        return frozenset(hook()) if callable(hook) else frozenset()

    def telemetry_source_globs(self, command: str) -> tuple[str, ...]:
        hook = getattr(self.plugin, "get_telemetry_source_globs", None)
        return tuple(hook(command)) if callable(hook) else ()

    def extra_provenance_paths(self, case_root: Path) -> tuple[Path, ...]:
        hook = getattr(self.plugin, "get_extra_provenance_paths", None)
        return tuple(hook(case_root)) if callable(hook) else ()

    def artifact_value_reader(self, artifact_format: str):
        hook = getattr(self.plugin, "get_artifact_value_reader", None)
        return hook(artifact_format) if callable(hook) else None


@dataclass(frozen=True)
class _CaseProvenanceAdapter:
    plugin: "SolverPlugin"

    def required_inputs(
        self,
        case_root: Path,
        resolved_case: dict[str, Any],
        selected_start_time: str,
    ) -> tuple[ResolvedInput, ...]:
        hook = getattr(self.plugin, "get_required_inputs", None)
        if callable(hook):
            return tuple(hook(case_root, resolved_case, selected_start_time))
        return ()

    def generated_output_globs(
        self,
        case_root: Path,
        resolved_case: dict[str, Any],
        selected_start_time: str,
    ) -> tuple[str, ...]:
        hook = getattr(self.plugin, "get_generated_output_globs", None)
        if callable(hook):
            return tuple(hook(case_root, resolved_case, selected_start_time))
        return ()


@dataclass(frozen=True)
class _ReportCatalogAdapter:
    plugin: "SolverPlugin"

    def reports(self) -> tuple["ReportDefinition", ...]:
        hook = getattr(self.plugin, "get_report_catalog", None)
        if callable(hook):
            return tuple(hook())
        from .compatibility import legacy_report_catalog

        return legacy_report_catalog(self.plugin)


@dataclass(frozen=True)
class _NamedCatalogsAdapter:
    plugin: "SolverPlugin"

    def catalogs(self) -> dict[str, Any]:
        hook = getattr(self.plugin, "get_named_catalogs", None)
        if callable(hook):
            return dict(hook())
        from .compatibility import legacy_named_catalogs

        return legacy_named_catalogs(self.plugin)


@dataclass(frozen=True)
class PluginCapabilities:
    """Focused internal view over the unchanged public plugin object."""

    tutorials: TutorialCatalogCapability
    dictionaries: DictionaryCatalogCapability
    manifest: CapabilityManifestCapability
    configuration_validator: ConfigurationValidatorCapability
    run_semantic_validator: RunSemanticValidatorCapability
    artifacts: ArtifactPredictorCapability
    run_document_configuration: RunDocumentConfigurationCapability
    cxx_mapping: CxxMappingCapability
    mesh_diagnostic_policy: MeshDiagnosticPolicyCapability
    case_compatibility: CaseCompatibilityCapability
    sweep_materializer: SweepMaterializerCapability
    command_authorization: CommandAuthorizationCapability
    case_introspection: CaseIntrospectionCapability
    case_files: CaseFileContractCapability
    override_schema: OverrideSchemaCapability
    runtime_evidence: RuntimeEvidenceCapability
    case_provenance: CaseProvenanceCapability
    report_catalog: ReportCatalogCapability
    named_catalogs: NamedCatalogsCapability


def adapt_plugin_capabilities(plugin: "SolverPlugin") -> PluginCapabilities:
    """Build a behavior-preserving capability bundle for ``plugin``."""

    return PluginCapabilities(
        tutorials=_TutorialCatalogAdapter(plugin),
        dictionaries=_DictionaryCatalogAdapter(plugin),
        manifest=_CapabilityManifestAdapter(plugin),
        configuration_validator=_ConfigurationValidatorAdapter(plugin),
        run_semantic_validator=_RunSemanticValidatorAdapter(plugin),
        artifacts=_ArtifactPredictorAdapter(plugin),
        run_document_configuration=_RunDocumentConfigurationAdapter(plugin),
        cxx_mapping=_CxxMappingAdapter(plugin),
        mesh_diagnostic_policy=_MeshDiagnosticPolicyAdapter(plugin),
        case_compatibility=_CaseCompatibilityAdapter(plugin),
        sweep_materializer=_SweepMaterializerAdapter(plugin),
        command_authorization=_CommandAuthorizationAdapter(plugin),
        case_introspection=_CaseIntrospectionAdapter(plugin),
        case_files=_CaseFileContractAdapter(plugin),
        override_schema=_OverrideSchemaAdapter(plugin),
        runtime_evidence=_RuntimeEvidenceAdapter(plugin),
        case_provenance=_CaseProvenanceAdapter(plugin),
        report_catalog=_ReportCatalogAdapter(plugin),
        named_catalogs=_NamedCatalogsAdapter(plugin),
    )
