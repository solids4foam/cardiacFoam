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
    from .runtime.models import DataArtifact, TutorialSpec
    from ..planning_types import StrictDiagnostic


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
    )
