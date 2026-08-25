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
    from openfoam_driver.core.planning_types import StrictDiagnostic
    from openfoam_driver.core.report_catalog import ReportDefinition
    from openfoam_driver.core.specs.apply_overrides import OverrideScope, RegenerationScope


@dataclass(frozen=True)
class ConfigurationValidationRequest:
    """Input to :class:`ConfigurationValidatorCapability`: the resolved spec
    whose configuration the plugin should judge."""

    spec: "TutorialSpec"


@dataclass(frozen=True)
class RunSemanticValidationRequest:
    """Input to :class:`RunSemanticValidatorCapability`: the loose run-context
    mapping assembled at execution time, not a resolved ``TutorialSpec``."""

    context: dict[str, Any]


@dataclass(frozen=True)
class ArtifactPredictionRequest:
    """Input to :class:`ArtifactPredictorCapability`: the case to inspect and
    the spec it was built from. May name a case that does not exist yet."""

    case_root: Path
    spec: "TutorialSpec"


@dataclass(frozen=True)
class RunDocumentConfigurationRequest:
    """Input to :class:`RunDocumentConfigurationCapability`: the spec whose
    plugin-owned RunDocument ``config`` object is to be built."""

    spec: "TutorialSpec"


@dataclass(frozen=True)
class CaseCompatibilityRequest:
    """Input to :class:`CaseCompatibilityCapability`: the case folder to judge
    by filesystem evidence, before any dictionary is parsed."""

    case_root: Path


@dataclass(frozen=True)
class SweepRoutingRequest:
    """Input to :class:`SweepMaterializerCapability` ``route``: the sweep's
    static base values plus one expanded axis combination from core."""

    base: dict[str, Any]
    resolved_axis_values: dict[str, Any]


@dataclass(frozen=True)
class SweepMaterializationRequest:
    """Input to :class:`SweepMaterializerCapability` ``materialize``: where to
    write, and the routed values ``route`` produced for this case."""

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
    """The tutorials this plugin registers, and how to display them.

    ``catalog`` returns the plugin's registry keyed by tutorial name -- the
    entry names ``driverFoam`` accepts. ``displays`` returns the presentation
    metadata ``describe`` renders. Both are required v1 members, so there is
    no fallback: a plugin that registers no tutorials returns empty rather
    than omitting the member.

    :adapts: get_tutorial_catalog, get_tutorial_displays
    :consumed-by: openfoam_driver/core/runtime/registry.py, openfoam_driver/plugins/cardiacfoam/dict_builder.py
    :fallback: none
    :status: mandatory
    """

    def catalog(self) -> dict[str, Any]: ...
    def displays(self) -> tuple[Any, ...]: ...


class DictionaryCatalogCapability(Protocol):
    """The plugin's dictionary vocabulary, in three shapes for three callers.

    ``entries`` is the flat tuple of ``DictEntry`` values; ``catalog`` is the
    same data as a queryable ``DictionaryCatalog`` (``entries_for(document)``);
    ``groups`` buckets entries by the plugin's own document names. Core does
    not know those names -- ``electroProperties`` is cardiac vocabulary, and a
    solids4foam plugin would say ``solidProperties`` instead.

    All three are required v1 members with no fallback. This is the seam that
    keeps dictionary *syntax* knowledge (core's) apart from dictionary
    *meaning* (the plugin's).

    :adapts: get_dict_entries, get_dict_groups, get_dictionary_catalog
    :consumed-by: openfoam_driver/dict_entries.py, openfoam_driver/plugins/cardiacfoam/sweep.py, openfoam_driver/core/specs/apply_overrides.py, openfoam_driver/core/specs/dict_builder.py, openfoam_driver/core/specs/validation.py, openfoam_driver/core/strict_planning.py
    :fallback: none
    :status: mandatory
    """

    def entries(self) -> tuple[Any, ...]: ...
    def catalog(self) -> Any: ...
    def groups(self) -> dict[str, tuple[Any, ...]]: ...


class CapabilityManifestCapability(Protocol):
    """The plugin's self-description of what it can model.

    Deliberately untyped at the core boundary: ``CapabilityManifest`` in
    ``plugin_interface`` is an empty Protocol, because the axes a plugin
    advertises are its own (cardiacFoam declares ionic models and solvers; a
    different plugin would declare something else entirely). Core namespaces
    and serialises it for ``describe`` without interpreting it.

    :adapts: get_capabilities
    :consumed-by: openfoam_driver/dict_entries.py, openfoam_driver/core/introspection.py, openfoam_driver/core/strict_planning.py
    :fallback: none
    :status: mandatory
    """

    def manifest(self) -> Any: ...


class ConfigurationValidatorCapability(Protocol):
    """Plan-time validation of a resolved tutorial spec, in the plugin's terms.

    Returns ``StrictDiagnostic`` values rather than raising, so the strict
    planner can report every problem in one pass instead of stopping at the
    first -- the property that lets an agent self-heal a case in one edit
    round. An empty tuple means "nothing this plugin can object to", never
    "not checked".

    :adapts: validate_configuration
    :consumed-by: openfoam_driver/core/strict_planning.py
    :fallback: none
    :status: mandatory
    """

    def validate(
        self, request: ConfigurationValidationRequest,
    ) -> tuple["StrictDiagnostic", ...]: ...


class RunSemanticValidatorCapability(Protocol):
    """Validation of a run's semantics, as opposed to its configuration.

    Distinct from :class:`ConfigurationValidatorCapability`: that one judges a
    ``TutorialSpec``, this one judges a looser run context dictionary at
    execution time. Required v1 member, no fallback.

    :adapts: validate_run_semantics
    :consumed-by: openfoam_driver/core/specs/validation.py
    :fallback: none
    :status: mandatory
    """

    def validate(self, request: RunSemanticValidationRequest) -> tuple[Any, ...]: ...


class ArtifactPredictorCapability(Protocol):
    """What files this case will produce, predicted before it runs.

    The prediction is a contract: the strict runner compares it against what
    actually appeared and fails with ``missing_expected_artifacts`` on a
    mismatch, which is how a silently-not-writing solver gets caught.

    That makes the predictor's honesty load-bearing. It must never raise --
    agents call it against partly-mutated cases -- and it must distinguish
    "I could not determine the exports" from "the exports are declared
    empty". Conflating those two through a falsy empty tuple is a real defect
    this code has already had.

    :adapts: predict_data_artifacts
    :consumed-by: openfoam_driver/core/runtime/artifacts.py
    :fallback: none
    :status: mandatory
    """

    def predict(self, request: ArtifactPredictionRequest) -> tuple["DataArtifact", ...]: ...


class RunDocumentConfigurationCapability(Protocol):
    """The plugin's half of a RunDocument: its ``config`` object and schema.

    ``schemas/run-document.json`` declares ``config`` as an open object
    (``additionalProperties: true``) with no fixed key set, so the whole
    vocabulary inside it belongs to the plugin. ``build`` produces the object
    and any diagnostics; ``schema`` produces the JSON Schema core validates it
    against dynamically, turning a plugin's own rules into structured
    diagnostics an agent can act on.

    The fallback returns an empty config for a non-cardiac plugin. It used to
    return the cardiac phase vocabulary
    (``anatomy``/``physics``/``stimulus``/``solver``) to every plugin --
    exactly the fixed phases RunDocument v3 removed from core.

    :adapts: build_run_document_config, get_run_document_config_schema
    :consumed-by: openfoam_driver/core/runtime/run_document_adapter.py, openfoam_driver/core/runtime/run_document_exec.py
    :fallback: legacy_run_document_config, legacy_run_document_config_schema
    :status: optional
    """

    def build(
        self, request: RunDocumentConfigurationRequest,
    ) -> tuple[dict[str, dict[str, Any]], tuple["StrictDiagnostic", ...]]: ...
    def schema(self) -> dict[str, Any]: ...


class CxxMappingCapability(Protocol):
    """The plugin's declarative profile: case-file rules and C++ provenance.

    Sourced from the plugin's ``plugin.yaml`` via ``get_profile()``. Named for
    the C++ source mapping it carries (which solver sources back which
    dictionary keys, used for provenance fingerprinting), but the same profile
    also backs :class:`CaseFileContractCapability`.

    :adapts: get_profile
    :consumed-by: openfoam_driver/core/strict_planning.py
    :fallback: none
    :status: mandatory
    """

    def profile(self) -> Any: ...


class MeshDiagnosticPolicyCapability(Protocol):
    """Plugin-owned exemptions from, and additions to, core's mesh diagnostics.

    Core classifies the physical scale of every polyMesh region and warns when
    it looks wrong. Two plugin-specific escapes exist: a case may be
    deliberately non-dimensional (a manufactured-solution verifier, a
    single-cell model) and should not be judged against SI expectations at
    all; and a plugin may own point sets that are not polyMesh regions and
    that core therefore cannot check (cardiacFoam's ``constant/purkinjeGraph*``
    conduction trees).

    ``is_nondimensional`` falls back to ``False`` for a non-cardiac plugin,
    which keeps the diagnostics on -- the conservative direction, since the
    failure mode of a wrong exemption is silence.

    :adapts: get_mesh_geometry_diagnostics, is_nondimensional_case
    :consumed-by: openfoam_driver/core/strict_planning.py
    :fallback: legacy_nondimensional_case
    :status: optional
    """

    def is_nondimensional(self, spec: "TutorialSpec") -> bool: ...
    def extra_geometry_diagnostics(self, case_root: Path) -> tuple[Any, ...]: ...


class CaseCompatibilityCapability(Protocol):
    """Whether a case folder on disk belongs to this plugin, and whether it
    can run without driver-owned workflow metadata.

    Both questions are answered from filesystem evidence alone, before any
    dictionary is parsed, so both are necessarily plugin-specific: cardiacFoam
    recognises its own cases by ``constant/electroProperties*``, which means
    nothing to any other solver.

    Consulted only when core cannot answer from plugin-neutral evidence first
    (an executable ``Allrun``). The fallback returns ``False`` for a
    non-cardiac plugin. Before it was gated it returned the *cardiac* answer
    whichever plugin was loaded, so a case carrying an ``electroProperties``
    file was claimed by a plugin that had never heard of it.

    :adapts: has_case_marker, is_case_runnable_without_workflow
    :consumed-by: openfoam_driver/core/runtime/registry.py
    :fallback: legacy_case_marker, legacy_case_runnable_without_workflow
    :status: optional
    """

    def has_case_marker(self, request: CaseCompatibilityRequest) -> bool: ...
    def is_runnable_without_workflow(self, request: CaseCompatibilityRequest) -> bool: ...


class SweepMaterializerCapability(Protocol):
    """How one resolved sweep-axis combination becomes a runnable case.

    Core owns sweep *expansion* -- ``sweep_expansion.py`` computes the
    cross-product or zip of axes, validates lengths, and caps case counts
    without knowing what any axis means. This capability owns what a resolved
    combination *is*: which of the plugin's dictionaries and keys each axis
    lands in, and how the case is written.

    The split between the two methods is pure/impure, not two kinds of sweep.
    ``route`` is a total function from axis values to a routed mapping and
    must not touch the filesystem, which is what makes ``sweep-plan``
    non-destructive; ``materialize`` does every write.

    Uniquely among these capabilities, the fallback cannot be neutral. An
    empty routing would silently yield a case that is not the one the sweep
    asked for, so a plugin without these hooks is refused by name instead.
    Before it was gated, the cardiac materializer ran for any plugin --
    writing an ``Allrun`` invoking the ``cardiacFoam`` binary.

    :adapts: materialize_sweep_case, route_sweep_case_values
    :consumed-by: openfoam_driver/sweep_materialize.py, openfoam_driver/sweep_routing.py
    :fallback: legacy_materialize_sweep_case, legacy_route_sweep_case
    :status: optional
    """

    def route(self, request: SweepRoutingRequest, *, driver_context: Any) -> dict[str, Any]: ...
    def materialize(self, request: SweepMaterializationRequest) -> None: ...


class CommandAuthorizationCapability(Protocol):
    """What the active plugin authorizes a workflow step to invoke.

    ``solver_commands`` and ``auxiliary_commands`` are both authorized, but
    only ``solver_commands`` names binaries that produce a run's artifacts;
    core's artifact-producer heuristic must consult that one alone.

    :adapts: get_auxiliary_commands, get_solver_commands, get_utility_manifests, get_utility_roots
    :consumed-by: openfoam_driver/core/runtime/artifacts.py, openfoam_driver/core/runtime/workflow.py, openfoam_driver/core/strict_planning.py
    :fallback: legacy_auxiliary_commands, legacy_solver_commands, legacy_utility_manifests, legacy_utility_roots
    :status: optional
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

    :adapts: get_samplable_fields, resolve_case_models
    :consumed-by: openfoam_driver/core/capability_manifest.py, openfoam_driver/core/runtime/provenance_inputs.py
    :fallback: legacy_resolve_case_models, legacy_samplable_fields
    :status: optional
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

    :adapts: get_profile, get_config_resolution_description
    :consumed-by: openfoam_driver/core/runtime/strict_audit.py, openfoam_driver/core/tutorial_contracts.py
    :fallback: legacy_describe_config_resolution
    :status: mixed
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

    :adapts: get_dict_entry_catalog, get_override_schema
    :consumed-by: openfoam_driver/core/introspection.py
    :fallback: legacy_dict_entry_catalog, legacy_override_schema
    :status: optional
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

    :adapts: get_artifact_value_reader, get_extra_provenance_paths, get_solve_step_commands, get_telemetry_source_globs
    :consumed-by: openfoam_driver/core/runtime/provenance_inputs.py
    :fallback: none
    :status: optional
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
    capability -- deliberately **not** a mandatory ``SolverPlugin``
    member, so existing v2 third-party plugins keep loading. The adapter's
    fallback returns empty for both, which under the resolution precedence
    (a DAG step's ``consumes``, then a plugin's ``required_inputs``, then
    ``generated_output_globs``, then: unknown files are ``required_input``)
    means "everything unknown is a required input" -- the safe default for
    a plugin that declares nothing.

    :adapts: get_generated_output_globs, get_required_inputs
    :consumed-by: openfoam_driver/core/runtime/provenance_inputs.py
    :fallback: none
    :status: optional
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

    ``report_catalog`` (:mod:`openfoam_driver.core.report_catalog`) owns the
    solver-neutral machinery -- ``ReportDefinition``, the ``applicable_when``
    predicate evaluator, the JSON record shape -- but the *catalog itself*
    (which reports exist, e.g. "Vm field" or "activation map") is
    solver-specific data. Not a mandatory ``SolverPlugin`` member, so
    existing v2 third-party plugins keep loading; the fallback
    (``legacy_report_catalog``) is cardiac-shaped only for the built-in
    cardiac plugin and empty for everyone else -- the honest answer for a
    plugin that declares no reports, matching the pattern already used by
    :class:`CaseProvenanceCapability`.

    :adapts: get_report_catalog
    :consumed-by: scripts/export-report-catalog.py
    :fallback: legacy_report_catalog
    :status: optional
    """

    def reports(self) -> tuple["ReportDefinition", ...]: ...


class NamedCatalogsCapability(Protocol):
    """The plugin's own named catalogs, namespaced generically.

    ``catalogs`` returns a mapping from plugin-chosen catalog name to
    plugin-chosen catalog content (e.g. the cardiac plugin's
    ``ionic_model_catalog``/``active_tension_catalog``) -- core imposes no
    key set, it only namespaces the whole mapping under
    ``describe_entry``'s ``plugin_catalogs`` key and serializes it. Not a
    mandatory ``SolverPlugin`` member, so existing v2 third-party plugins
    keep loading; the fallback (``legacy_named_catalogs``) is cardiac-shaped
    only for the built-in cardiac plugin and empty for everyone else,
    matching the pattern already used by :class:`ReportCatalogCapability`.

    :adapts: get_named_catalogs
    :consumed-by: openfoam_driver/core/introspection.py
    :fallback: legacy_named_catalogs
    :status: optional
    """

    def catalogs(self) -> dict[str, Any]: ...


class OverrideScopeCapability(Protocol):
    """Plugin-declared ``$TOKEN.`` override scopes for the agent-facing
    ``step --strict --apply`` path (:mod:`openfoam_driver.core.specs.apply_overrides`).

    Generalizes what was previously a single hardcoded cardiac scope
    (``$ELECTRO_MODEL_COEFFS`` -> ``constant/electroProperties``): core no
    longer assumes there is exactly one scope, or that it lives at that one
    path. Not a mandatory ``SolverPlugin`` member, so existing v2
    third-party plugins keep loading; the fallback (``legacy_override_scopes``)
    declares the cardiac plugin's one scope and an empty tuple for everyone
    else, matching the pattern already used by
    :class:`ReportCatalogCapability`/:class:`NamedCatalogsCapability`.

    :adapts: get_override_scopes
    :consumed-by: openfoam_driver/core/specs/apply_overrides.py
    :fallback: legacy_override_scopes
    :status: optional
    """

    def scopes(self) -> tuple["OverrideScope", ...]: ...


class DictRegenerationCapability(Protocol):
    """Plugin-declared bare "selector" overrides that must REGENERATE a
    dict file rather than key-patch it, for the agent-facing
    ``step --strict --apply`` path (:mod:`openfoam_driver.core.specs.apply_overrides`).

    A sibling of :class:`OverrideScopeCapability`: that one covers
    ``$TOKEN.``-scoped leaves that patch in place; this one covers bare
    selectors (e.g. cardiacFoam's ``myocardiumSolver``) whose value change
    restructures the file -- renames a sub-block, changes which sibling
    keys are legal -- so a single key/value/scope patch cannot express it.
    Not a mandatory ``SolverPlugin`` member, so existing v2 third-party
    plugins keep loading; the fallback (``legacy_dict_regeneration_scopes``)
    declares the cardiac plugin's one scope and an empty tuple for everyone
    else, matching :class:`OverrideScopeCapability`.

    :adapts: get_regeneration_scopes
    :consumed-by: openfoam_driver/core/specs/apply_overrides.py
    :fallback: legacy_dict_regeneration_scopes
    :status: optional
    """

    def scopes(self) -> tuple["RegenerationScope", ...]: ...


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

        return legacy_run_document_config(self.plugin, request.spec)

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

        return legacy_nondimensional_case(self.plugin, spec)

    def extra_geometry_diagnostics(self, case_root: Path) -> tuple[Any, ...]:
        """Plugin-owned plan-time geometry checks core cannot express.

        Core classifies the scale of every polyMesh region; a plugin may own
        further point sets in the case that are not mesh regions (cardiacFoam's
        ``constant/purkinjeGraph*`` conduction trees, for example). A plugin
        that declares no such check contributes nothing -- there is no legacy
        fallback here, because "no extra checks" is the correct answer for a
        plugin that never had any.
        """
        hook = getattr(self.plugin, "get_mesh_geometry_diagnostics", None)
        if callable(hook):
            return tuple(hook(case_root))
        return ()


@dataclass(frozen=True)
class _CaseCompatibilityAdapter:
    plugin: "SolverPlugin"

    def has_case_marker(self, request: CaseCompatibilityRequest) -> bool:
        hook = getattr(self.plugin, "has_case_marker", None)
        if callable(hook):
            return bool(hook(request.case_root))
        from .compatibility import legacy_case_marker

        return legacy_case_marker(self.plugin, request.case_root)

    def is_runnable_without_workflow(self, request: CaseCompatibilityRequest) -> bool:
        hook = getattr(self.plugin, "is_case_runnable_without_workflow", None)
        if callable(hook):
            return bool(hook(request.case_root))
        from .compatibility import legacy_case_runnable_without_workflow

        return legacy_case_runnable_without_workflow(self.plugin, request.case_root)


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
            self.plugin,
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

        legacy_materialize_sweep_case(
            self.plugin, case_dir=request.case_dir, routed=request.routed
        )


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
class _OverrideScopeAdapter:
    plugin: "SolverPlugin"

    def scopes(self) -> tuple["OverrideScope", ...]:
        hook = getattr(self.plugin, "get_override_scopes", None)
        if callable(hook):
            return tuple(hook())
        from .compatibility import legacy_override_scopes

        return legacy_override_scopes(self.plugin)


@dataclass(frozen=True)
class _DictRegenerationAdapter:
    plugin: "SolverPlugin"

    def scopes(self) -> tuple["RegenerationScope", ...]:
        hook = getattr(self.plugin, "get_regeneration_scopes", None)
        if callable(hook):
            return tuple(hook())
        from .compatibility import legacy_dict_regeneration_scopes

        return legacy_dict_regeneration_scopes(self.plugin)


@dataclass(frozen=True)
class PluginCapabilities:
    """Core's focused, internal view over one loaded plugin.

    **Direction matters.** This is not an authoring surface. A plugin author
    implements :class:`~openfoam_driver.core.plugin_interface.SolverPlugin`
    and optionally ``SolverPluginOptionalHooks``; this bundle is what *core*
    holds to consult that plugin, pointing the other way. Nothing here is
    implemented by a plugin.

    Its purpose is to stop core reaching through ``DriverContext.plugin``
    directly: each field is a narrow seam over one concern, so a core module
    can depend on the one capability it needs instead of the whole plugin.
    ``test_plugin_dependency_boundary.py`` enforces that -- no production
    module may use ``driver_context.plugin.``.

    **Reading a capability.** Every capability Protocol below carries prose
    explaining why the seam exists, then four structured fields:

    ``:adapts:``
        the plugin member(s) the adapter calls, or ``none``
    ``:consumed-by:``
        core modules that really touch this capability (subset, not exhaustive)
    ``:fallback:``
        the ``compatibility.py`` function used when the hook is absent
    ``:status:``
        ``mandatory`` (called unconditionally), ``optional`` (probed via
        ``getattr`` and degraded), or ``mixed``

    Those fields are the single source of the "Plugin capability seams" table
    in ``ARCHITECTURE.md``, rendered by
    ``scripts/export-capability-seams.py`` and kept honest by
    ``test_capability_seam_documentation.py`` -- which checks that every
    ``:adapts:`` names a real plugin member and every ``:fallback:`` a real
    compatibility function, so a stale reference fails rather than rots.

    **What a missing optional hook means.** The named fallback runs, returning
    cardiac data only for ``org.cardiacfoam`` and a neutral value for every
    other plugin. Two cannot be neutral: a plugin without the sweep hooks is
    refused by name rather than swept by another plugin's writer.
    """

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
    override_scopes: OverrideScopeCapability
    dict_regeneration: DictRegenerationCapability


def adapt_plugin_capabilities(plugin: "SolverPlugin") -> PluginCapabilities:
    """Wrap one loaded plugin in the capability bundle core consumes.

    Called once per :class:`~openfoam_driver.core.plugin_interface.DriverContext`
    and cheap: every adapter is a frozen dataclass holding the plugin, and no
    plugin member is called here. Optional hooks are probed lazily at each
    call site, so a plugin missing one is adapted successfully and degrades
    only when that capability is actually used.
    """

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
        override_scopes=_OverrideScopeAdapter(plugin),
        dict_regeneration=_DictRegenerationAdapter(plugin),
    )
