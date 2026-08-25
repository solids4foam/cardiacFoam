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
#     plugin_interface
#
# Description
#     Solver-agnostic plugin contract for driverFOAM.
#     Defines the boundary between the generic OpenFOAM execution engine and
#     any domain-specific solver plugin (e.g. cardiacFoam, shallowWaterFoam).
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

"""Solver-agnostic plugin contract for driverFOAM.

Two Protocol classes define what a solver plugin must implement:

- :class:`SolverPlugin` — the plugin contract; 27 required members,
  enforced in full by :func:`validate_plugin`.
- :class:`SolverPluginOptionalHooks` — 14 probe-based optional hooks that
  unlock additional capabilities (sweeps, mesh diagnostics, report catalogs,
  override scopes, …).  **Read this class** to discover all extension points
  before deciding your plugin is complete.

Use :func:`driver_context` or :func:`load_plugin_context` to create a
validated, immutable :class:`DriverContext` for each public operation.
Use :func:`generic_openfoam_context` for the built-in no-domain stub and
:func:`default_driver_context` only at compatibility boundaries.

To build a new plugin start from ``core/generic_plugin.py`` and follow
``.agents/skills/driverfoam-plugin-builder/SKILL.md``.
"""


import re
from dataclasses import dataclass
from functools import cached_property
from importlib import import_module
from typing import Any, Protocol, TYPE_CHECKING, runtime_checkable

if TYPE_CHECKING:
    from openfoam_driver.core.plugin_capabilities import PluginCapabilities, RuntimeDependency
    from openfoam_driver.core.contracts.dictionary import DictEntry
    from openfoam_driver.core.runtime.models import TutorialSpec, CaseConfig, DataArtifact
    from openfoam_driver.core.planning_types import StrictDiagnostic
    from openfoam_driver.core.tutorials_display import TutorialDisplay
    from openfoam_driver.core.plugin_capabilities import ResolvedInput
    from openfoam_driver.core.report_catalog import ReportDefinition
    from openfoam_driver.core.specs.apply_overrides import OverrideScope, RegenerationScope
    from pathlib import Path


class CapabilityManifest(Protocol):
    """Protocol for a solver's capability manifest."""
    # This can be expanded based on the solver's specific domain (e.g. models, physics)
    pass


@runtime_checkable
class SolverPlugin(Protocol):
    """
    The strict contract that any OpenFOAM solver must implement 
    to be orchestrated by driverFOAM. 
    
    This interface creates a clean boundary between the generic OpenFOAM execution 
    engine and the domain-specific solver logic (e.g., cardiacFoam, fireFoam).
    """
    
    @property
    def plugin_name(self) -> str:
        """Name of the solver plugin (e.g., 'cardiacFoam')."""
        ...

    @property
    def plugin_id(self) -> str:
        """Stable machine identifier, independent of the display name."""
        ...

    @property
    def plugin_version(self) -> str:
        """Version of the plugin semantics used to construct a plan."""
        ...

    @property
    def plugin_api_version(self) -> str:
        """Version of the driverFOAM plugin contract implemented by this plugin."""
        ...

    # -- Command authorization -----------------------------------------------
    def get_solver_commands(self) -> frozenset[str]:
        """Binaries that produce a run's artifacts. Core's artifact-producer
        heuristic consults this set alone, never the auxiliary one."""
        ...

    def get_auxiliary_commands(self) -> frozenset[str]:
        """Additionally authorized binaries that produce no artifacts of their
        own -- meshers, decomposers, reconstructors."""
        ...

    def get_utility_manifests(self) -> dict[str, Any]:
        """Per-utility declarations of what each pre/post-solve utility
        consumes and produces, so workflow steps can be checked before they
        run."""
        ...

    def get_utility_roots(self) -> tuple["Path", ...]:
        """Directories holding this plugin's utility sources, for provenance
        fingerprinting."""
        ...

    # -- Case introspection ---------------------------------------------------
    def resolve_case_models(self, case_root: "Path") -> dict[str, Any]:
        """Best-effort read of a case's on-disk model selections. Must never
        raise: agents call it against partly-written cases."""
        ...

    def get_samplable_fields(self, resolved: dict[str, Any]) -> dict[str, tuple[str, ...]]:
        """Fields the resolved model exposes for sampling by function objects,
        keyed by region."""
        ...

    # -- Configuration vocabulary --------------------------------------------
    def get_override_schema(
        self, tutorial_name: str, make_spec_info: dict[str, Any],
    ) -> dict[str, Any]:
        """Machine-readable description of the ``--config`` JSON an agent may
        write for this tutorial, including a worked example."""
        ...

    def get_run_document_config_schema(self) -> dict[str, Any]:
        """JSON Schema for this plugin's RunDocument ``config`` object. Core
        validates against it dynamically and reports structured diagnostics,
        which is what lets an agent repair its own document."""
        ...

    def get_dict_entry_catalog(self) -> dict[str, Any]:
        """The plugin's dictionary entries arranged by its own document names,
        unserialized -- core owns serialization, the plugin owns vocabulary."""
        ...

    # -- Runtime evidence ------------------------------------------------------
    def get_solve_step_commands(self) -> frozenset[str]:
        """Which commands count as the solve step, for telemetry attribution."""
        ...

    def get_telemetry_source_globs(self, command: str) -> tuple[str, ...]:
        """Where a given command writes the logs telemetry is parsed from."""
        ...

    def get_extra_provenance_paths(self, case_root: "Path") -> tuple["RuntimeDependency", ...]:
        """Run-time dependencies outside the case tree: the solver binary, a
        linked library, a case-local shared object.

        Returns ``RuntimeDependency`` rather than bare paths so that "required
        but not found" is expressible. A tuple of paths can only omit, and
        omission reads as "nothing to check" -- the gap that let a rebuilt
        solver replay a resumed run's numbers as fresh."""
        ...

    def get_artifact_value_reader(self, artifact_format: str) -> Any | None:
        """Reader for a plugin-specific artifact format, or ``None`` if this
        plugin cannot read that format."""
        ...

    def get_profile(self):
        """Return declarative case/C++ provenance metadata for this plugin."""
        ...
        
    def get_dict_entries(self) -> tuple[DictEntry, ...]:
        """
        Return the catalog of all solver-specific dictionary entries.
        Agents will use this to introspect capabilities deterministically.
        """
        ...

    def get_dictionary_catalog(self):
        """Return dictionary entries partitioned by plugin-owned document name."""
        ...

    def get_dict_groups(self) -> dict[str, tuple[DictEntry, ...]]:
        """
        Return the dictionary entries organized by logical group.
        """
        ...

    def get_capabilities(self) -> CapabilityManifest:
        """
        Return the capabilities of the solver (e.g., supported physics, 
        models, regions).
        """
        ...

    def get_tutorial_catalog(self) -> dict:
        """
        Return the tutorial specs provided by this solver.
        """
        ...

    def get_tutorial_displays(self) -> tuple[TutorialDisplay, ...]:
        """
        Return the UI display cards for the registered tutorials.
        """
        ...

    def validate_configuration(self, spec: TutorialSpec) -> tuple[StrictDiagnostic, ...]:
        """
        Solver-specific validation logic that goes beyond simple DictEntry constraints.
        Returns a tuple of diagnostics (errors/warnings).
        """
        ...

    def validate_run_semantics(self, context: dict[str, Any]) -> tuple[Any, ...]:
        """Return solver-specific validation errors for a flattened config."""
        ...

    def predict_data_artifacts(self, case_root: Path, spec: TutorialSpec) -> tuple[DataArtifact, ...]:
        """
        Predict the domain-specific artifacts (like ECGs or Purkinje VTK files) 
        that this solver expects to produce.
        """
        ...

@dataclass(frozen=True)
class PluginIdentity:
    """Stable description of the plugin semantics attached to an operation."""

    id: str
    version: str
    api_version: str
    source: str
    capability_digest: str

    def to_json(self) -> dict[str, str]:
        """Serialise this identity for provenance records and ``describe``."""
        return {
            "id": self.id,
            "version": self.version,
            "api_version": self.api_version,
            "source": self.source,
            "capability_digest": self.capability_digest,
        }


@dataclass(frozen=True)
class DriverContext:
    """Per-operation dependency bundle for solver-specific behaviour.

    A context is deliberately immutable and must be passed through planning,
    discovery, and execution.  It replaces the former process-global active
    plugin, which allowed one CLI invocation or test to change another one's
    solver semantics.
    """

    plugin: SolverPlugin
    identity: PluginIdentity

    @cached_property
    def capabilities(self) -> "PluginCapabilities":
        """Return the focused internal view without changing dataclass fields."""
        from .plugin_capabilities import adapt_plugin_capabilities

        return adapt_plugin_capabilities(self.plugin)


@runtime_checkable
class SolverPluginOptionalHooks(Protocol):
    """Optional hooks a plugin MAY implement. Documentation, not enforcement.

    Every member here is probed with ``getattr`` by an adapter in
    :mod:`openfoam_driver.core.plugin_capabilities`. None is listed in
    ``_REQUIRED_PLUGIN_MEMBERS`` or ``_REQUIRED_V2_MEMBERS``, so this class is
    inert at load time: ``validate_plugin`` never consults it, and declaring or
    omitting any of these changes no plugin's loading behaviour.

    **Why this class exists.** Until it did, these fourteen hooks appeared
    nowhere in the plugin contract. They were reachable only by reading the
    private ``_*Adapter`` bodies, so a plugin author reading this file could
    not discover that the extension points existed at all -- while *not*
    implementing one silently routed them into a compatibility fallback.

    **Not implementing a hook is a real choice, not a no-op.** When the hook
    is absent, the adapter falls back to
    :mod:`openfoam_driver.core.compatibility`, whose ``legacy_*`` functions
    return cardiac data for the built-in cardiac plugin and a neutral value
    for everyone else. Two of them cannot be neutral and refuse instead:
    a plugin that does not implement ``route_sweep_case_values`` and
    ``materialize_sweep_case`` cannot be swept, and will be told so by name.

    Hooks are grouped by the capability they back; see that capability's
    docstring in ``plugin_capabilities.py`` for the full contract.
    """

    # -- CaseCompatibilityCapability -----------------------------------------
    def has_case_marker(self, case_root: "Path") -> bool:
        """Whether this case folder belongs to this plugin, by filesystem
        evidence alone. Absent -> ``False`` for non-cardiac plugins."""
        ...

    def is_case_runnable_without_workflow(self, case_root: "Path") -> bool:
        """Whether a case with no driver-owned workflow metadata and no
        ``Allrun`` is still runnable.

        Absent -> ``False``; core then relies on an executable ``Allrun``.
        """
        ...

    # -- RunDocumentConfigurationCapability ----------------------------------
    def build_run_document_config(
        self, spec: "TutorialSpec",
    ) -> tuple[dict[str, dict[str, Any]], tuple["StrictDiagnostic", ...]]:
        """Build this plugin's RunDocument ``config`` object and any
        diagnostics. Core imposes no key set (``schemas/run-document.json``
        declares ``config`` open). Absent -> ``({}, ())``."""
        ...

    # -- MeshDiagnosticPolicyCapability --------------------------------------
    def is_nondimensional_case(self, spec: "TutorialSpec") -> bool:
        """Whether SI mesh-scale diagnostics should be skipped for this case.
        Absent -> ``False``, keeping the diagnostics on."""
        ...

    def get_mesh_geometry_diagnostics(self, case_root: "Path") -> tuple[Any, ...]:
        """Plan-time geometry checks over plugin-owned point sets that are not
        polyMesh regions. Absent -> ``()``; there is no fallback, because "no
        extra checks" is correct for a plugin that has none."""
        ...

    # -- SweepMaterializerCapability -----------------------------------------
    def route_sweep_case_values(
        self,
        *,
        base: dict[str, Any],
        resolved_axis_values: dict[str, Any],
        driver_context: Any,
    ) -> dict[str, Any]:
        """Map one resolved sweep-axis combination onto this plugin's own
        case vocabulary. Must be pure -- no writes; ``materialize_sweep_case``
        does those. Absent -> sweeps are refused by name."""
        ...

    def materialize_sweep_case(self, *, case_dir: "Path", routed: dict[str, Any]) -> None:
        """Write one routed sweep case to disk. Absent -> sweeps are refused
        by name rather than materialized by another plugin's writer."""
        ...

    # -- CaseFileContractCapability ------------------------------------------
    def get_config_resolution_description(self) -> str:
        """One human-readable sentence naming which files resolve into a valid
        RunDocument config. Absent -> a plugin-neutral sentence."""
        ...

    # -- CaseProvenanceCapability --------------------------------------------
    def get_required_inputs(
        self,
        case_root: "Path",
        resolved_case: dict[str, Any],
        selected_start_time: str,
    ) -> tuple["ResolvedInput", ...]:
        """Already-resolved input paths this case reads. Resolved, not globs:
        field names are dictionary-configurable and locations resolve by a
        backward ``Time::findInstance`` search. Absent -> ``()``, which under
        the resolution precedence means "every unknown file is a required
        input" -- safe, but coarse."""
        ...

    def get_generated_output_globs(
        self,
        case_root: "Path",
        resolved_case: dict[str, Any],
        selected_start_time: str,
    ) -> tuple[str, ...]:
        """Globs for files this case generates rather than consumes. Globs are
        fine here: generated diagnostics have fixed names. Absent -> ``()``."""
        ...

    # -- ReportCatalogCapability ---------------------------------------------
    def get_report_catalog(self) -> tuple["ReportDefinition", ...]:
        """Post-run reports this plugin offers. Core owns the machinery; the
        catalog itself is plugin data. Absent -> ``()``."""
        ...

    # -- NamedCatalogsCapability ---------------------------------------------
    def get_named_catalogs(self) -> dict[str, Any]:
        """Plugin-chosen catalogs, namespaced under ``plugin_catalogs`` in
        ``describe``. Core imposes no key set. Absent -> ``{}``."""
        ...

    # -- OverrideScopeCapability / DictRegenerationCapability ----------------
    def get_override_scopes(self) -> tuple["OverrideScope", ...]:
        """``$TOKEN.``-scoped override targets that patch a dict in place.
        Absent -> ``()``."""
        ...

    def get_regeneration_scopes(self) -> tuple["RegenerationScope", ...]:
        """Bare selector overrides whose value change REGENERATES a dict file
        rather than patching it -- renaming sub-blocks or changing which
        sibling keys are legal. Absent -> ``()``."""
        ...


# The single plugin contract version this core can drive. Anything else is
# refused before any plugin catalog code runs.
SUPPORTED_PLUGIN_API_VERSIONS: frozenset[str] = frozenset({"2"})


_REQUIRED_PLUGIN_MEMBERS = (
    "plugin_name",
    "plugin_id",
    "plugin_version",
    "plugin_api_version",
    "get_profile",
    "get_dict_entries",
    "get_dictionary_catalog",
    "get_dict_groups",
    "get_capabilities",
    "get_tutorial_catalog",
    "get_tutorial_displays",
    "validate_configuration",
    "validate_run_semantics",
    "predict_data_artifacts",
    "get_solver_commands",
    "get_auxiliary_commands",
    "get_utility_manifests",
    "get_utility_roots",
    "resolve_case_models",
    "get_samplable_fields",
    "get_override_schema",
    "get_run_document_config_schema",
    "get_dict_entry_catalog",
    "get_solve_step_commands",
    "get_telemetry_source_globs",
    "get_extra_provenance_paths",
    "get_artifact_value_reader",
)

_PLUGIN_ID_RE = re.compile(r"[a-z0-9](?:[a-z0-9.-]*[a-z0-9])?")


def validate_plugin(plugin: Any) -> SolverPlugin:
    """Reject malformed plugin objects before they enter a driver context.

    This is an interface guard, not a sandbox: an imported plugin is trusted
    in-process Python code.  Methods are checked structurally here; their
    returned values are validated by their owning core consumers.
    """

    missing = [name for name in _REQUIRED_PLUGIN_MEMBERS if not hasattr(plugin, name)]
    if missing:
        raise TypeError(
            "SolverPlugin is missing required members: " + ", ".join(sorted(missing))
        )
    for name in ("plugin_name", "plugin_id", "plugin_version", "plugin_api_version"):
        value = getattr(plugin, name)
        if not isinstance(value, str) or not value.strip():
            raise TypeError(f"SolverPlugin.{name} must be a non-empty string")
    # Refused here -- before the callable checks, before get_profile(), and
    # before get_dict_entries() -- so an unsupported plugin's catalog code
    # never executes.
    if plugin.plugin_api_version not in SUPPORTED_PLUGIN_API_VERSIONS:
        raise TypeError(
            f"SolverPlugin.plugin_api_version {plugin.plugin_api_version!r} is "
            "not supported; this driverFOAM core drives "
            f"{sorted(SUPPORTED_PLUGIN_API_VERSIONS)}"
        )
    if not _PLUGIN_ID_RE.fullmatch(plugin.plugin_id):
        raise TypeError(
            "SolverPlugin.plugin_id must use lowercase letters, digits, dots, "
            "or hyphens and cannot start or end with punctuation"
        )
    # A member's presence is not enough -- it must be callable. Collected as
    # one batch rather than raised on the first miss, so a partial
    # implementation is reported completely instead of one name at a time.
    non_callable = [
        name for name in _REQUIRED_PLUGIN_MEMBERS
        if name not in ("plugin_name", "plugin_id", "plugin_version", "plugin_api_version")
        and not callable(getattr(plugin, name, None))
    ]
    if non_callable:
        raise TypeError(
            "SolverPlugin does not implement the plugin contract; missing: "
            + ", ".join(sorted(non_callable))
        )
    return plugin


def driver_context(plugin: SolverPlugin, *, source: str) -> DriverContext:
    """Create a validated immutable context for one public operation."""

    checked = validate_plugin(plugin)
    profile = checked.get_profile()
    if profile.plugin_id != checked.plugin_id:
        raise TypeError("SolverPlugin profile id does not match plugin_id")
    if profile.api_version != checked.plugin_api_version:
        raise TypeError("SolverPlugin profile API version does not match plugin_api_version")
    from .contracts.dictionary import DictEntry

    entries = tuple(checked.get_dict_entries())
    invalid_entries = [
        entry for entry in entries
        if not isinstance(entry, DictEntry) or not entry.driver_path.strip()
    ]
    if invalid_entries:
        raise TypeError("SolverPlugin.get_dict_entries() must return DictEntry values with paths")
    paths = [entry.driver_path for entry in entries]
    duplicates = sorted({path for path in paths if paths.count(path) > 1})
    if duplicates:
        raise TypeError(
            "SolverPlugin dictionary catalog has duplicate paths: "
            + ", ".join(duplicates)
        )
    return DriverContext(
        plugin=checked,
        identity=PluginIdentity(
            id=checked.plugin_id,
            version=checked.plugin_version,
            api_version=checked.plugin_api_version,
            source=source,
            capability_digest=profile.digest,
        ),
    )


def load_plugin_context(target: str) -> DriverContext:
    """Load a plugin by discovered id, or by trusted ``module:Class`` import.

    A colon always means the trusted local-development import form, which the
    CLI labels unsafe. Without a colon the argument names an installed plugin
    from the ``driverfoam.plugins`` entry-point group. Neither form is
    sandboxed: loading a plugin executes its Python code.
    """

    if ":" not in target:
        from .plugin_discovery import load_discovered_plugin

        return load_discovered_plugin(target)

    try:
        module_path, class_name = target.split(":", maxsplit=1)
        if not module_path or not class_name:
            raise ValueError
    except ValueError as exc:
        raise ValueError("Plugin target must use the form 'module.path:ClassName'") from exc
    module = import_module(module_path)
    plugin_class = getattr(module, class_name)
    return driver_context(plugin_class(), source=f"trusted-import:{target}")


def default_driver_context() -> DriverContext:
    """Return a fresh compatibility context for the built-in cardiac plugin.

    This function exists at public compatibility boundaries only.  Core
    internals must receive a :class:`DriverContext` explicitly and must not
    retain it in module state.
    """

    from .compatibility import legacy_default_driver_context

    return legacy_default_driver_context()


def generic_openfoam_context() -> DriverContext:
    """Return the built-in context with no solver-specific semantics."""

    from .generic_plugin import GenericOpenFOAMPlugin

    return driver_context(GenericOpenFOAMPlugin(), source="built-in:generic-openfoam")
