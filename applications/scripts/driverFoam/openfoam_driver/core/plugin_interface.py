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
#     Defines the contract for expanding driverFOAM to other solvers.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

from __future__ import annotations

import re
from dataclasses import dataclass
from importlib import import_module
from typing import Any, Protocol, TYPE_CHECKING, runtime_checkable

if TYPE_CHECKING:
    from openfoam_driver.core.contracts.dictionary import DictEntry
    from openfoam_driver.core.runtime.models import TutorialSpec, CaseConfig, DataArtifact
    from openfoam_driver.planning_types import StrictDiagnostic
    from openfoam_driver.tutorials_display import TutorialDisplay
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

    def get_profile(self):
        """Return declarative case/C++ provenance metadata for this plugin."""
        ...
        
    def get_dict_entries(self) -> tuple[DictEntry, ...]:
        """
        Return the catalog of all solver-specific dictionary entries.
        Agents will use this to introspect capabilities deterministically.
        """
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


_REQUIRED_PLUGIN_MEMBERS = (
    "plugin_name",
    "plugin_id",
    "plugin_version",
    "plugin_api_version",
    "get_profile",
    "get_dict_entries",
    "get_dict_groups",
    "get_capabilities",
    "get_tutorial_catalog",
    "get_tutorial_displays",
    "validate_configuration",
    "validate_run_semantics",
    "predict_data_artifacts",
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
    if not _PLUGIN_ID_RE.fullmatch(plugin.plugin_id):
        raise TypeError(
            "SolverPlugin.plugin_id must use lowercase letters, digits, dots, "
            "or hyphens and cannot start or end with punctuation"
        )
    for name in (
        "get_dict_entries",
        "get_dict_groups",
        "get_capabilities",
        "get_tutorial_catalog",
        "get_tutorial_displays",
        "validate_configuration",
        "validate_run_semantics",
        "predict_data_artifacts",
        "get_profile",
    ):
        if not callable(getattr(plugin, name)):
            raise TypeError(f"SolverPlugin.{name} must be callable")
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
    """Load a trusted development plugin specified as ``module:Class``.

    This compatibility loader intentionally does not pretend to sandbox the
    import.  The CLI labels this form as an unsafe local-development option;
    installed plugin discovery can replace it in a later migration.
    """

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

    return load_plugin_context(
        "openfoam_driver.plugins.cardiacfoam_plugin:CardiacFoamPlugin"
    )


def generic_openfoam_context() -> DriverContext:
    """Return the built-in context with no solver-specific semantics."""

    from .generic_plugin import GenericOpenFOAMPlugin

    return driver_context(GenericOpenFOAMPlugin(), source="built-in:generic-openfoam")
