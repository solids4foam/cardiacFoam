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

from typing import Protocol, TYPE_CHECKING

if TYPE_CHECKING:
    from openfoam_driver.dict_entries import DictEntry
    from openfoam_driver.core.runtime.models import TutorialSpec, CaseConfig, DataArtifact
    from openfoam_driver.planning_types import StrictDiagnostic
    from openfoam_driver.tutorials_display import TutorialDisplay
    from pathlib import Path


class CapabilityManifest(Protocol):
    """Protocol for a solver's capability manifest."""
    # This can be expanded based on the solver's specific domain (e.g. models, physics)
    pass


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

    def predict_data_artifacts(self, case_root: Path, spec: TutorialSpec) -> tuple[DataArtifact, ...]:
        """
        Predict the domain-specific artifacts (like ECGs or Purkinje VTK files) 
        that this solver expects to produce.
        """
        ...

_ACTIVE_PLUGIN: SolverPlugin | None = None

def set_active_plugin(plugin: SolverPlugin) -> None:
    """Inject the active solver plugin for this session."""
    global _ACTIVE_PLUGIN
    _ACTIVE_PLUGIN = plugin

def get_active_plugin() -> SolverPlugin:
    """
    Retrieve the globally active solver plugin. 
    Falls back to CardiacFoamPlugin if none is set to maintain backward compatibility.
    """
    global _ACTIVE_PLUGIN
    if _ACTIVE_PLUGIN is None:
        from openfoam_driver.plugins.cardiacfoam_plugin import CardiacFoamPlugin
        _ACTIVE_PLUGIN = CardiacFoamPlugin()
    return _ACTIVE_PLUGIN
