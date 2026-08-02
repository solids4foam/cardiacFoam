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
from openfoam_driver.core.plugin_interface import SolverPlugin, CapabilityManifest

# Note: now imported from the local plugin catalog instead of dict_entries
from openfoam_driver.plugins.cardiacfoam.dict_entries_catalog import ELECTRO_PROPERTY_ENTRY_GROUPS, HETEROGENEITY_MODELS
from openfoam_driver.dict_entries import PHYSICS_PROPERTY_ENTRIES
from openfoam_driver.plugins.cardiacfoam.active_tension_catalog import ACTIVE_TENSION_MODEL_CATALOG
from openfoam_driver.plugins.cardiacfoam.ionic_model_catalog import IONIC_MODEL_CATALOG
from openfoam_driver.capability_manifest import build_capability_manifest
from openfoam_driver.plugins.cardiacfoam.solver_coupling import SOLVER_COMPATIBILITY_RULES
from openfoam_driver.core.runtime.registry import list_tutorials
from openfoam_driver.planning_types import StrictDiagnostic, diagnostic

if TYPE_CHECKING:
    from openfoam_driver.dict_entries import DictEntry
    from openfoam_driver.core.runtime.models import TutorialSpec, CaseConfig


class CardiacFoamPlugin:
    """
    Plugin implementation for cardiacFoam.
    Provides domain-specific dictionaries, tutorials, and capabilities 
    to the generic driverFOAM engine.
    """
    
    @property
    def plugin_name(self) -> str:
        return "cardiacFoam"
        
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
        # Note: A complete implementation would also aggregate other cardiac-specific dicts
        return tuple(entries)

    def get_capabilities(self) -> CapabilityManifest:
        """
        Return the cardiacFoam capabilities (models, solvers, etc.).
        """
        manifest = build_capability_manifest(ionic_model_catalog=IONIC_MODEL_CATALOG, active_tension_model_catalog=ACTIVE_TENSION_MODEL_CATALOG)
        manifest["heterogeneity_models"] = HETEROGENEITY_MODELS
        manifest["ionic_models"] = IONIC_MODEL_CATALOG
        manifest["active_tension_models"] = ACTIVE_TENSION_MODEL_CATALOG
        manifest["solver_compatibility_rules"] = SOLVER_COMPATIBILITY_RULES
        return manifest

    def get_tutorials(self) -> tuple[TutorialSpec, ...]:
        """
        Return the list of tutorials provided by cardiacFoam.
        """
        return tuple(list_tutorials())

    def validate_configuration(self, config: CaseConfig) -> tuple[StrictDiagnostic, ...]:
        """
        Perform any cardiacFoam-specific validation (e.g., solver coupling rules).
        """
        diagnostics = []
        
        # Example validation: Check if a selected ionic model is compatible with the solver
        params = config.params
        
        # This is a stub for where complex, domain-specific validation logic goes.
        # e.g. checking SOLVER_COMPATIBILITY_RULES from solver_coupling.py
        
        return tuple(diagnostics)


# Ensure CardiacFoamPlugin satisfies the SolverPlugin protocol
def _check_protocol() -> None:
    plugin: SolverPlugin = CardiacFoamPlugin()
