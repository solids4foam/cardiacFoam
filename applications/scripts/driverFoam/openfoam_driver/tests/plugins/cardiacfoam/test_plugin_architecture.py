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
#     test_plugin_architecture
#
# Description
#     Verifies the robustness of the driverFOAM plugin contract.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

import pytest

from openfoam_driver.core.plugin_interface import SolverPlugin
from openfoam_driver.plugins.cardiacfoam_plugin import CardiacFoamPlugin
from openfoam_driver.core.runtime.models import CaseConfig

def test_cardiacfoam_plugin_satisfies_protocol():
    """
    Test that the CardiacFoamPlugin implements the SolverPlugin protocol 
    correctly and returns the expected structured data types.
    """
    plugin: SolverPlugin = CardiacFoamPlugin()
    
    assert plugin.plugin_name == "cardiacFoam"
    
    # 1. Test Dict Entries Contract
    entries = plugin.get_dict_entries()
    assert isinstance(entries, tuple), "get_dict_entries must return a tuple"
    assert len(entries) > 0, "Plugin must expose dictionary entries"
    assert hasattr(entries[0], "driver_path"), "Entries must be DictEntry objects"
    
    # 3. Test Tutorials Contract
    catalog = plugin.get_tutorial_catalog()
    assert isinstance(catalog, dict), "get_tutorial_catalog must return a dict"
    displays = plugin.get_tutorial_displays()
    assert isinstance(displays, tuple), "get_tutorial_displays must return a tuple"
    assert "registered_tutorials" in catalog, "Catalog must contain registered_tutorials"
    assert "spec_factories" in catalog, "Catalog must contain spec_factories"
    
    
    # 2. Test Capabilities Contract
    capabilities = plugin.get_capabilities()
    assert capabilities is not None
    
    
    # 4. Test Validation Contract
    from unittest.mock import MagicMock
    from openfoam_driver.core.runtime.models import TutorialSpec
    dummy_spec = MagicMock(spec=TutorialSpec)
    dummy_spec.case_root = "/tmp/dummy"
    diagnostics = plugin.validate_configuration(dummy_spec)


def test_mock_solver_plugin():
    """
    Simulate an entirely new solver (e.g., heatTransferFoam) injecting its
    own plugin, ensuring the contract is agnostic to the physics domain.
    """
    class MockHeatTransferPlugin:
        @property
        def plugin_name(self) -> str:
            return "heatTransferFoam"
            
        def get_dict_entries(self):
            return ()
            
        def get_capabilities(self):
            class MockCapabilities:
                pass
            return MockCapabilities()
            
        def get_tutorials(self):
            return ()
            
        def validate_configuration(self, config):
            return ()
            
    plugin: SolverPlugin = MockHeatTransferPlugin()
    assert plugin.plugin_name == "heatTransferFoam"
    assert plugin.get_dict_entries() == ()
