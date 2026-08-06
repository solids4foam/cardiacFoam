import json
from pathlib import Path

import pytest
from openfoam_driver.tests.conftest import skip_without_monorepo
pytestmark = skip_without_monorepo

REPO_ROOT = Path(__file__).resolve().parent.parent.parent.parent.parent.parent.parent


def test_paperI_dx_json_files_are_valid():
    """Verify the paperI explicit/implicit run configurations are valid JSON
    and contain the expected 'niederer2012' configuration object."""
    convergence_dir = REPO_ROOT / "tutorials" / "NiedererEtAl2011" / "NiedererEtAl2011verification" / "setup" / "studies" / "cartesianConvergence"
    
    expected_files = [
        "paperI_dx05_implicit.json",
        "paperI_dx02_implicit.json",
        "paperI_dx01_implicit.json"
    ]
    
    for filename in expected_files:
        json_path = convergence_dir / filename
        assert json_path.is_file(), f"Missing required Paper I config: {json_path}"
        
        with open(json_path) as f:
            data = json.load(f)
            
        assert "niederer2012" in data, f"Config {filename} missing 'niederer2012' root key"
        config = data["niederer2012"]
        assert "solvers" in config
        assert "dx_values" in config
        assert "dt_values" in config
        assert "output_dir_name" in config


def test_paperI_sweep_hex_convergence_is_valid():
    """Verify the bathBidomain sweep spec used in Paper I is valid JSON."""
    sweep_spec_path = REPO_ROOT / "tutorials" / "manufacturedSolutions" / "bathBidomain" / "setup" / "studies" / "cartesianConvergence" / "sweep_hex_convergence.json"
    
    assert sweep_spec_path.is_file(), f"Missing required sweep spec: {sweep_spec_path}"
    
    with open(sweep_spec_path) as f:
        data = json.load(f)
        
    assert "base" in data
    assert "sweep" in data
    assert "mode" in data["sweep"]
    assert data["sweep"]["mode"] == "zip"
