"""Required case files come from the plugin profile, not a core constant."""

from __future__ import annotations

from openfoam_driver.core.plugin_interface import (
    default_driver_context,
    generic_openfoam_context,
)


def test_cardiac_required_files_include_its_dictionaries() -> None:
    contract = default_driver_context().capabilities.case_files
    required = contract.required_files()
    assert "constant/electroProperties" in required
    assert "constant/physicsProperties" in required
    assert "system/controlDict" in required


def test_generic_plugin_requires_no_solver_dictionaries() -> None:
    contract = generic_openfoam_context().capabilities.case_files
    required = contract.required_files()
    assert "constant/electroProperties" not in required
    assert "constant/physicsProperties" not in required


def test_conditional_files_are_separated_from_required() -> None:
    contract = default_driver_context().capabilities.case_files
    conditional = contract.conditional_files()
    assert "system/blockMeshDict" in conditional
    assert "system/blockMeshDict" not in contract.required_files()
