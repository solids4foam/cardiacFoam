"""Samplable fields and case-model resolution are plugin knowledge."""

from __future__ import annotations

from pathlib import Path

from openfoam_driver.core.plugin_interface import (
    default_driver_context,
    generic_openfoam_context,
)


def _cardiac_case(root: Path) -> Path:
    (root / "constant").mkdir(parents=True)
    (root / "constant" / "electroProperties").write_text(
        "myocardiumSolver monodomainSolver;\n"
    )
    return root


def test_generic_plugin_exposes_no_cardiac_fields(tmp_path: Path) -> None:
    introspection = generic_openfoam_context().capabilities.case_introspection
    resolved = introspection.resolve_case_models(_cardiac_case(tmp_path))
    fields = introspection.samplable_fields(resolved)
    flat = [name for names in fields.values() for name in names]
    assert "Vm" not in flat
    assert "phiE" not in flat


def test_cardiac_plugin_exposes_its_fixed_fields(tmp_path: Path) -> None:
    introspection = default_driver_context().capabilities.case_introspection
    resolved = introspection.resolve_case_models(_cardiac_case(tmp_path))
    assert resolved["solver"] == "monodomainSolver"
    electro = introspection.samplable_fields(resolved)["electro"]
    assert "Vm" in electro
    assert "activationTime" in electro


def test_missing_case_file_resolves_to_none_without_raising(tmp_path: Path) -> None:
    introspection = default_driver_context().capabilities.case_introspection
    resolved = introspection.resolve_case_models(tmp_path)
    assert resolved == {"solver": None, "ionic_model": None, "active_tension": None}


def test_no_active_tension_means_no_solid_region(tmp_path: Path) -> None:
    """A spatial EP solver alone must not imply a mechanics region."""
    introspection = default_driver_context().capabilities.case_introspection
    resolved = introspection.resolve_case_models(_cardiac_case(tmp_path))
    assert introspection.samplable_fields(resolved)["solid"] == ()


def test_deprecated_tuple_shim_still_returns_three_values(tmp_path: Path) -> None:
    from openfoam_driver.core.capability_manifest import resolve_case_models

    assert resolve_case_models("/nonexistent/case") == (None, None, None)
    assert resolve_case_models(_cardiac_case(tmp_path))[0] == "monodomainSolver"
