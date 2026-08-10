"""cardiacFoam mutation hook for the legacy generic-case specification."""

from __future__ import annotations

from pathlib import Path

from openfoam_driver.specs.common import (
    apply_electro_property_overrides,
    apply_physics_property_overrides,
)


def apply_case_mutation(
    case_root: Path,
    case,
    *,
    electro_properties_relpath: Path,
    physics_properties_relpath: Path,
) -> None:
    apply_electro_property_overrides(
        case_root / electro_properties_relpath,
        case.params.get("electro_property_overrides"),
    )
    apply_physics_property_overrides(
        case_root / physics_properties_relpath,
        case.params.get("physics_property_overrides"),
    )
