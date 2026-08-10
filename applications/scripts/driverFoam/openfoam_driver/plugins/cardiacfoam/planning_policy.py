"""cardiacFoam-specific strict-planning policy decisions."""

from __future__ import annotations

from pathlib import Path

from openfoam_driver.specs.common import (
    detect_myocardium_solver_name,
    detect_verification_model_type,
)


def is_nondimensional_case(spec) -> bool:
    electro_path = Path(spec.case_root) / "constant" / "electroProperties"
    if not electro_path.exists():
        return False
    try:
        return (
            detect_myocardium_solver_name(electro_path) == "singleCellSolver"
            or detect_verification_model_type(electro_path) is not None
        )
    except Exception:
        return False
