"""Build the legacy RunDocument-v2 cardiac configuration from case files.

This is the unchanged electroProperties/physicsProperties parser formerly
embedded in core.  RunDocument v2 remains cardiac-shaped during Plan 1; Plan 2
may replace this capability with a solver-neutral v3 configuration contract.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

from foamlib import FoamFile

from openfoam_driver.core.planning_types import StrictDiagnostic, diagnostic
from openfoam_driver.specs.dict_builder import populate_values
from openfoam_driver.plugins.cardiacfoam.dict_builder import (
    build_electro_properties,
    build_physics_properties,
    parse_electro_properties,
    resolve_context,
    select_applicable_entries,
)
from openfoam_driver.specs.validation import primary_phase, slot_key


def _read_physics_type(path: Path) -> str | None:
    if not path.exists():
        return None
    try:
        return str(FoamFile(path)["type"])
    except (KeyError, ValueError):
        return None


def build_config(spec) -> tuple[dict[str, dict[str, Any]], tuple[StrictDiagnostic, ...]]:
    diagnostics: list[StrictDiagnostic] = []
    config: dict[str, dict[str, Any]] = {
        "anatomy": {},
        "physics": {},
        "stimulus": {},
        "solver": {},
    }
    case_root = Path(spec.case_root)
    generic_case = bool(spec.metadata.get("generic_case")) if spec.metadata else False
    if generic_case:
        return config, ()

    electro_path = case_root / "constant" / "electroProperties"
    physics_path = case_root / "constant" / "physicsProperties"
    physics_type = _read_physics_type(physics_path)
    if physics_type is None:
        diagnostics.append(diagnostic(
            "error",
            "missing_physics_properties",
            f"Could not read physicsProperties type from {physics_path}",
            source=str(physics_path),
            field="type",
        ))
    else:
        config["physics"]["type"] = physics_type
        try:
            build_physics_properties({"type": physics_type})
        except Exception as exc:
            diagnostics.append(diagnostic(
                "error", "invalid_physics_properties", str(exc), source=str(physics_path),
            ))

    if not electro_path.exists():
        diagnostics.append(diagnostic(
            "error",
            "missing_electro_properties",
            f"Missing electroProperties at {electro_path}",
            source=str(electro_path),
        ))
    else:
        try:
            parsed = parse_electro_properties(electro_path)
            selectors = parsed["selectors"]
            overrides = parsed.get("overrides", {})
            try:
                build_electro_properties(selectors, overrides=overrides or None)
            except Exception as exc:
                diagnostics.append(diagnostic(
                    "error", "invalid_electro_properties", str(exc), source=str(electro_path),
                ))
            context = resolve_context(selectors, overrides=overrides or None)
            applicable_entries = select_applicable_entries(context)
            populated = populate_values(applicable_entries, context)
            for entry_obj in applicable_entries:
                key = slot_key(entry_obj.driver_path)
                if entry_obj.dynamic_path and key not in context:
                    continue
                if key not in populated:
                    continue
                phase = primary_phase(entry_obj) or "physics"
                config[phase][key] = populated[key]
        except Exception as exc:
            diagnostics.append(diagnostic(
                "error", "unparseable_electro_properties", str(exc), source=str(electro_path),
            ))

    return config, tuple(diagnostics)
