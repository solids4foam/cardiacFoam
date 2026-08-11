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
#     capability_manifest
#
# Description
#     Assembles the machine-readable surface of what the driver will accept:
#     the allowed workflow-command set and the field names each solver can
#     sample. Sourced only from the single owners (the core-neutral workflow
#     allowlist plus the calling plugin's commands, utility manifests, and
#     ionic / active-tension catalogs) so the manifest cannot drift from the
#     enforcers that actually gate execution.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

from __future__ import annotations

from pathlib import Path
from typing import Any, Iterable

from .core.runtime.workflow import (
    CASE_SCRIPT_COMMANDS,
    CORE_NEUTRAL_COMMANDS,
)

# Fixed fields the solvers expose regardless of the ionic / active-tension model
# (documented under "Function objects" in AGENT_GUIDE.md). Model-specific field
# names come from the catalogs below.
_ELECTRO_SOLVER_FIELDS = ("Vm", "activationTime", "Iion", "phiE", "phiI")
_SOLID_SOLVER_FIELDS = ("Ta", "lambda")


def resolve_case_models(
    case_root: str | Path,
) -> tuple[str | None, str | None, str | None]:
    """Best-effort ``(solver, ionic_model, active_tension)`` from a case's
    ``constant/electroProperties``. Any of the three may be ``None`` when the
    file is absent or the entry is not declared; never raises. Used by the
    discovery paths (``describe_entry``, ``strict_plan``) so the capability
    manifest can name the resolved model's fields."""
    from .specs.common import (
        detect_active_tension_model_name,
        detect_ionic_model_name,
        detect_myocardium_solver_name,
    )

    electro_path = Path(case_root) / "constant" / "electroProperties"
    if not electro_path.exists():
        return None, None, None
    solver = ionic = active_tension = None
    try:
        solver = detect_myocardium_solver_name(electro_path)
    except (OSError, KeyError):
        solver = None
    try:
        ionic = detect_ionic_model_name(electro_path)
    except (OSError, KeyError):
        ionic = None
    try:
        active_tension = detect_active_tension_model_name(electro_path)
    except (OSError, KeyError):
        active_tension = None
    return solver, ionic, active_tension


def _utility_commands(utility_manifests: dict[str, Any]) -> dict[str, list[str]]:
    """The plugin utility commands the allowlist accepts, keyed to what they
    produce. Mirrors the acceptance rule in ``validate_workflow_commands``
    (only utilities that declare ``produces`` are accepted)."""

    return {
        command: [produce.artifact_id for produce in manifest.produces]
        for command, manifest in utility_manifests.items()
        if manifest.produces
    }


def build_capability_manifest(
    *,
    resolved_solver: str | None = None,
    resolved_ionic_model: str | None = None,
    resolved_active_tension: str | None = None,
    ionic_model_catalog: dict | None = None,
    active_tension_model_catalog: dict | None = None,
    plugin_commands: Iterable[str] = (),
    utility_manifests: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Return the driver's accept-surface as a plain JSON-able dict.

    ``plugin_commands`` and ``utility_manifests`` are supplied by the calling
    plugin (the same injection style as the model catalogs) so that core names
    no solver here; together with :data:`CORE_NEUTRAL_COMMANDS` they reproduce
    exactly what ``validate_workflow_commands`` accepts for that plugin.

    ``allowed_commands`` names exactly what a workflow DAG step may invoke;
    ``samplable_fields`` names the fields a function object may sample for the
    resolved model, split by region. Unknown model names are ignored (the
    manifest degrades to the fixed solver fields) rather than raising, so this
    is always safe to call during discovery.
    """

    electro = set(_ELECTRO_SOLVER_FIELDS)
    ionic_entry = (
        (ionic_model_catalog or {}).get(resolved_ionic_model) if resolved_ionic_model else None
    )
    if ionic_entry is not None:
        electro.update(ionic_entry.states)
        electro.update(ionic_entry.algebraic)
        electro.update(ionic_entry.recommended_exports)

    solid: set[str] = set()
    # A spatial active-tension model is positive evidence of electromechanical
    # coupling. A spatial EP solver alone does not imply a mechanics region.
    has_solid_region = (
        resolved_active_tension is not None
        and resolved_solver is not None
        and resolved_solver != "singleCellSolver"
    )
    if has_solid_region:
        solid.update(_SOLID_SOLVER_FIELDS)
    at_entry = (
        (active_tension_model_catalog or {}).get(resolved_active_tension)
        if resolved_active_tension
        else None
    )
    if has_solid_region and at_entry is not None:
        solid.update(at_entry.states)
        solid.update(at_entry.algebraic)

    return {
        "allowed_commands": {
            "core": sorted(set(CORE_NEUTRAL_COMMANDS) | set(plugin_commands)),
            "case_scripts": sorted(CASE_SCRIPT_COMMANDS),
            "utilities": _utility_commands(utility_manifests or {}),
            "installed_openfoam_apps_note": (
                "When OpenFOAM is sourced, any executable under $FOAM_APPBIN or "
                "$FOAM_USER_APPBIN is also accepted (core apps + your compiled "
                "utilities). Unsourced, only core + case_scripts + utilities apply."
            ),
        },
        "samplable_fields": {
            "electro": sorted(electro),
            "solid": sorted(solid),
            "note": (
                "Function objects are OpenFOAM's; these are the field NAMES this "
                "solver exposes. Sampling a name not listed here is silently "
                "dropped by the solver (strict planning emits unknown_sampled_field "
                "warnings for such names)."
            ),
        },
    }
