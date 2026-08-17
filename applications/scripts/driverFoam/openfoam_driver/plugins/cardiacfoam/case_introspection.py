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
#     case_introspection
#
# Description
#     Cardiac-specific case-model resolution and the field names each
#     resolved model exposes for sampling. Reads a case's
#     constant/electroProperties and unions the fixed solver fields with the
#     resolved ionic / active-tension catalog entries. Owned by the plugin so
#     core (capability_manifest.py) holds no solver knowledge.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

from __future__ import annotations

from pathlib import Path

# Fixed fields the cardiac solvers expose regardless of the ionic /
# active-tension model. Model-specific names come from the catalogs.
_ELECTRO_SOLVER_FIELDS = ("Vm", "activationTime", "Iion", "phiE", "phiI")

# Active-tension / electromechanics coupling fields -- this repo's own
# catalog, always available regardless of build configuration.
_SOLID_SOLVER_FIELDS = ("Ta", "lambda")

# NOT YET ADDED: base solids4foam mechanics fields (D, DD, sigmaHyd, ...).
# tutorials/manufacturedSolutions/monodomainTotalLagrangianEM verifies these
# are the real, correct names (system/solid/fvSolution solves "D|DD|sigmaHyd";
# its own postprocessing checks {"Vm", "D", "lambda", "Ta"} with L1/L2/Linf
# error norms on D) -- so this is a known-good list, not a guess. Deliberately
# withheld from _SOLID_SOLVER_FIELDS until solids4foam is a build
# configuration this repo can actually run everywhere this catalog is
# consulted (some builds set FORCE_LIGHTWEIGHT_PHYSICSMODEL=1 and never
# link solids4foam in at all -- see buildAndTest.yml), so that "samplable"
# never claims a field a given build genuinely cannot produce.


def resolve_case_models(case_root: str | Path) -> dict[str, str | None]:
    """Best-effort resolution from ``constant/electroProperties``. Never raises;
    any of the three values may be ``None`` when the file or entry is absent."""
    from openfoam_driver.plugins.cardiacfoam.detection import (
        detect_active_tension_model_name,
        detect_ionic_model_name,
        detect_myocardium_solver_name,
    )

    electro_path = Path(case_root) / "constant" / "electroProperties"
    if not electro_path.exists():
        return {"solver": None, "ionic_model": None, "active_tension": None}
    resolved: dict[str, str | None] = {}
    for key, detect in (
        ("solver", detect_myocardium_solver_name),
        ("ionic_model", detect_ionic_model_name),
        ("active_tension", detect_active_tension_model_name),
    ):
        try:
            resolved[key] = detect(electro_path)
        except (OSError, KeyError):
            resolved[key] = None
    return resolved


def samplable_fields(resolved: dict[str, str | None]) -> dict[str, tuple[str, ...]]:
    """Field names the resolved cardiac model exposes, by region."""
    from openfoam_driver.plugins.cardiacfoam.active_tension_catalog import (
        ACTIVE_TENSION_MODEL_CATALOG,
    )
    from openfoam_driver.plugins.cardiacfoam.ionic_model_catalog import (
        IONIC_MODEL_CATALOG,
    )

    electro = set(_ELECTRO_SOLVER_FIELDS)
    ionic_entry = IONIC_MODEL_CATALOG.get(resolved.get("ionic_model") or "")
    if ionic_entry is not None:
        electro.update(ionic_entry.states)
        electro.update(ionic_entry.algebraic)
        electro.update(ionic_entry.recommended_exports)

    solid: set[str] = set()
    # A spatial active-tension model is positive evidence of electromechanical
    # coupling. A spatial EP solver alone does not imply a mechanics region.
    # This is the only detection signal available today -- every solid-region
    # case in this repo also declares an active-tension model. It is not a
    # guarantee: a hypothetical passive-only mechanics case (no active
    # contraction) would have a genuine solid region this check would miss.
    active_tension = resolved.get("active_tension")
    solver = resolved.get("solver")
    has_solid_region = (
        active_tension is not None
        and solver is not None
        and solver != "singleCellSolver"
    )
    if has_solid_region:
        solid.update(_SOLID_SOLVER_FIELDS)
        at_entry = ACTIVE_TENSION_MODEL_CATALOG.get(active_tension or "")
        if at_entry is not None:
            solid.update(at_entry.states)
            solid.update(at_entry.algebraic)

    return {"electro": tuple(sorted(electro)), "solid": tuple(sorted(solid))}
