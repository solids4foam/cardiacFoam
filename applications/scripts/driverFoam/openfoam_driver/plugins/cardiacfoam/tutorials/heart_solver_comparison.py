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
#     heart_solver_comparison
#
# Description
#     Compares whole solver stacks (eikonal / monodomain / monodomain-eikonal
#     / bidomain) against ONE shared real heart anatomy (mesh + Purkinje
#     graph + material fields), held fixed in tutorials/heartSim3D-1D/
#     eikonalHeart -- reused in place per the owner's own instruction,
#     not duplicated. Each variant's config lives as a complete template
#     under that case's setup/solverVariants/<variant>/ (electroProperties,
#     fvSchemes, fvSolution, controlDict), since the differences between
#     solver stacks are whole structural blocks (different fields solved,
#     different coupler, different numerics), not small deltas -- the
#     apply_electro_property_overrides patch mechanism used elsewhere this
#     session cannot add or remove whole sub-dictionaries, only tweak
#     existing keys. Selector-driven synthesis (dict_builder.
#     build_electro_properties) was considered but rejected: nothing in
#     this codebase's own test suite has ever exercised it for a combined
#     myocardium + conductionNetworkDomains + domainCouplings + ecgDomains
#     configuration, an untested combination not worth risking here.
#
#     This is the catalog/registry the owner asked for: HEART_SOLVER_VARIANTS
#     names the four studies; each maps to a fixed template directory whose
#     four files apply_case swaps into constant/ and system/ verbatim.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

from __future__ import annotations

import shutil
import subprocess
from functools import partial
from pathlib import Path
from typing import Final

from openfoam_driver.core.runtime.models import CaseConfig, TutorialSpec
from openfoam_driver.core.runtime.parallel_execution import solve_steps
from openfoam_driver.specs.paths import resolve_spec_paths

# The owner's four named studies -- each key is a directory name under
# <case_root>/setup/solverVariants/. Adding a fifth study means adding a
# fifth template directory (four files: electroProperties, fvSchemes,
# fvSolution, controlDict) and a fifth entry here; nothing else changes.
HEART_SOLVER_VARIANTS: Final[tuple[str, ...]] = (
    "eikonal",
    "monodomain",
    "monodomain-eikonal",
    "bidomain",
)

_TEMPLATE_FILES: Final[tuple[str, ...]] = (
    "electroProperties",
    "fvSchemes",
    "fvSolution",
    "controlDict",
)

_DEFAULT_CASE_DIR_NAME: Final[str] = "heartSim3D-1D/eikonalHeart"


def _validate_variant(solver_variant: str) -> None:
    if solver_variant not in HEART_SOLVER_VARIANTS:
        known = ", ".join(HEART_SOLVER_VARIANTS)
        raise ValueError(f"solver_variant must be one of: {known}; got {solver_variant!r}")


def _template_dir(case_root: Path, solver_variant: str) -> Path:
    return case_root / "setup" / "solverVariants" / solver_variant


def _build_cases(*, solver_variant: str) -> list[CaseConfig]:
    return [CaseConfig(case_id=solver_variant, params={"solver_variant": solver_variant})]


def _apply_case(case_root: Path, case: CaseConfig) -> None:
    solver_variant = str(case.params["solver_variant"])
    template_dir = _template_dir(case_root, solver_variant)
    shutil.copy(template_dir / "electroProperties", case_root / "constant" / "electroProperties")
    for name in ("fvSchemes", "fvSolution", "controlDict"):
        shutil.copy(template_dir / name, case_root / "system" / name)


def _run_case(
    case_root: Path,
    setup_root: Path,
    case: CaseConfig,
    *,
    run_in_parallel: bool,
) -> None:
    # Structural TutorialSpec requirement only -- run --strict/sweep-run
    # execute the workflow_dag below exclusively, never this function.
    del setup_root, case
    if run_in_parallel:
        subprocess.run(["decomposePar"], cwd=case_root, check=True)
        subprocess.run(["mpirun", "-np", "6", "cardiacFoam", "-parallel"], cwd=case_root, check=True)
        subprocess.run(["reconstructPar"], cwd=case_root, check=True)
    else:
        subprocess.run(["cardiacFoam"], cwd=case_root, check=True)


def make_spec(
    *,
    tutorials_root: Path | None = None,
    case_dir_name: str = _DEFAULT_CASE_DIR_NAME,
    setup_dir_name: str | None = "setup",
    output_dir_name: str | None = None,
    solver_variant: str = "eikonal",
    run_in_parallel: bool = False,
) -> TutorialSpec:
    _validate_variant(solver_variant)

    case_root, setup_root, output_dir = resolve_spec_paths(
        tutorials_root=tutorials_root,
        case_dir_name=case_dir_name,
        setup_dir_name=setup_dir_name,
        output_dir_name=output_dir_name,
        default_output_dir_name=f"sweepRun_{solver_variant}",
    )

    steps, _final_id = solve_steps(
        solve_id="solve",
        solve_command="cardiacFoam",
        depends_on=[],
        run_in_parallel=run_in_parallel,
        case_root=case_root,
    )

    return TutorialSpec(
        name=f"heartSolverComparison-{solver_variant}",
        case_root=case_root,
        setup_root=setup_root,
        output_dir=output_dir,
        build_cases=partial(_build_cases, solver_variant=solver_variant),
        apply_case=_apply_case,
        run_case=partial(_run_case, run_in_parallel=run_in_parallel),
        metadata={
            "notes": "Solver-stack comparison over one shared real heart anatomy.",
            "workflow_dag": {"steps": steps},
            "solver_variant": solver_variant,
            "run_in_parallel": run_in_parallel,
        },
    )
