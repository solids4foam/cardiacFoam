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
#     solver_coupling
#
# Description
#     Validates and enforces compatibility rules between solver pairs.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

"""Cross-domain solver coupling rules.

Defines which `myocardiumSolver × purkinjeSolver × required_coupler`
combinations are physically valid in cardiacFoam, and explains why the
invalid ones are forbidden. Lives outside `ionic_model_catalog.py` because
the rules are about *solver* compatibility, not about any individual ionic
model — keeping them here makes it clearer to readers (and agents) which
concern owns the data.

Consumer today: `introspection.py` exposes the table to LLM agents through
the describe-tutorial JSON payload so they can reject incompatible solver
combinations before launching a run.
"""
from __future__ import annotations

from typing import Final


SOLVER_COMPATIBILITY_RULES: Final[tuple[dict, ...]] = (
    {
        "myocardium_solver": "monodomainSolver",
        "purkinje_solver": "monodomain1DSolver",
        "required_coupler": "reactionDiffusionPvjCoupler",
        "valid": True,
    },
    {
        "myocardium_solver": "eikonalSolver",
        "purkinje_solver": "eikonalSolver1D",
        "required_coupler": "eikonalPvjCoupler",
        "valid": True,
    },
    {
        "myocardium_solver": "eikonalSolver",
        "purkinje_solver": "restitutionEikonalSolver1D",
        "required_coupler": "eikonalPvjCoupler",
        "valid": True,
    },
    {
        "myocardium_solver": "monodomainSolver",
        "purkinje_solver": "eikonalSolver",
        "required_coupler": None,
        "valid": False,
        "reason": "Incompatible physics: reaction-diffusion myocardium cannot couple to eikonal Purkinje",
    },
    {
        "myocardium_solver": "monodomainSolver",
        "purkinje_solver": "eikonalSolver1D",
        "required_coupler": "eikonalMonodomainPvjCoupler",
        "valid": True,
    },
    {
        "myocardium_solver": "monodomainSolver",
        "purkinje_solver": "restitutionEikonalSolver1D",
        "required_coupler": "eikonalMonodomainPvjCoupler",
        "valid": True,
    },
    {
        "myocardium_solver": "eikonalSolver",
        "purkinje_solver": "monodomain1DSolver",
        "required_coupler": None,
        "valid": False,
        "reason": "Incompatible physics: eikonal myocardium cannot couple to reaction-diffusion Purkinje",
    },
    {
        "myocardium_solver": "bidomainSolver",
        "purkinje_solver": "monodomain1DSolver",
        "required_coupler": "reactionDiffusionPvjCoupler",
        "valid": True,
    },
    {
        "myocardium_solver": "bidomainSolver",
        "purkinje_solver": "*",
        "required_coupler": None,
        "valid": False,
        "reason": "bidomainSolver only supports reaction-diffusion coupling",
    },
    {
        "myocardium_solver": "singleCellSolver",
        "purkinje_solver": "*",
        "required_coupler": None,
        "valid": False,
        "reason": "singleCellSolver has no PDE domain; Purkinje coupling not applicable",
    },
)
