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
        "purkinje_solver": "eikonalSolver",
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
        "myocardium_solver": "eikonalSolver",
        "purkinje_solver": "monodomain1DSolver",
        "required_coupler": None,
        "valid": False,
        "reason": "Incompatible physics: eikonal myocardium cannot couple to reaction-diffusion Purkinje",
    },
    {
        "myocardium_solver": "bidomainSolver",
        "purkinje_solver": "*",
        "required_coupler": None,
        "valid": False,
        "reason": "bidomainSolver does not support Purkinje network coupling",
    },
    {
        "myocardium_solver": "singleCellSolver",
        "purkinje_solver": "*",
        "required_coupler": None,
        "valid": False,
        "reason": "singleCellSolver has no PDE domain; Purkinje coupling not applicable",
    },
)
