#----------------------------------------------------------------------------#
# Module
#     linear_solves
#
# Description
#     The largest iteration count and final residual of each linear solve in
#     the log.cardiacFoam files of a run copy. Imported by summarize.py of the
#     bidomain studies and of the bidomain sphericalShell studies.
#----------------------------------------------------------------------------#

from __future__ import annotations

import re
from pathlib import Path

SOLVE_RE = re.compile(
    r"Solving for (\w+), Initial residual = \S+, "
    r"Final residual = ([0-9.eE+-]+), No Iterations (\d+)"
)


def read_solves(run: Path) -> dict[str, tuple[int, float]]:
    """Largest iteration count and final residual of each solved field."""
    solves: dict[str, tuple[int, float]] = {}
    for log in sorted(run.glob("log.cardiacFoam*")):
        for field, residual, iterations in SOLVE_RE.findall(log.read_text()):
            its, res = solves.get(field, (0, 0.0))
            solves[field] = (max(its, int(iterations)), max(res, float(residual)))
    return solves
