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
#     plugins.cardiacfoam.mesh_provisioning
#
# Description
#     Chooses and provisions a mesh for a cardiacFoam case built from scratch
#     by build_and_launch, based on the case's myocardiumSolver.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

"""electroModel.C requires a real fvMesh regardless of solver
(`refCast<const fvMesh>(mesh())` at electroModel.C:344) -- even
singleCellSolver needs one, which is why the real singleCell tutorial ships a
static trivial 1-cell `constant/polyMesh/`. Neither `build_and_launch` nor
`sweep_runner.materialize_case` provisioned any mesh for a from-scratch
`case_folder` before this module existed.

Two provisioning strategies, chosen by `myocardiumSolver`:
  - `singleCellSolver` has no real spatial geometry: copy the bundled static
    1-cell polyMesh fixture directly into `constant/polyMesh/`.
  - `monodomainSolver`/`bidomainSolver`/`eikonalSolver` need real geometry:
    write the generic default `system/blockMeshDict` that
    `specs/mesh_provisioning.py` renders (a small slab, "walls" patch) and let
    `blockMesh` generate the mesh at run time. This is a sane generic default,
    not a scientifically tuned geometry for any specific tutorial -- callers
    that need particular dimensions should still author their own
    blockMeshDict (this only fills the gap for a case built purely from
    selectors/overrides, which have no geometry concept at all today).

Which solver needs which strategy is cardiacFoam vocabulary, so it lives here;
the solver-neutral halves -- rendering the default `blockMeshDict` and the
`dx`-to-cell-count arithmetic -- stay in `specs/mesh_provisioning.py`.
"""

from __future__ import annotations

import shutil
from pathlib import Path

from openfoam_driver.specs.mesh_provisioning import default_block_mesh_dict_text

_FIXTURES_DIR = Path(__file__).parent / "fixtures"
_SINGLE_CELL_POLYMESH_DIR = _FIXTURES_DIR / "single_cell_polymesh"

MESHLESS_SOLVERS = frozenset({"singleCellSolver"})
BLOCK_MESH_SOLVERS = frozenset({"monodomainSolver", "bidomainSolver", "eikonalSolver"})

_POLYMESH_FILES = ("points", "faces", "owner", "neighbour", "boundary")


def provision_mesh(
    *, case_dir: Path, myocardium_solver: str, dx_m: float | None = None,
) -> bool:
    """Provision whatever mesh `myocardium_solver` needs under `case_dir`.

    Returns True if the case now needs a `blockMesh` run before solving
    (a `system/blockMeshDict` was written or already exists), False if a
    concrete mesh was copied directly (or the solver needs no mesh at all,
    which no current solver does).

    Unlike `build_and_launch`'s other generated files, a mesh is never
    clobbered on a repeat call regardless of that call's own `overwrite`
    flag: re-materializing a case (e.g. a retried sweep-run case) should not
    wipe out an already-valid mesh, and a hand-authored custom
    blockMeshDict/polyMesh must never be silently replaced by this generic
    default.

    `dx_m` (metres) only means something for the generic default
    `blockMeshDict` (`BLOCK_MESH_SOLVERS`) -- it is meaningless for
    `MESHLESS_SOLVERS` (no spatial geometry at all) and rejected outright
    rather than silently having no effect, and it has no bearing on real
    anatomical meshes imported via `vtkUnstructuredToFoam`, which this
    function never touches.
    """
    if myocardium_solver in MESHLESS_SOLVERS:
        if dx_m is not None:
            raise ValueError(
                f"dx has no effect for myocardiumSolver={myocardium_solver!r} "
                "(no spatial mesh -- it has no geometry for dx to resolve)."
            )
        poly_mesh_dir = case_dir / "constant" / "polyMesh"
        already_present = all((poly_mesh_dir / name).exists() for name in _POLYMESH_FILES)
        if not already_present:
            poly_mesh_dir.mkdir(parents=True, exist_ok=True)
            for name in _POLYMESH_FILES:
                shutil.copyfile(_SINGLE_CELL_POLYMESH_DIR / name, poly_mesh_dir / name)
        return False

    if myocardium_solver in BLOCK_MESH_SOLVERS:
        block_mesh_dict = case_dir / "system" / "blockMeshDict"
        if not block_mesh_dict.exists():
            block_mesh_dict.parent.mkdir(parents=True, exist_ok=True)
            block_mesh_dict.write_text(default_block_mesh_dict_text(dx_m=dx_m))
        return True

    # Unknown/future solver: leave mesh provisioning to the caller.
    return False
