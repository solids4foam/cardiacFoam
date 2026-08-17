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
#     test_solver_mesh_provisioning
#
# Description
#     Tests the myocardiumSolver-keyed mesh provisioning strategy for
#     from-scratch case_folder cases.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

import re

import pytest

from openfoam_driver.plugins.cardiacfoam.mesh_provisioning import provision_mesh
from openfoam_driver.specs.mesh_provisioning import default_block_mesh_dict_text


def _cell_counts(text: str) -> tuple[int, int, int]:
    match = re.search(r"hex \([^)]*\)\s*\((\d+)\s+(\d+)\s+(\d+)\)", text)
    assert match is not None, text
    return tuple(int(g) for g in match.groups())


def test_provision_mesh_spatial_solver_honours_dx(tmp_path):
    case_dir = tmp_path / "case"
    provision_mesh(case_dir=case_dir, myocardium_solver="monodomainSolver", dx_m=0.0004)
    text = (case_dir / "system" / "blockMeshDict").read_text()
    default_cells = _cell_counts(default_block_mesh_dict_text())
    assert _cell_counts(text)[0] > default_cells[0]


def test_provision_mesh_rejects_dx_for_meshless_solver(tmp_path):
    # singleCellSolver has no spatial geometry at all -- dx would silently
    # have zero effect, same silent-no-op failure mode this whole fix pass
    # exists to close off.
    case_dir = tmp_path / "case"
    with pytest.raises(ValueError, match="dx"):
        provision_mesh(case_dir=case_dir, myocardium_solver="singleCellSolver", dx_m=0.0004)


def test_provision_mesh_copies_the_bundled_single_cell_polymesh(tmp_path):
    # The 1-cell polyMesh fixture ships with the plugin, so this also pins
    # that the packaged fixture path still resolves.
    case_dir = tmp_path / "case"
    needs_block_mesh = provision_mesh(
        case_dir=case_dir, myocardium_solver="singleCellSolver",
    )
    assert needs_block_mesh is False
    poly_mesh = case_dir / "constant" / "polyMesh"
    for name in ("points", "faces", "owner", "neighbour", "boundary"):
        assert (poly_mesh / name).is_file(), name


def test_provision_mesh_leaves_an_unknown_solver_to_the_caller(tmp_path):
    case_dir = tmp_path / "case"
    assert provision_mesh(case_dir=case_dir, myocardium_solver="futureSolver") is False
    assert not (case_dir / "system" / "blockMeshDict").exists()
    assert not (case_dir / "constant" / "polyMesh").exists()
