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
#     test_mesh_provisioning
#
# Description
#     Tests default mesh provisioning for from-scratch case_folder cases.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

import re

import pytest

from openfoam_driver.specs.mesh_provisioning import (
    cell_counts_from_dx,
    default_block_mesh_dict_text,
    provision_mesh,
)


def _cell_counts(text: str) -> tuple[int, int, int]:
    match = re.search(r"hex \([^)]*\)\s*\((\d+)\s+(\d+)\s+(\d+)\)", text)
    assert match is not None, text
    return tuple(int(g) for g in match.groups())


def test_cell_counts_from_dx_divides_exactly():
    assert cell_counts_from_dx(0.001, (0.002, 0.002, 0.002)) == (2, 2, 2)
    assert cell_counts_from_dx(0.0004, (0.002, 0.002, 0.002)) == (5, 5, 5)


def test_cell_counts_from_dx_rejects_non_exact_division():
    # Deliberately no silent rounding: dx that doesn't fit the domain is a
    # caller error, not something to approximate quietly (see
    # project_driverfoam_sweep_bugs_found memory -- this mirrors
    # niederer_2012.py's own established rigor for the same problem).
    with pytest.raises(ValueError, match="does not evenly divide"):
        cell_counts_from_dx(0.0003, (0.002, 0.002, 0.002))


def test_cell_counts_from_dx_rejects_non_positive_dx():
    with pytest.raises(ValueError, match="dx must be positive"):
        cell_counts_from_dx(0.0, (0.002, 0.002, 0.002))


def test_default_block_mesh_dict_has_a_fixed_default_cell_count_with_no_dx():
    text = default_block_mesh_dict_text()
    assert _cell_counts(text) == (4, 4, 4)


def test_dx_controls_cell_count_finer_mesh_for_smaller_dx():
    coarse = default_block_mesh_dict_text(dx_m=0.001)
    fine = default_block_mesh_dict_text(dx_m=0.0004)
    coarse_cells = _cell_counts(coarse)
    fine_cells = _cell_counts(fine)
    assert coarse_cells == (2, 2, 2)
    assert fine_cells == (5, 5, 5)
    assert fine_cells[0] > coarse_cells[0]


def test_dx_that_does_not_evenly_divide_default_slab_raises():
    with pytest.raises(ValueError, match="does not evenly divide"):
        default_block_mesh_dict_text(dx_m=0.0003)


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
