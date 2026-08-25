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
#     test_sweep_materialize
#
# Description
#     Tests sweep case materialization via build_and_launch.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

import pytest

from openfoam_driver.sweep_materialize import materialize_case


def test_materialize_case_writes_dict_files_and_allrun_only(tmp_path):
    case_dir = tmp_path / "TNNP_1e-06"
    materialize_case(
        case_dir=case_dir,
        routed={
            "electro_selectors": {"myocardiumSolver": "singleCellSolver", "tissue": "epicardialCells", "ionicModel": "TNNP"},
            "physics_selectors": {"type": "electroModel"},
            "electro_overrides": {},
            "physics_overrides": {},
            "delta_t": 1e-6,
            "end_time": None,
        },
    )
    assert (case_dir / "constant" / "electroProperties").exists()
    assert (case_dir / "constant" / "physicsProperties").exists()
    assert (case_dir / "system" / "controlDict").exists()
    assert "1e-06" in (case_dir / "system" / "controlDict").read_text()

    allrun = case_dir / "Allrun"
    assert allrun.exists()
    assert allrun.stat().st_mode & 0o111  # executable
    assert not (case_dir / "workflow_contract.json").exists()


def test_materialize_case_two_cases_do_not_collide(tmp_path):
    materialize_case(
        case_dir=tmp_path / "caseA",
        routed={"electro_selectors": {"myocardiumSolver": "singleCellSolver", "tissue": "epicardialCells", "ionicModel": "TNNP"},
                "physics_selectors": {"type": "electroModel"}, "electro_overrides": {}, "physics_overrides": {},
                "delta_t": None, "end_time": None},
    )
    materialize_case(
        case_dir=tmp_path / "caseB",
        routed={"electro_selectors": {"myocardiumSolver": "singleCellSolver", "tissue": "epicardialCells", "ionicModel": "BuenoOrovio"},
                "physics_selectors": {"type": "electroModel"}, "electro_overrides": {}, "physics_overrides": {},
                "delta_t": None, "end_time": None},
    )
    a = (tmp_path / "caseA" / "constant" / "electroProperties").read_text()
    b = (tmp_path / "caseB" / "constant" / "electroProperties").read_text()
    assert "TNNP" in a and "BuenoOrovio" not in a
    assert "BuenoOrovio" in b and "TNNP" not in b


def test_materialize_case_runs_block_mesh_first_for_spatial_solver(tmp_path):
    # monodomainSolver needs a real fvMesh; provision_mesh writes a default
    # blockMeshDict for it (see mesh_provisioning.py), so the generated
    # Allrun must run blockMesh before cardiacFoam -- otherwise cardiacFoam
    # crashes with "Cannot find file points in polyMesh" (the exact failure
    # this fix addresses, see project_driverfoam_sweep_bugs_found memory).
    case_dir = tmp_path / "TNNP_monodomain"
    materialize_case(
        case_dir=case_dir,
        routed={
            "electro_selectors": {"myocardiumSolver": "monodomainSolver", "tissue": "epicardialCells", "ionicModel": "TNNP"},
            "physics_selectors": {"type": "electroModel"},
            "electro_overrides": {}, "physics_overrides": {},
            "delta_t": None, "end_time": None,
        },
    )
    assert (case_dir / "system" / "blockMeshDict").exists()
    allrun_text = (case_dir / "Allrun").read_text()
    assert allrun_text.index("blockMesh") < allrun_text.index("cardiacFoam")


def test_materialize_case_honours_dx_for_spatial_solver(tmp_path):
    from openfoam_driver.core.specs.mesh_provisioning import default_block_mesh_dict_text

    case_dir = tmp_path / "TNNP_monodomain_fine"
    materialize_case(
        case_dir=case_dir,
        routed={
            "electro_selectors": {"myocardiumSolver": "monodomainSolver", "tissue": "epicardialCells", "ionicModel": "TNNP"},
            "physics_selectors": {"type": "electroModel"},
            "electro_overrides": {}, "physics_overrides": {},
            "delta_t": None, "end_time": None, "dx": 0.0004,
        },
    )
    written = (case_dir / "system" / "blockMeshDict").read_text()
    assert written == default_block_mesh_dict_text(dx_m=0.0004)


def test_materialize_case_raises_on_invalid_combination(tmp_path):
    with pytest.raises(ValueError, match="tissue"):
        materialize_case(
            case_dir=tmp_path / "bad_case",
            routed={"electro_selectors": {"myocardiumSolver": "singleCellSolver", "tissue": "myocyte", "ionicModel": "TNNP"},
                    "physics_selectors": {"type": "electroModel"}, "electro_overrides": {}, "physics_overrides": {},
                    "delta_t": None, "end_time": None},
        )
