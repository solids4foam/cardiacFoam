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
#     test_manufactured_bidomain_tet
#
# Description
#     Confirms manufactured_bidomain.make_spec's wrapper threads the new
#     tet-mesh logic through the base
#     manufactured_monodomain_pseudo_ecg.make_spec (pass-through only, no new logic here).
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

from pathlib import Path
from unittest import mock

from openfoam_driver.plugins.cardiacfoam.tutorials.manufactured_bidomain import make_spec


def test_mesh_family_tet_reaches_workflow_dag(tmp_path):
    case_root = tmp_path / "manufacturedSolutions" / "bidomain"
    (case_root / "constant").mkdir(parents=True)
    (case_root / "system").mkdir(parents=True)

    spec = make_spec(
        tutorials_root=tmp_path,
        case_dir_name="manufacturedSolutions/bidomain",
        dimensions=["3D"],
        number_cells=[10],
        dt_values=[0.00892857],
        mesh_family="tet",
        run_in_parallel=False,
    )
    commands = [s["command"] for s in spec.metadata["workflow_dag"]["steps"]]
    assert commands == ["Allclean", "gmsh", "gmshToFoam", "checkMesh", "cardiacFoam"]


def test_hex_is_still_the_default(tmp_path):
    (tmp_path / "manufacturedSolutions" / "bidomain" / "constant").mkdir(parents=True)
    (tmp_path / "manufacturedSolutions" / "bidomain" / "system").mkdir(parents=True)

    spec = make_spec(
        tutorials_root=tmp_path,
        case_dir_name="manufacturedSolutions/bidomain",
        dimensions=["3D"],
        number_cells=[10],
        dt_values=[0.00892857],
        run_in_parallel=False,
    )
    commands = [s["command"] for s in spec.metadata["workflow_dag"]["steps"]]
    assert commands == ["blockMesh", "cardiacFoam"]


def test_convergence_axis_is_forwarded_to_the_shared_factory(tmp_path):
    (tmp_path / "manufacturedSolutions" / "bidomain" / "constant").mkdir(parents=True)
    (tmp_path / "manufacturedSolutions" / "bidomain" / "system").mkdir(parents=True)

    spec = make_spec(
        tutorials_root=tmp_path,
        case_dir_name="manufacturedSolutions/bidomain",
        dimensions=["1D"],
        number_cells=[640],
        dt_values=[0.001],
        convergence_axis="temporal",
        run_in_parallel=False,
    )

    assert spec.metadata["convergence_axis"] == "temporal"


def test_ode_tolerances_are_forwarded_to_the_shared_factory():
    with mock.patch(
        "openfoam_driver.plugins.cardiacfoam.tutorials.manufactured_bidomain.make_base_spec"
    ) as make_base_spec:
        make_spec(ode_abs_tolerance=1e-10, ode_rel_tolerance=1e-8)

    assert make_base_spec.call_args.kwargs["ode_abs_tolerance"] == 1e-10
    assert make_base_spec.call_args.kwargs["ode_rel_tolerance"] == 1e-8
