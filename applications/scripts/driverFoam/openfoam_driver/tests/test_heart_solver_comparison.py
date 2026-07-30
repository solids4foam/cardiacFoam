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
#     test_heart_solver_comparison
#
# Description
#     Tests heart_solver_comparison.make_spec: the solver-variant catalog,
#     apply_case's template-file swap, and workflow_dag construction
#     (including run_in_parallel via the shared solve_steps helper).
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

from pathlib import Path

import pytest

from openfoam_driver.specs.tutorials.heart_solver_comparison import (
    HEART_SOLVER_VARIANTS,
    make_spec,
)

_TEMPLATE_FILES = ("electroProperties", "fvSchemes", "fvSolution", "controlDict")


def _write_case(tmp_path: Path) -> Path:
    case_root = tmp_path / "heartSim3D-1D" / "eikonalHeart"
    (case_root / "constant").mkdir(parents=True)
    (case_root / "system").mkdir(parents=True)
    (case_root / "system" / "decomposeParDict").write_text(
        "FoamFile\n{\n}\nnumberOfSubdomains 6;\n"
    )
    for variant in HEART_SOLVER_VARIANTS:
        template_dir = case_root / "setup" / "solverVariants" / variant
        template_dir.mkdir(parents=True)
        for name in _TEMPLATE_FILES:
            (template_dir / name).write_text(f"// {variant} {name}\n")
    return case_root


def test_unknown_variant_is_rejected(tmp_path):
    _write_case(tmp_path)
    with pytest.raises(ValueError, match="solver_variant"):
        make_spec(
            tutorials_root=tmp_path,
            case_dir_name="heartSim3D-1D/eikonalHeart",
            solver_variant="not-a-real-variant",
        )


@pytest.mark.parametrize("variant", HEART_SOLVER_VARIANTS)
def test_apply_case_copies_the_matching_template(tmp_path, variant):
    case_root = _write_case(tmp_path)
    spec = make_spec(
        tutorials_root=tmp_path,
        case_dir_name="heartSim3D-1D/eikonalHeart",
        solver_variant=variant,
    )
    cases = spec.build_cases()
    assert len(cases) == 1
    assert cases[0].case_id == variant
    spec.apply_case(spec.case_root, cases[0])

    assert (case_root / "constant" / "electroProperties").read_text() == f"// {variant} electroProperties\n"
    assert (case_root / "system" / "fvSchemes").read_text() == f"// {variant} fvSchemes\n"
    assert (case_root / "system" / "fvSolution").read_text() == f"// {variant} fvSolution\n"
    assert (case_root / "system" / "controlDict").read_text() == f"// {variant} controlDict\n"


def test_serial_workflow_dag_is_a_single_solve_step(tmp_path):
    _write_case(tmp_path)
    spec = make_spec(
        tutorials_root=tmp_path,
        case_dir_name="heartSim3D-1D/eikonalHeart",
        solver_variant="eikonal",
        run_in_parallel=False,
    )
    steps = spec.metadata["workflow_dag"]["steps"]
    assert steps == [{"id": "solve", "command": "cardiacFoam", "depends_on": []}]


def test_parallel_workflow_dag_wraps_decompose_reconstruct(tmp_path):
    _write_case(tmp_path)
    spec = make_spec(
        tutorials_root=tmp_path,
        case_dir_name="heartSim3D-1D/eikonalHeart",
        solver_variant="bidomain",
        run_in_parallel=True,
    )
    steps = spec.metadata["workflow_dag"]["steps"]
    assert [s["id"] for s in steps] == ["decomposePar", "solve", "reconstructPar"]
    by_id = {s["id"]: s for s in steps}
    assert by_id["solve"]["command"] == "mpirun"
    assert by_id["solve"]["args"] == ["-np", "6", "cardiacFoam", "-parallel"]


def test_metadata_records_the_solver_variant(tmp_path):
    _write_case(tmp_path)
    spec = make_spec(
        tutorials_root=tmp_path,
        case_dir_name="heartSim3D-1D/eikonalHeart",
        solver_variant="monodomain-eikonal",
    )
    assert spec.metadata["solver_variant"] == "monodomain-eikonal"
