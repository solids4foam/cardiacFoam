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
#     test_manufactured_monodomain_pseudo_ecg_tet
#
# Description
#     Tests mesh_family="tet" support in manufactured_monodomain_pseudo_ecg.py: workflow_dag
#     branching, apply_case's render-only materialization (never invokes
#     gmsh), numerics_profile overlay selection, grad_scheme/phi_tolerance/
#     end_time application, and strict kwarg validation.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

from pathlib import Path
from unittest import mock

import pytest

from openfoam_driver.core.runtime.models import CaseConfig
from openfoam_driver.tests.conftest import assert_foam_entry
from openfoam_driver.plugins.cardiacfoam.tutorials.manufactured_monodomain_pseudo_ecg import (
    make_spec,
)

_CONTROL_DICT = """FoamFile
{
    version 2.0;
    format ascii;
    class dictionary;
    object controlDict;
}
application cardiacFoam;
startTime 0;
stopAt endTime;
endTime 0.055;
deltaT 1e-06;
"""

_BLOCK_MESH_DICT = """FoamFile
{
    version 2.0;
    format ascii;
    class dictionary;
    object blockMeshDict;
}
scale 1;
vertices
(
    (0 0 0) (1 0 0) (1 1 0) (0 1 0)
    (0 0 1) (1 0 1) (1 1 1) (0 1 1)
);
blocks
(
    hex (0 1 2 3 4 5 6 7) (10 10 10) simpleGrading (1 1 1)
);
"""

_HEX_FV_SCHEMES = """FoamFile
{
    version 2.0;
    format ascii;
    class dictionary;
    object fvSchemes;
}
gradSchemes
{
    default Gauss linear;
}
"""

_TET_FV_SCHEMES = """FoamFile
{
    version 2.0;
    format ascii;
    class dictionary;
    object fvSchemes;
}
gradSchemes
{
    default leastSquares;
}
"""

_HEX_FV_SOLUTION = """FoamFile
{
    version 2.0;
    format ascii;
    class dictionary;
    object fvSolution;
}
solvers
{
    "phiE|phiEFinal|phiI|phiIFinal"
    {
        solver PCG;
        tolerance 1e-15;
    }
}
"""

_TET_FV_SOLUTION = """FoamFile
{
    version 2.0;
    format ascii;
    class dictionary;
    object fvSolution;
}
solvers
{
    "phiE|phiEFinal|phiI|phiIFinal"
    {
        solver PCG;
        tolerance 1e-15;
    }
}
"""

_ELECTRO_PROPERTIES = """FoamFile
{
    version 2.0;
    format ascii;
    class dictionary;
    object electroProperties;
}
myocardiumSolver monodomainSolver;
monodomainSolverCoeffs
{
    ionicModel monodomainFDAManufactured;
    conductivity [-1 -3 3 0 0 2 0] (0.111453302 0 0 0.121585420 0 0.030396355);
    dimension "3D";
    solutionAlgorithm implicit;
    verificationModel
    {
        type manufacturedFDAMonodomainVerifier;
    }
    ecgDomains
    {
        ECG
        {
            ecgSolver none;
            verificationModel
            {
                type manufacturedPseudoECGVerifier;
                enabled false;
                dimension "3D";
                referenceQuadratureOrder 1;
                checkQuadratureOrders (1);
            }
            electrodePositions
            {
                E1 (0 0 0);
                E2 (0 0 0);
                E3 (0 0 0);
                E4 (0 0 0);
                E5 (0 0 0);
            }
        }
    }
}
"""

_PHYSICS_PROPERTIES = """FoamFile
{
    version 2.0;
    format ascii;
    class dictionary;
    object physicsProperties;
}
type electroModel;
"""

_GEO_TEMPLATE = "lc = __LC__;\nBox(1) = {0, 0, 0, 1, 1, 1};\n"

_DECOMPOSE_PAR_DICT = """FoamFile
{
    version 2.0;
    format ascii;
    class dictionary;
    object decomposeParDict;
}
numberOfSubdomains 6;
method scotch;
"""


def _write_case(tutorials_root: Path, *, case_dir_name: str = "manufacturedSolutions/bidomain") -> Path:
    case_root = tutorials_root / case_dir_name
    (case_root / "constant").mkdir(parents=True)
    (case_root / "system").mkdir(parents=True)
    (case_root / "setup" / "studies" / "tetConvergence").mkdir(parents=True)

    (case_root / "constant" / "electroProperties").write_text(_ELECTRO_PROPERTIES)
    (case_root / "constant" / "physicsProperties").write_text(_PHYSICS_PROPERTIES)
    (case_root / "system" / "controlDict").write_text(_CONTROL_DICT)
    (case_root / "system" / "blockMeshDict.3D").write_text(_BLOCK_MESH_DICT)
    (case_root / "system" / "fvSchemes").write_text(_HEX_FV_SCHEMES)
    (case_root / "system" / "fvSolution").write_text(_HEX_FV_SOLUTION)
    (case_root / "system" / "decomposeParDict").write_text(_DECOMPOSE_PAR_DICT)
    (case_root / "setup" / "studies" / "tetConvergence" / "box.geo.template").write_text(_GEO_TEMPLATE)
    (case_root / "setup" / "studies" / "tetConvergence" / "fvSchemes").write_text(_TET_FV_SCHEMES)
    (case_root / "setup" / "studies" / "tetConvergence" / "fvSolution").write_text(_TET_FV_SOLUTION)
    return case_root


def _call_make_spec(tmp_path, **overrides):
    kwargs = {
        "tutorials_root": tmp_path,
        "case_dir_name": "manufacturedSolutions/bidomain",
        "dimensions": ["3D"],
        "number_cells": [10],
        "dt_values": [0.00892857],
        "ecg_enabled": False,
        # Defaults to False here so mesh_family-branching tests aren't also
        # coupled to run_in_parallel behavior; the run_in_parallel tests
        # override this explicitly.
        "run_in_parallel": False,
        **overrides,
    }
    return make_spec(**kwargs)


def _make_case(tmp_path, **overrides):
    """Write the fixture case and call make_spec against it, for tests that
    only need the resulting spec (not case_root itself)."""
    _write_case(tmp_path)
    return _call_make_spec(tmp_path, **overrides)


# --- Validation -------------------------------------------------------------

def test_invalid_mesh_family_is_rejected(tmp_path):
    with pytest.raises(ValueError, match="mesh_family"):
        _make_case(tmp_path, mesh_family="quad")


def test_tet_requires_3d_only(tmp_path):
    with pytest.raises(ValueError, match="3D"):
        _make_case(tmp_path, mesh_family="tet", dimensions=["1D"])


def test_hex_mesh_family_is_default(tmp_path):
    spec = _make_case(tmp_path)
    steps = spec.metadata["workflow_dag"]["steps"]
    assert [s["id"] for s in steps] == ["mesh", "solve"]


# --- workflow_dag branching --------------------------------------------------

def test_tet_workflow_dag_has_gmsh_pipeline_not_blockmesh(tmp_path):
    spec = _make_case(tmp_path, mesh_family="tet")
    steps = spec.metadata["workflow_dag"]["steps"]
    ids = [s["id"] for s in steps]
    commands = [s["command"] for s in steps]
    assert "blockMesh" not in commands
    assert ids == ["clean", "gmsh", "gmshToFoam", "checkMesh", "solve"]
    assert commands == ["Allclean", "gmsh", "gmshToFoam", "checkMesh", "cardiacFoam"]


def test_tet_workflow_dag_gmsh_steps_have_correct_args(tmp_path):
    # Found via a real (non-mocked) sweep-run: gmsh with no args launches its
    # GUI and hangs indefinitely instead of meshing anything. Each step needs
    # the same args the real bash script passes: gmsh -3 <geo> -o box.msh
    # -format msh2; gmshToFoam box.msh. checkMesh/Allclean/cardiacFoam take
    # none (matching the hex "solve" step's own no-args precedent).
    spec = _make_case(tmp_path, mesh_family="tet")
    by_id = {s["id"]: s for s in spec.metadata["workflow_dag"]["steps"]}
    assert by_id["gmsh"]["args"] == [
        "-3", "setup/studies/tetConvergence/box.geo", "-o", "box.msh", "-format", "msh2",
    ]
    assert by_id["gmshToFoam"]["args"] == ["box.msh"]


def test_tet_workflow_dag_step_dependencies_are_sequential(tmp_path):
    spec = _make_case(tmp_path, mesh_family="tet")
    by_id = {s["id"]: s for s in spec.metadata["workflow_dag"]["steps"]}
    assert by_id["clean"]["depends_on"] == []
    assert by_id["gmsh"]["depends_on"] == ["clean"]
    assert by_id["gmshToFoam"]["depends_on"] == ["gmsh"]
    assert by_id["checkMesh"]["depends_on"] == ["gmshToFoam"]
    assert by_id["solve"]["depends_on"] == ["checkMesh"]


# --- run_in_parallel: decomposePar/mpirun/reconstructPar wrapping -----------

def test_hex_run_in_parallel_wraps_solve_with_decompose_reconstruct(tmp_path):
    spec = _make_case(tmp_path, run_in_parallel=True)
    steps = spec.metadata["workflow_dag"]["steps"]
    ids = [s["id"] for s in steps]
    assert ids == ["mesh", "decomposePar", "solve", "reconstructPar"]
    by_id = {s["id"]: s for s in steps}
    assert by_id["decomposePar"]["command"] == "decomposePar"
    assert by_id["decomposePar"]["depends_on"] == ["mesh"]
    assert by_id["solve"]["command"] == "mpirun"
    assert by_id["solve"]["args"] == ["-np", "6", "cardiacFoam", "-parallel"]
    assert by_id["solve"]["depends_on"] == ["decomposePar"]
    assert by_id["reconstructPar"]["command"] == "reconstructPar"
    assert by_id["reconstructPar"]["depends_on"] == ["solve"]


def test_hex_run_in_parallel_false_is_unaffected(tmp_path):
    spec = _make_case(tmp_path, run_in_parallel=False)
    steps = spec.metadata["workflow_dag"]["steps"]
    assert [s["id"] for s in steps] == ["mesh", "solve"]
    by_id = {s["id"]: s for s in steps}
    assert by_id["solve"]["command"] == "cardiacFoam"


def test_tet_run_in_parallel_wraps_solve_after_check_mesh(tmp_path):
    spec = _make_case(tmp_path, mesh_family="tet", run_in_parallel=True)
    steps = spec.metadata["workflow_dag"]["steps"]
    ids = [s["id"] for s in steps]
    assert ids == [
        "clean", "gmsh", "gmshToFoam", "checkMesh",
        "decomposePar", "solve", "reconstructPar",
    ]
    by_id = {s["id"]: s for s in steps}
    assert by_id["decomposePar"]["depends_on"] == ["checkMesh"]
    assert by_id["solve"]["command"] == "mpirun"
    assert by_id["solve"]["args"] == ["-np", "6", "cardiacFoam", "-parallel"]


# --- apply_case: render-only, never executes gmsh ---------------------------

def test_apply_case_renders_geo_but_never_calls_gmsh(tmp_path):
    case_root = _write_case(tmp_path)
    spec = _call_make_spec(tmp_path, mesh_family="tet", numerics_profile="bidomain_tet")
    cases = spec.build_cases()
    assert len(cases) == 1

    with mock.patch("openfoam_driver.plugins.cardiacfoam.tutorials.manufactured_monodomain_pseudo_ecg.subprocess") as mock_subprocess:
        spec.apply_case(spec.case_root, cases[0])

    mock_subprocess.run.assert_not_called()
    geo_text = (case_root / "setup" / "studies" / "tetConvergence" / "box.geo").read_text()
    assert "__LC__" not in geo_text
    assert "lc = 0.1;" in geo_text


def test_apply_case_can_select_an_optimised_geo_template(tmp_path):
    case_root = _write_case(tmp_path)
    optimised_template = case_root / "setup" / "studies" / "tetConvergence" / "box.geo.template.optimised"
    optimised_template.write_text(
        _GEO_TEMPLATE + "Mesh.Algorithm3D = 4;\nMesh.OptimizeNetgen = 1;\n"
    )
    spec = _call_make_spec(
        tmp_path,
        mesh_family="tet",
        numerics_profile="bidomain_tet",
        tet_geo_template_relpath="setup/studies/tetConvergence/box.geo.template.optimised",
    )

    spec.apply_case(spec.case_root, spec.build_cases()[0])

    geo_text = (case_root / "setup" / "studies" / "tetConvergence" / "box.geo").read_text()
    assert "lc = 0.1;" in geo_text
    assert "Mesh.Algorithm3D = 4;" in geo_text
    assert spec.metadata["tet_geo_template_relpath"].endswith("box.geo.template.optimised")





def test_conductivity_shorthand_updates_monodomain_tensor(tmp_path):
    case_root = _write_case(tmp_path)
    tensor = "[-1 -3 3 0 0 2 0] (0.1 -0.01 -0.02 0.12 -0.013 0.04)"
    spec = _call_make_spec(tmp_path, conductivity=tensor)

    spec.apply_case(spec.case_root, spec.build_cases()[0])

    assert_foam_entry(
        case_root / "constant" / "electroProperties",
        "conductivity",
        tensor,
        scope="monodomainSolverCoeffs",
    )


def test_apply_case_hex_does_not_render_geo(tmp_path):
    case_root = _write_case(tmp_path)
    spec = _call_make_spec(tmp_path)
    cases = spec.build_cases()
    spec.apply_case(spec.case_root, cases[0])
    assert not (case_root / "setup" / "mesh" / "tet" / "box.geo").exists()


# --- numerics_profile: case-specific overlay sets ---------------------------

def test_bidomain_tet_profile_installs_fvschemes_only(tmp_path):
    case_root = _write_case(tmp_path)
    spec = _call_make_spec(tmp_path, mesh_family="tet", numerics_profile="bidomain_tet")
    cases = spec.build_cases()
    spec.apply_case(spec.case_root, cases[0])
    assert "leastSquares" in (case_root / "system" / "fvSchemes").read_text()


def test_monodomain_tet_profile_installs_fvschemes_and_fvsolution(tmp_path):
    case_root = _write_case(tmp_path)
    spec = _call_make_spec(tmp_path, mesh_family="tet", numerics_profile="monodomain_tet")
    cases = spec.build_cases()
    spec.apply_case(spec.case_root, cases[0])
    assert "leastSquares" in (case_root / "system" / "fvSchemes").read_text()
    assert (case_root / "system" / "fvSolution").exists()
    assert "1e-15" in (case_root / "system" / "fvSolution").read_text()


# --- grad_scheme: validated enum, correct OpenFOAM token --------------------

def test_grad_scheme_gauss_linear_writes_literal_openfoam_token(tmp_path):
    case_root = _write_case(tmp_path)
    spec = _call_make_spec(
        tmp_path, mesh_family="tet", numerics_profile="bidomain_tet", grad_scheme="gauss_linear",
    )
    cases = spec.build_cases()
    spec.apply_case(spec.case_root, cases[0])
    text = (case_root / "system" / "fvSchemes").read_text()
    assert "Gauss linear;" in text
    assert "GaussLinear" not in text


def test_grad_scheme_least_squares_writes_literal_openfoam_token(tmp_path):
    case_root = _write_case(tmp_path)
    spec = _call_make_spec(
        tmp_path, mesh_family="tet", numerics_profile="bidomain_tet", grad_scheme="least_squares",
    )
    cases = spec.build_cases()
    spec.apply_case(spec.case_root, cases[0])
    text = (case_root / "system" / "fvSchemes").read_text()
    assert "leastSquares;" in text


def test_unknown_grad_scheme_name_is_rejected(tmp_path):
    with pytest.raises(ValueError, match="grad_scheme"):
        _make_case(tmp_path, mesh_family="tet", numerics_profile="bidomain_tet", grad_scheme="GaussLinear")


# --- phi_tolerance -----------------------------------------------------------

def test_phi_tolerance_changes_only_phie_phii_block(tmp_path):
    case_root = _write_case(tmp_path)
    spec = _call_make_spec(
        tmp_path, mesh_family="tet", numerics_profile="monodomain_tet", phi_tolerance=1e-6,
    )
    cases = spec.build_cases()
    spec.apply_case(spec.case_root, cases[0])
    assert_foam_entry(
        case_root / "system" / "fvSolution",
        "tolerance",
        "1e-6",
        scope=("solvers", '"phiE|phiEFinal|phiI|phiIFinal"'),
    )


def test_non_positive_phi_tolerance_is_rejected(tmp_path):
    _write_case(tmp_path)
    with pytest.raises(ValueError, match="phi_tolerance"):
        _call_make_spec(tmp_path, mesh_family="tet", numerics_profile="monodomain_tet", phi_tolerance=0)
    with pytest.raises(ValueError, match="phi_tolerance"):
        _call_make_spec(tmp_path, mesh_family="tet", numerics_profile="monodomain_tet", phi_tolerance=-1e-6)


# --- ECG + tet interaction ---------------------------------------------------

def test_tet_mesh_family_works_with_ecg_enabled(tmp_path):
    # monodomainPseudoECG's real tet study runs with ECG active (it reports
    # ecg_L2/ecg_Linf alongside mono_L2/mono_Linf for every case, confirmed
    # from run_mono_tet.sh) -- ECG_ENABLED defaults to True in this shared
    # module, so mesh_family="tet" must not implicitly assume ECG is off.
    # ECG overrides touch electroProperties; mesh_family/grad_scheme/
    # phi_tolerance touch fvSchemes/fvSolution/controlDict -- disjoint files,
    # but verified here rather than assumed.
    case_root = _write_case(tmp_path)
    spec = _call_make_spec(
        tmp_path, mesh_family="tet", numerics_profile="monodomain_tet",
        grad_scheme="least_squares", ecg_enabled=True,
    )
    cases = spec.build_cases()
    spec.apply_case(spec.case_root, cases[0])

    electro = case_root / "constant" / "electroProperties"
    ecg = ("monodomainSolverCoeffs", "ecgDomains", "ECG")
    assert_foam_entry(electro, "ecgSolver", "pseudoECG", scope=ecg)
    assert_foam_entry(electro, "enabled", "yes", scope=ecg + ("verificationModel",))
    assert_foam_entry(
        electro, "E3", "(1.2 0.23 0.61)", scope=ecg + ("electrodePositions",)
    )
    scheme_text = (case_root / "system" / "fvSchemes").read_text()
    assert "leastSquares;" in scheme_text


def test_ecg_disabled_removes_block_from_a_reused_entry_case(tmp_path):
    case_root = _write_case(tmp_path)
    enabled = _call_make_spec(tmp_path, ecg_enabled=True)
    enabled.apply_case(enabled.case_root, enabled.build_cases()[0])

    disabled = _call_make_spec(tmp_path, ecg_enabled=False)
    disabled.apply_case(disabled.case_root, disabled.build_cases()[0])

    electro_text = (case_root / "constant" / "electroProperties").read_text()
    assert "ecgDomains" not in electro_text
    assert "ecgSolver" not in electro_text


# --- end_time -------------------------------------------------------------

def test_end_time_overrides_control_dict(tmp_path):
    case_root = _write_case(tmp_path)
    spec = _call_make_spec(tmp_path, mesh_family="tet", numerics_profile="bidomain_tet", end_time=0.2)
    cases = spec.build_cases()
    spec.apply_case(spec.case_root, cases[0])
    assert_foam_entry(case_root / "system" / "controlDict", "endTime", "0.2")
