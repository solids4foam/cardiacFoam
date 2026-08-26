import inspect
from pathlib import Path
from unittest import mock

import pytest

from openfoam_driver.plugins.cardiacfoam.tutorials.manufactured_eikonal_ecg import make_spec


_ELECTRO_PROPERTIES = """FoamFile
{
    version 2.0;
    format ascii;
    class dictionary;
    object electroProperties;
}
myocardiumSolver eikonalSolver;
eikonalSolverCoeffs
{
    conductivity [ -1 -3 3 0 0 2 0 ] (1 0 0 1 0 1);
    eikonalAdvectionDiffusionApproach true;
    verificationModel
    {
        type manufacturedEikonalVerifier;
    }
    ecgDomains
    {
        ECG
        {
            ecgSolver none;
            verificationModel
            {
                type manufacturedEikonalECGVerifier;
                enabled false;
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

_FVS = """FoamFile
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

_FVSOL = """FoamFile
{
    version 2.0;
    format ascii;
    class dictionary;
    object fvSolution;
}
solvers
{
}
"""


def _write_case(tutorials_root: Path) -> Path:
    case_root = tutorials_root / "manufacturedSolutions" / "eikonalECG"
    (case_root / "constant").mkdir(parents=True, exist_ok=True)
    (case_root / "system").mkdir(parents=True, exist_ok=True)
    (case_root / "setup" / "studies" / "tetConvergence").mkdir(parents=True, exist_ok=True)

    (case_root / "constant" / "electroProperties").write_text(_ELECTRO_PROPERTIES)
    (case_root / "constant" / "physicsProperties").write_text(_PHYSICS_PROPERTIES)
    (case_root / "system" / "blockMeshDict.3D").write_text("blocks ((10 10 10));\n")
    (case_root / "system" / "fvSchemes").write_text(_FVS)
    (case_root / "system" / "fvSolution").write_text(_FVSOL)
    (case_root / "setup" / "studies" / "tetConvergence" / "box.geo.template").write_text(
        "lc = __LC__;\nBox(1) = {0, 0, 0, 1, 1, 1};\n"
    )
    (case_root / "setup" / "studies" / "tetConvergence" / "fvSolution").write_text(
        _FVSOL.replace("solvers", "tetSolvers")
    )
    return case_root


def _make_spec(tmp_path, **overrides):
    _write_case(tmp_path)
    kwargs = {
        "tutorials_root": tmp_path,
        "case_dir_name": "manufacturedSolutions/eikonalECG",
        "dimensions": ["3D"],
        "number_cells": [10],
        # These tests are about mesh_family/gmsh-pipeline branching, not
        # parallel execution -- keep them decoupled from run_in_parallel's
        # decomposeParDict requirement (see parallel_execution.py).
        "run_in_parallel": False,
        **overrides,
    }
    return make_spec(**kwargs)


def test_eikonal_factory_does_not_advertise_a_nonexistent_convergence_axis():
    assert "convergence_axis" not in inspect.signature(make_spec).parameters


def test_tet_workflow_dag_matches_gmsh_pipeline(tmp_path):
    spec = _make_spec(tmp_path, mesh_family="tet")
    steps = spec.metadata["workflow_dag"]["steps"]
    assert [step["id"] for step in steps] == [
        "clean",
        "gmsh",
        "gmshToFoam",
        "checkMesh",
        "solve",
    ]
    assert [step["command"] for step in steps] == [
        "Allclean",
        "gmsh",
        "gmshToFoam",
        "checkMesh",
        "cardiacFoam",
    ]
    assert steps[1]["args"] == [
        "-3",
        "setup/studies/tetConvergence/box.geo",
        "-o",
        "box.msh",
        "-format",
        "msh2",
    ]


def test_tet_workflow_dag_appends_gradient_reconstruction_after_solve(tmp_path):
    spec = _make_spec(tmp_path, mesh_family="tet", gradient_reconstruction=True)
    steps = spec.metadata["workflow_dag"]["steps"]
    assert [step["id"] for step in steps] == [
        "clean",
        "gmsh",
        "gmshToFoam",
        "checkMesh",
        "solve",
        "gradientReconstructionOrder",
    ]
    assert steps[-1]["command"] == "gradientReconstructionOrder"
    assert steps[-1]["depends_on"] == ["solve"]


def test_tet_workflow_dag_appends_write_cell_centres_after_solve(tmp_path):
    spec = _make_spec(tmp_path, mesh_family="tet", error_localisation_analysis=True)
    steps = spec.metadata["workflow_dag"]["steps"]
    assert [step["id"] for step in steps] == [
        "clean",
        "gmsh",
        "gmshToFoam",
        "checkMesh",
        "solve",
        "writeCellCentres",
    ]
    assert steps[-1]["command"] == "postProcess"
    assert steps[-1]["args"] == ["-func", "writeCellCentres", "-latestTime"]
    assert steps[-1]["depends_on"] == ["solve"]


def test_tet_apply_case_renders_geo_installs_overlay_and_grad_scheme(tmp_path):
    case_root = _write_case(tmp_path)
    spec = make_spec(
        tutorials_root=tmp_path,
        case_dir_name="manufacturedSolutions/eikonalECG",
        dimensions=["3D"],
        number_cells=[10],
        mesh_family="tet",
        numerics_profile="eikonal_tet",
        grad_scheme="least_squares",
        run_in_parallel=False,
    )

    with mock.patch("openfoam_driver.plugins.cardiacfoam.tutorials.manufactured_eikonal_ecg.subprocess") as mock_subprocess:
        spec.apply_case(spec.case_root, spec.build_cases()[0])

    mock_subprocess.run.assert_not_called()
    assert "lc = 0.1;" in (case_root / "setup" / "studies" / "tetConvergence" / "box.geo").read_text()
    assert "tetSolvers" in (case_root / "system" / "fvSolution").read_text()
    assert "leastSquares;" in (case_root / "system" / "fvSchemes").read_text()


def test_tet_apply_case_forwards_conductivity_and_advection_approach(tmp_path):
    case_root = _write_case(tmp_path)
    conductivity = "[ -1 -3 3 0 0 2 0 ] (0.111 0 0 0.122 0 0.030)"
    spec = make_spec(
        tutorials_root=tmp_path,
        case_dir_name="manufacturedSolutions/eikonalECG",
        dimensions=["3D"],
        number_cells=[10],
        mesh_family="tet",
        run_in_parallel=False,
        conductivity=conductivity,
        eikonal_advection_diffusion_approach="false",
    )

    spec.apply_case(spec.case_root, spec.build_cases()[0])

    properties = (case_root / "constant" / "electroProperties").read_text()
    assert conductivity in properties
    assert "eikonalAdvectionDiffusionApproach    false;" in properties


def test_tet_validation_rejects_invalid_options(tmp_path):
    with pytest.raises(ValueError, match="3D"):
        _make_spec(tmp_path, mesh_family="tet", dimensions=["2D"])
    with pytest.raises(ValueError, match="numerics_profile"):
        _make_spec(tmp_path, mesh_family="tet", numerics_profile="bidomain_tet")
    with pytest.raises(ValueError, match="grad_scheme"):
        _make_spec(tmp_path, mesh_family="tet", grad_scheme="GaussLinear")
