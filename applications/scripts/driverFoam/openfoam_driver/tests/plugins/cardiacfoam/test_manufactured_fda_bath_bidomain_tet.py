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
#     test_manufactured_fda_bath_bidomain_tet
#
# Description
#     Tests bath-bidomain tet workflow_dag/materialization and the explicit
#     bath predictor/corrector setting.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

from pathlib import Path

import pytest

from openfoam_driver.plugins.cardiacfoam.tutorials.manufactured_fda_bath_bidomain import make_spec


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
blocks
(
    hex (0 1 2 3 4 5 6 7) (10 10 10) simpleGrading (1 1 1)
    hex (0 1 2 3 4 5 6 7) (10 10 10) simpleGrading (1 1 1)
    hex (0 1 2 3 4 5 6 7) (10 10 10) simpleGrading (1 1 1)
);
"""

_HEX_ELECTRO_PROPERTIES = """FoamFile
{
    version 2.0;
    format ascii;
    class dictionary;
    object electroProperties;
}
myocardiumSolver bidomainSolver;
bidomainSolverCoeffs
{
    cellZone myocardium;
    bathPredictorCorrector true;
    ionicModel bathBidomainFDAManufactured;
    verificationModel
    {
        type manufacturedFDABathBidomainVerifier;
        groundElectrode yes;
    }
    manufacturedBidomain
    {
        groundElectrode yes;
    }
    dimension "1D";
    solutionAlgorithm implicit;
    bathPotentialDomain
    {
        bathCellZones (bath);
        bathConductivityField bodyAndOrgansConductivity;
    }
}
"""

_TET_ELECTRO_PROPERTIES = _HEX_ELECTRO_PROPERTIES.replace(
    'dimension "1D";',
    'dimension "3D";',
).replace(
    "bathConductivityField bodyAndOrgansConductivity;",
    """bathConductivityField bodyAndOrgansConductivity;
        interfaceConductivityInterpolation distanceWeightedHarmonic;""",
)

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

_TET_FV_SCHEMES = _HEX_FV_SCHEMES.replace("default Gauss linear;", "default leastSquares;")

_FV_SOLUTION = """FoamFile
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

_PHYSICS_PROPERTIES = """FoamFile
{
    version 2.0;
    format ascii;
    class dictionary;
    object physicsProperties;
}
type electroModel;
"""

_GEO_TEMPLATE = "lc = __LC__; // __LC__ in comment should be ignored\n"


def _write_case(tutorials_root: Path) -> Path:
    case_root = tutorials_root / "manufacturedSolutions" / "bathBidomain"
    (case_root / "constant").mkdir(parents=True)
    (case_root / "system").mkdir(parents=True)
    (case_root / "setup" / "mesh" / "tet").mkdir(parents=True)

    (case_root / "constant" / "electroProperties").write_text(_HEX_ELECTRO_PROPERTIES)
    (case_root / "constant" / "physicsProperties").write_text(_PHYSICS_PROPERTIES)
    (case_root / "system" / "controlDict").write_text(_CONTROL_DICT)
    (case_root / "system" / "blockMeshDict.3D").write_text(_BLOCK_MESH_DICT)
    (case_root / "system" / "fvSchemes").write_text(_HEX_FV_SCHEMES)
    (case_root / "system" / "fvSolution").write_text(_FV_SOLUTION)
    (case_root / "setup" / "mesh" / "tet" / "electroProperties").write_text(
        _TET_ELECTRO_PROPERTIES
    )
    (case_root / "setup" / "mesh" / "tet" / "fvSchemes").write_text(_TET_FV_SCHEMES)
    (case_root / "setup" / "mesh" / "tet" / "three_domain_box.geo.template").write_text(
        _GEO_TEMPLATE
    )
    return case_root


def _make_spec(tmp_path, **overrides):
    _write_case(tmp_path)
    kwargs = {
        "tutorials_root": tmp_path,
        "case_dir_name": "manufacturedSolutions/bathBidomain",
        "dimensions": ["3D"],
        "number_cells": [10],
        "dt_values": [0.00892857],
        # These tests are about mesh_family/gmsh-pipeline branching and dict
        # overrides, not parallel execution -- keep them decoupled from
        # run_in_parallel's decomposeParDict requirement (parallel_execution.py).
        "run_in_parallel": False,
        **overrides,
    }
    return make_spec(**kwargs)


def test_tet_workflow_dag_uses_three_domain_gmsh_pipeline(tmp_path):
    spec = _make_spec(tmp_path, mesh_family="tet")
    steps = spec.metadata["workflow_dag"]["steps"]
    assert [s["id"] for s in steps] == [
        "clean", "gmsh", "gmshToFoam", "checkMesh", "setConductivity", "solve",
        "interfaceMetrics",
    ]
    assert [s["command"] for s in steps] == [
        "Allclean", "gmsh", "gmshToFoam", "checkMesh",
        "setTorsoOrganConductivityField", "cardiacFoam",
        "bathBidomainInterfaceMetrics",
    ]
    by_id = {s["id"]: s for s in steps}
    assert by_id["gmsh"]["args"] == [
        "-3", "setup/mesh/tet/three_domain_box.geo", "-o",
        "three_domain_box.msh", "-format", "msh2",
    ]
    assert by_id["gmshToFoam"]["args"] == ["three_domain_box.msh"]
    # Reads the reconstructed final-time solution, not something the live
    # (potentially parallel-decomposed) verifier can compute itself -- see
    # applications/utilities/bathBidomainInterfaceMetrics: it does
    # heart/bath fvMeshSubset + interface-face analysis on "a reconstructed
    # serial mesh", which only exists once the solve step has fully exited.
    assert by_id["interfaceMetrics"]["args"] == ["-latestTime"]
    assert by_id["interfaceMetrics"]["depends_on"] == ["solve"]


def test_tet_unclaimed_artifacts_are_credited_to_the_solve_step(tmp_path):
    """The last *artifact-producing solver* step owns unclaimed artifacts.

    ``bathBidomainInterfaceMetrics`` is authorized to run but is a
    post-processing utility, not a solver, so it must not be credited with the
    run's outputs -- ``produces`` is enforced post-step
    (``workflow_runner._run_step``'s ``missing_artifacts`` check), so crediting
    the metrics utility would make a silent solver blame ``interfaceMetrics``.
    """
    from openfoam_driver.core.plugin_interface import default_driver_context
    from openfoam_driver.core.runtime.models import DataArtifact
    from openfoam_driver.core.runtime.workflow import normalize_workflow_dag

    spec = _make_spec(tmp_path, mesh_family="tet")
    artifact = DataArtifact(
        artifact_id="vm_field",
        path_pattern="{time}/Vm",
        format="openfoam_field",
    )
    dag, diagnostics = normalize_workflow_dag(
        spec.metadata["workflow_dag"],
        expected_artifacts=(artifact,),
        driver_context=default_driver_context(),
    )
    assert [d for d in diagnostics if d.level == "error"] == []
    producers = {
        step["id"] for step in dag["steps"] if "vm_field" in step.get("produces", ())
    }
    assert producers == {"solve"}


def test_hex_workflow_dag_is_still_blockmesh_toposet_pipeline(tmp_path):
    spec = _make_spec(tmp_path)
    assert [s["command"] for s in spec.metadata["workflow_dag"]["steps"]] == [
        "Allclean", "blockMesh", "topoSet", "setTorsoOrganConductivityField", "cardiacFoam",
    ]


def test_tet_rejects_unsupported_options(tmp_path):
    _write_case(tmp_path)
    base = {
        "tutorials_root": tmp_path,
        "case_dir_name": "manufacturedSolutions/bathBidomain",
        "dimensions": ["3D"],
        "number_cells": [10],
        "dt_values": [0.00892857],
    }
    with pytest.raises(ValueError, match="mesh_family"):
        make_spec(**base, mesh_family="quad")
    with pytest.raises(ValueError, match="3D"):
        make_spec(**{**base, "dimensions": ["2D"]}, mesh_family="tet")
    with pytest.raises(ValueError, match="numerics_profile"):
        make_spec(**base, mesh_family="tet", numerics_profile="bidomain_tet")
    with pytest.raises(ValueError, match="grad_scheme"):
        make_spec(**base, mesh_family="tet", grad_scheme="GaussLinear")
    with pytest.raises(ValueError, match="phi_tolerance"):
        make_spec(**base, mesh_family="tet", phi_tolerance=0)
