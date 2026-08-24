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
#     test_template_contract
#
# Description
#     Tests template contract logic and specification contracts.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

from __future__ import annotations

import unittest
from pathlib import Path

from openfoam_driver.dict_entries import all_documented_driver_paths
from openfoam_driver.tests.conftest import monorepo_root, skip_without_monorepo


def _template_path() -> Path:
    """Resolve the canonical electroProperties template.

    Tier 1: monorepo layout — ``tutorials/template/constant/electroProperties``
        found somewhere in the ancestor tree.
    Tier 2: bundled package fixture — shipped inside
        ``openfoam_driver/specs/fixtures/template/constant/electroProperties``
        for standalone / CI installs.
    """
    current = Path(__file__).resolve()
    # Tier 1: look for the live monorepo template
    for parent in current.parents:
        candidate = parent / "tutorials" / "template" / "constant" / "electroProperties"
        if candidate.exists():
            return candidate
    # Tier 2: bundled fixture co-located with the package
    # __file__ = …/tests/plugins/cardiacfoam/test_template_contract.py
    # parents[0]=cardiacfoam/ parents[1]=plugins/ parents[2]=tests/ parents[3]=openfoam_driver/
    fixture = current.parents[3] / "specs" / "fixtures" / "template" / "constant" / "electroProperties"
    if fixture.exists():
        return fixture
    raise FileNotFoundError(
        "Cannot locate the electroProperties template. "
        "Expected either the monorepo tutorials/template/ tree or the bundled fixture at "
        f"{fixture}"
    )


def _read(path: Path) -> str:
    return path.read_text()


class TestTemplateAndSchemaContract(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.template_path = _template_path()

    def test_template_uses_code_backed_selector_keys(self) -> None:
        template = _read(self.template_path)

        self.assertIn("myocardiumSolver", template)
        self.assertIn("monodomainSolverCoeffs", template)
        self.assertIn("singleCellSolverCoeffs", template)
        self.assertIn("ecgSolver", template)
        self.assertIn("electrodePositions", template)
        self.assertIn("conductivityIntracellular", template)
        self.assertIn("conductivityExtracellular", template)
        self.assertIn("monodomain1DSolver", template)
        self.assertIn("phiERefPoint", template)
        self.assertIn("purkinjeGraphModelCoeffs", template)

        # Bath-coupled ECG support keys (canonical C++ key set; see spec §3.1).
        self.assertIn("bathPotentialDomain", template)
        self.assertIn("extracellularPotentialDomain", template)
        self.assertIn("bathCellZones", template)
        self.assertIn("bathConductivityField", template)
        self.assertIn("torsoECG", template)
        self.assertIn("groundPatches", template)
        self.assertIn("surfaceCurrentPatches", template)

        self.assertNotIn("ecgDomainModel", template)
        self.assertNotIn("pseudoECGElectroCoeffs", template)
        self.assertNotIn("monoDomainElectroCoeffs", template)
        self.assertNotIn("phiEReferenceCell", template)
        self.assertNotIn("purkinjeNetworkModelCoeffs", template)
        self.assertNotIn("bidomainBathECG", template)

    def test_template_truth_markers_match_core_cpp_contracts(self) -> None:
        import pytest
        repo_root = monorepo_root
        if repo_root is None or not (repo_root / "src").exists():
            pytest.skip(
                "Requires the monorepo src/ tree — skipping C++ contract cross-check in standalone install."
            )
        electro_model_core = _read(
            repo_root / "src" / "electroModels" / "core" / "electroModel.C"
        )
        bidomain_solver = _read(
            repo_root
            / "src"
            / "electroModels"
            / "myocardiumModels"
            / "bidomainSolver"
            / "bidomainSolver.C"
        )
        ecg_solver = _read(
            repo_root
            / "src"
            / "electroModels"
            / "electroDomains"
            / "ecgDomain"
            / "ecgSolver.C"
        )
        ecg_domain = _read(
            repo_root
            / "src"
            / "electroModels"
            / "electroDomains"
            / "ecgDomain"
            / "ecgDomain.C"
        )
        system_builder = _read(
            repo_root
            / "src"
            / "electroModels"
            / "core"
            / "system"
            / "electrophysicsSystemBuilder.C"
        )
        electro_coupler = _read(
            repo_root
            / "src"
            / "electroModels"
            / "electroCouplers"
            / "electroDomainCoupler.C"
        )
        conduction_domain_selector = _read(
            repo_root
            / "src"
            / "electroModels"
            / "electroDomains"
            / "conductionSystemDomain"
            / "conductionSystemDomain.C"
        )
        monodomain_runtime = _read(
            repo_root
            / "src"
            / "electroModels"
            / "myocardiumModels"
            / "monodomainSolver"
            / "monodomainSolver.H"
        )
        bidomain_runtime = _read(
            repo_root
            / "src"
            / "electroModels"
            / "myocardiumModels"
            / "bidomainSolver"
            / "bidomainSolver.H"
        )
        single_cell_runtime = _read(
            repo_root
            / "src"
            / "electroModels"
            / "myocardiumModels"
            / "singleCellSolver"
            / "singleCellSolver.H"
        )
        # In the new architecture, eikonalSolver is no longer a standalone class
        # under myocardiumModels/.  It is registered as a named alias of
        # electrophysiologyModel via addNamedToRunTimeSelectionTable in
        # electrophysiologyModel.C, which then dispatches to
        # EikonalMyocardiumDomain.  Read the registration file instead.
        eikonal_runtime = _read(
            repo_root
            / "src"
            / "electroModels"
            / "core"
            / "electrophysiologyModel"
            / "electrophysiologyModel.C"
        )
        conduction_runtime = _read(
            repo_root
            / "src"
            / "electroModels"
            / "conductionSystemModels"
            / "monodomain1DSolver"
            / "monodomain1DSolver.H"
        )
        pseudo_ecg_runtime = _read(
            repo_root
            / "src"
            / "electroModels"
            / "ecgModels"
            / "pseudoECGSolver"
            / "pseudoECGSolver.H"
        )

        self.assertIn('lookup("myocardiumSolver")', electro_model_core)
        self.assertIn('"conductivityIntracellular"', bidomain_solver)
        self.assertIn('"conductivityExtracellular"', bidomain_solver)
        self.assertIn('"phiERefPoint"', bidomain_solver)
        self.assertIn('"ecgSolver"', ecg_solver)
        self.assertIn('"electrodePositions"', ecg_domain)
        self.assertIn("torsoECG requires myocardiumSolver", system_builder)
        self.assertIn("bidomainSolver because it samples the", system_builder)
        self.assertIn('"electroDomainCoupler"', electro_coupler)
        self.assertIn('"purkinjeGraphModel"', conduction_domain_selector)
        self.assertIn('dict.get<word>("graphFile")', conduction_domain_selector)
        self.assertIn('dict.subDict("rootStimulus")', conduction_domain_selector)

        self.assertIn('OverrideTypeName("monodomainSolver")', monodomain_runtime)
        self.assertIn('OverrideTypeName("bidomainSolver")', bidomain_runtime)
        self.assertIn('OverrideTypeName("singleCellSolver")', single_cell_runtime)
        # eikonalSolver is registered as a named alias of electrophysiologyModel
        # via addNamedToRunTimeSelectionTable, not via OverrideTypeName on a
        # separate class (the old standalone eikonalSolver class was removed).
        self.assertIn("addNamedToRunTimeSelectionTable", eikonal_runtime)
        self.assertIn("eikonalSolver", eikonal_runtime)
        self.assertIn('OverrideTypeName("monodomain1DSolver")', conduction_runtime)
        self.assertIn('OverrideTypeName("eikonalSolver1D")', _read(
            repo_root
            / "src"
            / "electroModels"
            / "conductionSystemModels"
            / "eikonalSolver1D"
            / "eikonalSolver1D.H"
        ))
        self.assertIn('OverrideTypeName("pseudoECG")', pseudo_ecg_runtime)

    def test_driver_schema_paths_follow_template_truth_family(self) -> None:
        documented = set(all_documented_driver_paths())

        expected = {
            "myocardiumSolver",
            "$ELECTRO_MODEL_COEFFS.solutionAlgorithm",
            "$ELECTRO_MODEL_COEFFS.verificationModel.type",
            "$ELECTRO_MODEL_COEFFS.bathPotentialDomain.bathCellZones",
            "$ELECTRO_MODEL_COEFFS.ecgDomains.<name>.ecgSolver",
            "$ELECTRO_MODEL_COEFFS.ecgDomains.<name>.electrodePositions.<electrode>",
            "$ELECTRO_MODEL_COEFFS.conductivityIntracellular",
            "$ELECTRO_MODEL_COEFFS.conductivityExtracellular",
            "$ELECTRO_MODEL_COEFFS.conductionNetworkDomains.<name>.conductionSystemDomain",
            "$ELECTRO_MODEL_COEFFS.conductionNetworkDomains.<name>.purkinjeGraphModelCoeffs.graphFile",
            "$ELECTRO_MODEL_COEFFS.domainCouplings.<name>.electroDomainCoupler",
            "$ELECTRO_MODEL_COEFFS.phiERefPoint",
        }
        self.assertTrue(expected.issubset(documented))

        forbidden_fragments = (
            "ecgDomainModel",
            "pseudoECGElectroCoeffs",
            ".electrodes.",
            "monoDomainElectroCoeffs",
            "singleCellElectroCoeffs",
            "phiEReferenceCell",
            "purkinjeNetworkModelCoeffs",
            ".pvjNodes",
            ".pvjLocations",
            "bidomainBath",
        )
        self.assertFalse(
            [path for path in documented if any(fragment in path for fragment in forbidden_fragments)]
        )

class TestMakeSpecDirectRun(unittest.TestCase):
    """make_spec records the solver command it was given.

    Formerly also covered _run_direct vs _run_case dispatch; both were
    removed with TutorialSpec.run_case, which nothing ever invoked. The
    live pre-solve-then-solver ordering is covered by
    test_dict_builder.py::TestBuildAndLaunchDirectRun.
    """

    def test_spec_accepts_solver_command(self) -> None:
        import tempfile
        from openfoam_driver.plugins.cardiacfoam.tutorials.generic_case import make_spec
        with tempfile.TemporaryDirectory() as d:
            case_dir = Path(d) / "mycase"
            case_dir.mkdir()
            spec = make_spec(
                tutorials_root=Path(d),
                case_dir_name="mycase",
                solver_command="cardiacFoam",
            )
            self.assertIsNotNone(spec)

    def test_spec_accepts_pre_solve_commands(self) -> None:
        import tempfile
        from openfoam_driver.plugins.cardiacfoam.tutorials.generic_case import make_spec
        with tempfile.TemporaryDirectory() as d:
            case_dir = Path(d) / "mycase"
            case_dir.mkdir()
            spec = make_spec(
                tutorials_root=Path(d),
                case_dir_name="mycase",
                solver_command="cardiacFoam",
                pre_solve_commands=["vtkUnstructuredToFoam", "setTorsoOrganConductivityField"],
            )
            self.assertIsNotNone(spec)

    def test_metadata_records_solver_command(self) -> None:
        import tempfile
        from openfoam_driver.plugins.cardiacfoam.tutorials.generic_case import make_spec
        with tempfile.TemporaryDirectory() as d:
            case_dir = Path(d) / "mycase"
            case_dir.mkdir()
            spec = make_spec(
                tutorials_root=Path(d),
                case_dir_name="mycase",
                solver_command="cardiacFoam",
                pre_solve_commands=["vtkUnstructuredToFoam"],
            )
            self.assertEqual(spec.metadata["solver_command"], "cardiacFoam")
            self.assertEqual(spec.metadata["pre_solve_commands"], ["vtkUnstructuredToFoam"])


if __name__ == "__main__":
    unittest.main()
