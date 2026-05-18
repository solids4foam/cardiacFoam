from __future__ import annotations

import unittest
from pathlib import Path

from openfoam_driver.dict_entries import all_documented_driver_paths


def _repo_root_from_test() -> Path:
    current = Path(__file__).resolve()
    for parent in current.parents:
        if (parent / "tutorials").exists() and (parent / "applications").exists():
            return parent
    raise RuntimeError("Could not locate repository root from test path")


def _read(path: Path) -> str:
    return path.read_text()


class TestTemplateAndSchemaContract(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.repo_root = _repo_root_from_test()
        cls.template_path = cls.repo_root / "tutorials" / "electroProperties.template"

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
        self.assertIn("potentialDomain", template)
        self.assertIn("extracellularPotentialDomain", template)
        self.assertIn("bathCellZones", template)
        self.assertIn("bathConductivityField", template)
        self.assertIn("bathECGProbe", template)
        self.assertIn("groundPatches", template)
        self.assertIn("surfaceCurrentPatches", template)

        self.assertNotIn("ecgDomainModel", template)
        self.assertNotIn("pseudoECGElectroCoeffs", template)
        self.assertNotIn("monoDomainElectroCoeffs", template)
        self.assertNotIn("phiEReferenceCell", template)
        self.assertNotIn("purkinjeNetworkModelCoeffs", template)
        self.assertNotIn("bidomainBathECG", template)

    def test_template_truth_markers_match_core_cpp_contracts(self) -> None:
        electro_model_core = _read(
            self.repo_root / "src" / "electroModels" / "core" / "electroModel.C"
        )
        bidomain_solver = _read(
            self.repo_root
            / "src"
            / "electroModels"
            / "myocardiumModels"
            / "bidomainSolver"
            / "bidomainSolver.C"
        )
        ecg_solver = _read(
            self.repo_root
            / "src"
            / "electroModels"
            / "electroDomains"
            / "ecgDomain"
            / "ecgSolver.C"
        )
        ecg_domain = _read(
            self.repo_root
            / "src"
            / "electroModels"
            / "electroDomains"
            / "ecgDomain"
            / "ecgDomain.C"
        )
        system_builder = _read(
            self.repo_root
            / "src"
            / "electroModels"
            / "core"
            / "system"
            / "electrophysicsSystemBuilder.C"
        )
        electro_coupler = _read(
            self.repo_root
            / "src"
            / "electroModels"
            / "electroCouplers"
            / "electroDomainCoupler.C"
        )
        conduction_domain_selector = _read(
            self.repo_root
            / "src"
            / "electroModels"
            / "electroDomains"
            / "conductionSystemDomain"
            / "conductionSystemDomain.C"
        )
        monodomain_runtime = _read(
            self.repo_root
            / "src"
            / "electroModels"
            / "myocardiumModels"
            / "monodomainSolver"
            / "monodomainSolver.H"
        )
        bidomain_runtime = _read(
            self.repo_root
            / "src"
            / "electroModels"
            / "myocardiumModels"
            / "bidomainSolver"
            / "bidomainSolver.H"
        )
        single_cell_runtime = _read(
            self.repo_root
            / "src"
            / "electroModels"
            / "myocardiumModels"
            / "singleCellSolver"
            / "singleCellSolver.H"
        )
        eikonal_runtime = _read(
            self.repo_root
            / "src"
            / "electroModels"
            / "myocardiumModels"
            / "eikonalSolver"
            / "eikonalSolver.H"
        )
        conduction_runtime = _read(
            self.repo_root
            / "src"
            / "electroModels"
            / "conductionSystemModels"
            / "monodomain1DSolver"
            / "monodomain1DSolver.H"
        )
        pseudo_ecg_runtime = _read(
            self.repo_root
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
        self.assertIn("bathECGProbe requires myocardiumSolver", system_builder)
        self.assertIn("bidomainSolver because it samples the", system_builder)
        self.assertIn('"electroDomainCoupler"', electro_coupler)
        self.assertIn('"purkinjeGraphModel"', conduction_domain_selector)
        self.assertIn('dict.get<word>("graphFile")', conduction_domain_selector)
        self.assertIn('dict.subDict("rootStimulus")', conduction_domain_selector)

        self.assertIn('OverrideTypeName("monodomainSolver")', monodomain_runtime)
        self.assertIn('OverrideTypeName("bidomainSolver")', bidomain_runtime)
        self.assertIn('OverrideTypeName("singleCellSolver")', single_cell_runtime)
        self.assertIn('OverrideTypeName("eikonalSolver")', eikonal_runtime)
        self.assertIn('OverrideTypeName("monodomain1DSolver")', conduction_runtime)
        self.assertIn('OverrideTypeName("eikonalSolver1D")', _read(
            self.repo_root
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

if __name__ == "__main__":
    unittest.main()
