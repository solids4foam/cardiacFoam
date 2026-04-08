from __future__ import annotations

import unittest
from pathlib import Path

from openfoam_driver.dict_entries import all_documented_driver_paths
from openfoam_driver.specs.common import detect_electro_coeffs_scope


def _repo_root_from_test() -> Path:
    current = Path(__file__).resolve()
    for parent in current.parents:
        if (parent / "tutorials").exists() and (parent / "applications").exists():
            return parent
    raise RuntimeError("Could not locate repository root from test path")


def _read(path: Path) -> str:
    return path.read_text()


def _classify_electro_properties(text: str) -> str:
    has_template_truth = "myocardiumSolver" in text
    has_solver_scopes = any(
        token in text
        for token in (
            "monodomainSolverCoeffs",
            "bidomainSolverCoeffs",
            "singleCellSolverCoeffs",
            "eikonalSolverCoeffs",
        )
    )
    has_legacy_family = any(
        token in text
        for token in (
            "electroModel",
            "monoDomainElectroCoeffs",
            "singleCellElectroCoeffs",
            "ecgDomainModel",
            "pseudoECGElectroCoeffs",
            "electrodes",
        )
    )
    has_mixed_scope = "electroModelCoeffs" in text

    if has_template_truth and has_solver_scopes and not has_legacy_family and not has_mixed_scope:
        return "template_truth"
    if has_legacy_family and not has_template_truth:
        return "legacy_family"
    if has_template_truth or has_legacy_family or has_mixed_scope:
        return "mixed"
    return "unclassified"


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

        self.assertNotIn("ecgDomainModel", template)
        self.assertNotIn("pseudoECGElectroCoeffs", template)
        self.assertNotIn("monoDomainElectroCoeffs", template)

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
        electro_coupler = _read(
            self.repo_root
            / "src"
            / "electroModels"
            / "electroCouplers"
            / "electroDomainCoupler.C"
        )
        conduction_domain = _read(
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
        self.assertIn('"ecgSolver"', ecg_solver)
        self.assertIn('"electrodePositions"', ecg_domain)
        self.assertIn('"electroDomainCoupler"', electro_coupler)
        self.assertIn('subDict("purkinjeNetworkModelCoeffs").get<scalar>("cm")', conduction_domain)

        self.assertIn('OverrideTypeName("monodomainSolver")', monodomain_runtime)
        self.assertIn('OverrideTypeName("bidomainSolver")', bidomain_runtime)
        self.assertIn('OverrideTypeName("singleCellSolver")', single_cell_runtime)
        self.assertIn('OverrideTypeName("eikonalSolver")', eikonal_runtime)
        self.assertIn('OverrideTypeName("monodomain1DSolver")', conduction_runtime)
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
            "$ELECTRO_MODEL_COEFFS.domainCouplings.<name>.electroDomainCoupler",
        }
        self.assertTrue(expected.issubset(documented))

        forbidden_fragments = (
            "ecgDomainModel",
            "pseudoECGElectroCoeffs",
            ".electrodes.",
            "monoDomainElectroCoeffs",
            "singleCellElectroCoeffs",
        )
        self.assertFalse(
            [path for path in documented if any(fragment in path for fragment in forbidden_fragments)]
        )


class TestTutorialElectroPropertiesAudit(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.repo_root = _repo_root_from_test()
        cls.tutorials_root = cls.repo_root / "tutorials"

    def test_detect_electro_coeffs_scope_tracks_template_truth_cases(self) -> None:
        expected = {
            "tutorials/NiedererEtAl2012/constant/electroProperties": "monodomainSolverCoeffs",
            "tutorials/NiedererEtAl2012/referenceTest/constant/electroProperties": "monodomainSolverCoeffs",
            "tutorials/manufacturedFDA/constant/electroProperties": "monodomainSolverCoeffs",
            "tutorials/restitutionCurves_s1s2Protocol/constant/electroProperties": "singleCellSolverCoeffs",
            "tutorials/singleCell/constant/electroProperties": "singleCellSolverCoeffs",
            "tutorials/singleCell/referenceTest/constant/electroProperties": "singleCellSolverCoeffs",
        }

        for relpath, scope in expected.items():
            with self.subTest(path=relpath):
                self.assertEqual(
                    detect_electro_coeffs_scope(self.repo_root / relpath),
                    scope,
                )

    def test_all_live_tutorial_electroproperties_files_follow_template_truth(self) -> None:
        discovered = sorted(
            path
            for path in self.tutorials_root.rglob("electroProperties*")
            if path.parent.name == "constant"
        )
        actual = {
            str(path.relative_to(self.repo_root)): _classify_electro_properties(_read(path))
            for path in discovered
        }
        self.assertEqual(
            set(actual.values()),
            {"template_truth"},
        )

    def test_obsolete_with_default_values_files_are_absent(self) -> None:
        obsolete = list(self.tutorials_root.rglob("electroProperties.withDefaultValues"))
        self.assertEqual(obsolete, [])


if __name__ == "__main__":
    unittest.main()
