from __future__ import annotations

import io
import json
import tempfile
import unittest
from contextlib import redirect_stdout
from pathlib import Path

from openfoam_driver.cli import main
from openfoam_driver.introspection import describe_tutorial


def _repo_root_from_test() -> Path:
    current = Path(__file__).resolve()
    for parent in current.parents:
        if (parent / "tutorials").exists() and (parent / "applications").exists():
            return parent
    raise RuntimeError("Could not locate repository root from test path")


class TestIntrospection(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.repo_root = _repo_root_from_test()
        cls.tutorials_root = cls.repo_root / "tutorials"

    def test_describe_tutorial_reports_registered_spec_schema(self) -> None:
        payload = describe_tutorial(
            "singleCell",
            overrides={"tutorials_root": self.tutorials_root},
        )

        self.assertEqual(payload["resolution"], "registered")
        self.assertEqual(payload["resolved_name"], "singleCell")
        self.assertIn("singleCell", payload["registered_tutorials"])
        self.assertIn("ionic_models", payload["make_spec"]["parameters"])
        self.assertIn("gui_schema", payload)
        self.assertEqual(payload["gui_schema"]["routes"][0]["path"], "/tutorials")
        self.assertIn(
            "$ELECTRO_MODEL_COEFFS.singleCellStimulus.stim_period_S1",
            {
                item["driver_path"]
                for item in payload["dict_entries"]["electroProperties"]["single_cell_stimulus"]
            },
        )
        solution_algorithm = next(
            item
            for item in payload["dict_entries"]["electroProperties"]["common_model_coeffs"]
            if item["driver_path"] == "$ELECTRO_MODEL_COEFFS.solutionAlgorithm"
        )
        self.assertEqual(solution_algorithm["value_kind"], "enum")
        self.assertIn("implicit", solution_algorithm["enum_values"])
        self.assertGreater(payload["spec"]["cases"]["count"], 0)
        self.assertIn("tutorial_contract", payload)
        self.assertEqual(payload["tutorial_contract"]["name"], "singleCell")
        self.assertIn(
            "constant/electroProperties",
            payload["tutorial_contract"]["core_required_files"],
        )
        self.assertIn(
            "system/fvSchemes",
            payload["tutorial_contract"]["solver_required_files"],
        )
        self.assertIn(
            "system/fvSolution",
            payload["tutorial_contract"]["solver_required_files"],
        )
        self.assertIn("Allrun", payload["tutorial_contract"]["conditional_files"])
        self.assertIn("README.md", payload["tutorial_contract"]["conditional_files"])
        self.assertIn("ionicModel", payload["tutorial_contract"]["case_parameters"])
        self.assertIn(
            "regressionTests/singleCell",
            payload["tutorial_contract"]["reference_cases"],
        )
        self.assertIn("launch", payload)
        self.assertEqual(payload["launch"]["sim"]["action"], "sim")
        self.assertTrue(payload["launch"]["all"]["manifest_path"].endswith("run_manifest.json"))

    def test_describe_tutorial_reports_generic_case_folder_resolution(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            tutorials_root = Path(temp_dir)
            case_root = tutorials_root / "randomCase"
            (case_root / "constant").mkdir(parents=True, exist_ok=True)
            (case_root / "constant" / "electroProperties").write_text(
                "\n".join(
                    [
                        "myocardiumSolver singleCellSolver;",
                        "",
                        "singleCellSolverCoeffs",
                        "{",
                        "    ionicModel BuenoOrovio;",
                        "}",
                        "",
                    ]
                )
            )
            (case_root / "constant" / "physicsProperties").write_text("type electroModel;\n")

            payload = describe_tutorial(
                "randomCase",
                overrides={"tutorials_root": tutorials_root},
            )

            self.assertEqual(payload["resolution"], "case_folder")
            self.assertEqual(payload["resolved_name"], "randomCase")
            self.assertIn("randomCase", payload["available_tutorials"])
            self.assertEqual(
                payload["make_spec"]["callable"],
                "openfoam_driver.specs.tutorials.generic_case.make_spec",
            )
            self.assertEqual(payload["spec"]["cases"]["count"], 1)
            self.assertEqual(payload["tutorial_contract"]["name"], "randomCase")
            self.assertEqual(
                payload["tutorial_contract"]["core_required_files"],
                ["constant/electroProperties", "constant/physicsProperties"],
            )
            self.assertEqual(payload["tutorial_contract"]["solver_required_files"], [])
            self.assertEqual(payload["tutorial_contract"]["reference_cases"], [])

    def test_cli_describe_prints_json_payload(self) -> None:
        stream = io.StringIO()
        with redirect_stdout(stream):
            exit_code = main(
                [
                    "describe",
                    "--tutorial",
                    "singleCell",
                    "--tutorials-root",
                    str(self.tutorials_root),
                ]
            )

        self.assertEqual(exit_code, 0)
        payload = json.loads(stream.getvalue())
        self.assertEqual(payload["resolved_name"], "singleCell")
        self.assertIn("common_override_keys", payload)
        self.assertIn("launch", payload)

    def test_heartsimtemplate_exposes_machine_readable_authoring_contract(self) -> None:
        payload = describe_tutorial(
            "HeartSimTemplate",
            overrides={"tutorials_root": self.tutorials_root},
        )

        contract = payload["tutorial_contract"]
        self.assertEqual(contract["name"], "HeartSimTemplate")
        self.assertEqual(contract["resolution"], "case_folder")
        self.assertEqual(contract["authoring_contract_path"], "workflow_contract.json")
        self.assertIsNotNone(contract["authoring_contract"])
        self.assertEqual(
            contract["authoring_contract"]["tutorial_family"],
            "HeartSimTemplate",
        )
        self.assertEqual(
            contract["authoring_contract"]["workflow_templates"][0]["workflow_id"],
            "monodomain_purkinje_pseudo_ecg",
        )


if __name__ == "__main__":
    unittest.main()
