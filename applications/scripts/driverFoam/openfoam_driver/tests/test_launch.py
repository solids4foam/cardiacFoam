from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from openfoam_driver.launch import describe_launch, describe_launch_matrix


class TestLaunchDescription(unittest.TestCase):
    def test_describe_launch_reports_manifest_and_command_for_registered_tutorial(self) -> None:
        repo_root = Path(__file__).resolve()
        for parent in repo_root.parents:
            if (parent / "tutorials").exists() and (parent / "applications").exists():
                repo_root = parent
                break
        else:
            self.fail("Could not locate repository root from test path")

        tutorials_root = repo_root / "tutorials"
        payload = describe_launch(
            "sim",
            "singleCell",
            tutorials_root=tutorials_root,
            continue_on_error=True,
        )

        self.assertEqual(payload["action"], "sim")
        self.assertEqual(payload["resolved_name"], "singleCell")
        self.assertTrue(payload["manifest_path"].endswith("run_manifest.json"))
        self.assertIn("--continue-on-error", payload["command"])
        self.assertIn("singleCell", payload["command_display"])

    def test_describe_launch_matrix_handles_case_folder_resolution(self) -> None:
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

            payload = describe_launch_matrix(
                "randomCase",
                tutorials_root=tutorials_root,
            )

            self.assertEqual(payload["sim"]["resolution"], "case_folder")
            self.assertEqual(payload["all"]["resolved_name"], "randomCase")
            self.assertTrue(payload["post"]["manifest_path"].endswith("run_manifest.json"))


if __name__ == "__main__":
    unittest.main()
