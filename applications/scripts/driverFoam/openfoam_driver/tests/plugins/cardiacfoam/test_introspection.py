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
#     test_introspection
#
# Description
#     Tests introspection logic and specification contracts.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

from __future__ import annotations

import io
import json
import tempfile
import unittest
from contextlib import redirect_stdout
from pathlib import Path

from openfoam_driver.cli import main
from openfoam_driver.introspection import describe_tutorial
from openfoam_driver.tests.conftest import monorepo_root, skip_without_monorepo


@skip_without_monorepo

class TestIntrospection(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.repo_root = monorepo_root
        cls.tutorials_root = monorepo_root / "tutorials"  # type: ignore[operator]

    def test_describe_tutorial_reports_registered_spec_schema(self) -> None:
        payload = describe_tutorial(
            "singleCell",
            overrides={"tutorials_root": self.tutorials_root},
        )

        self.assertEqual(payload["resolution"], "registered")
        self.assertEqual(payload["entry_kind"], "registered_tutorial")
        self.assertEqual(payload["resolved_name"], "singleCell")
        self.assertEqual(payload["entry"]["entry_name"], "singleCell")
        self.assertIn("singleCell", payload["registered_tutorials"])
        self.assertIn("ionic_models", payload["make_spec"]["parameters"])
        self.assertIn("entry_catalog", payload)
        self.assertIsInstance(payload["entry_catalog"], list)
        self.assertIn("workflow_catalog", payload)
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
        self.assertTrue(
            payload["spec"]["case_root"].endswith(
                "tutorials/electrophysiologyProtocols/singleCell"
            )
        )
        self.assertTrue(Path(payload["spec"]["case_root"]).exists())
        self.assertTrue(Path(payload["spec"]["setup_root"]).exists())
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
        # strict_launch is the canonical way to actually execute this entry.
        self.assertIn("strict_launch", payload)
        self.assertEqual(payload["strict_launch"]["action"], "run")
        self.assertIn("--strict", payload["strict_launch"]["command"])
        self.assertIn("singleCell", payload["strict_launch"]["command"])
        self.assertTrue(payload["strict_launch"]["case_root"])
        self.assertTrue(payload["strict_launch"]["output_dir"])

    def test_describe_tutorial_ionic_model_catalog_present_and_complete(self) -> None:
        payload = describe_tutorial(
            "singleCell",
            overrides={"tutorials_root": self.tutorials_root},
        )

        self.assertIn("ionic_model_catalog", payload)
        catalog = payload["ionic_model_catalog"]
        self.assertEqual(catalog["schema_version"], "1.0")

        # ionic_models
        self.assertIn("ionic_models", catalog)
        self.assertIn("TNNP", catalog["ionic_models"])
        tnnp = catalog["ionic_models"]["TNNP"]
        for key in ("states", "algebraic", "recommended_exports", "compatible_tissues",
                    "native_tissue_labels", "approximate_tissue_labels",
                    "compatible_solvers", "species", "cardiac_region", "model_type"):
            self.assertIn(key, tnnp, msg=f"Missing key {key!r} on TNNP")
        # tuples must be serialised to lists
        self.assertIsInstance(tnnp["states"], list)
        self.assertGreater(len(tnnp["states"]), 0)

        # solver_compatibility
        self.assertIn("solver_compatibility", catalog)
        self.assertIsInstance(catalog["solver_compatibility"], list)
        self.assertGreater(len(catalog["solver_compatibility"]), 0)

        self.assertIn("active_tension_catalog", payload)
        at_catalog = payload["active_tension_catalog"]
        self.assertEqual(at_catalog["schema_version"], "1.0")
        self.assertIn("active_tension_models", at_catalog)
        self.assertGreater(len(at_catalog["active_tension_models"]), 0)

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
            self.assertEqual(payload["entry_kind"], "case_folder")
            self.assertEqual(payload["resolved_name"], "randomCase")
            self.assertIn("randomCase", payload["available_tutorials"])
            self.assertEqual(
                payload["make_spec"]["callable"],
                "openfoam_driver.core.runtime.generic_case.make_spec",
            )
            self.assertEqual(payload["spec"]["cases"]["count"], 1)
            self.assertEqual(payload["tutorial_contract"]["name"], "randomCase")
            self.assertEqual(
                payload["tutorial_contract"]["core_required_files"],
                ["constant/electroProperties", "constant/physicsProperties"],
            )
            self.assertEqual(payload["tutorial_contract"]["solver_required_files"], [])
            self.assertEqual(payload["tutorial_contract"]["reference_cases"], [])

    def test_describe_tutorial_contract_file_lists_preserve_profile_order(self) -> None:
        # Regression test: describe_tutorial_contract's three public file
        # lists (core_required_files, solver_required_files,
        # conditional_files) must preserve the active plugin profile's
        # declaration order verbatim, not an alphabetically sorted order.
        # A prior fix wrapped these in sorted() to make an unrelated,
        # order-sensitive assertion pass; that silently reordered
        # conditional_files in the public payload. Every file referenced
        # below exists on disk so _existing_relpaths does not filter
        # anything out, and the assertions use assertEqual (not assertIn)
        # so a reordering fails this test.
        with tempfile.TemporaryDirectory() as temp_dir:
            tutorials_root = Path(temp_dir)
            case_root = tutorials_root / "orderedCase"
            (case_root / "constant").mkdir(parents=True, exist_ok=True)
            (case_root / "system").mkdir(parents=True, exist_ok=True)

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
            (case_root / "constant" / "physicsProperties").write_text(
                "type electroModel;\n"
            )
            for relpath in (
                "system/controlDict",
                "system/fvSchemes",
                "system/fvSolution",
                "system/decomposeParDict",
                "system/blockMeshDict",
            ):
                (case_root / relpath).write_text("// placeholder\n")
            for relpath in ("Allrun", "Allclean", "runRegressionTest.sh"):
                path = case_root / relpath
                path.write_text("#!/bin/sh\n")
                path.chmod(0o755)
            (case_root / "README.md").write_text("# orderedCase\n")

            payload = describe_tutorial(
                "orderedCase",
                overrides={"tutorials_root": tutorials_root},
            )

            contract = payload["tutorial_contract"]
            self.assertEqual(
                contract["core_required_files"],
                ["constant/electroProperties", "constant/physicsProperties"],
            )
            self.assertEqual(
                contract["solver_required_files"],
                ["system/controlDict", "system/fvSchemes", "system/fvSolution"],
            )
            self.assertEqual(
                contract["conditional_files"],
                [
                    "system/decomposeParDict",
                    "system/blockMeshDict",
                    "Allrun",
                    "Allclean",
                    "README.md",
                    "runRegressionTest.sh",
                ],
            )

    def test_cli_describe_prints_json_payload(self) -> None:
        stream = io.StringIO()
        with redirect_stdout(stream):
            exit_code = main(
                [
                    "describe",
                    "--entry",
                    "singleCell",
                    "--tutorials-root",
                    str(self.tutorials_root),
                ]
            )

        self.assertEqual(exit_code, 0)
        payload = json.loads(stream.getvalue())
        self.assertEqual(payload["resolved_name"], "singleCell")
        self.assertIn("common_override_keys", payload)
        self.assertIn("strict_launch", payload)
        self.assertEqual(payload["entry_kind"], "registered_tutorial")

    def test_cli_describe_requires_entry_name(self) -> None:
        stream = io.StringIO()
        with self.assertRaises(SystemExit):
            with redirect_stdout(stream):
                main(["describe", "--tutorials-root", str(self.tutorials_root)])

if __name__ == "__main__":
    unittest.main()
