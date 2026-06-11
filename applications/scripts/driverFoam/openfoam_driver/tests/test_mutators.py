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
#     test_mutators
#
# Description
#     Tests mutators logic and specification contracts.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from openfoam_driver.core.runtime.mutators import (
    ensure_foam_dict,
    remove_foam_dict,
    update_foam_entry,
)
from openfoam_driver.specs.common import (
    apply_electro_property_overrides,
    apply_physics_property_overrides,
    detect_electro_coeffs_scope,
    ensure_electro_property_dict,
    normalize_entry_overrides,
    remove_electro_property_dict,
)


class TestScopedMutators(unittest.TestCase):
    def test_updates_only_within_scope(self) -> None:
        text = "\n".join(
            [
                "myocardiumSolver monodomainSolver;",
                "",
                "monodomainSolverCoeffs",
                "{",
                "    ionicModel TNNP;",
                "}",
                "",
                "singleCellSolverCoeffs",
                "{",
                "    ionicModel BuenoOrovio;",
                "}",
                "",
            ]
        )

        with tempfile.TemporaryDirectory() as temp_dir:
            path = Path(temp_dir) / "electroProperties"
            path.write_text(text)

            update_foam_entry(
                path,
                "ionicModel",
                "Gaur",
                scope="singleCellSolverCoeffs",
            )

            updated = path.read_text()
            self.assertIn("ionicModel TNNP;", updated)
            self.assertIn("ionicModel    Gaur;", updated)

    def test_nested_scope_path(self) -> None:
        text = "\n".join(
            [
                "outer",
                "{",
                "    inner",
                "    {",
                "        target 1;",
                "    }",
                "}",
                "",
            ]
        )

        with tempfile.TemporaryDirectory() as temp_dir:
            path = Path(temp_dir) / "dict"
            path.write_text(text)

            update_foam_entry(path, "target", 2, scope=("outer", "inner"))
            self.assertIn("target    2;", path.read_text())

    def test_missing_scope_raises(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            path = Path(temp_dir) / "dict"
            path.write_text("a { b 1; }\n")

            with self.assertRaises(KeyError):
                update_foam_entry(path, "b", 2, scope="missing")

    def test_python_parser_fails_on_c_style_comments(self) -> None:
        text = "\n".join(
            [
                "someDict",
                "{",
                "    /* This is a block comment with a brace { inside it */",
                "    value 1;",
                "}",
                "",
            ]
        )

        with tempfile.TemporaryDirectory() as temp_dir:
            path = Path(temp_dir) / "dict"
            path.write_text(text)

            # The current update_foam_entry uses brace counting, so the extra {
            # inside the block comment throws off the parser, causing it to
            # incorrectly raise a KeyError for unbalanced braces.
            with self.assertRaises(KeyError):
                update_foam_entry(path, "value", 2, scope="someDict")

    def test_remove_foam_dict_removes_nested_dictionary(self) -> None:
        text = "\n".join(
            [
                "outer",
                "{",
                "    keep 1;",
                "    removeMe",
                "    {",
                "        nested",
                "        {",
                "            value 1;",
                "        }",
                "    }",
                "    after 2;",
                "}",
                "",
            ]
        )

        with tempfile.TemporaryDirectory() as temp_dir:
            path = Path(temp_dir) / "dict"
            path.write_text(text)

            remove_foam_dict(path, "removeMe", scope="outer")

            updated = path.read_text()
            self.assertIn("keep 1;", updated)
            self.assertIn("after 2;", updated)
            self.assertNotIn("removeMe", updated)
            self.assertNotIn("value 1;", updated)

    def test_single_cell_stimulus_updates_use_nested_scope(self) -> None:
        text = "\n".join(
            [
                "singleCellSolverCoeffs",
                "{",
                "    singleCellStimulus",
                "    {",
                "        stim_amplitude 0.4;",
                "        stim_period_S1 1000;",
                "        stim_period_S2 250;",
                "        nstim1 10;",
                "        nstim2 2;",
                "    }",
                "}",
                "",
            ]
        )

        with tempfile.TemporaryDirectory() as temp_dir:
            path = Path(temp_dir) / "electroProperties"
            path.write_text(text)

            apply_electro_property_overrides(
                path,
                {
                    "singleCellSolverCoeffs.singleCellStimulus.stim_amplitude": 0.8,
                    "singleCellSolverCoeffs.singleCellStimulus.stim_period_S1": 1200,
                    "singleCellSolverCoeffs.singleCellStimulus.stim_period_S2": 300,
                    "singleCellSolverCoeffs.singleCellStimulus.nstim1": 12,
                    "singleCellSolverCoeffs.singleCellStimulus.nstim2": 3,
                },
            )

            updated = path.read_text()
            self.assertIn("stim_amplitude    0.8;", updated)
            self.assertIn("stim_period_S1    1200;", updated)
            self.assertIn("stim_period_S2    300;", updated)
            self.assertIn("nstim1    12;", updated)
            self.assertIn("nstim2    3;", updated)

    def test_detect_electro_coeffs_scope(self) -> None:
        text = "\n".join(
            [
                "myocardiumSolver monodomainSolver;",
                "",
                "monodomainSolverCoeffs",
                "{",
                "    ionicModel TNNP;",
                "}",
                "",
            ]
        )

        with tempfile.TemporaryDirectory() as temp_dir:
            path = Path(temp_dir) / "electroProperties"
            path.write_text(text)
            self.assertEqual(detect_electro_coeffs_scope(path), "monodomainSolverCoeffs")

    def test_normalize_entry_overrides_supports_electro_scope_token(self) -> None:
        text = "\n".join(
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

        with tempfile.TemporaryDirectory() as temp_dir:
            path = Path(temp_dir) / "electroProperties"
            path.write_text(text)

            normalized = normalize_entry_overrides(
                {"$ELECTRO_MODEL_COEFFS.ionicModel": "Gaur"},
                electro_properties_path=path,
            )

            self.assertEqual(
                normalized,
                [{"key": "ionicModel", "value": "Gaur", "scope": ("singleCellSolverCoeffs",)}],
            )

    def test_apply_electro_property_overrides_handles_nested_paths(self) -> None:
        text = "\n".join(
            [
                "myocardiumSolver singleCellSolver;",
                "",
                "singleCellSolverCoeffs",
                "{",
                "    ionicModel BuenoOrovio;",
                "    singleCellStimulus",
                "    {",
                "        stim_period_S1 1000;",
                "    }",
                "}",
                "",
            ]
        )

        with tempfile.TemporaryDirectory() as temp_dir:
            path = Path(temp_dir) / "electroProperties"
            path.write_text(text)

            apply_electro_property_overrides(
                path,
                {
                    "$ELECTRO_MODEL_COEFFS.ionicModel": "Gaur",
                    "$ELECTRO_MODEL_COEFFS.singleCellStimulus.stim_period_S1": 750,
                },
            )

            updated = path.read_text()
            self.assertIn("ionicModel    Gaur;", updated)
            self.assertIn("stim_period_S1    750;", updated)

    def test_remove_electro_property_dict_supports_electro_scope_token(self) -> None:
        text = "\n".join(
            [
                "myocardiumSolver bidomainSolver;",
                "",
                "bidomainSolverCoeffs",
                "{",
                "    ecgDomains",
                "    {",
                "        ECG",
                "        {",
                "            ecgSolver torsoECG;",
                "        }",
                "    }",
                "    bathPotentialDomain",
                "    {",
                "        bathCellZones (bath);",
                "    }",
                "}",
                "",
            ]
        )

        with tempfile.TemporaryDirectory() as temp_dir:
            path = Path(temp_dir) / "electroProperties"
            path.write_text(text)

            remove_electro_property_dict(
                path,
                "ecgDomains",
                scope="$ELECTRO_MODEL_COEFFS",
            )

            updated = path.read_text()
            self.assertNotIn("ecgDomains", updated)
            self.assertIn("bathPotentialDomain", updated)

    def test_ensure_foam_dict_inserts_missing_dict_in_scope(self) -> None:
        text = "\n".join(
            [
                "root",
                "{",
                "    existing yes;",
                "}",
                "",
            ]
        )

        with tempfile.TemporaryDirectory() as temp_dir:
            path = Path(temp_dir) / "dict"
            path.write_text(text)

            inserted = ensure_foam_dict(
                path,
                "newBlock",
                "    newBlock\n    {\n        value 1;\n    }\n",
                scope="root",
            )
            inserted_again = ensure_foam_dict(
                path,
                "newBlock",
                "    newBlock\n    {\n        value 2;\n    }\n",
                scope="root",
            )

            updated = path.read_text()
            self.assertTrue(inserted)
            self.assertFalse(inserted_again)
            self.assertEqual(updated.count("newBlock"), 1)
            self.assertIn("value 1;", updated)

    def test_ensure_electro_property_dict_supports_electro_scope_token(self) -> None:
        text = "\n".join(
            [
                "myocardiumSolver bidomainSolver;",
                "",
                "bidomainSolverCoeffs",
                "{",
                "    bathPotentialDomain",
                "    {",
                "        bathCellZones (bath);",
                "    }",
                "}",
                "",
            ]
        )

        with tempfile.TemporaryDirectory() as temp_dir:
            path = Path(temp_dir) / "electroProperties"
            path.write_text(text)

            inserted = ensure_electro_property_dict(
                path,
                "ecgDomains",
                "    ecgDomains\n    {\n        ECG {}\n    }\n",
                scope="$ELECTRO_MODEL_COEFFS",
            )

            updated = path.read_text()
            self.assertTrue(inserted)
            self.assertIn("ecgDomains", updated)
            self.assertIn("bathPotentialDomain", updated)

    def test_apply_physics_property_overrides_updates_root_dictionary(self) -> None:
        text = "\n".join(
            [
                "type electroModel;",
                "",
            ]
        )

        with tempfile.TemporaryDirectory() as temp_dir:
            path = Path(temp_dir) / "physicsProperties"
            path.write_text(text)

            apply_physics_property_overrides(path, {"type": "electroMechanicalModel"})
            self.assertIn("type    electroMechanicalModel;", path.read_text())


if __name__ == "__main__":
    unittest.main()
