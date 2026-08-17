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
#     test_detection_and_overrides
#
# Description
#     Tests the cardiacFoam plugin's electroProperties detection helpers and
#     the electro/physics-property override appliers that sit on top of them,
#     including the plugin-local `$ELECTRO_MODEL_COEFFS` scope token. These
#     moved out of `tests/core/test_mutators.py` alongside their modules'
#     move out of `specs/`: the generic scope-resolution behaviour they lean
#     on is still covered there, this file covers only the cardiac layer.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from openfoam_driver.plugins.cardiacfoam.detection import detect_electro_coeffs_scope
from openfoam_driver.plugins.cardiacfoam.overrides import (
    apply_electro_property_overrides,
    apply_physics_property_overrides,
    ensure_electro_property_dict,
    normalize_entry_overrides,
    remove_electro_property_dict,
)
from openfoam_driver.tests.conftest import assert_foam_entry


class TestCardiacPropertyOverrides(unittest.TestCase):
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

            stimulus = ("singleCellSolverCoeffs", "singleCellStimulus")
            for key, expected in (
                ("stim_amplitude", "0.8"),
                ("stim_period_S1", "1200"),
                ("stim_period_S2", "300"),
                ("nstim1", "12"),
                ("nstim2", "3"),
            ):
                assert_foam_entry(path, key, expected, scope=stimulus)

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

            assert_foam_entry(
                path, "ionicModel", "Gaur", scope="singleCellSolverCoeffs"
            )
            assert_foam_entry(
                path,
                "stim_period_S1",
                "750",
                scope=("singleCellSolverCoeffs", "singleCellStimulus"),
            )

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
            assert_foam_entry(path, "type", "electroMechanicalModel")


if __name__ == "__main__":
    unittest.main()
