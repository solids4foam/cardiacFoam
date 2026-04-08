from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from openfoam_driver.core.runtime.mutators import update_foam_entry
from openfoam_driver.specs.common import (
    apply_electro_property_overrides,
    apply_physics_property_overrides,
    detect_electro_coeffs_scope,
    normalize_entry_overrides,
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
