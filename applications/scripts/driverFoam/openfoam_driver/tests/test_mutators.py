from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from openfoam_driver.core.runtime.mutators import update_foam_entry
from openfoam_driver.specs.common import (
    set_n_stim1,
    set_n_stim2,
    set_s1_period,
    set_s2_period,
    set_stimulus_amplitude,
)


class TestScopedMutators(unittest.TestCase):
    def test_updates_only_within_scope(self) -> None:
        text = "\n".join(
            [
                "electroModel monoDomainElectro;",
                "",
                "monoDomainElectroCoeffs",
                "{",
                "    ionicModel TNNP;",
                "}",
                "",
                "singleCellElectroCoeffs",
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
                scope="singleCellElectroCoeffs",
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

    def test_single_cell_stimulus_helpers_use_nested_scope(self) -> None:
        text = "\n".join(
            [
                "singleCellElectroCoeffs",
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

            set_stimulus_amplitude(path, 0.8, scope="singleCellElectroCoeffs")
            set_s1_period(path, 1200, scope="singleCellElectroCoeffs")
            set_s2_period(path, 300, scope="singleCellElectroCoeffs")
            set_n_stim1(path, 12, scope="singleCellElectroCoeffs")
            set_n_stim2(path, 3, scope="singleCellElectroCoeffs")

            updated = path.read_text()
            self.assertIn("stim_amplitude    0.8;", updated)
            self.assertIn("stim_period_S1    1200;", updated)
            self.assertIn("stim_period_S2    300;", updated)
            self.assertIn("nstim1    12;", updated)
            self.assertIn("nstim2    3;", updated)


if __name__ == "__main__":
    unittest.main()
