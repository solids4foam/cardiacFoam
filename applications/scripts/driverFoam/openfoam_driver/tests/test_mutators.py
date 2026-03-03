from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from openfoam_driver.core.runtime.mutators import update_foam_entry


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


if __name__ == "__main__":
    unittest.main()
