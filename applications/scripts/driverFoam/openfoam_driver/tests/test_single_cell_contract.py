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
#     test_single_cell_contract
#
# Description
#     Tests single cell contract logic and specification contracts.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

from __future__ import annotations

import ast
import unittest
from pathlib import Path

from openfoam_driver.core.runtime.registry import load_tutorial_spec


class TestSingleCellContract(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        current = Path(__file__).resolve()
        repo_root = None
        for parent in current.parents:
            if (parent / "tutorials").exists() and (parent / "applications").exists():
                repo_root = parent
                break
        if repo_root is None:
            raise RuntimeError("Could not locate repository root from test path")

        spec = load_tutorial_spec(
            "singleCell",
            overrides={"tutorials_root": repo_root / "tutorials"},
        )
        cls.module_path = spec.setup_root / "singleCellinteractivePlots.py"
        cls.tree = ast.parse(cls.module_path.read_text())

    def test_load_simulation_data_no_output_folder_dependency(self) -> None:
        target = None
        for node in self.tree.body:
            if isinstance(node, ast.FunctionDef) and node.name == "load_simulation_data":
                target = node
                break
        self.assertIsNotNone(target, "load_simulation_data function not found")

        arg_names = [arg.arg for arg in target.args.args]
        kwonly_names = [arg.arg for arg in target.args.kwonlyargs]
        self.assertIn("filename", arg_names)
        self.assertIn("base_folder", kwonly_names)

        names = {node.id for node in ast.walk(target) if isinstance(node, ast.Name)}
        self.assertNotIn("OUTPUT_FOLDER", names)

    def test_run_postprocessing_exposes_automation_kwargs(self) -> None:
        target = None
        for node in self.tree.body:
            if isinstance(node, ast.FunctionDef) and node.name == "run_postprocessing":
                target = node
                break
        self.assertIsNotNone(target, "run_postprocessing function not found")

        kwonly_names = [arg.arg for arg in target.args.kwonlyargs]
        self.assertIn("files", kwonly_names)
        self.assertIn("categories", kwonly_names)
        self.assertIn("rename_legends", kwonly_names)


if __name__ == "__main__":
    unittest.main()
