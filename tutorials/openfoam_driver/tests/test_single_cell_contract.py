from __future__ import annotations

import ast
import unittest
from pathlib import Path


class TestSingleCellContract(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.module_path = Path(__file__).resolve().parents[2] / "singleCell" / "setupSingleCell" / "singleCellinteractivePlots.py"
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
