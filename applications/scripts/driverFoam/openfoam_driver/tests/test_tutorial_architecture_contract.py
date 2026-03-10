from __future__ import annotations

import inspect
import unittest
from pathlib import Path

from openfoam_driver.core.defaults import niederer_2012 as niederer_defaults
from openfoam_driver.core.runtime.models import CaseConfig, TutorialSpec
from openfoam_driver.core.runtime.registry import list_tutorials, load_tutorial_spec
from openfoam_driver.specs.tutorials import (
    manufactured_fda,
    niederer_2012,
    restitution_curves,
    single_cell,
)


def _repo_root_from_test() -> Path:
    current = Path(__file__).resolve()
    for parent in current.parents:
        if (parent / "tutorials").exists() and (parent / "applications").exists():
            return parent
    raise RuntimeError("Could not locate repository root from test path")


class TestTutorialArchitectureContract(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.repo_root = _repo_root_from_test()
        cls.tutorials_root = cls.repo_root / "tutorials"
        cls.module_map = {
            "singleCell": single_cell,
            "niederer2012": niederer_2012,
            "manufacturedFDA": manufactured_fda,
            "restitutionCurves": restitution_curves,
        }

    def _load_spec(self, tutorial: str) -> TutorialSpec:
        return load_tutorial_spec(
            tutorial,
            overrides={"tutorials_root": self.tutorials_root},
        )

    def test_registry_contains_expected_tutorials(self) -> None:
        self.assertEqual(
            set(list_tutorials()),
            {"singleCell", "niederer2012", "manufacturedFDA", "restitutionCurves"},
        )

    def test_all_tutorial_specs_load(self) -> None:
        for tutorial in list_tutorials():
            with self.subTest(tutorial=tutorial):
                spec = self._load_spec(tutorial)
                self.assertIsInstance(spec, TutorialSpec)
                self.assertTrue(callable(spec.build_cases))
                self.assertTrue(callable(spec.apply_case))
                self.assertTrue(callable(spec.run_case))
                self.assertIsNotNone(spec.postprocess)
                self.assertIn("notes", spec.metadata)
                self.assertIn("run_script_relpath", spec.metadata)
                self.assertIn("postprocess_strict_artifacts", spec.metadata)
                self.assertEqual(spec.case_root.parent, self.tutorials_root)
                self.assertEqual(spec.setup_root.parent, spec.case_root)
                self.assertEqual(spec.output_dir.parent, spec.case_root)

    def test_case_sweeps_are_non_empty_and_unique(self) -> None:
        for tutorial in list_tutorials():
            with self.subTest(tutorial=tutorial):
                spec = self._load_spec(tutorial)
                cases = spec.build_cases()
                self.assertTrue(cases, "build_cases() should return at least one case")
                self.assertTrue(all(isinstance(case, CaseConfig) for case in cases))
                ids = [case.case_id for case in cases]
                self.assertEqual(len(ids), len(set(ids)), "case_id values must be unique")

    def test_make_spec_has_common_keyword_contract(self) -> None:
        required_kwonly = {
            "tutorials_root",
            "case_dir_name",
            "setup_dir_name",
            "output_dir_name",
            "run_script_relpath",
            "postprocess_strict_artifacts",
        }
        for tutorial, module in self.module_map.items():
            with self.subTest(tutorial=tutorial):
                signature = inspect.signature(module.make_spec)
                kwonly = {
                    name
                    for name, param in signature.parameters.items()
                    if param.kind is inspect.Parameter.KEYWORD_ONLY
                }
                missing = sorted(required_kwonly - kwonly)
                self.assertFalse(
                    missing,
                    f"make_spec() is missing required keyword-only params: {missing}",
                )

    def test_output_collection_strategy_is_explicit(self) -> None:
        for tutorial in list_tutorials():
            with self.subTest(tutorial=tutorial):
                spec = self._load_spec(tutorial)
                if spec.collect_outputs is not None:
                    self.assertTrue(callable(spec.collect_outputs))
                    continue

                run_case_keywords = getattr(spec.run_case, "keywords", {}) or {}
                self.assertTrue(
                    "output_dir" in run_case_keywords or "output_relpath" in run_case_keywords,
                    "Specs without collect_outputs must bind an output path in run_case",
                )

    def test_niederer_solver_default_tracks_defaults_module(self) -> None:
        signature = inspect.signature(niederer_2012.make_spec)
        self.assertEqual(
            signature.parameters["solvers"].default,
            niederer_defaults.SOLVERS,
        )


if __name__ == "__main__":
    unittest.main()
