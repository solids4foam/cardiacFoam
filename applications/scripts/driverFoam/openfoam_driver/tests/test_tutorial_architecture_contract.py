from __future__ import annotations

import inspect
import tempfile
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

    def test_manufactured_defaults_to_multid_implicit_cases(self) -> None:
        spec = self._load_spec("manufacturedFDA")
        cases = spec.build_cases()

        self.assertTrue(cases)
        self.assertEqual(
            {case.params["dimension"] for case in cases},
            {"1D", "2D", "3D"},
        )
        self.assertTrue(all(case.params["solver"] == "implicit" for case in cases))
        self.assertTrue(spec.run_case.keywords["run_in_parallel"])
        self.assertTrue(spec.metadata["ecg_enabled"])
        self.assertEqual(sorted({int(case.params["cells"]) for case in cases}), [10, 20, 40, 80])

    def test_manufactured_stage_output_falls_back_to_processor0_postprocessing(self) -> None:
        case = CaseConfig(
            case_id="manufactured",
            params={"dimension": "1D", "solver": "implicit", "cells": 20, "dt": 0.1},
        )

        with tempfile.TemporaryDirectory() as temp_dir:
            case_root = Path(temp_dir)
            source_dir = case_root / "processor0" / "postProcessing"
            source_dir.mkdir(parents=True, exist_ok=True)
            source_file = source_dir / "1D_20_cells_implicit.dat"
            source_file.write_text("parallel-output")

            staged = manufactured_fda._stage_case_output(case_root, case)

            self.assertEqual(staged, case_root / "postProcessing" / "1D_20_cells_implicit.dat")
            self.assertEqual(staged.read_text(), "parallel-output")

    def test_manufactured_collect_outputs_copies_archived_logs(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            case_root = root / "case"
            output_dir = root / "output"
            (case_root / "logs" / "caseA").mkdir(parents=True, exist_ok=True)
            (case_root / "logs" / "caseA" / "log.blockMesh").write_text("mesh-log")

            manufactured_fda._collect_outputs(case_root, output_dir)

            self.assertEqual(
                (output_dir / "logs" / "caseA" / "log.blockMesh").read_text(),
                "mesh-log",
            )

    def test_manufactured_collect_outputs_copies_dat_files_without_removing_case_archive(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            case_root = root / "case"
            output_dir = root / "output"
            case_root.mkdir(parents=True, exist_ok=True)
            output_dir.mkdir(parents=True, exist_ok=True)
            archived_dir = case_root / "postProcessing"
            archived_dir.mkdir(parents=True, exist_ok=True)
            archived = archived_dir / "1D_80_cells_implicit.dat"
            archived.write_text("archived-output")

            manufactured_fda._collect_outputs(case_root, output_dir)

            self.assertEqual((output_dir / archived.name).read_text(), "archived-output")
            self.assertTrue(archived.exists(), "collect_outputs should preserve postProcessing archives")

    def test_manufactured_collect_outputs_preserves_existing_output_when_case_archive_missing(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            case_root = root / "case"
            output_dir = root / "output"
            case_root.mkdir(parents=True, exist_ok=True)
            output_dir.mkdir(parents=True, exist_ok=True)
            existing = output_dir / "1D_80_cells_implicit.dat"
            existing.write_text("existing-output")

            manufactured_fda._collect_outputs(case_root, output_dir)

            self.assertEqual(existing.read_text(), "existing-output")

    def test_manufactured_collect_outputs_is_noop_when_archive_dir_is_output_dir(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            case_root = Path(temp_dir)
            output_dir = case_root / "postProcessing"
            output_dir.mkdir(parents=True, exist_ok=True)
            archived = output_dir / "1D_10_cells_implicit.dat"
            archived.write_text("archived-output")

            manufactured_fda._collect_outputs(case_root, output_dir)

            self.assertEqual(archived.read_text(), "archived-output")

    def test_manufactured_stages_ecg_outputs_with_case_prefix(self) -> None:
        case = CaseConfig(
            case_id="manufactured",
            params={"dimension": "3D", "solver": "implicit", "cells": 20, "dt": 0.1},
        )

        with tempfile.TemporaryDirectory() as temp_dir:
            case_root = Path(temp_dir)
            source_dir = case_root / "processor0" / "postProcessing"
            source_dir.mkdir(parents=True, exist_ok=True)
            (source_dir / "pseudoECG.dat").write_text("pseudo-output")
            (source_dir / "manufacturedPseudoECG.dat").write_text("manufactured-output")

            staged = manufactured_fda._stage_case_ecg_outputs(case_root, case)

            staged_names = {path.name for path in staged}
            self.assertIn("ECG_manufactured_pseudoECG.dat", staged_names)
            self.assertIn("ECG_manufactured_manufacturedPseudoECG.dat", staged_names)
            self.assertEqual(
                (case_root / "postProcessing" / "ECG_manufactured_pseudoECG.dat").read_text(),
                "pseudo-output",
            )

    def test_manufactured_moves_root_postprocessing_ecg_files_to_prefixed_archives(self) -> None:
        case = CaseConfig(
            case_id="manufactured",
            params={"dimension": "3D", "solver": "implicit", "cells": 20, "dt": 0.1},
        )

        with tempfile.TemporaryDirectory() as temp_dir:
            case_root = Path(temp_dir)
            source_dir = case_root / "postProcessing"
            source_dir.mkdir(parents=True, exist_ok=True)
            raw = source_dir / "pseudoECG.dat"
            raw.write_text("pseudo-output")

            staged = manufactured_fda._stage_case_ecg_outputs(case_root, case)

            self.assertIn(source_dir / "ECG_manufactured_pseudoECG.dat", staged)
            self.assertFalse(raw.exists(), "raw root postProcessing ECG file should be renamed away")

    def test_registry_can_load_non_registered_case_folder_via_generic_spec(self) -> None:
        spec = load_tutorial_spec(
            "ECG",
            overrides={"tutorials_root": self.tutorials_root},
        )

        self.assertEqual(spec.name, "ECG")
        self.assertEqual(len(spec.build_cases()), 1)
        self.assertIn("Generic case runner", spec.metadata["notes"])
        self.assertIsNone(spec.collect_outputs)
        self.assertIsNone(spec.postprocess)

    def test_generic_case_overrides_update_electro_and_physics_properties(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            tutorials_root = Path(temp_dir)
            case_root = tutorials_root / "randomCase"
            (case_root / "constant").mkdir(parents=True, exist_ok=True)

            (case_root / "constant" / "electroProperties").write_text(
                "\n".join(
                    [
                        "electroModel singleCellElectro;",
                        "",
                        "singleCellElectroCoeffs",
                        "{",
                        "    ionicModel BuenoOrovio;",
                        "}",
                        "",
                    ]
                )
            )
            (case_root / "constant" / "physicsProperties").write_text("type electroModel;\n")

            spec = load_tutorial_spec(
                "randomCase",
                overrides={
                    "tutorials_root": tutorials_root,
                    "electro_property_overrides": {
                        "$ELECTRO_MODEL_COEFFS.ionicModel": "Gaur",
                    },
                    "physics_property_overrides": {
                        "type": "electroMechanicalModel",
                    },
                },
            )

            case = spec.build_cases()[0]
            spec.apply_case(spec.case_root, case)

            electro_text = (case_root / "constant" / "electroProperties").read_text()
            physics_text = (case_root / "constant" / "physicsProperties").read_text()
            self.assertIn("ionicModel    Gaur;", electro_text)
            self.assertIn("type    electroMechanicalModel;", physics_text)

    def test_generic_case_supports_multi_case_sweeps(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            tutorials_root = Path(temp_dir)
            case_root = tutorials_root / "ECG"
            (case_root / "constant").mkdir(parents=True, exist_ok=True)
            (case_root / "constant" / "electroProperties").write_text(
                "\n".join(
                    [
                        "electroModel monoDomainElectro;",
                        "",
                        "monoDomainElectroCoeffs",
                        "{",
                        "    ionicModel BuenoOrovio;",
                        "    tissue epicardialCells;",
                        "}",
                        "",
                    ]
                )
            )
            (case_root / "constant" / "physicsProperties").write_text("type electroModel;\n")

            spec = load_tutorial_spec(
                "genericCase",
                overrides={
                    "tutorials_root": tutorials_root,
                    "case_dir_name": "ECG",
                    "electro_property_overrides": {
                        "$ELECTRO_MODEL_COEFFS.tissue": "epicardialCells",
                    },
                    "cases": [
                        {
                            "case_id": "epi",
                            "electro_property_overrides": {
                                "$ELECTRO_MODEL_COEFFS.tissue": "epicardialCells",
                            },
                        },
                        {
                            "case_id": "endo",
                            "electro_property_overrides": {
                                "$ELECTRO_MODEL_COEFFS.tissue": "endocardialCells",
                                "$ELECTRO_MODEL_COEFFS.ECG.pseudoECGElectroCoeffs.electrodes.V1": "(1 2 3)",
                            },
                        },
                    ],
                },
            )

            cases = spec.build_cases()
            self.assertEqual([case.case_id for case in cases], ["epi", "endo"])
            self.assertEqual(
                cases[1].params["electro_property_overrides"][
                    "$ELECTRO_MODEL_COEFFS.ECG.pseudoECGElectroCoeffs.electrodes.V1"
                ],
                "(1 2 3)",
            )

    def test_curated_spec_applies_case_overrides_then_user_overrides(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            tutorials_root = Path(temp_dir)
            case_root = tutorials_root / "singleCell"
            (case_root / "constant").mkdir(parents=True, exist_ok=True)

            (case_root / "constant" / "electroProperties").write_text(
                "\n".join(
                    [
                        "electroModel singleCellElectro;",
                        "",
                        "singleCellElectroCoeffs",
                        "{",
                        "    ionicModel BuenoOrovio;",
                        "    tissue myocyte;",
                        "    singleCellStimulus",
                        "    {",
                        "        stim_amplitude 0.4;",
                        "    }",
                        "}",
                        "",
                    ]
                )
            )
            (case_root / "constant" / "physicsProperties").write_text("type electroModel;\n")

            spec = single_cell.make_spec(
                tutorials_root=tutorials_root,
                case_dir_name="singleCell",
                ionic_models=["TNNP"],
                ionic_model_tissue_map={"TNNP": ["epicardialCells"]},
                stimulus_map={"TNNP": 60.0},
                electro_property_overrides={
                    "singleCellElectroCoeffs.ionicModel": "Gaur",
                },
                physics_property_overrides={
                    "type": "electroMechanicalModel",
                },
            )

            case = spec.build_cases()[0]
            spec.apply_case(spec.case_root, case)

            electro_text = (case_root / "constant" / "electroProperties").read_text()
            physics_text = (case_root / "constant" / "physicsProperties").read_text()

            self.assertIn("ionicModel    Gaur;", electro_text)
            self.assertIn("tissue    epicardialCells;", electro_text)
            self.assertIn("stim_amplitude    60.0;", electro_text)
            self.assertIn("type    electroMechanicalModel;", physics_text)


if __name__ == "__main__":
    unittest.main()
