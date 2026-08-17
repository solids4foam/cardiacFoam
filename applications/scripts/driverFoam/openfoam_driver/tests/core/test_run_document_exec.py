"""Tests for the RunDocument execution adapter.

Happy-path execution (a valid config that passes validate_run + a DAG that
actually runs) is proven end-to-end in test_cli_run_document.py. These unit
tests cover loading, v1 migration, and the diagnostic behaviours that make a
document non-executable — none of which need a catalog-valid config.
"""
from __future__ import annotations

import json
import os
import tempfile
import unittest
from pathlib import Path
from unittest import mock

from openfoam_driver.core.runtime.run_document_exec import (
    build_execution_inputs,
    load_run_document,
)
from openfoam_driver.core.runtime.run_model import RunDocument


def _empty_config() -> dict:
    return {"anatomy": {}, "physics": {}, "stimulus": {}, "solver": {}}


def _make_runnable_case(root: Path) -> Path:
    """Create a minimal but _case_is_runnable-passing OpenFOAM case dir."""
    case = root / "case"
    (case / "constant").mkdir(parents=True)
    (case / "system").mkdir()
    for rel in (
        "constant/electroProperties",
        "constant/physicsProperties",
        "system/controlDict",
        "system/fvSchemes",
        "system/fvSolution",
    ):
        (case / rel).write_text("\n")
    return case


def _minimal_doc(**overrides) -> RunDocument:
    base = dict(
        id="d1",
        name="doc-one",
        status="planned",
        config=_empty_config(),
        launch={"caseRoot": "/tmp/case", "outputDir": "/tmp/case/output"},
        workflowDag={
            "schema_version": "1",
            "step_status_values": ["pending", "running", "completed", "failed", "skipped"],
            "steps": [
                {"id": "solve", "command": "Allrun", "args": [], "cwd": ".",
                 "depends_on": [], "produces": [], "consumes": [],
                 "retry_policy": {}, "command_display": "Allrun"},
            ],
        },
        expectedArtifacts=[],
    )
    base.update(overrides)
    return RunDocument(**base)


class TestLoadRunDocument(unittest.TestCase):
    def test_loads_a_v3_document(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / "run.json"
            path.write_text(json.dumps(_minimal_doc().to_json()))
            doc = load_run_document(path)
            self.assertEqual(doc.name, "doc-one")
            self.assertEqual(doc.version, "3")

    def test_migrates_a_v1_document(self) -> None:
        v1 = {
            "version": "1",
            "id": "old",
            "name": "legacy",
            "status": "draft",
            "config": _empty_config(),
        }
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / "run.json"
            path.write_text(json.dumps(v1))
            doc = load_run_document(path)
            self.assertEqual(doc.version, "3")
            self.assertEqual(doc.name, "legacy")

    def test_v2_document_is_rejected(self) -> None:
        v2 = {
            "version": "2",
            "id": "old",
            "name": "archived",
            "status": "draft",
            "config": _empty_config(),
        }
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / "run.json"
            path.write_text(json.dumps(v2))
            with self.assertRaises(ValueError):
                load_run_document(path)

    def test_non_object_json_raises(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / "run.json"
            path.write_text("[]")
            with self.assertRaises(ValueError):
                load_run_document(path)


class TestBuildExecutionInputsDiagnostics(unittest.TestCase):
    def test_missing_workflow_dag_is_not_executable(self) -> None:
        doc = _minimal_doc(workflowDag=None)
        inputs, diagnostics = build_execution_inputs(doc)
        self.assertIsNone(inputs)
        codes = {d["code"] for d in diagnostics}
        self.assertIn("missing_workflow_dag", codes)

    def test_missing_launch_paths_are_not_executable(self) -> None:
        doc = _minimal_doc(launch={})
        inputs, diagnostics = build_execution_inputs(doc)
        self.assertIsNone(inputs)
        codes = {d["code"] for d in diagnostics}
        self.assertIn("missing_case_root", codes)
        self.assertIn("missing_output_dir", codes)

    def test_unknown_command_in_dag_is_rejected(self) -> None:
        doc = _minimal_doc(workflowDag={
            "steps": [{"id": "s", "command": "rm", "args": [], "cwd": ".",
                       "depends_on": [], "produces": [], "consumes": []}],
        })
        inputs, diagnostics = build_execution_inputs(doc)
        self.assertIsNone(inputs)
        codes = {d["code"] for d in diagnostics}
        self.assertIn("unknown_workflow_command", codes)

    def test_invalid_expected_artifact_is_reported(self) -> None:
        doc = _minimal_doc(expectedArtifacts=[
            {"artifact_id": "x", "path_pattern": "a/{bogus}.dat", "format": "csv_probe"},
        ])
        inputs, diagnostics = build_execution_inputs(doc)
        self.assertIsNone(inputs)
        codes = {d["code"] for d in diagnostics}
        self.assertIn("invalid_expected_artifact", codes)

    def test_non_dict_launch_does_not_raise(self) -> None:
        doc = _minimal_doc(launch="not-a-dict")
        inputs, diagnostics = build_execution_inputs(doc)  # must not raise
        self.assertIsNone(inputs)
        codes = {d["code"] for d in diagnostics}
        self.assertIn("invalid_launch", codes)

    def test_non_iterable_expected_artifacts_does_not_raise(self) -> None:
        doc = _minimal_doc(expectedArtifacts=42)
        inputs, diagnostics = build_execution_inputs(doc)  # must not raise
        self.assertIsNone(inputs)
        codes = {d["code"] for d in diagnostics}
        self.assertIn("invalid_expected_artifacts", codes)

    def test_malformed_workflow_state_is_reported(self) -> None:
        doc = _minimal_doc(workflowState={"bogus": "shape"})
        inputs, diagnostics = build_execution_inputs(doc)  # must not raise
        self.assertIsNone(inputs)
        codes = {d["code"] for d in diagnostics}
        self.assertIn("invalid_workflow_state", codes)

    def test_every_diagnostic_has_required_keys(self) -> None:
        doc = _minimal_doc(workflowDag=None, launch={})
        _inputs, diagnostics = build_execution_inputs(doc)
        for d in diagnostics:
            self.assertIn("level", d)
            self.assertIn("code", d)
            self.assertIn("message", d)
            self.assertIn("field", d)


class TestCaseRootValidation(unittest.TestCase):
    def test_nonexistent_case_root_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            missing = Path(temp) / "does-not-exist"
            doc = _minimal_doc(launch={
                "caseRoot": str(missing),
                "outputDir": str(missing / "out"),
            })
            inputs, diagnostics = build_execution_inputs(doc)
            self.assertIsNone(inputs)
            codes = {d["code"] for d in diagnostics}
            self.assertIn("case_root_missing", codes)

    def test_non_runnable_case_root_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            empty = Path(temp) / "empty"
            empty.mkdir()
            doc = _minimal_doc(launch={
                "caseRoot": str(empty),
                "outputDir": str(empty / "out"),
            })
            inputs, diagnostics = build_execution_inputs(doc)
            self.assertIsNone(inputs)
            codes = {d["code"] for d in diagnostics}
            self.assertIn("case_root_not_a_runnable_case", codes)

    def test_canonical_paths_are_resolved_and_stored(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            case = _make_runnable_case(Path(temp))
            doc = _minimal_doc(launch={
                "caseRoot": str(case),
                "outputDir": "output",  # relative -> resolves under caseRoot
            })
            with mock.patch(
                "openfoam_driver.core.runtime.run_document_exec.validate_run",
                return_value=[],
            ):
                inputs, diagnostics = build_execution_inputs(doc)
            self.assertIsNotNone(inputs, diagnostics)
            self.assertEqual(inputs.case_root, case.resolve())
            self.assertEqual(inputs.output_dir, (case.resolve() / "output"))


class TestAllowedRunsRoot(unittest.TestCase):
    def test_case_root_outside_allowed_root_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            allowed = root / "allowed"
            allowed.mkdir()
            case = _make_runnable_case(root)  # under root, NOT under allowed
            doc = _minimal_doc(launch={
                "caseRoot": str(case),
                "outputDir": str(case / "output"),
            })
            with mock.patch.dict(os.environ, {"DRIVERFOAM_ALLOWED_RUNS_ROOT": str(allowed)}), \
                 mock.patch(
                     "openfoam_driver.core.runtime.run_document_exec.validate_run",
                     return_value=[],
                 ):
                inputs, diagnostics = build_execution_inputs(doc)
            self.assertIsNone(inputs)
            codes = {d["code"] for d in diagnostics}
            self.assertIn("case_root_outside_allowed_root", codes)

    def test_output_dir_outside_allowed_root_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            allowed = root / "allowed"
            allowed.mkdir()
            case = _make_runnable_case(allowed)  # case under allowed
            outside_out = root / "elsewhere"      # outputDir NOT under allowed
            doc = _minimal_doc(launch={
                "caseRoot": str(case),
                "outputDir": str(outside_out),
            })
            with mock.patch.dict(os.environ, {"DRIVERFOAM_ALLOWED_RUNS_ROOT": str(allowed)}), \
                 mock.patch(
                     "openfoam_driver.core.runtime.run_document_exec.validate_run",
                     return_value=[],
                 ):
                inputs, diagnostics = build_execution_inputs(doc)
            self.assertIsNone(inputs)
            codes = {d["code"] for d in diagnostics}
            self.assertIn("output_dir_outside_allowed_root", codes)

    def test_both_inside_allowed_root_is_accepted(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            allowed = Path(temp)  # everything under the temp dir
            case = _make_runnable_case(allowed)
            doc = _minimal_doc(launch={
                "caseRoot": str(case),
                "outputDir": str(case / "output"),
            })
            with mock.patch.dict(os.environ, {"DRIVERFOAM_ALLOWED_RUNS_ROOT": str(allowed)}), \
                 mock.patch(
                     "openfoam_driver.core.runtime.run_document_exec.validate_run",
                     return_value=[],
                 ):
                inputs, diagnostics = build_execution_inputs(doc)
            self.assertIsNotNone(inputs, diagnostics)

    def test_unset_allowed_root_permits_separate_output_dir(self) -> None:
        # No DRIVERFOAM_ALLOWED_RUNS_ROOT: an absolute outputDir outside the case
        # is allowed (matches resolve_spec_paths separate-results-dir layout).
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            case = _make_runnable_case(root)
            separate_out = root / "results"
            doc = _minimal_doc(launch={
                "caseRoot": str(case),
                "outputDir": str(separate_out),
            })
            # patch.dict (clear=False) snapshots + restores os.environ on exit,
            # so popping the var here is safe and does not disturb PATH etc.
            with mock.patch.dict(os.environ, {}), \
                 mock.patch(
                     "openfoam_driver.core.runtime.run_document_exec.validate_run",
                     return_value=[],
                 ):
                os.environ.pop("DRIVERFOAM_ALLOWED_RUNS_ROOT", None)
                inputs, diagnostics = build_execution_inputs(doc)
            self.assertIsNotNone(inputs, diagnostics)
            self.assertEqual(inputs.output_dir, separate_out.resolve())


if __name__ == "__main__":
    unittest.main()
