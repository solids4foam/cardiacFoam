"""Tests for the RunDocument execution adapter (Task 3).

Happy-path execution (a valid config that passes validate_run + a DAG that
actually runs) is proven end-to-end in test_cli_run_document.py. These unit
tests cover loading, v1 migration, and the diagnostic behaviours that make a
document non-executable — none of which need a catalog-valid config.
"""
from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

from openfoam_driver.core.runtime.run_document_exec import (
    build_execution_inputs,
    load_run_document,
)
from openfoam_driver.core.runtime.run_model import RunDocument


def _empty_config() -> dict:
    return {"anatomy": {}, "physics": {}, "stimulus": {}, "solver": {}}


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
    def test_loads_a_v2_document(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / "run.json"
            path.write_text(json.dumps(_minimal_doc().to_json()))
            doc = load_run_document(path)
            self.assertEqual(doc.name, "doc-one")
            self.assertEqual(doc.version, "2")

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
            self.assertEqual(doc.version, "2")
            self.assertEqual(doc.name, "legacy")

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


if __name__ == "__main__":
    unittest.main()
