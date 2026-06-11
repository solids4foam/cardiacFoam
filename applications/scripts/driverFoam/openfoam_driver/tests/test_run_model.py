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
#     test_run_model
#
# Description
#     Tests run model logic and specification contracts.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

"""Tests for the Run document JSON Schema and ``RunDocument`` dataclass.

The JSON Schema in ``schemas/run-document.json`` is the single source of truth
for its shape, and ``RunDocument`` is the Python model used by validation
helpers and catalog exporters.
"""

from __future__ import annotations

import json
from pathlib import Path

import jsonschema
import pytest

from openfoam_driver.core.runtime.run_model import RunDocument

SCHEMA_PATH = (
    Path(__file__).resolve().parents[2] / "schemas" / "run-document.json"
)


@pytest.fixture
def schema():
    return json.loads(SCHEMA_PATH.read_text())


def _valid_run_dict():
    return {
        "version": "2",
        "id": "run-0001",
        "name": "demo",
        "createdAt": "2026-04-20T10:00:00Z",
        "lastModified": "2026-04-20T10:00:00Z",
        "status": "draft",
        "config": {
            "anatomy": {},
            "physics": {},
            "stimulus": {},
            "solver": {},
        },
        "validation": {},
        "resolvedEntry": None,
        "workflowDag": None,
        "workflowState": None,
        "launch": None,
        "expectedArtifacts": [],
        "terminalStatusValues": ["completed", "failed"],
    }


def test_schema_validates_minimal_valid_run(schema):
    jsonschema.validate(_valid_run_dict(), schema)  # does not raise


def test_schema_rejects_unknown_status(schema):
    bad = _valid_run_dict()
    bad["status"] = "nonsense"
    with pytest.raises(jsonschema.ValidationError):
        jsonschema.validate(bad, schema)


def test_run_document_round_trip():
    doc = RunDocument.from_json(_valid_run_dict())
    back = doc.to_json()
    assert back["id"] == "run-0001"
    assert back["config"]["anatomy"] == {}
    assert back["status"] == "draft"
    assert back["version"] == "2"


def test_schema_accepts_valid_heterogeneity_block(schema):
    doc = _valid_run_dict()
    doc["config"]["physics"] = {
        "myocardiumSolver": "monodomainSolver",
        "ionicModel": "BuenoOrovio",
        "tissue": "epicardialCells",
        "ionicHeterogeneity.field": "t",
        "ionicHeterogeneity.mode": "transmuralBands",
        "ionicHeterogeneity.endoMInterface": "0.3",
        "ionicHeterogeneity.mEpiInterface": "0.7",
        "ionicHeterogeneity.transitionMode": "blend",
        "ionicHeterogeneity.smoothing": "smoothstep",
    }
    jsonschema.validate(doc, schema)  # does not raise


def test_schema_rejects_unknown_heterogeneity_mode(schema):
    doc = _valid_run_dict()
    doc["config"]["physics"] = {"ionicHeterogeneity.mode": "bogusMode"}
    with pytest.raises(jsonschema.ValidationError):
        jsonschema.validate(doc, schema)


def test_schema_rejects_unknown_tissue(schema):
    doc = _valid_run_dict()
    doc["config"]["physics"] = {"tissue": "notATissue"}
    with pytest.raises(jsonschema.ValidationError):
        jsonschema.validate(doc, schema)


def test_schema_still_allows_unlisted_physics_keys(schema):
    # additionalProperties stays open: the dict-catalog has far more keys
    # than the schema enumerates.
    doc = _valid_run_dict()
    doc["config"]["physics"] = {"someUncataloguedKey": "x"}
    jsonschema.validate(doc, schema)  # does not raise


def test_schema_accepts_normalized_workflow_dag(schema):
    doc = _valid_run_dict()
    doc["workflowDag"] = {
        "schema_version": "1",
        "step_status_values": ["pending", "running", "completed", "failed", "skipped"],
        "steps": [
            {
                "id": "solve",
                "command": "cardiacFoam",
                "args": [],
                "cwd": ".",
                "depends_on": [],
                "produces": ["single_cell_trace"],
                "consumes": [],
                "retry_policy": {"max_attempts": 1},
                "command_display": "cardiacFoam",
            }
        ],
    }
    jsonschema.validate(doc, schema)  # does not raise


def test_schema_rejects_raw_workflow_step_shape(schema):
    doc = _valid_run_dict()
    doc["workflowDag"] = {
        "schema_version": "1",
        "step_status_values": ["pending", "running", "completed", "failed", "skipped"],
        "steps": [{"id": "solve", "command": "cardiacFoam"}],
    }
    with pytest.raises(jsonschema.ValidationError):
        jsonschema.validate(doc, schema)


def test_schema_accepts_initial_workflow_state(schema):
    doc = _valid_run_dict()
    doc["workflowState"] = {
        "status": "pending",
        "current_step_id": "solve",
        "completed_steps": [],
        "failed_step_id": None,
        "steps": [
            {
                "step_id": "solve",
                "status": "pending",
                "attempt": 0,
                "command": "cardiacFoam",
                "args": [],
                "cwd": ".",
                "started_at": None,
                "finished_at": None,
                "exit_code": None,
                "stdout_log": None,
                "stderr_log": None,
                "produced_artifacts": [],
                "diagnostics": [],
            }
        ],
    }
    jsonschema.validate(doc, schema)  # does not raise


def test_schema_rejects_unknown_workflow_state_status(schema):
    doc = _valid_run_dict()
    doc["workflowState"] = {
        "status": "waiting",
        "current_step_id": None,
        "completed_steps": [],
        "failed_step_id": None,
        "steps": [],
    }
    with pytest.raises(jsonschema.ValidationError):
        jsonschema.validate(doc, schema)


def test_run_document_rejects_implicit_v1_from_json():
    old = _valid_run_dict()
    old["version"] = "1"
    with pytest.raises(ValueError, match="migrate_v1"):
        RunDocument.from_json(old)


def test_run_document_migrates_v1_explicitly():
    old = _valid_run_dict()
    old["version"] = "1"
    old.pop("resolvedEntry")
    old.pop("workflowDag")
    old.pop("workflowState")
    old.pop("launch")
    old.pop("expectedArtifacts")
    old.pop("terminalStatusValues")

    doc = RunDocument.migrate_v1(old)
    payload = doc.to_json()
    assert payload["version"] == "2"
    assert payload["config"] == old["config"]
    assert payload["resolvedEntry"] is None
    assert payload["workflowState"] is None
