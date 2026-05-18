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
        "version": "1",
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
