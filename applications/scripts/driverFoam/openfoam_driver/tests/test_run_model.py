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
