"""Stable semantic snapshots for every registered cardiac tutorial.

The fixture stores hashes instead of duplicating large OpenFOAM dictionaries
and capability reports.  Each hash is computed from canonical JSON after
repository roots and interpreter paths are sanitized.  Volatile timestamps
and environment-preflight diagnostics are deliberately outside this contract.
"""

from __future__ import annotations

import hashlib
import json
import sys
from copy import deepcopy
from pathlib import Path
from typing import Any

import pytest

from openfoam_driver.core.plugin_interface import default_driver_context
from openfoam_driver.core.runtime.registry import load_entry_spec
from openfoam_driver.core.strict_planning import strict_plan
from openfoam_driver.tests.conftest import monorepo_root, skip_without_monorepo

pytestmark = skip_without_monorepo

_FIXTURE = Path(__file__).parents[1] / "fixtures" / "cardiac_tutorial_characterization.json"
_DICTIONARY_PATHS = (
    "constant/electroProperties",
    "constant/physicsProperties",
    "system/controlDict",
    "system/fvSchemes",
    "system/fvSolution",
)


def _sanitize(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(key): _sanitize(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_sanitize(item) for item in value]
    if isinstance(value, str):
        # sys.executable FIRST: a project-local virtualenv lives inside
        # monorepo_root, so replacing the repo prefix first mangles the
        # interpreter path and the exact-match replace below then finds
        # nothing. Always substitute the more specific token first.
        value = value.replace(sys.executable, "<python>")
        value = value.replace(str(monorepo_root), "<repo>")
        return value
    return value


def _digest(value: Any) -> str:
    payload = json.dumps(
        _sanitize(value),
        ensure_ascii=False,
        separators=(",", ":"),
        sort_keys=True,
    ).encode()
    return "sha256:" + hashlib.sha256(payload).hexdigest()


def _characterize(name: str) -> dict[str, Any]:
    context = default_driver_context()
    spec = load_entry_spec(name, driver_context=context)
    report = strict_plan(name, driver_context=context).to_json()

    run_document = deepcopy(report["run_document"])
    run_document["createdAt"] = "<volatile>"
    run_document["lastModified"] = "<volatile>"
    # Environment availability and executable age are host properties.  Keep
    # the semantic validation status while snapshotting strict diagnostics in
    # their dedicated, environment-independent section below.
    run_document["validation"] = {
        "status": run_document["validation"]["status"],
        "diagnostics": "<environment-independent-section>",
    }

    diagnostics = {
        key: report[key]
        for key in (
            "validation_diagnostics",
            "workflow_diagnostics",
            "catalog_coverage_errors",
            "artifact_diagnostics",
            "mesh_geometry_diagnostics",
            "function_object_diagnostics",
        )
    }
    dictionary_contents = {
        relative: (spec.case_root / relative).read_text()
        for relative in _DICTIONARY_PATHS
        if (spec.case_root / relative).is_file()
    }

    sections = {
        "entry_resolution": report["resolved_entry"],
        "workflow_dag": report["workflow_dag"],
        "strict_diagnostics": diagnostics,
        "run_document_v2": run_document,
        "capability_report": report["capability_manifest"],
        "artifacts": report["expected_artifacts"],
        "dictionary_contents": dictionary_contents,
    }
    return {
        "case_count": len(spec.build_cases()),
        "dictionary_paths": sorted(dictionary_contents),
        "section_digests": {
            section: _digest(payload) for section, payload in sections.items()
        },
    }


def _current_characterization() -> dict[str, Any]:
    context = default_driver_context()
    registered = context.capabilities.tutorials.catalog()["registered_tutorials"]
    return {name: _characterize(name) for name in registered}


def test_all_registered_cardiac_tutorials_match_characterization_fixture() -> None:
    expected = json.loads(_FIXTURE.read_text())
    actual = _current_characterization()

    assert set(actual) == set(expected["tutorials"])
    assert actual == expected["tutorials"]

