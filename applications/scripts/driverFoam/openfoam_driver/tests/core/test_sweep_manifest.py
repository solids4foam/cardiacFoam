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
#     test_sweep_manifest
#
# Description
#     Tests sweep manifest schema, hashing, and atomic I/O.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

import json
import tempfile
from pathlib import Path

from openfoam_driver.core.runtime.sweep_manifest import (
    CaseManifestEntry,
    SweepManifest,
    compute_spec_hash,
    compute_override_hash,
    write_manifest,
    read_manifest,
)


def test_compute_spec_hash_is_stable_for_equivalent_dicts():
    a = {"sweep": {"mode": "cross_product", "independent": {"x": [1, 2]}}}
    b = {"sweep": {"independent": {"x": [1, 2]}, "mode": "cross_product"}}
    assert compute_spec_hash(a) == compute_spec_hash(b)


def test_compute_spec_hash_differs_for_different_specs():
    a = {"sweep": {"mode": "cross_product", "independent": {"x": [1, 2]}}}
    b = {"sweep": {"mode": "cross_product", "independent": {"x": [1, 3]}}}
    assert compute_spec_hash(a) != compute_spec_hash(b)


def test_compute_override_hash_differs_by_value():
    assert compute_override_hash({"deltaT": 1e-6}) != compute_override_hash({"deltaT": 2e-6})


def test_write_then_read_manifest_round_trips(tmp_path):
    manifest = SweepManifest(
        schema_version="1.0",
        sweep_spec_hash="sha256:abc",
        created_at="2026-07-01T12:00:00Z",
        updated_at="2026-07-01T12:00:00Z",
        cases=[
            CaseManifestEntry(
                case_id="TNNP_1e-06",
                resolved_axis_values={"deltaT": 1e-6},
                override_hash="sha256:def",
                run_document_path="TNNP_1e-06/run_document.json",
                workflow_state_path="TNNP_1e-06/postProcessing/workflow_state.json",
                status="pending",
                outcome="fresh",
                started_at=None,
                updated_at="2026-07-01T12:00:00Z",
            )
        ],
    )
    manifest_path = tmp_path / "sweep_manifest.json"
    write_manifest(manifest_path, manifest)
    loaded = read_manifest(manifest_path)
    assert loaded.sweep_spec_hash == "sha256:abc"
    assert loaded.cases[0].case_id == "TNNP_1e-06"
    assert loaded.cases[0].resolved_axis_values == {"deltaT": 1e-6}


def test_case_manifest_entry_defaults_case_record_path_to_empty_string():
    entry = CaseManifestEntry(
        case_id="TNNP_1e-06",
        resolved_axis_values={"deltaT": 1e-6},
        override_hash="sha256:def",
        run_document_path="TNNP_1e-06/run_document.json",
        workflow_state_path="TNNP_1e-06/postProcessing/workflow_state.json",
        status="pending",
        outcome="fresh",
        started_at=None,
        updated_at="2026-07-01T12:00:00Z",
    )
    assert entry.case_record_path == ""


def test_write_then_read_manifest_round_trips_case_record_path(tmp_path):
    manifest = SweepManifest(
        schema_version="1.0",
        sweep_spec_hash="sha256:abc",
        created_at="2026-08-25T12:00:00Z",
        updated_at="2026-08-25T12:00:00Z",
        cases=[
            CaseManifestEntry(
                case_id="TNNP_1e-06",
                resolved_axis_values={"deltaT": 1e-6},
                override_hash="sha256:def",
                run_document_path="TNNP_1e-06/run_document.json",
                workflow_state_path="TNNP_1e-06/postProcessing/workflow_state.json",
                status="completed",
                outcome="fresh",
                started_at="2026-08-25T12:00:00Z",
                updated_at="2026-08-25T12:01:00Z",
                case_record_path="TNNP_1e-06/case_record.json",
            )
        ],
    )
    manifest_path = tmp_path / "sweep_manifest.json"
    write_manifest(manifest_path, manifest)
    loaded = read_manifest(manifest_path)
    assert loaded.cases[0].case_record_path == "TNNP_1e-06/case_record.json"


def test_read_manifest_defaults_case_record_path_for_pre_existing_manifest(tmp_path):
    # A manifest written before this field existed must still load -- a
    # resumed sweep from an older driverFOAM version cannot be forced to
    # regenerate everything just because one new field appeared.
    manifest_path = tmp_path / "sweep_manifest.json"
    manifest_path.write_text(json.dumps({
        "schema_version": "1.0",
        "sweep_spec_hash": "sha256:abc",
        "created_at": "t0",
        "updated_at": "t0",
        "cases": [{
            "case_id": "TNNP_1e-06",
            "resolved_axis_values": {"deltaT": 1e-6},
            "override_hash": "sha256:def",
            "run_document_path": "TNNP_1e-06/run_document.json",
            "workflow_state_path": "TNNP_1e-06/postProcessing/workflow_state.json",
            "status": "completed",
            "outcome": "fresh",
            "started_at": "t0",
            "updated_at": "t0",
        }],
    }))
    loaded = read_manifest(manifest_path)
    assert loaded.cases[0].case_record_path == ""


def test_write_manifest_uses_atomic_replace(tmp_path):
    manifest_path = tmp_path / "sweep_manifest.json"
    manifest = SweepManifest(
        schema_version="1.0", sweep_spec_hash="sha256:abc",
        created_at="t0", updated_at="t0", cases=[],
    )
    write_manifest(manifest_path, manifest)
    assert manifest_path.exists()
    assert not manifest_path.with_name(manifest_path.name + ".tmp").exists()
    data = json.loads(manifest_path.read_text())
    assert data["schema_version"] == "1.0"
