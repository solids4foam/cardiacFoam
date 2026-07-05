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
#     test_sweep_runner
#
# Description
#     Tests sweep_plan/sweep_run orchestration.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

import json
import subprocess
from pathlib import Path
from unittest import mock

import pytest

from openfoam_driver.core.runtime.sweep_runner import sweep_plan
from openfoam_driver.sweep_expansion import SweepValidationError


def _write_spec(path: Path, models=("TNNP", "BuenoOrovio")):
    spec = {
        "base": {
            "electro_selectors": {"myocardiumSolver": "singleCellSolver", "tissue": "epicardialCells"},
            "physics_selectors": {"type": "electroModel"},
        },
        "sweep": {
            "mode": "cross_product",
            "independent": {"ionicModel": list(models)},
            "dependent": [{"name": "caseId", "derive": "case_id_template", "of": ["ionicModel"]}],
        },
    }
    path.write_text(json.dumps(spec))
    return spec


def test_sweep_plan_materializes_and_audits_each_case_for_real(tmp_path):
    spec_path = tmp_path / "sweep.json"
    _write_spec(spec_path)
    output_dir = tmp_path / "out"

    result = sweep_plan(spec_path, output_dir=output_dir)

    assert result["case_count"] == 2
    assert {c["case_id"] for c in result["cases"]} == {"TNNP", "BuenoOrovio"}
    for case in result["cases"]:
        assert case["plan"]["status"] == "ok"
        assert (output_dir / case["case_id"] / "constant" / "electroProperties").exists()


def test_sweep_plan_refuses_over_cap_without_expanding(tmp_path):
    spec = {
        "base": {"electro_selectors": {"myocardiumSolver": "singleCellSolver", "tissue": "epicardialCells"},
                 "physics_selectors": {"type": "electroModel"}},
        "sweep": {"mode": "cross_product", "independent": {"a": list(range(20)), "b": list(range(20))}, "dependent": []},
    }
    spec_path = tmp_path / "sweep.json"
    spec_path.write_text(json.dumps(spec))

    with mock.patch("openfoam_driver.core.runtime.sweep_runner.materialize_case") as mock_materialize:
        with pytest.raises(SweepValidationError):
            sweep_plan(spec_path, output_dir=tmp_path / "out")
    mock_materialize.assert_not_called()


def test_sweep_plan_records_materialization_failure_and_continues(tmp_path):
    spec = {
        "base": {
            "electro_selectors": {"myocardiumSolver": "singleCellSolver", "tissue": "epicardialCells"},
            "physics_selectors": {"type": "electroModel"},
        },
        "sweep": {
            "mode": "cross_product",
            "independent": {"ionicModel": ["TNNP", "NotARealModel"]},
            "dependent": [{"name": "caseId", "derive": "case_id_template", "of": ["ionicModel"]}],
        },
    }
    spec_path = tmp_path / "sweep.json"
    spec_path.write_text(json.dumps(spec))

    result = sweep_plan(spec_path, output_dir=tmp_path / "out")

    by_id = {case["case_id"]: case for case in result["cases"]}
    assert by_id["TNNP"]["status"] == "ok"
    assert by_id["NotARealModel"]["status"] == "failed"
    assert "materialization_error" in by_id["NotARealModel"]


def test_sweep_plan_records_unrecognized_axis_as_per_case_failure(tmp_path):
    # route_case_values now raises SweepValidationError for an unrecognized
    # axis like "bogusAxis" (see sweep_routing.py fix). That per-case error
    # must be caught and recorded like any other materialization failure,
    # not propagate uncaught and crash the whole sweep_plan call -- a caller
    # sweeping N cases with one bad axis should still see a clean per-case
    # report, the same as an invalid ionicModel does today.
    spec = {
        "base": {
            "electro_selectors": {"myocardiumSolver": "singleCellSolver", "tissue": "myocyte"},
            "physics_selectors": {"type": "electroModel"},
        },
        "sweep": {
            "mode": "cross_product",
            "independent": {"bogusAxis": [0.5, 0.2]},
            "dependent": [{"name": "caseId", "derive": "case_id_template", "of": ["bogusAxis"]}],
        },
    }
    spec_path = tmp_path / "sweep.json"
    spec_path.write_text(json.dumps(spec))

    result = sweep_plan(spec_path, output_dir=tmp_path / "out")

    assert result["case_count"] == 2
    for case in result["cases"]:
        assert case["status"] == "failed"
        assert "bogusAxis" in case["materialization_error"]


def test_sweep_run_writes_run_documents_and_continues_past_failure(tmp_path):
    spec_path = tmp_path / "sweep.json"
    _write_spec(spec_path)
    output_dir = tmp_path / "out"

    call_log = []
    real_subprocess_run = subprocess.run  # captured before patching, so real internal
    # subprocess calls materialize_case makes (e.g. foamDictionary, when it's on PATH)
    # still execute for real instead of being swallowed by this fake.

    def fake_subprocess_run(cmd, **kwargs):
        if "--run-document" not in cmd:
            return real_subprocess_run(cmd, **kwargs)
        call_log.append(cmd)
        run_doc_path = Path(cmd[cmd.index("--run-document") + 1])
        run_doc = json.loads(run_doc_path.read_text())
        workflow_state_path = Path(run_doc["launch"]["outputDir"]) / "workflow_state.json"
        workflow_state_path.parent.mkdir(parents=True, exist_ok=True)
        failed = "BuenoOrovio" in str(run_doc_path.parent)
        state = {"status": "failed" if failed else "completed"}
        workflow_state_path.write_text(json.dumps(state))
        return mock.Mock(returncode=1 if failed else 0, stdout="", stderr="")

    with mock.patch("openfoam_driver.core.runtime.sweep_runner.subprocess.run", side_effect=fake_subprocess_run):
        from openfoam_driver.core.runtime.sweep_runner import sweep_run
        result = sweep_run(spec_path, output_dir=output_dir)

    assert len(call_log) == 2
    assert (output_dir / "TNNP" / "run_document.json").exists()
    assert (output_dir / "BuenoOrovio" / "run_document.json").exists()

    manifest_path = output_dir / "sweep_manifest.json"
    assert manifest_path.exists()
    from openfoam_driver.core.runtime.sweep_manifest import read_manifest
    manifest = read_manifest(manifest_path)
    statuses = {c.case_id: c.status for c in manifest.cases}
    assert statuses["TNNP"] == "completed"
    assert statuses["BuenoOrovio"] == "failed"
    state_paths = {c.case_id: c.workflow_state_path for c in manifest.cases}
    assert state_paths["TNNP"] == "TNNP/postProcessing/workflow_state.json"
    assert result["failed_count"] == 1
    assert result["completed_count"] == 1


def test_sweep_run_records_unrecognized_axis_as_per_case_failure(tmp_path):
    # Mirrors test_sweep_plan_records_unrecognized_axis_as_per_case_failure:
    # route_case_values's SweepValidationError must be caught per-case inside
    # sweep_run's loop too, not crash the whole call.
    spec = {
        "base": {
            "electro_selectors": {"myocardiumSolver": "singleCellSolver", "tissue": "myocyte"},
            "physics_selectors": {"type": "electroModel"},
        },
        "sweep": {
            "mode": "cross_product",
            "independent": {"bogusAxis": [0.5, 0.2]},
            "dependent": [{"name": "caseId", "derive": "case_id_template", "of": ["bogusAxis"]}],
        },
    }
    spec_path = tmp_path / "sweep.json"
    spec_path.write_text(json.dumps(spec))

    from openfoam_driver.core.runtime.sweep_runner import sweep_run
    with mock.patch("openfoam_driver.core.runtime.sweep_runner.subprocess.run") as mock_run:
        result = sweep_run(spec_path, output_dir=tmp_path / "out")

    mock_run.assert_not_called()
    assert result["failed_count"] == 2
    assert result["completed_count"] == 0
    for case in result["cases"]:
        assert case["status"] == "failed"
        assert "bogusAxis" in case["materialization_error"]


def test_sweep_run_refuses_over_cap_without_expanding(tmp_path):
    spec = {
        "base": {"electro_selectors": {"myocardiumSolver": "singleCellSolver", "tissue": "epicardialCells"},
                 "physics_selectors": {"type": "electroModel"}},
        "sweep": {"mode": "cross_product", "independent": {"a": list(range(20)), "b": list(range(20))}, "dependent": []},
    }
    spec_path = tmp_path / "sweep.json"
    spec_path.write_text(json.dumps(spec))

    with mock.patch("openfoam_driver.core.runtime.sweep_runner.materialize_case") as mock_materialize, \
         mock.patch("openfoam_driver.core.runtime.sweep_runner.subprocess.run") as mock_run:
        from openfoam_driver.core.runtime.sweep_runner import sweep_run
        with pytest.raises(SweepValidationError):
            sweep_run(spec_path, output_dir=tmp_path / "out")
    mock_materialize.assert_not_called()
    mock_run.assert_not_called()


def test_sweep_run_accepts_over_cap_with_explicit_override(tmp_path):
    spec = {
        "base": {"electro_selectors": {"myocardiumSolver": "singleCellSolver", "tissue": "epicardialCells"},
                 "physics_selectors": {"type": "electroModel"}},
        "sweep": {"mode": "cross_product", "independent": {"ionicModel": ["TNNP"] * 250}, "dependent": []},
    }
    spec_path = tmp_path / "sweep.json"
    spec_path.write_text(json.dumps(spec))
    output_dir = tmp_path / "out"

    def fake_materialize(*, case_dir, routed):
        case_dir.mkdir(parents=True, exist_ok=True)

    fake_report = mock.Mock()
    fake_report.status = "ok"
    fake_report.to_json.return_value = {
        "status": "ok",
        "run_document": {
            "version": "2",
            "launch": {"outputDir": str(output_dir / "dummy" / "postProcessing")},
        },
    }

    with mock.patch("openfoam_driver.core.runtime.sweep_runner.materialize_case", side_effect=fake_materialize), \
         mock.patch("openfoam_driver.core.runtime.sweep_runner.strict_plan", return_value=fake_report), \
         mock.patch("openfoam_driver.core.runtime.sweep_runner.subprocess.run"):
        from openfoam_driver.core.runtime.sweep_runner import sweep_run
        result = sweep_run(spec_path, output_dir=output_dir, max_cases=300)
    assert result["case_count"] == 250


def test_resume_skips_terminal_completed_case(tmp_path):
    spec_path = tmp_path / "sweep.json"
    _write_spec(spec_path, models=("TNNP",))
    spec = json.loads(spec_path.read_text())
    output_dir = tmp_path / "out"
    output_dir.mkdir()

    from openfoam_driver.core.runtime.sweep_manifest import (
        CaseManifestEntry, SweepManifest, compute_spec_hash, write_manifest,
    )
    case_dir = output_dir / "TNNP"
    state_dir = case_dir / "postProcessing"
    state_dir.mkdir(parents=True)
    (state_dir / "workflow_state.json").write_text('{"status": "completed"}')
    manifest = SweepManifest(
        schema_version="1.0", sweep_spec_hash=compute_spec_hash(spec),
        created_at="t0", updated_at="t0",
        cases=[CaseManifestEntry(
            case_id="TNNP", resolved_axis_values={"ionicModel": "TNNP"},
            override_hash="sha256:x", run_document_path="TNNP/run_document.json",
            workflow_state_path="TNNP/postProcessing/workflow_state.json",
            status="completed", outcome="fresh", started_at="t0", updated_at="t0",
        )],
    )
    write_manifest(output_dir / "sweep_manifest.json", manifest)

    with mock.patch("openfoam_driver.core.runtime.sweep_runner.materialize_case") as mock_materialize, \
         mock.patch("openfoam_driver.core.runtime.sweep_runner.subprocess.run") as mock_run:
        from openfoam_driver.core.runtime.sweep_runner import sweep_run
        result = sweep_run(spec_path, output_dir=output_dir)

    mock_materialize.assert_not_called()
    mock_run.assert_not_called()
    assert result["skipped_count"] == 1


def test_resume_leaves_terminal_failed_alone_without_retry_flag(tmp_path):
    spec_path = tmp_path / "sweep.json"
    _write_spec(spec_path, models=("TNNP",))
    spec = json.loads(spec_path.read_text())
    output_dir = tmp_path / "out"
    output_dir.mkdir()

    from openfoam_driver.core.runtime.sweep_manifest import (
        CaseManifestEntry, SweepManifest, compute_spec_hash, write_manifest,
    )
    case_dir = output_dir / "TNNP"
    state_dir = case_dir / "postProcessing"
    state_dir.mkdir(parents=True)
    (state_dir / "workflow_state.json").write_text('{"status": "failed"}')
    manifest = SweepManifest(
        schema_version="1.0", sweep_spec_hash=compute_spec_hash(spec),
        created_at="t0", updated_at="t0",
        cases=[CaseManifestEntry(
            case_id="TNNP", resolved_axis_values={"ionicModel": "TNNP"},
            override_hash="sha256:x", run_document_path="TNNP/run_document.json",
            workflow_state_path="TNNP/postProcessing/workflow_state.json",
            status="failed", outcome="fresh", started_at="t0", updated_at="t0",
        )],
    )
    write_manifest(output_dir / "sweep_manifest.json", manifest)

    with mock.patch("openfoam_driver.core.runtime.sweep_runner.materialize_case") as mock_materialize, \
         mock.patch("openfoam_driver.core.runtime.sweep_runner.subprocess.run") as mock_run:
        from openfoam_driver.core.runtime.sweep_runner import sweep_run
        result = sweep_run(spec_path, output_dir=output_dir, retry_failed=False)

    mock_materialize.assert_not_called()
    mock_run.assert_not_called()
    assert result["failed_count"] == 1


def test_resume_retries_terminal_failed_case_with_retry_flag(tmp_path):
    spec_path = tmp_path / "sweep.json"
    _write_spec(spec_path, models=("TNNP",))
    spec = json.loads(spec_path.read_text())
    output_dir = tmp_path / "out"
    output_dir.mkdir()

    from openfoam_driver.core.runtime.sweep_manifest import (
        CaseManifestEntry, SweepManifest, compute_spec_hash, write_manifest,
    )
    case_dir = output_dir / "TNNP"
    state_dir = case_dir / "postProcessing"
    state_dir.mkdir(parents=True)
    (state_dir / "workflow_state.json").write_text('{"status": "failed"}')
    manifest = SweepManifest(
        schema_version="1.0", sweep_spec_hash=compute_spec_hash(spec),
        created_at="t0", updated_at="t0",
        cases=[CaseManifestEntry(
            case_id="TNNP", resolved_axis_values={"ionicModel": "TNNP"},
            override_hash="sha256:x", run_document_path="TNNP/run_document.json",
            workflow_state_path="TNNP/postProcessing/workflow_state.json",
            status="failed", outcome="fresh", started_at="t0", updated_at="t0",
        )],
    )
    write_manifest(output_dir / "sweep_manifest.json", manifest)

    fake_report = mock.Mock()
    fake_report.status = "ok"
    fake_report.to_json.return_value = {
        "status": "ok",
        "run_document": {"version": "2", "launch": {"outputDir": str(state_dir)}},
    }

    def fake_subprocess_run(cmd, **kwargs):
        state_dir.mkdir(parents=True, exist_ok=True)
        (state_dir / "workflow_state.json").write_text('{"status": "completed"}')
        return mock.Mock(returncode=0, stdout="", stderr="")

    with mock.patch("openfoam_driver.core.runtime.sweep_runner.materialize_case") as mock_materialize, \
         mock.patch("openfoam_driver.core.runtime.sweep_runner.strict_plan", return_value=fake_report), \
         mock.patch("openfoam_driver.core.runtime.sweep_runner.subprocess.run", side_effect=fake_subprocess_run) as mock_run:
        from openfoam_driver.core.runtime.sweep_runner import sweep_run
        result = sweep_run(spec_path, output_dir=output_dir, retry_failed=True)

    mock_materialize.assert_called_once()
    mock_run.assert_called_once()
    assert result["completed_count"] == 1
    assert result["failed_count"] == 0
    by_id = {case["case_id"]: case for case in result["cases"]}
    assert by_id["TNNP"]["outcome"] == "retried"
    assert by_id["TNNP"]["status"] == "completed"


def test_sweep_run_case_timeout_marks_failed_and_continues(tmp_path):
    # A case whose run subprocess exceeds case_timeout_s must be recorded as a
    # per-case failure (not crash the whole sweep), and the timeout must be
    # passed through to subprocess.run.
    spec_path = tmp_path / "sweep.json"
    _write_spec(spec_path, models=("TNNP",))
    output_dir = tmp_path / "out"

    def fake_materialize(*, case_dir, routed):
        case_dir.mkdir(parents=True, exist_ok=True)

    fake_report = mock.Mock()
    fake_report.status = "ok"
    fake_report.to_json.return_value = {
        "status": "ok",
        "run_document": {
            "version": "2",
            "launch": {"outputDir": str(output_dir / "TNNP" / "postProcessing")},
        },
    }

    seen_kwargs = {}

    def fake_subprocess_run(cmd, **kwargs):
        seen_kwargs.update(kwargs)
        raise subprocess.TimeoutExpired(cmd, kwargs.get("timeout"))

    with mock.patch("openfoam_driver.core.runtime.sweep_runner.materialize_case", side_effect=fake_materialize), \
         mock.patch("openfoam_driver.core.runtime.sweep_runner.strict_plan", return_value=fake_report), \
         mock.patch("openfoam_driver.core.runtime.sweep_runner.subprocess.run", side_effect=fake_subprocess_run):
        from openfoam_driver.core.runtime.sweep_runner import sweep_run
        result = sweep_run(spec_path, output_dir=output_dir, case_timeout_s=0.01)

    assert seen_kwargs.get("timeout") == 0.01
    assert result["failed_count"] == 1
    assert result["completed_count"] == 0
    by_id = {case["case_id"]: case for case in result["cases"]}
    assert by_id["TNNP"]["status"] == "failed"
    assert "timeout" in by_id["TNNP"]["timeout_error"].lower()
    # sweep stayed resumable: manifest was still written
    assert (output_dir / "sweep_manifest.json").exists()


def test_spec_hash_mismatch_is_refused(tmp_path):
    spec_path = tmp_path / "sweep.json"
    _write_spec(spec_path, models=("TNNP",))
    output_dir = tmp_path / "out"
    output_dir.mkdir()
    from openfoam_driver.core.runtime.sweep_manifest import SweepManifest, write_manifest
    write_manifest(
        output_dir / "sweep_manifest.json",
        SweepManifest(schema_version="1.0", sweep_spec_hash="sha256:stale", created_at="t0", updated_at="t0", cases=[]),
    )
    from openfoam_driver.core.runtime.sweep_runner import sweep_run
    with pytest.raises(SweepValidationError, match="hash|spec changed"):
        sweep_run(spec_path, output_dir=output_dir)
