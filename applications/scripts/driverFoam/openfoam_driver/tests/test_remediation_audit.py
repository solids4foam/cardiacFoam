import json
from pathlib import Path
from openfoam_driver.core.runtime.remediation_audit import append_remediation_record


def test_append_writes_one_jsonl_record(tmp_path):
    append_remediation_record(
        tmp_path, step_id="solve", attempt=2,
        applied_overrides=[{"driver_path": "deltaT", "value": "0.0001"}],
        resulting_status="failed",
    )
    path = tmp_path / "remediation_history.jsonl"
    lines = path.read_text().splitlines()
    assert len(lines) == 1
    rec = json.loads(lines[0])
    assert rec["step_id"] == "solve"
    assert rec["attempt"] == 2
    assert rec["resulting_status"] == "failed"
    assert rec["applied_overrides"][0]["driver_path"] == "deltaT"
    assert "timestamp" in rec


def test_append_is_additive(tmp_path):
    for _ in range(3):
        append_remediation_record(
            tmp_path, step_id="solve", attempt=1, applied_overrides=[],
            resulting_status="ok",
        )
    path = tmp_path / "remediation_history.jsonl"
    assert len(path.read_text().splitlines()) == 3


def test_append_never_raises_on_bad_dir(tmp_path):
    # A non-existent nested output dir must not crash the rerun.
    append_remediation_record(
        tmp_path / "does" / "not" / "exist", step_id="s", attempt=1,
        applied_overrides=[], resulting_status="ok",
    )  # should silently no-op, not raise
