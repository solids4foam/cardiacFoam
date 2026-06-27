from openfoam_driver.core.runtime.remediation import (
    RemediationHint,
    build_candidate_remediations,
)


def _fc(diagnostics, stdout_tail="", stderr_tail=""):
    """Minimal failure_context dict, matching build_failure_context output."""
    return {
        "step_id": "solve",
        "attempt": 1,
        "exit_code": 1,
        "diagnostics": list(diagnostics),
        "stdout_log": None,
        "stderr_log": None,
        "stdout_tail": stdout_tail,
        "stderr_tail": stderr_tail,
        "stdout_truncated": False,
        "stderr_truncated": False,
    }


def test_missing_artifacts_yields_advisory_hint():
    hints = build_candidate_remediations(
        _fc([{"level": "error", "code": "missing_artifacts", "message": "x", "field": "solve"}])
    )
    assert len(hints) == 1
    h = hints[0]
    assert isinstance(h, RemediationHint)
    assert h.diagnostic_code == "missing_artifacts"
    assert h.driver_path == ""        # advisory: no mutation
    assert h.change == ""
    assert h.rationale                 # non-empty
    assert h.source == "static"


def test_exec_error_yields_advisory_hint():
    hints = build_candidate_remediations(
        _fc([{"level": "error", "code": "workflow_step_exec_error", "message": "x", "field": "solve"}])
    )
    assert [h.diagnostic_code for h in hints] == ["workflow_step_exec_error"]
    assert hints[0].driver_path == ""


def test_to_json_round_trips_fields():
    h = RemediationHint(
        diagnostic_code="missing_artifacts",
        driver_path="",
        change="",
        rationale="check producer",
        source="static",
        confidence="low",
    )
    assert h.to_json() == {
        "diagnostic_code": "missing_artifacts",
        "driver_path": "",
        "change": "",
        "rationale": "check producer",
        "source": "static",
        "confidence": "low",
    }






