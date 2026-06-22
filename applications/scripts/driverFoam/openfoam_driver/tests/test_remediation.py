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


from openfoam_driver.core.runtime.remediation import interpret_log_signatures


def test_log_layer_suggests_halving_deltat_on_recognized_signature():
    hints = interpret_log_signatures(
        _fc([], stderr_tail="--> FOAM FATAL ERROR:\nMaximum number of iterations exceeded")
    )
    assert len(hints) == 1
    h = hints[0]
    assert h.driver_path == "deltaT"
    assert h.change == "halve"
    assert h.source == "log_signature"


def test_log_layer_returns_empty_on_unrecognized_tail():
    # No recognized divergence marker -> the driver recognizes, it does not guess.
    hints = interpret_log_signatures(_fc([], stdout_tail="some unrecognized solver output"))
    assert hints == ()


def test_ladder_prefers_static_and_skips_log_layer():
    # missing_artifacts is a coded (static) failure: the log layer must NOT add a deltaT hint.
    hints = build_candidate_remediations(
        _fc(
            [{"level": "error", "code": "missing_artifacts", "message": "x", "field": "solve"}],
            stderr_tail="FOAM FATAL ERROR",
        )
    )
    assert [h.source for h in hints] == ["static"]
    assert all(h.driver_path != "deltaT" for h in hints)


def test_ladder_falls_back_to_log_layer_when_diagnostics_empty():
    # Bare nonzero divergence (no diagnostic code) + recognized signature -> log layer fires.
    hints = build_candidate_remediations(_fc([], stderr_tail="FOAM FATAL ERROR: singularity"))
    assert [h.source for h in hints] == ["log_signature"]
    assert hints[0].driver_path == "deltaT"


def test_ladder_skips_log_layer_for_coded_failure_without_static_hint():
    # workflow_step_timeout is coded but has no static hint; the log layer must NOT fire,
    # so a timeout never gets a (counter-productive) "halve deltaT" suggestion even when its
    # tail happens to contain a divergence-looking marker.
    hints = build_candidate_remediations(
        _fc(
            [{"level": "error", "code": "workflow_step_timeout", "message": "x", "field": "solve"}],
            stderr_tail="FOAM FATAL ERROR",
        )
    )
    assert hints == ()
