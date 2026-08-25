"""P2.1: one launch predicate, truth-table tested across every state category
(structural, environment, warning, execution) the roadmap requires covered.

The real states exercised here come from what the CLI already computes:
- ``StrictPlanReport.status`` ("ok"/"failed") -- structural/semantic plan
  validity, derived from ``plan_diagnostics`` only (strict_planning.py).
- ``StrictPlanReport.environment_diagnostics`` -- execution-readiness
  diagnostics (level "error" or "warning"), excluded from ``.status`` by
  design so a sampled-field/environment warning can never fail a plan.
- ``WorkflowRunState``/``WorkflowStepState.status`` -- post-execution
  outcome ("completed" is the only success terminal state; see
  ``cli._terminal_status_label``).
"""
from __future__ import annotations

import pytest

from openfoam_driver.core.runtime.launch_readiness import (
    LaunchReadiness,
    is_execution_successful,
    is_launchable,
)
from openfoam_driver.core.planning_types import StrictDiagnostic


def _env_diagnostics(*levels: str) -> tuple[StrictDiagnostic, ...]:
    return tuple(
        StrictDiagnostic(level=level, code=f"fake_{level}", message="fake diagnostic")
        for level in levels
    )


@pytest.mark.parametrize(
    "plan_status, environment_diagnostics, expected",
    [
        ("ok", (), True),
        ("ok", _env_diagnostics("warning"), True),        # warnings alone never block launch
        ("ok", _env_diagnostics("error"), False),          # environment errors block launch
        ("failed", (), False),                              # structural errors block launch
        ("failed", _env_diagnostics("error"), False),      # both bad
    ],
)
def test_is_launchable_truth_table(plan_status, environment_diagnostics, expected) -> None:
    readiness = is_launchable(
        plan_status=plan_status,
        environment_diagnostics=environment_diagnostics,
    )
    assert isinstance(readiness, LaunchReadiness)
    assert readiness.launchable == expected
    assert readiness.structural_ok == (plan_status == "ok")
    assert readiness.environment_ok == (
        not any(d.level == "error" for d in environment_diagnostics)
    )


def test_is_launchable_has_warnings_reflects_warning_level_diagnostics() -> None:
    readiness = is_launchable(
        plan_status="ok",
        environment_diagnostics=_env_diagnostics("warning"),
    )
    assert readiness.has_warnings is True
    assert readiness.launchable is True

    readiness_no_warnings = is_launchable(plan_status="ok", environment_diagnostics=())
    assert readiness_no_warnings.has_warnings is False


def test_is_launchable_blocking_reason_prefers_structural_over_environment() -> None:
    # Structural failure should be reported even when environment also has
    # errors -- the CLI's --entry gate short-circuits before the environment
    # is even checked, so structural is the "first" failure in practice.
    readiness = is_launchable(
        plan_status="failed",
        environment_diagnostics=_env_diagnostics("error"),
    )
    assert readiness.launchable is False
    assert readiness.blocking_reason is not None
    assert "structur" in readiness.blocking_reason.lower()


def test_is_launchable_blocking_reason_environment_only() -> None:
    readiness = is_launchable(
        plan_status="ok",
        environment_diagnostics=_env_diagnostics("error"),
    )
    assert readiness.launchable is False
    assert readiness.blocking_reason is not None
    assert "environment" in readiness.blocking_reason.lower()


def test_is_launchable_no_blocking_reason_when_launchable() -> None:
    readiness = is_launchable(plan_status="ok", environment_diagnostics=())
    assert readiness.launchable is True
    assert readiness.blocking_reason is None


@pytest.mark.parametrize(
    "workflow_status, expected",
    [
        ("completed", True),
        ("pending", False),
        ("running", False),
        ("failed", False),
        ("skipped", False),
    ],
)
def test_is_execution_successful_matches_strict_success_contract(workflow_status, expected) -> None:
    # Mirrors WORKFLOW_STATE_STATUS_VALUES in core/runtime/workflow.py: only
    # "completed" is a success terminal state at this boundary.
    assert is_execution_successful(workflow_status) is expected
