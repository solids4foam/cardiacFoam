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
#     test_failure_classification
#
# Description
#     Failure classification tests.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

from openfoam_driver.core.runtime.failure_classification import classify_failure
from openfoam_driver.core.runtime.workflow_state import WorkflowStepState


def _failed_step(*codes):
    return WorkflowStepState(
        step_id="solve",
        status="failed",
        attempt=1,
        command="cardiacFoam",
        args=(),
        cwd=".",
        exit_code=1,
        diagnostics=tuple({"level": "error", "code": c, "message": "x"} for c in codes),
    )


def test_timeout_is_retryable():
    assert classify_failure(_failed_step("workflow_step_timeout")) == "retryable"


def test_missing_artifacts_is_fatal():
    assert classify_failure(_failed_step("missing_artifacts")) == "fatal"


def test_exec_error_is_fatal():
    assert classify_failure(_failed_step("workflow_step_exec_error")) == "fatal"


def test_generic_failure_no_code_is_fatal():
    assert classify_failure(_failed_step()) == "fatal"


def test_override_forces_fatal_code_to_retryable():
    result = classify_failure(
        _failed_step("missing_artifacts"),
        overrides={"missing_artifacts": "retryable"},
    )
    assert result == "retryable"


def test_override_forces_retryable_code_to_fatal():
    result = classify_failure(
        _failed_step("workflow_step_timeout"),
        overrides={"workflow_step_timeout": "fatal"},
    )
    assert result == "fatal"
