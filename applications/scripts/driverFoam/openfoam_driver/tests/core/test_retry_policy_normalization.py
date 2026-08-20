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
#     test_retry_policy_normalization
#
# Description
#     retry_policy normalization tests.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

from openfoam_driver.core.runtime.workflow import normalize_workflow_dag


_ABSENT = object()


def _normalize(retry_policy_marker, **step_extra):
    step = {"id": "solve", "command": "cardiacFoam"}
    if retry_policy_marker is not _ABSENT:
        step["retry_policy"] = retry_policy_marker
    step.update(step_extra)
    # No expected_artifacts here, so the context only has to be passed, not to
    # authorize anything -- these tests are about retry_policy normalization.
    dag, diagnostics = normalize_workflow_dag({"steps": [step]}, driver_context=None)
    return dag, diagnostics


def _codes(diagnostics):
    return {d.code for d in diagnostics}


def test_absent_retry_policy_normalizes_to_empty_dict():
    dag, diagnostics = _normalize(_ABSENT)
    assert dag["steps"][0]["retry_policy"] == {}
    assert "invalid_workflow_field" not in _codes(diagnostics)


def test_valid_retry_policy_preserved():
    dag, diagnostics = _normalize({"max_attempts": 3, "backoff_seconds": 2})
    assert dag["steps"][0]["retry_policy"] == {"max_attempts": 3, "backoff_seconds": 2}
    assert "invalid_workflow_field" not in _codes(diagnostics)


def test_non_dict_retry_policy_is_error_and_resets():
    dag, diagnostics = _normalize("nope")
    assert dag["steps"][0]["retry_policy"] == {}
    assert "invalid_workflow_field" in _codes(diagnostics)


def test_bad_max_attempts_dropped_with_diagnostic():
    for bad in ("3", 0, -1, True):
        dag, diagnostics = _normalize({"max_attempts": bad})
        assert "max_attempts" not in dag["steps"][0]["retry_policy"], bad
        assert "invalid_workflow_field" in _codes(diagnostics), bad


def test_bad_backoff_seconds_dropped_with_diagnostic():
    for bad in (-1, True, "1"):
        dag, diagnostics = _normalize({"backoff_seconds": bad})
        assert "backoff_seconds" not in dag["steps"][0]["retry_policy"], bad
        assert "invalid_workflow_field" in _codes(diagnostics), bad


def test_good_backoff_seconds_float_preserved():
    dag, diagnostics = _normalize({"backoff_seconds": 1.5})
    assert dag["steps"][0]["retry_policy"]["backoff_seconds"] == 1.5
    assert "invalid_workflow_field" not in _codes(diagnostics)


def test_zero_backoff_seconds_accepted():
    dag, diagnostics = _normalize({"backoff_seconds": 0})
    assert dag["steps"][0]["retry_policy"]["backoff_seconds"] == 0
    assert "invalid_workflow_field" not in _codes(diagnostics)


def test_missing_dag_message_names_contract_files():
    _, diagnostics = normalize_workflow_dag(None, driver_context=None)
    missing = [d for d in diagnostics if d.code == "missing_workflow_dag"]
    assert len(missing) == 1
    assert "Allrun" in missing[0].message
