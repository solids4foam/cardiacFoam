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
#     test_sweep_expansion
#
# Description
#     Tests sweep axis expansion (cross product / zip).
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

import pytest
from openfoam_driver.sweep_expansion import expand_sweep, compute_case_count, SweepValidationError


def test_cross_product_produces_all_combinations():
    spec = {
        "base": {},
        "sweep": {
            "mode": "cross_product",
            "independent": {"ionicModel": ["TNNP", "BuenoOrovio"], "deltaT": [1e-6, 2e-6]},
            "dependent": [],
        },
    }
    cases = expand_sweep(spec)
    assert len(cases) == 4
    combos = {(c.resolved_axis_values["ionicModel"], c.resolved_axis_values["deltaT"]) for c in cases}
    assert combos == {
        ("TNNP", 1e-6), ("TNNP", 2e-6),
        ("BuenoOrovio", 1e-6), ("BuenoOrovio", 2e-6),
    }


def test_zip_mode_pairs_by_position():
    spec = {
        "base": {},
        "sweep": {
            "mode": "zip",
            "independent": {"deltaT": [1e-6, 2e-6], "endTime": [0.1, 0.2]},
            "dependent": [],
        },
    }
    cases = expand_sweep(spec)
    assert len(cases) == 2
    pairs = [(c.resolved_axis_values["deltaT"], c.resolved_axis_values["endTime"]) for c in cases]
    assert pairs == [(1e-6, 0.1), (2e-6, 0.2)]


def test_zip_mode_rejects_mismatched_lengths():
    spec = {
        "base": {},
        "sweep": {
            "mode": "zip",
            "independent": {"deltaT": [1e-6, 2e-6], "endTime": [0.1, 0.2, 0.3]},
            "dependent": [],
        },
    }
    with pytest.raises(SweepValidationError, match="deltaT.*endTime|endTime.*deltaT"):
        expand_sweep(spec)


def test_compute_case_count_cross_product():
    spec = {"sweep": {"mode": "cross_product", "independent": {"a": [1, 2, 3], "b": [1, 2]}, "dependent": []}}
    assert compute_case_count(spec) == 6


def test_compute_case_count_zip():
    spec = {"sweep": {"mode": "zip", "independent": {"a": [1, 2, 3], "b": [4, 5, 6]}, "dependent": []}}
    assert compute_case_count(spec) == 3


def test_unknown_mode_is_rejected_before_expansion():
    spec = {"sweep": {"mode": "pairwise", "independent": {"a": [1]}, "dependent": []}}
    with pytest.raises(SweepValidationError, match="mode"):
        compute_case_count(spec)


def test_independent_axis_values_must_be_lists():
    spec = {"sweep": {"mode": "cross_product", "independent": {"a": "not-a-list"}, "dependent": []}}
    with pytest.raises(SweepValidationError, match="independent.*a"):
        compute_case_count(spec)


def test_empty_axis_is_rejected():
    spec = {"sweep": {"mode": "cross_product", "independent": {"a": []}, "dependent": []}}
    with pytest.raises(SweepValidationError, match="empty"):
        compute_case_count(spec)


from openfoam_driver.sweep_expansion import DEFAULT_MAX_CASES, check_case_count_cap


def test_check_case_count_cap_passes_within_default():
    spec = {"sweep": {"mode": "cross_product", "independent": {"a": [1, 2]}, "dependent": []}}
    check_case_count_cap(spec)  # must not raise


def test_check_case_count_cap_rejects_over_default():
    spec = {
        "sweep": {
            "mode": "cross_product",
            "independent": {"a": list(range(20)), "b": list(range(20))},
            "dependent": [],
        }
    }
    assert compute_case_count(spec) > DEFAULT_MAX_CASES
    with pytest.raises(SweepValidationError, match="max_cases"):
        check_case_count_cap(spec)


def test_check_case_count_cap_accepts_explicit_override():
    spec = {
        "sweep": {
            "mode": "cross_product",
            "independent": {"a": list(range(20)), "b": list(range(20))},
            "dependent": [],
        }
    }
    check_case_count_cap(spec, max_cases=500)  # must not raise
