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
from openfoam_driver.core.sweep.sweep_expansion import expand_sweep, compute_case_count, SweepValidationError


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


from openfoam_driver.core.sweep.sweep_expansion import DEFAULT_MAX_CASES, check_case_count_cap


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


def _fake_lookup(catalog):
    def lookup(name):
        if name not in catalog:
            raise SweepValidationError(f"Unknown derivation '{name}'")
        return catalog[name]
    return lookup


def test_dependent_derivation_folds_into_resolved_values():
    spec = {
        "base": {},
        "sweep": {
            "mode": "cross_product",
            "independent": {"ionicModel": ["TNNP"], "deltaT": [1e-6]},
            "dependent": [{"name": "caseId", "derive": "join", "of": ["ionicModel", "deltaT"]}],
        },
    }
    lookup = _fake_lookup({"join": lambda values: {"caseId": "-".join(str(v) for v in values.values())}})
    cases = expand_sweep(spec, get_derivation=lookup)
    assert cases[0].resolved_axis_values["caseId"] == "TNNP-1e-06"


def test_later_derivation_can_read_earlier_derivation_output():
    spec = {
        "base": {},
        "sweep": {
            "mode": "cross_product",
            "independent": {"deltaT": [1e-6]},
            "dependent": [
                {"name": "doubled", "derive": "double", "of": ["deltaT"]},
                {"name": "quadrupled", "derive": "double2", "of": ["doubled"]},
            ],
        },
    }

    def _double_dt(values):
        return {"doubled": values["deltaT"] * 2}

    def _double_doubled(values):
        return {"quadrupled": values["doubled"] * 2}

    lookup = _fake_lookup({"double": _double_dt, "double2": _double_doubled})
    cases = expand_sweep(spec, get_derivation=lookup)
    assert cases[0].resolved_axis_values["quadrupled"] == pytest.approx(4e-6)


def test_forward_reference_is_rejected():
    spec = {
        "base": {},
        "sweep": {
            "mode": "cross_product",
            "independent": {"deltaT": [1e-6]},
            "dependent": [
                {"name": "a", "derive": "noop", "of": ["b"]},  # "b" declared after "a"
                {"name": "b", "derive": "noop", "of": ["deltaT"]},
            ],
        },
    }
    lookup = _fake_lookup({"noop": lambda values: {}})
    with pytest.raises(SweepValidationError, match="forward|unknown|not yet"):
        expand_sweep(spec, get_derivation=lookup)


def test_unknown_of_reference_is_rejected():
    spec = {
        "base": {},
        "sweep": {
            "mode": "cross_product",
            "independent": {"deltaT": [1e-6]},
            "dependent": [{"name": "a", "derive": "noop", "of": ["notAnAxis"]}],
        },
    }
    lookup = _fake_lookup({"noop": lambda values: {}})
    with pytest.raises(SweepValidationError, match="notAnAxis"):
        expand_sweep(spec, get_derivation=lookup)


def test_derivation_key_collision_is_rejected():
    spec = {
        "base": {},
        "sweep": {
            "mode": "cross_product",
            "independent": {"deltaT": [1e-6]},
            "dependent": [{"name": "deltaT", "derive": "noop", "of": ["deltaT"]}],  # collides
        },
    }
    lookup = _fake_lookup({"noop": lambda values: {"deltaT": 999}})
    with pytest.raises(SweepValidationError, match="collis"):
        expand_sweep(spec, get_derivation=lookup)


def test_dependent_entries_require_lookup_function():
    spec = {
        "base": {},
        "sweep": {
            "mode": "cross_product",
            "independent": {"deltaT": [1e-6]},
            "dependent": [{"name": "caseId", "derive": "join", "of": ["deltaT"]}],
        },
    }
    with pytest.raises(SweepValidationError, match="get_derivation|lookup"):
        expand_sweep(spec)


def test_declared_case_id_becomes_case_id():
    spec = {
        "base": {},
        "sweep": {
            "mode": "cross_product",
            "independent": {"ionicModel": ["TNNP"], "deltaT": [1e-6]},
            "dependent": [{"name": "caseId", "derive": "join", "of": ["ionicModel", "deltaT"]}],
        },
    }
    lookup = _fake_lookup({"join": lambda values: {"caseId": "-".join(str(v) for v in values.values())}})
    cases = expand_sweep(spec, get_derivation=lookup)
    assert cases[0].case_id == "TNNP-1e-06"


def test_no_case_id_falls_back_to_ordinal():
    spec = {
        "base": {},
        "sweep": {"mode": "cross_product", "independent": {"ionicModel": ["TNNP", "BuenoOrovio"]}, "dependent": []},
    }
    cases = expand_sweep(spec)
    assert [c.case_id for c in cases] == ["case_0001", "case_0002"]


def test_non_unique_case_id_is_rejected():
    spec = {
        "base": {},
        "sweep": {
            "mode": "cross_product",
            # varies both ionicModel and deltaT, but caseId only encodes ionicModel -> collision
            "independent": {"ionicModel": ["TNNP", "TNNP"], "deltaT": [1e-6, 2e-6]},
            "dependent": [{"name": "caseId", "derive": "join_model_only", "of": ["ionicModel"]}],
        },
    }
    lookup = _fake_lookup({"join_model_only": lambda values: {"caseId": str(values["ionicModel"])}})
    with pytest.raises(SweepValidationError, match="not unique|collis|duplicate"):
        expand_sweep(spec, get_derivation=lookup)


def test_path_unsafe_case_id_is_rejected_before_materialization():
    spec = {
        "base": {},
        "sweep": {
            "mode": "cross_product",
            "independent": {"ionicModel": ["TNNP"]},
            "dependent": [{"name": "caseId", "derive": "bad_label", "of": ["ionicModel"]}],
        },
    }
    lookup = _fake_lookup({"bad_label": lambda values: {"caseId": "../bad"}})
    with pytest.raises(SweepValidationError, match="path-safe|caseId"):
        expand_sweep(spec, get_derivation=lookup)
