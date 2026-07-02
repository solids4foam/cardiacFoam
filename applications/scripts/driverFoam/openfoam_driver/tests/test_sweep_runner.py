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
