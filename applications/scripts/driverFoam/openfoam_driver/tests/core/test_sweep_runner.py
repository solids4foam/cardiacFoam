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

from openfoam_driver.core.runtime.sweep_runner import (
    _stage_entry_case,
    sweep_plan,
    sweep_run,
)
from openfoam_driver.core.sweep.sweep_expansion import SweepValidationError


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


def _write_entry_spec(path, entry="niederer2012", values=(0.5, 0.2)):
    spec = {
        "base": {"entry": entry},
        "sweep": {
            "mode": "cross_product",
            "independent": {"dx_values": [[v] for v in values]},
            "dependent": [{"name": "caseId", "derive": "output_dir_name_template", "of": ["dx_values"]}],
        },
    }
    path.write_text(json.dumps(spec))
    return spec


def test_entry_case_staging_keeps_authored_case_clean(tmp_path):
    source = tmp_path / "tutorials" / "case"
    source.mkdir(parents=True)
    (source / "system").mkdir()
    (source / "system" / "controlDict").write_text("endTime 0.2;\n")
    (source / "0").mkdir()
    (source / "0" / "Vm").write_text("initial field")
    (source / "postProcessing").mkdir()
    (source / "postProcessing" / "old.dat").write_text("stale")
    (source / "processor0").mkdir()
    (source / "processor0" / "old").write_text("stale")
    (source / "0.2").mkdir()
    (source / "0.2" / "Vm").write_text("stale")
    (source / "workflow_state.json").write_text("{}")
    generated_case = source / "gauss_linear_40_manufacturedVerifier"
    (generated_case / "workflow_logs").mkdir(parents=True)
    (generated_case / "system").mkdir()
    (generated_case / "system" / "controlDict").write_text("generated")

    staged = tmp_path / "scratch" / "case_0001"
    _stage_entry_case(source, staged)

    assert (staged / "system" / "controlDict").read_text() == "endTime 0.2;\n"
    assert (staged / "0" / "Vm").exists()
    assert not (staged / "postProcessing").exists()
    assert not (staged / "processor0").exists()
    assert not (staged / "0.2").exists()
    assert not (staged / "workflow_state.json").exists()
    assert not (staged / generated_case.name).exists()
    assert (source / "postProcessing" / "old.dat").exists()


def test_sweep_plan_entry_mode_materializes_via_apply_case_and_audits(tmp_path):
    # Entry-based sweeps target an existing registered tutorial whose
    # apply_case()/build_cases() mutate its own shared case_root in place
    # (confirmed empirically for niederer2012 -- it is not a from-scratch
    # case_folder). sweep_plan must call spec.build_cases() + spec.apply_case()
    # directly instead of materialize_case()/build_and_launch, then audit via
    # strict_plan with the same routed overrides.
    spec_path = tmp_path / "sweep.json"
    _write_entry_spec(spec_path)

    fake_case_config = mock.Mock(case_id="implicit_TNNP_DX0.5")
    fake_spec = mock.Mock()
    fake_spec.case_root = tmp_path / "case_root"
    fake_spec.build_cases.return_value = [fake_case_config]

    fake_report = mock.Mock()
    fake_report.status = "ok"
    fake_report.to_json.return_value = {
        "status": "ok",
        "run_document": {"version": "3", "launch": {"outputDir": str(tmp_path / "out")}},
    }

    with mock.patch("openfoam_driver.core.runtime.sweep_runner.load_entry_spec", return_value=fake_spec) as mock_load, \
         mock.patch("openfoam_driver.core.runtime.sweep_runner.strict_plan", return_value=fake_report) as mock_strict_plan, \
         mock.patch("openfoam_driver.core.runtime.sweep_runner.materialize_case") as mock_materialize:
        result = sweep_plan(spec_path, output_dir=tmp_path / "out")

    mock_materialize.assert_not_called()
    assert mock_load.call_count == 2
    for call in mock_load.call_args_list:
        args, kwargs = call
        assert args[0] == "niederer2012"
        assert "dx_values" in kwargs["overrides"]
        assert "caseId" not in kwargs["overrides"]
    fake_spec.apply_case.assert_has_calls(
        [mock.call(fake_spec.case_root, fake_case_config)] * 2
    )
    assert mock_strict_plan.call_count == 2
    assert result["case_count"] == 2
    for case in result["cases"]:
        assert case["status"] == "ok"


def test_sweep_plan_entry_mode_rejects_axis_combination_resolving_to_multiple_cases(tmp_path):
    # sweep-run's per-axis-combination model assumes exactly one case per
    # resolved combination (see route_entry_case_values docstring); a
    # combination that still fans out inside the tutorial's own build_cases()
    # (e.g. missing a constraining kwarg like "solvers") must fail loudly as
    # a per-case error, not silently apply_case() only the first of several.
    spec_path = tmp_path / "sweep.json"
    _write_entry_spec(spec_path, values=(0.5,))

    fake_spec = mock.Mock()
    fake_spec.case_root = tmp_path / "case_root"
    fake_spec.build_cases.return_value = [mock.Mock(), mock.Mock()]

    with mock.patch("openfoam_driver.core.runtime.sweep_runner.load_entry_spec", return_value=fake_spec):
        result = sweep_plan(spec_path, output_dir=tmp_path / "out")

    fake_spec.apply_case.assert_not_called()
    assert result["cases"][0]["status"] == "failed"
    assert "2 cases" in result["cases"][0]["materialization_error"]


def test_sweep_run_entry_mode_executes_run_document_sequentially(tmp_path):
    # Because apply_case mutates the tutorial's shared case_root in place,
    # entry-mode sweep-run must process cases strictly one at a time (never
    # in parallel) -- already guaranteed by sweep_run's plain synchronous
    # for-loop, verified here by asserting apply_case/subprocess.run calls
    # happen in resolved-case order.
    spec_path = tmp_path / "sweep.json"
    _write_entry_spec(spec_path)
    output_dir = tmp_path / "out"

    call_order = []
    fake_case_config = mock.Mock(case_id="implicit_TNNP")
    fake_spec = mock.Mock()
    fake_spec.case_root = tmp_path / "case_root"
    fake_spec.build_cases.return_value = [fake_case_config]
    fake_spec.apply_case.side_effect = lambda *a, **k: call_order.append("apply_case")

    fake_report = mock.Mock()
    fake_report.status = "ok"

    def fake_to_json():
        # Realistic entry-mode path: the tutorial's own case_root/output_dir_name
        # tree, which is NOT a subdirectory of the sweep's own --output-dir --
        # found via a real (non-mocked) sweep-run: relative_to(output_dir)
        # raised ValueError because these are two unrelated directory trees.
        state_dir = fake_spec.case_root / f"state_{len(call_order)}"
        return {"status": "ok", "run_document": {"version": "3", "launch": {"caseRoot": str(fake_spec.case_root), "outputDir": str(state_dir)}}}
    fake_report.to_json.side_effect = fake_to_json

    def fake_subprocess_run(cmd, **kwargs):
        call_order.append("run")
        run_doc_path = Path(cmd[cmd.index("--run-document") + 1])
        run_doc = json.loads(run_doc_path.read_text())
        workflow_state_path = Path(run_doc["launch"]["outputDir"]) / "workflow_state.json"
        workflow_state_path.parent.mkdir(parents=True, exist_ok=True)
        workflow_state_path.write_text(json.dumps({"status": "completed"}))
        return mock.Mock(returncode=0, stdout="", stderr="")

    with mock.patch("openfoam_driver.core.runtime.sweep_runner.load_entry_spec", return_value=fake_spec), \
         mock.patch("openfoam_driver.core.runtime.sweep_runner.strict_plan", return_value=fake_report), \
         mock.patch("openfoam_driver.core.runtime.sweep_runner.subprocess.run", side_effect=fake_subprocess_run):
        result = sweep_run(spec_path, output_dir=output_dir)

    assert call_order == ["apply_case", "run", "apply_case", "run"]
    assert result["completed_count"] == 2
    assert result["failed_count"] == 0
    assert result["postprocess"]["status"] == "stub"


def test_sweep_run_writes_case_record_json_for_every_case(tmp_path):
    spec_path = tmp_path / "sweep.json"
    _write_spec(spec_path)
    output_dir = tmp_path / "out"

    result = sweep_run(spec_path, output_dir=output_dir)

    assert result["completed_count"] == 2
    for case_id in ("TNNP", "BuenoOrovio"):
        record_path = output_dir / case_id / "case_record.json"
        assert record_path.is_file()
        record = json.loads(record_path.read_text())
        assert record["case_id"] == case_id
        assert record["status"] == "completed"


def test_sweep_run_writes_case_record_json_even_when_a_case_fails(tmp_path):
    # case_record.json must exist for every case regardless of whether the
    # whole sweep succeeded -- an agent diagnosing a partially-failed sweep
    # needs the successful cases' records just as much as a clean sweep does.
    spec_path = tmp_path / "sweep.json"
    spec = {
        "base": {
            "electro_selectors": {"myocardiumSolver": "singleCellSolver", "tissue": "epicardialCells"},
            "physics_selectors": {"type": "electroModel"},
        },
        "sweep": {
            "mode": "cross_product",
            "independent": {"ionicModel": ["TNNP", "not_a_real_model"]},
            "dependent": [{"name": "caseId", "derive": "case_id_template", "of": ["ionicModel"]}],
        },
    }
    spec_path.write_text(json.dumps(spec))
    output_dir = tmp_path / "out"

    result = sweep_run(spec_path, output_dir=output_dir)

    assert result["failed_count"] >= 1
    assert (output_dir / "TNNP" / "case_record.json").is_file()


def test_sweep_run_archives_each_case_postprocessing_output_when_configured(tmp_path):
    # base.archive_dir_name opts an entry-mode sweep into the generic
    # snapshot/diff collection (output_collection.py): real bug this
    # reproduces -- hex workflow_dags have no "clean" step, so
    # case_root/postProcessing/ persists and accumulates across sequential
    # cases sharing one case_root. Each case's own new/changed file must land
    # inside that case's own output_dir_name folder (workflow_state_path's
    # parent, the same directory workflow_state.json lives in) under
    # <archive_dir_name>/, distinctly, without needing the tutorial's own
    # bespoke staging code or a separate cache location.
    spec_path = tmp_path / "sweep.json"
    spec = {
        "base": {"entry": "niederer2012", "archive_dir_name": "sweepCases"},
        "sweep": {
            "mode": "cross_product",
            "independent": {"dx_values": [[0.5], [0.2]]},
            "dependent": [{"name": "caseId", "derive": "case_id_template", "of": ["dx_values"]}],
        },
    }
    spec_path.write_text(json.dumps(spec))
    output_dir = tmp_path / "out"
    case_root = tmp_path / "case_root"
    (case_root / "postProcessing").mkdir(parents=True)

    call_order = []
    case_output_dirs: dict[int, Path] = {}
    fake_case_config = mock.Mock(case_id="dx0.5")
    fake_spec = mock.Mock()
    fake_spec.case_root = case_root
    fake_spec.build_cases.return_value = [fake_case_config]
    fake_spec.apply_case.side_effect = lambda *a, **k: call_order.append("apply_case")

    fake_report = mock.Mock()
    fake_report.status = "ok"

    def fake_to_json():
        state_dir = case_root / f"state_{len(call_order)}"
        return {"status": "ok", "run_document": {"version": "3", "launch": {"caseRoot": str(case_root), "outputDir": str(state_dir)}}}
    fake_report.to_json.side_effect = fake_to_json

    def fake_subprocess_run(cmd, **kwargs):
        call_order.append("run")
        n = len([c for c in call_order if c == "run"])
        run_doc_path = Path(cmd[cmd.index("--run-document") + 1])
        run_doc = json.loads(run_doc_path.read_text())
        output_dir_for_case = Path(run_doc["launch"]["outputDir"])
        case_output_dirs[n] = output_dir_for_case
        workflow_state_path = output_dir_for_case / "workflow_state.json"
        workflow_state_path.parent.mkdir(parents=True, exist_ok=True)
        workflow_state_path.write_text(json.dumps({"status": "completed"}))
        # Simulate the solver writing this case's own deterministically-named
        # output into the SHARED case_root's postProcessing/ dir.
        (case_root / "postProcessing" / f"case_{n}.dat").write_text(f"result {n}")
        return mock.Mock(returncode=0, stdout="", stderr="")

    with mock.patch("openfoam_driver.core.runtime.sweep_runner.load_entry_spec", return_value=fake_spec), \
         mock.patch("openfoam_driver.core.runtime.sweep_runner.strict_plan", return_value=fake_report), \
         mock.patch("openfoam_driver.core.runtime.sweep_runner.subprocess.run", side_effect=fake_subprocess_run):
        result = sweep_run(spec_path, output_dir=output_dir)

    assert result["completed_count"] == 2
    # Each case's archived output lands inside that case's own output_dir --
    # the same directory workflow_state.json lives in -- not a separate
    # shared cache keyed by case_id.
    assert (case_output_dirs[1] / "sweepCases" / "case_1.dat").read_text() == "result 1"
    assert (case_output_dirs[2] / "sweepCases" / "case_2.dat").read_text() == "result 2"
    assert case_output_dirs[1] != case_output_dirs[2]


def test_sweep_run_archives_each_case_postprocessing_output_by_default(tmp_path):
    # Same setup as test_sweep_run_archives_each_case_postprocessing_output_when_configured,
    # but the spec does NOT supply base.archive_dir_name -- archival must still
    # happen, using the built-in default name, not be skipped entirely.
    spec_path = tmp_path / "sweep.json"
    spec = {
        "base": {"entry": "niederer2012"},
        "sweep": {
            "mode": "cross_product",
            "independent": {"dx_values": [[0.5], [0.2]]},
            "dependent": [{"name": "caseId", "derive": "case_id_template", "of": ["dx_values"]}],
        },
    }
    spec_path.write_text(json.dumps(spec))
    output_dir = tmp_path / "out"
    case_root = tmp_path / "case_root"
    (case_root / "postProcessing").mkdir(parents=True)

    call_order = []
    case_output_dirs: dict[int, Path] = {}
    fake_case_config = mock.Mock(case_id="dx0.5")
    fake_spec = mock.Mock()
    fake_spec.case_root = case_root
    fake_spec.build_cases.return_value = [fake_case_config]
    fake_spec.apply_case.side_effect = lambda *a, **k: call_order.append("apply_case")

    fake_report = mock.Mock()
    fake_report.status = "ok"

    def fake_to_json():
        state_dir = case_root / f"state_{len(call_order)}"
        return {"status": "ok", "run_document": {"version": "3", "launch": {"caseRoot": str(case_root), "outputDir": str(state_dir)}}}
    fake_report.to_json.side_effect = fake_to_json

    def fake_subprocess_run(cmd, **kwargs):
        call_order.append("run")
        n = len([c for c in call_order if c == "run"])
        run_doc_path = Path(cmd[cmd.index("--run-document") + 1])
        run_doc = json.loads(run_doc_path.read_text())
        output_dir_for_case = Path(run_doc["launch"]["outputDir"])
        case_output_dirs[n] = output_dir_for_case
        workflow_state_path = output_dir_for_case / "workflow_state.json"
        workflow_state_path.parent.mkdir(parents=True, exist_ok=True)
        workflow_state_path.write_text(json.dumps({"status": "completed"}))
        (case_root / "postProcessing" / f"case_{n}.dat").write_text(f"result {n}")
        return mock.Mock(returncode=0, stdout="", stderr="")

    with mock.patch("openfoam_driver.core.runtime.sweep_runner.load_entry_spec", return_value=fake_spec), \
         mock.patch("openfoam_driver.core.runtime.sweep_runner.strict_plan", return_value=fake_report), \
         mock.patch("openfoam_driver.core.runtime.sweep_runner.subprocess.run", side_effect=fake_subprocess_run):
        result = sweep_run(spec_path, output_dir=output_dir)

    assert result["completed_count"] == 2
    assert (case_output_dirs[1] / "collectedOutput" / "case_1.dat").read_text() == "result 1"
    assert (case_output_dirs[2] / "collectedOutput" / "case_2.dat").read_text() == "result 2"


def test_sweep_run_archives_nothing_for_generic_case_folder_sweeps(tmp_path):
    # archive_dir_name only applies to entry mode (see sweep_runner.sweep_run's
    # archive_dir_name assignment) -- a generic/case-folder sweep must not
    # gain a spurious "collectedOutput" subfolder just because the default
    # changed from None to a string.
    spec_path = tmp_path / "sweep.json"
    _write_spec(spec_path)
    output_dir = tmp_path / "out"

    result = sweep_run(spec_path, output_dir=output_dir)

    assert result["completed_count"] == 2
    for case_id in ("TNNP", "BuenoOrovio"):
        assert not (output_dir / case_id / "collectedOutput").exists()


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
    run_doc = json.loads((output_dir / "TNNP" / "run_document.json").read_text())
    step = run_doc["workflowDag"]["steps"][0]
    assert step["id"] == "run"
    assert step["command"] == "Allrun"
    assert step["depends_on"] == []
    state_step = run_doc["workflowState"]["steps"][0]
    assert state_step["step_id"] == "run"
    assert state_step["status"] == "pending"
    assert state_step["attempt"] == 0

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
    assert result["postprocess"]["status"] == "skipped"


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
            "version": "3",
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


def test_fresh_reruns_case_reported_as_completed_and_wipes_stray_files(tmp_path):
    # Mirrors test_resume_skips_terminal_completed_case, but with fresh=True:
    # this reproduces the 2026-08-05 incident (a case directory whose
    # workflow_state.json says "completed" from a previous session was
    # silently reported as fresh) and asserts --fresh actually reruns it and
    # wipes the whole output_dir, not just the state file.
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
    stray_path = output_dir / "stray_from_previous_run.txt"
    stray_path.write_text("should be wiped by --fresh")
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

    fake_report = mock.Mock()
    fake_report.status = "ok"
    fake_report.to_json.return_value = {
        "status": "ok",
        "run_document": {"version": "3", "launch": {"outputDir": str(state_dir)}},
    }

    def fake_subprocess_run(cmd, **kwargs):
        state_dir.mkdir(parents=True, exist_ok=True)
        (state_dir / "workflow_state.json").write_text('{"status": "completed"}')
        return mock.Mock(returncode=0, stdout="", stderr="")

    with mock.patch("openfoam_driver.core.runtime.sweep_runner.materialize_case") as mock_materialize, \
         mock.patch("openfoam_driver.core.runtime.sweep_runner.strict_plan", return_value=fake_report), \
         mock.patch("openfoam_driver.core.runtime.sweep_runner.subprocess.run", side_effect=fake_subprocess_run) as mock_run:
        from openfoam_driver.core.runtime.sweep_runner import sweep_run
        result = sweep_run(spec_path, output_dir=output_dir, fresh=True)

    mock_materialize.assert_called_once()
    mock_run.assert_called_once()
    assert result["completed_count"] == 1
    assert not stray_path.exists()
    by_id = {case["case_id"]: case for case in result["cases"]}
    assert by_id["TNNP"]["outcome"] != "skipped"


def test_fresh_defaults_to_false_and_preserves_resume_behavior(tmp_path):
    # Regression guard: omitting fresh (or passing fresh=False) must keep the
    # existing skip-if-completed behavior exactly as test_resume_skips_
    # terminal_completed_case already verifies -- this just re-asserts it
    # with fresh explicitly passed as False, at the new call signature.
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
        result = sweep_run(spec_path, output_dir=output_dir, fresh=False)

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
        "run_document": {"version": "3", "launch": {"outputDir": str(state_dir)}},
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
            "version": "3",
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
