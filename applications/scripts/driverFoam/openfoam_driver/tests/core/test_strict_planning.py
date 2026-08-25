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
#     test_strict_planning
#
# Description
#     Tests strict planning logic and specification contracts.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

from __future__ import annotations

import json
import tempfile
from contextlib import redirect_stdout
from io import StringIO
from pathlib import Path
from unittest import mock

import pytest
from openfoam_driver.tests.conftest import skip_without_monorepo
pytestmark = skip_without_monorepo

from openfoam_driver.core import strict_planning
from openfoam_driver.cli import main
from openfoam_driver.scripts._dict_keys_scanner import (
    compute_dict_key_drift,
    strict_dict_key_report,
)
from openfoam_driver.plugins.cardiacfoam_plugin import CardiacFoamPlugin
from types import SimpleNamespace

from openfoam_driver.core.runtime.models import CaseConfig, TutorialSpec
from openfoam_driver.core.strict_planning import (
    StrictPlanReport,
    _is_nondimensional_entry,
    _mesh_geometry_diagnostics,
    strict_plan,
)


REPO_ROOT = Path(__file__).resolve().parents[6]
CARDIAC_PLUGIN = CardiacFoamPlugin()
CARDIAC_MAPPING = CARDIAC_PLUGIN.get_profile().cxx_mapping


_FOAM_HEADER = (
    "FoamFile\n{\n version 2.0;\n format ascii;\n"
    ' arch "LSB;label=32;scalar=64";\n class vectorField;\n'
    " object points;\n}\n"
)


def _write_unit_mesh(case_root: Path) -> None:
    """A [0,1] unit-domain mesh: max_dim == 1.0."""
    pm = case_root / "constant" / "polyMesh"
    pm.mkdir(parents=True)
    pm.joinpath("points").write_text(_FOAM_HEADER + "\n2\n(\n(0 0 0)\n(1 1 1)\n)\n")


def _spec_with_workflow(case_root: Path, *, steps: list[dict]) -> TutorialSpec:
    return TutorialSpec(
        name=case_root.name,
        case_root=case_root,
        setup_root=case_root,
        output_dir=case_root / "postProcessing",
        build_cases=lambda: [CaseConfig(case_id="default", params={})],
        apply_case=lambda *_args, **_kwargs: None,
        metadata={
            "entry_name": case_root.name,
            "entry_kind": "case_folder",
            "entry_path": case_root.name,
            "source_type": "filesystem_case",
            "workflow_family": None,
            "workflow_dag": {"steps": steps},
        },
    )


def test_report_has_mesh_geometry_field() -> None:
    report = StrictPlanReport(status="ok", entry="x", resolved_entry={})
    payload = report.to_json()
    assert "mesh_geometry_diagnostics" in payload
    assert payload["mesh_geometry_diagnostics"] == []
    assert payload["readiness_score"] == {}
    assert payload["simulation_audit"] == []


def test_mesh_adapter_flags_non_si(tmp_path: Path) -> None:
    pm = tmp_path / "constant" / "polyMesh"
    pm.mkdir(parents=True)
    pm.joinpath("points").write_text(
        _FOAM_HEADER + "\n2\n(\n(0 0 0)\n(50 50 50)\n)\n"
    )
    diags = _mesh_geometry_diagnostics(tmp_path)
    codes = {d.code for d in diags}
    assert "mesh_not_si" in codes
    assert all(d.source == "mesh_geometry" for d in diags)


def test_mesh_gate_skipped_by_env(tmp_path: Path, monkeypatch) -> None:
    monkeypatch.setenv("SKIP_MESH_DIAGNOSTICS", "1")
    assert _mesh_geometry_diagnostics(tmp_path) == ()


def test_exempt_short_circuits_unit_domain(tmp_path: Path) -> None:
    # A [0,1] mesh would classify "mm", but an exempt case must not be flagged.
    _write_unit_mesh(tmp_path)
    assert _mesh_geometry_diagnostics(tmp_path, exempt=False) != ()  # baseline
    assert _mesh_geometry_diagnostics(tmp_path, exempt=True) == ()


def test_manufactured_entry_is_nondimensional(tmp_path: Path) -> None:
    spec = SimpleNamespace(
        case_root=str(tmp_path),
        metadata={"entry_name": "manufacturedBidomain"},
    )
    assert _is_nondimensional_entry(spec) is True


def test_plain_entry_is_dimensional(tmp_path: Path) -> None:
    spec = SimpleNamespace(
        case_root=str(tmp_path),
        metadata={"entry_name": "singleCell", "workflow_family": "tutorial"},
    )
    assert _is_nondimensional_entry(spec) is False


def test_strict_plan_succeeds_for_single_cell() -> None:
    report = strict_plan("singleCell", openfoam_bashrc="/no/such/openfoam/bashrc")
    payload = report.to_json()

    assert payload["status"] == "ok"
    assert payload["readiness_score"]["score"] == 100
    assert payload["readiness_score"]["status"] == "ready"
    assert {
        item["stage"] for item in payload["simulation_audit"]
    } == {
        "simulation_generation",
        "case_preparation_files",
        "dictionary_resolution",
        "workflow_preparation",
        "artifact_prediction",
        "environment_preflight",
        "mesh_geometry",
    }
    generation_audit = next(
        item for item in payload["simulation_audit"]
        if item["stage"] == "simulation_generation"
    )
    assert generation_audit["evidence"]["case_count"] >= 1
    assert payload["resolved_entry"]["entry_kind"] == "registered_tutorial"
    assert payload["expected_artifacts"]
    assert payload["run_document"]["version"] == "3"
    assert payload["run_document"]["validation"]["status"] == "ok"
    assert payload["workflow_diagnostics"] == []
    assert payload["workflow_dag"]["schema_version"] == "1"
    assert payload["workflow_dag"]["step_status_values"] == [
        "pending",
        "running",
        "completed",
        "failed",
        "skipped",
    ]
    solve_step = payload["workflow_dag"]["steps"][0]
    assert solve_step["command"] == "cardiacFoam"
    assert solve_step["args"] == []
    assert solve_step["cwd"] == "."
    assert solve_step["retry_policy"] == {}
    assert {artifact["artifact_id"] for artifact in payload["expected_artifacts"]} <= set(
        solve_step["produces"]
    )
    assert payload["run_document"]["workflowDag"] == payload["workflow_dag"]
    assert payload["workflow_state"]["status"] == "pending"
    assert payload["workflow_state"]["current_step_id"] == "solve"
    assert payload["workflow_state"]["completed_steps"] == []
    assert payload["workflow_state"]["failed_step_id"] is None
    assert payload["workflow_state"]["steps"][0]["step_id"] == "solve"
    assert payload["workflow_state"]["steps"][0]["status"] == "pending"
    assert payload["workflow_state"]["steps"][0]["attempt"] == 0
    assert payload["workflow_state"]["steps"][0]["command"] == "cardiacFoam"
    assert payload["workflow_state"]["steps"][0]["args"] == []
    assert payload["workflow_state"]["steps"][0]["cwd"] == "."
    assert payload["workflow_state"]["steps"][0]["exit_code"] is None
    assert payload["workflow_state"]["steps"][0]["stdout_log"] is None
    assert payload["workflow_state"]["steps"][0]["stderr_log"] is None
    assert payload["run_document"]["workflowState"] == payload["workflow_state"]


def test_strict_plan_succeeds_for_manufactured_tutorial() -> None:
    report = strict_plan("manufacturedBidomain")
    payload = report.to_json()

    assert payload["status"] == "ok"
    assert payload["workflow_dag"]["steps"]
    assert payload["workflow_state"]["current_step_id"] == "mesh"
    # Step count is not asserted here -- run_in_parallel defaults to True and
    # wraps solve with decomposePar/reconstructPar (see parallel_execution.py),
    # so the exact count is an implementation detail of the real committed
    # decomposeParDict, not something this test should hardcode.
    assert [step["status"] for step in payload["workflow_state"]["steps"]] == (
        ["pending"] * len(payload["workflow_dag"]["steps"])
    )
    assert {
        step["step_id"] for step in payload["workflow_state"]["steps"]
    } == {
        step["id"] for step in payload["workflow_dag"]["steps"]
    }
    assert any(
        artifact["artifact_id"] == "verification_error_summary"
        for artifact in payload["expected_artifacts"]
    )


def test_cli_plan_strict_prints_json_and_returns_zero() -> None:
    out = StringIO()
    with redirect_stdout(out):
        code = main(["plan", "--strict", "--entry", "singleCell"])

    payload = json.loads(out.getvalue())
    assert code == 0
    assert payload["status"] == "ok"
    assert payload["launch"]["command"]


def test_strict_plan_status_ignores_environment_only_errors(monkeypatch) -> None:
    monkeypatch.delenv("SKIP_ENV_DIAGNOSTICS", raising=False)
    monkeypatch.delenv("WM_PROJECT_DIR", raising=False)
    monkeypatch.setattr(
        strict_planning.shutil,
        "which",
        lambda name, *_, **__: f"/usr/bin/{name}" if name == "cardiacFoam" else None,
    )

    report = strict_plan("singleCell", openfoam_bashrc="/no/such/openfoam/bashrc")
    payload = report.to_json()

    assert payload["status"] == "ok"
    assert payload["run_document"]["status"] == "planned"
    assert payload["run_document"]["validation"]["status"] == "ok"
    assert payload["readiness_score"]["status"] == "blocked"
    assert "environment_preflight" in payload["readiness_score"]["blocked_stages"]
    assert any(
        item["code"] == "missing_openfoam_env"
        for item in payload["environment_diagnostics"]
    )


def test_cli_run_strict_refuses_environment_errors_before_execution(monkeypatch) -> None:
    monkeypatch.delenv("SKIP_ENV_DIAGNOSTICS", raising=False)
    monkeypatch.delenv("WM_PROJECT_DIR", raising=False)
    monkeypatch.setattr(
        strict_planning.shutil,
        "which",
        lambda name, *_, **__: f"/usr/bin/{name}" if name == "cardiacFoam" else None,
    )

    out = StringIO()
    with redirect_stdout(out):
        code = main([
            "run",
            "--strict",
            "--entry",
            "singleCell",
            "--openfoam-bashrc",
            "/no/such/openfoam/bashrc",
        ])

    payload = json.loads(out.getvalue())
    assert code == 1
    assert payload["status"] == "failed"
    assert payload["error"] == "Execution environment preflight failed."
    assert any(
        item["code"] == "missing_openfoam_env"
        for item in payload["environment_diagnostics"]
    )


def test_strict_plan_fails_on_unknown_workflow_command() -> None:
    with tempfile.TemporaryDirectory() as temp_dir:
        tutorials_root = Path(temp_dir)
        case_root = tutorials_root / "badCase"
        (case_root / "constant").mkdir(parents=True)
        (case_root / "system").mkdir()
        (case_root / "constant" / "physicsProperties").write_text("type electroModel;\n")
        (case_root / "constant" / "electroProperties").write_text(
            "myocardiumSolver singleCellSolver;\n"
            "singleCellSolverCoeffs\n"
            "{\n"
            "    ionicModel AlievPanfilov;\n"
            "    tissue myocyte;\n"
            "    solutionAlgorithm explicit;\n"
            "}\n"
        )
        for name in ("controlDict", "fvSchemes", "fvSolution"):
            (case_root / "system" / name).write_text("\n")
        with mock.patch.object(
            strict_planning,
            "load_entry_spec",
            return_value=_spec_with_workflow(
                case_root,
                steps=[{"id": "unknown", "command": "notARealUtility", "depends_on": []}],
            ),
        ):
            report = strict_plan(
                "badCase",
                overrides={"tutorials_root": str(tutorials_root)},
            )

    payload = report.to_json()
    assert payload["status"] == "failed"
    assert any(
        item["code"] == "unknown_workflow_command"
        for item in payload["artifact_diagnostics"]
    )


def test_strict_plan_fails_on_unknown_workflow_dependency() -> None:
    with tempfile.TemporaryDirectory() as temp_dir:
        tutorials_root = Path(temp_dir)
        case_root = tutorials_root / "badDependency"
        (case_root / "constant").mkdir(parents=True)
        (case_root / "system").mkdir()
        (case_root / "constant" / "physicsProperties").write_text("type electroModel;\n")
        (case_root / "constant" / "electroProperties").write_text(
            "myocardiumSolver singleCellSolver;\n"
            "singleCellSolverCoeffs\n"
            "{\n"
            "    ionicModel AlievPanfilov;\n"
            "    tissue myocyte;\n"
            "    solutionAlgorithm explicit;\n"
            "}\n"
        )
        for name in ("controlDict", "fvSchemes", "fvSolution"):
            (case_root / "system" / name).write_text("\n")
        with mock.patch.object(
            strict_planning,
            "load_entry_spec",
            return_value=_spec_with_workflow(
                case_root,
                steps=[{"id": "solve", "command": "cardiacFoam", "depends_on": ["mesh"]}],
            ),
        ):
            report = strict_plan(
                "badDependency",
                overrides={"tutorials_root": str(tutorials_root)},
            )

    payload = report.to_json()
    assert payload["status"] == "failed"
    assert payload["readiness_score"]["status"] == "blocked"
    assert "workflow_preparation" in payload["readiness_score"]["blocked_stages"]
    workflow_audit = next(
        item for item in payload["simulation_audit"]
        if item["stage"] == "workflow_preparation"
    )
    assert workflow_audit["points"] == 0
    assert payload["run_document"]["status"] == "failed"
    assert payload["run_document"]["validation"]["status"] == "failed"
    assert any(
        item["code"] == "unknown_workflow_dependency"
        for item in payload["workflow_diagnostics"]
    )


def test_strict_plan_fails_when_artifact_prediction_is_empty() -> None:
    with tempfile.TemporaryDirectory() as temp_dir:
        tutorials_root = Path(temp_dir)
        case_root = tutorials_root / "missingArtifacts"
        (case_root / "constant").mkdir(parents=True)
        (case_root / "system").mkdir()
        (case_root / "constant" / "physicsProperties").write_text("type electroModel;\n")
        (case_root / "constant" / "electroProperties").write_text(
            "myocardiumSolver futureSolver;\n"
            "futureSolverCoeffs\n"
            "{\n"
            "    ionicModel AlievPanfilov;\n"
            "}\n"
        )
        for name in ("controlDict", "fvSchemes", "fvSolution"):
            (case_root / "system" / name).write_text("\n")

        report = strict_plan(
            "missingArtifacts",
            overrides={"tutorials_root": str(tutorials_root)},
        )

    payload = report.to_json()
    assert payload["status"] == "failed"
    assert any(
        item["code"] == "empty_artifact_prediction"
        for item in payload["artifact_diagnostics"]
    )


def test_strict_dict_key_scanner_allowlist_is_current() -> None:
    assert CARDIAC_MAPPING is not None
    report = strict_dict_key_report(
        REPO_ROOT / "src",
        allowlist_path=CARDIAC_MAPPING.allowlist_path,
        entries=CARDIAC_PLUGIN.get_dict_entries(),
    )
    assert report.status == "ok"
    assert report.to_json()["unused_allowlist"] == []


def test_strict_dict_key_scanner_fails_on_unallowlisted_key() -> None:
    with tempfile.TemporaryDirectory() as temp_dir:
        src_root = Path(temp_dir) / "src"
        src_root.mkdir()
        (src_root / "reader.C").write_text(
            'void read(const Foam::dictionary& dict) { dict.lookup("unlistedStrictKey"); }\n'
        )
        drift = compute_dict_key_drift(
            src_root,
            entries=CARDIAC_PLUGIN.get_dict_entries(),
        )
        allowlist_path = Path(temp_dir) / "allowlist.json"
        allowlist_path.write_text(json.dumps({
            "unmatched_cxx_reads": sorted(drift["unmatched_cxx_reads"] - {"unlistedStrictKey"}),
            "stale_paths": sorted(drift["stale_paths"]),
            "unmatched_subdicts": sorted(drift["unmatched_subdicts"]),
        }))

        report = strict_dict_key_report(
            src_root,
            allowlist_path=allowlist_path,
            entries=CARDIAC_PLUGIN.get_dict_entries(),
        )

    payload = report.to_json()
    assert payload["status"] == "failed"
    assert payload["unmatched_cxx_reads"] == ["unlistedStrictKey"]


def test_batched_ionic_model_does_not_require_optional_batched_keys():
    """batchedIntegrator/batchedSubsteps default in C++, so a batched case
    that omits them must still plan cleanly.

    Both are read only via lookupOrDefault -- batchedIonicModel.H:197,200 and
    batchedActiveTensionModel.C:46,48 (defaults 1 and "euler"). The catalog
    nonetheless marked them required_when the ionic model is batched, which
    rejected monodomain1DCableCV: it selects TWorldcompactBatched and sets
    neither key, which is legal.
    """
    from openfoam_driver.core.plugin_interface import default_driver_context
    from openfoam_driver.core import strict_planning as sp

    context = default_driver_context()
    report = sp.strict_plan(
        "cable1DCVConvergence", driver_context=context
    ).to_json()
    errors = [
        d for d in report["run_document"]["validation"].get("diagnostics", [])
        if d.get("level") == "error"
    ]
    assert errors == [], f"unexpected validation errors: {errors}"


def test_electromechanics_is_advertised_as_not_working_while_it_is_not():
    """Keep the agent-facing warning and reality in sync.

    Electromechanics is a deliberately deferred gap: the EM entry lays its
    dicts out per region (constant/electro/electroProperties) while the
    planner looks for constant/electroProperties, so it fails strict
    planning. Agents were finding that failure and trying to "fix" it.

    This asserts both halves. If EM is ever made to work, this test fails --
    which is the point: the display summary and AGENT_GUIDE warning must be
    removed in the same change, not left behind telling agents to stay away
    from something that now works.
    """
    from openfoam_driver.core.plugin_interface import default_driver_context
    from openfoam_driver.core import strict_planning as sp
    from openfoam_driver.plugins.cardiacfoam.tutorials.display import TUTORIALS

    entry = "manufacturedMonodomainTotalLagrangianEM"

    report = sp.strict_plan(entry, driver_context=default_driver_context()).to_json()
    errors = [
        d for d in report["run_document"]["validation"].get("diagnostics", [])
        if d.get("level") == "error"
    ]
    assert errors, (
        f"{entry} now plans cleanly. Electromechanics apparently works: drop "
        "the NOT CURRENTLY WORKING warning from tutorials/display.py and the "
        "electromechanics note from AGENT_GUIDE.md, then delete this test."
    )

    display = next(d for d in TUTORIALS if d.id == entry)
    haystack = f"{display.title} {display.summary}".lower()
    assert "not currently working" in haystack, (
        f"{entry} fails strict planning but its display does not say so; an "
        "agent will pick it and then try to repair the planner."
    )


def test_absent_stimulus_block_is_not_invented_from_defaults():
    """A case with no stimulus must not come back paced.

    stimulusIO.C:149-155 returns a no-op protocol when a case has no
    singleCellStimulus sub-dict at all; the FatalError at :159-176 only
    guards a block that exists and is incomplete. So "no stimulus" is legal.

    The catalog marked the whole family required whenever
    myocardiumSolver==singleCellSolver, and the builder satisfies a
    required-but-absent key by writing its typical_value -- so dropping the
    block yielded stim_amplitude 60 and nstim1 3, turning a quiescent run
    into a paced one.
    """
    from openfoam_driver.plugins.cardiacfoam.dict_builder import (
        build_electro_properties,
        parse_electro_properties,
    )
    from openfoam_driver.core.specs.common import tutorials_root_default

    committed = (
        tutorials_root_default()
        / "electrophysiologyProtocols/singleCell/constant/electroProperties"
    )
    parsed = parse_electro_properties(committed)
    without_stimulus = {
        k: v for k, v in parsed["overrides"].items()
        if "singleCellStimulus" not in k
    }

    text = build_electro_properties(parsed["selectors"], overrides=without_stimulus)

    invented = [
        line.strip() for line in text.splitlines()
        if any(k in line for k in ("stim_start", "stim_duration",
                                   "stim_amplitude", "stim_period", "nstim"))
    ]
    assert not invented, (
        "builder invented a stimulus the case did not ask for:\n  "
        + "\n  ".join(invented)
    )


def test_dictionary_resolution_audit_text_is_plugin_neutral_for_non_cardiac_plugin(
    tmp_path: Path,
) -> None:
    """P2.7: the dictionary_resolution audit stage's success text must come
    from the active plugin, not a core-hardcoded cardiac sentence. A
    non-cardiac plugin must not see "electroProperties"/"physicsProperties"
    in its own audit text."""
    from openfoam_driver.core.runtime.strict_audit import _build_simulation_audit
    from openfoam_driver.core.plugin_interface import driver_context
    from openfoam_driver.tests.plugins.minimal_plugin import MinimalOpenFOAMPlugin

    context = driver_context(MinimalOpenFOAMPlugin(), source="test:minimal")
    spec = SimpleNamespace(
        case_root=tmp_path,
        metadata={},  # not a generic_case, exercises the plugin-sourced branch
        build_cases=lambda: [],
    )

    audit_items, _generation_diagnostics, _readiness = _build_simulation_audit(
        spec=spec,
        driver_context=context,
        workflow_dag=None,
        artifacts=(),
        validation_diagnostics=(),
        workflow_diagnostics=(),
        artifact_diagnostics=(),
        environment_diagnostics=(),
        mesh_geometry_diagnostics=(),
    )

    resolution_item = next(
        item for item in audit_items if item.stage == "dictionary_resolution"
    )
    assert "electroProperties" not in resolution_item.summary
    assert "physicsProperties" not in resolution_item.summary
