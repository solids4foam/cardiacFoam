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
# Script
#     cli
#
# Description
#     Provides command-line interface for execution pipelines.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

from __future__ import annotations

import argparse
import json
from dataclasses import asdict, dataclass
from pathlib import Path

from .core.runtime.failure_context import build_failure_context
from .core.runtime.launch_readiness import is_execution_successful, is_launchable
from .core.runtime.openfoam_environment import load_openfoam_environment
from .core.runtime.remediation import build_candidate_remediations
from .core.runtime.remediation_audit import append_remediation_record
from .core.runtime.workflow_runner import run_workflow_step, _step_state_by_id
from .core.runtime.workflow_orchestrator import run_workflow
from .core.runtime.workflow_state import workflow_state_from_json
from .core.runtime.postprocess_phase import build_standalone_case_record, run_postprocess_phase, write_case_record
from .core.runtime.registry import ENTRY_KIND_VALUES, list_tutorials
from .core.runtime.sweep_runner import _stage_entry_case, sweep_plan, sweep_run
from openfoam_driver.core.introspection import describe_entry
from openfoam_driver.core.specs.common import default_setup_dir_name
from openfoam_driver.core.specs.paths import (
    default_sweep_output_dir,
    driverfoam_scratch_root,
    repo_root_default,
)
from openfoam_driver.core.specs.apply_overrides import validate_overrides, apply_overrides, OverrideError
from openfoam_driver.core.strict_planning import (
    StrictDiagnostic,
    _environment_diagnostics,
    _utility_produces_by_command,
    strict_plan,
)
from .core.runtime.run_document_exec import build_execution_inputs, load_run_document, _allowed_runs_root
from .core.runtime.fresh import ensure_fresh_output_dir


@dataclass(frozen=True)
class _ExecutionContext:
    entry_label: str
    workflow_dag: dict
    planned_state: object
    case_root: Path
    output_dir: Path
    expected_artifacts: tuple
    setup_root: Path | None = None
    environment_diagnostics: tuple[StrictDiagnostic, ...] = ()
    execution_env: dict[str, str] | None = None
    source_path: str | None = None


def _step_payload(
    *,
    status: str,
    entry: str,
    step: str,
    workflow_state_path: Path,
    workflow_state: dict,
    exit_code: int | None = None,
    stdout_log: str | None = None,
    stderr_log: str | None = None,
    error: str | None = None,
) -> dict:
    payload = {
        "status": status,
        "entry": entry,
        "step": step,
        "exit_code": exit_code,
        "stdout_log": stdout_log,
        "stderr_log": stderr_log,
        "workflow_state_path": str(workflow_state_path),
        "workflow_state": workflow_state,
    }
    if error is not None:
        payload["error"] = error
    return payload


def _terminal_status_label(workflow_status: str) -> str:
    """Map a WorkflowStepState/WorkflowRunState status to the CLI's ok/failed label.

    Delegates to ``is_execution_successful``, the shared post-execution
    predicate, so step and run cannot drift from each other. Decisions
    derive from status, never from a subprocess exit code.
    """
    return "ok" if is_execution_successful(workflow_status) else "failed"


def _refuse_environment_errors(context: _ExecutionContext, *, action: str) -> int | None:
    # Structural validity was already established when this context was
    # built (_context_from_entry / _context_from_run_document both refuse to
    # hand back a context otherwise), so only the environment half of
    # is_launchable is relevant at this dispatch-time gate.
    readiness = is_launchable(
        plan_status="ok",
        environment_diagnostics=context.environment_diagnostics,
    )
    if readiness.environment_ok:
        return None
    payload = {
        "status": "failed",
        "entry": context.entry_label,
        "action": action,
        "error": "Execution environment preflight failed.",
        "environment_diagnostics": [asdict(diagnostic) for diagnostic in context.environment_diagnostics],
    }
    if context.source_path is not None:
        payload["run_document"] = context.source_path
    print(json.dumps(payload, indent=2))
    return 1


def _attach_failure_context(payload: dict, state, step_id: str | None, *, tail_lines: int) -> None:
    """Attach a failure_context bundle when the named step state is failed.

    Mutates ``payload`` in place. Shared by step and run so the two paths
    surface failures identically. Never persisted into workflow_state.json.
    """
    if step_id is None:
        return
    step_state = _step_state_by_id(state, step_id)
    if step_state.status == "failed":
        fc = build_failure_context(step_state, max_lines=tail_lines)
        fc["candidate_remediations"] = [
            hint.to_json() for hint in build_candidate_remediations(fc)
        ]
        payload["failure_context"] = fc


def _execute_step(
    *,
    entry_label: str,
    step_id: str,
    workflow_dag: dict,
    planned_state,
    case_root: Path,
    output_dir: Path,
    expected_artifacts,
    tail_lines: int,
    execution_env: dict[str, str] | None = None,
    apply_overrides_path: str | None = None,
) -> int:
    """Run one workflow step, print the JSON payload, return the exit code.

    Shared by the --entry (strict_plan) path and the --run-document path.
    Resumes from an existing workflow_state.json under output_dir when present.
    """
    state_path = output_dir / "workflow_state.json"
    workflow_state = planned_state
    if state_path.exists():
        try:
            workflow_state = workflow_state_from_json(json.loads(state_path.read_text()))
        except Exception as exc:
            print(json.dumps({
                "status": "failed",
                "entry": entry_label,
                "step": step_id,
                "error": f"Could not read existing workflow state: {exc}",
                "workflow_state_path": str(state_path),
            }, indent=2))
            return 1
    overrides = None
    if apply_overrides_path is not None:
        try:
            overrides = json.loads(Path(apply_overrides_path).read_text())
            validate_overrides(overrides)
            apply_overrides(overrides, case_root=case_root)
        except (OSError, ValueError, OverrideError) as exc:
            print(json.dumps({
                "status": "failed",
                "entry": entry_label,
                "step": step_id,
                "error": f"--apply rejected: {exc}",
            }, indent=2))
            return 1
    try:
        result = run_workflow_step(
            workflow_dag,
            workflow_state,
            step_id,
            case_root=case_root,
            log_dir=output_dir / "workflow_logs",
            state_path=state_path,
            expected_artifacts=expected_artifacts,
            env=execution_env,
        )
    except Exception as exc:
        if overrides is not None:
            try:
                _attempt = _step_state_by_id(workflow_state, step_id).attempt
            except Exception:
                _attempt = 0
            append_remediation_record(
                output_dir, step_id=step_id, attempt=_attempt,
                applied_overrides=overrides, resulting_status="rerun_error",
            )
        print(json.dumps({
            "status": "failed",
            "entry": entry_label,
            "step": step_id,
            "error": str(exc),
            "workflow_state": workflow_state.to_json(),
        }, indent=2))
        return 1
    step_state = _step_state_by_id(result.state, step_id)
    status = _terminal_status_label(step_state.status)
    payload = _step_payload(
        status=status,
        entry=entry_label,
        step=step_id,
        workflow_state_path=state_path,
        workflow_state=result.state.to_json(),
        exit_code=result.exit_code,
        stdout_log=result.stdout_log,
        stderr_log=result.stderr_log,
    )
    _attach_failure_context(payload, result.state, step_id, tail_lines=tail_lines)
    if overrides is not None:
        append_remediation_record(
            output_dir,
            step_id=step_id,
            attempt=step_state.attempt,
            applied_overrides=overrides,
            resulting_status=status,
        )
    print(json.dumps(payload, indent=2))
    return 0 if status == "ok" else 1


def _reconciliation_payload(case_root: Path, expected_artifacts) -> dict:
    """Reconcile predicted artifacts against what actually landed on disk.

    Describes only: which predicted artifacts resolved to which files, their
    size and sha256. It makes no judgement about the values inside them --
    interpretation is the caller's job.
    """
    from .core.runtime.reconciler import reconcile_artifacts

    return reconcile_artifacts(case_root, expected_artifacts or ()).to_json()


def _execute_run(
    *,
    entry_label: str,
    workflow_dag: dict,
    planned_state,
    case_root: Path,
    output_dir: Path,
    expected_artifacts,
    tail_lines: int,
    setup_root: Path | None = None,
    execution_env: dict[str, str] | None = None,
    max_total_attempts: int | None = None,
) -> int:
    """Run a workflow to completion, print the JSON payload, return the exit code.

    Shared by the --entry (strict_plan) path and the --run-document path.
    Refuses to auto-resume a terminally-failed saved state (use action=step).
    """
    state_path = output_dir / "workflow_state.json"
    workflow_state = planned_state
    if state_path.exists():
        try:
            workflow_state = workflow_state_from_json(json.loads(state_path.read_text()))
        except Exception as exc:
            print(json.dumps({
                "status": "failed",
                "entry": entry_label,
                "error": f"Could not read existing workflow state: {exc}",
                "workflow_state_path": str(state_path),
            }, indent=2))
            return 1
    if workflow_state.status == "failed":
        print(json.dumps({
            "status": "failed",
            "entry": entry_label,
            "error": "workflow_state is failed; use action=step to rerun a failed step explicitly",
            "workflow_state_path": str(state_path),
            "workflow_state": workflow_state.to_json(),
        }, indent=2))
        return 1
    try:
        outcome = run_workflow(
            workflow_dag,
            workflow_state,
            case_root=case_root,
            output_dir=output_dir,
            expected_artifacts=expected_artifacts,
            state_path=state_path,
            env=execution_env,
            max_total_attempts=max_total_attempts,
        )
    except Exception as exc:
        try:
            error_state = workflow_state_from_json(json.loads(state_path.read_text()))
        except Exception:
            error_state = workflow_state
        print(json.dumps(_step_payload(
            status="failed",
            entry=entry_label,
            step=error_state.current_step_id or workflow_state.current_step_id,
            workflow_state_path=state_path,
            workflow_state=error_state.to_json(),
            error=str(exc),
        ), indent=2))
        return 1
    workflow_state = outcome.state
    results = list(outcome.steps)
    status = _terminal_status_label(workflow_state.status)
    payload = {
        "status": status,
        "entry": entry_label,
        "steps": results,
        "workflow_state_path": str(state_path),
        "workflow_state": workflow_state.to_json(),
    }
    if workflow_state.status == "pending" and workflow_state.current_step_id is None:
        payload["error"] = "workflow_state is pending but has no current_step_id"
    payload["artifact_reconciliation"] = _reconciliation_payload(
        case_root, expected_artifacts
    )
    if status == "ok":
        payload["postprocess"] = run_postprocess_phase(
            entry=entry_label, output_dir=output_dir,
        ).to_json()
    else:
        payload["postprocess"] = {
            "status": "skipped",
            "message": f"workflow did not complete (status={status}); postprocess not run",
        }
    case_record = build_standalone_case_record(
        entry=entry_label, case_root=case_root, setup_root=setup_root, output_dir=output_dir,
    )
    case_record_path = output_dir / "case_record.json"
    write_case_record(case_record_path, case_record)
    payload["case_record_path"] = str(case_record_path)
    _attach_failure_context(payload, workflow_state, workflow_state.failed_step_id, tail_lines=tail_lines)
    print(json.dumps(payload, indent=2))
    return 0 if status == "ok" else 1


def _context_from_run_document(args, driver_context) -> _ExecutionContext | None:
    """Load + validate an agent-authored RunDocument into executor inputs."""
    try:
        run_doc = load_run_document(args.run_document)
    except Exception as exc:
        print(json.dumps({
            "status": "failed",
            "error": f"Could not load run document: {exc}",
            "run_document": args.run_document,
        }, indent=2))
        return None
    if run_doc.plugin is not None:
        planned = run_doc.plugin
        selected = driver_context.identity.to_json()
        mismatched = [
            key for key in ("id", "version", "api_version", "capability_digest")
            if planned.get(key) != selected.get(key)
        ]
        if mismatched:
            print(json.dumps({
                "status": "failed",
                "run_document": args.run_document,
                "error": (
                    "RunDocument plugin does not match the selected plugin: "
                    + ", ".join(mismatched)
                ),
                "planned_plugin": planned,
                "selected_plugin": selected,
            }, indent=2))
            return None
    inputs, diagnostics = build_execution_inputs(
        run_doc,
        utility_produces=_utility_produces_by_command(driver_context),
        driver_context=driver_context,
    )
    if inputs is None:
        print(json.dumps({
            "status": "failed",
            "run_document": args.run_document,
            "diagnostics": list(diagnostics),
        }, indent=2))
        return None
    execution_env = load_openfoam_environment(
        explicit_bashrc=args.openfoam_bashrc,
        driver_context=driver_context,
    ).env
    setup_root_raw = (run_doc.launch or {}).get("setupRoot")
    return _ExecutionContext(
        entry_label=run_doc.name,
        workflow_dag=inputs.workflow_dag,
        planned_state=inputs.workflow_state,
        case_root=inputs.case_root,
        output_dir=inputs.output_dir,
        expected_artifacts=inputs.expected_artifacts,
        setup_root=Path(setup_root_raw) if setup_root_raw else None,
        environment_diagnostics=_environment_diagnostics(
            inputs.workflow_dag,
            openfoam_bashrc=args.openfoam_bashrc,
            driver_context=driver_context,
        ),
        execution_env=execution_env,
        source_path=args.run_document,
    )


def _context_from_entry(
    *,
    selected_entry: str,
    entry_kind: str | None,
    overrides: dict | None,
    config_path: str | None,
    openfoam_bashrc: str | None,
    driver_context,
    stage_for_execution: bool = False,
    fresh: bool = False,
) -> tuple[_ExecutionContext | None, int]:
    report = strict_plan(
        selected_entry,
        entry_kind=entry_kind,
        overrides=overrides,
        config_path=config_path,
        openfoam_bashrc=openfoam_bashrc,
        driver_context=driver_context,
    )
    readiness = is_launchable(
        plan_status=report.status,
        environment_diagnostics=report.environment_diagnostics,
    )
    if not readiness.structural_ok:
        print(json.dumps(report.to_json(), indent=2))
        return None, 1
    if report.workflow_dag is None or report.workflow_state is None:
        print(json.dumps({
            "status": "failed",
            "error": "strict plan did not produce workflow_dag and workflow_state",
        }, indent=2))
        return None, 1
    source_case_root = Path(report.launch["case_root"]).resolve()
    # Test/standalone case folders may intentionally live in a caller-owned
    # temporary root. The repository-template protection applies to registered
    # cases in this checkout; generic external cases retain their explicit
    # launch root contract.
    if stage_for_execution and source_case_root.is_relative_to(repo_root_default()):
        safe_entry = "".join(
            char if char.isalnum() or char in {"-", "_", "."} else "_"
            for char in selected_entry
        ).strip("._") or "entry"
        staged_case_root = driverfoam_scratch_root() / "runs" / safe_entry
        if fresh or not staged_case_root.exists():
            _stage_entry_case(source_case_root, staged_case_root)
        staged_overrides = dict(overrides or {})
        staged_overrides["tutorials_root"] = str(staged_case_root.parent)
        staged_overrides["case_dir_name"] = staged_case_root.name
        report = strict_plan(
            selected_entry,
            entry_kind=entry_kind,
            overrides=staged_overrides,
            config_path=config_path,
            openfoam_bashrc=openfoam_bashrc,
            driver_context=driver_context,
        )
        readiness = is_launchable(
            plan_status=report.status,
            environment_diagnostics=report.environment_diagnostics,
        )
        if not readiness.structural_ok:
            print(json.dumps(report.to_json(), indent=2))
            return None, 1
    execution_env = load_openfoam_environment(
        explicit_bashrc=openfoam_bashrc,
        driver_context=driver_context,
    ).env
    return (
        _ExecutionContext(
            entry_label=selected_entry,
            workflow_dag=report.workflow_dag,
            planned_state=report.workflow_state,
            case_root=Path(report.launch["case_root"]),
            output_dir=Path(report.launch["output_dir"]),
            expected_artifacts=report.expected_artifacts,
            setup_root=Path(report.launch["setup_root"]),
            environment_diagnostics=report.environment_diagnostics,
            execution_env=execution_env,
        ),
        0,
    )


def _dispatch_context(args, context: _ExecutionContext) -> int:
    blocked = _refuse_environment_errors(context, action=args.action)
    if blocked is not None:
        return blocked
    if args.fresh:
        fresh_error = ensure_fresh_output_dir(
            context.output_dir, fresh=True, allowed_root=_allowed_runs_root(),
        )
        if fresh_error is not None:
            print(json.dumps({
                "status": "failed",
                "entry": context.entry_label,
                "action": args.action,
                "error": fresh_error,
            }, indent=2))
            return 1
    if args.action == "step":
        return _execute_step(
            entry_label=context.entry_label,
            step_id=args.step,
            workflow_dag=context.workflow_dag,
            planned_state=context.planned_state,
            case_root=context.case_root,
            output_dir=context.output_dir,
            expected_artifacts=context.expected_artifacts,
            tail_lines=args.tail_lines,
            execution_env=context.execution_env,
            apply_overrides_path=args.apply,
        )
    return _execute_run(
        entry_label=context.entry_label,
        workflow_dag=context.workflow_dag,
        planned_state=context.planned_state,
        case_root=context.case_root,
        output_dir=context.output_dir,
        expected_artifacts=context.expected_artifacts,
        setup_root=context.setup_root,
        tail_lines=args.tail_lines,
        execution_env=context.execution_env,
        max_total_attempts=args.max_total_attempts,
    )


def _run_document_dispatch(args, driver_context) -> int:
    """Execute a validated RunDocument through the shared strict executor."""
    context = _context_from_run_document(args, driver_context)
    if context is None:
        return 1
    return _dispatch_context(args, context)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Generic OpenFOAM tutorial automation driver")
    parser.add_argument(
        "action",
        choices=[
            "describe", "plan", "step", "run", "sweep-plan", "sweep-run",
        ],
        help="Pipeline stage to execute",
    )
    parser.add_argument(
        "--plugin",
        help=(
            "Plugin to drive: an installed plugin id from the "
            "'driverfoam.plugins' entry-point group, a trusted "
            "local-development import target (module.path:PluginClass), or "
            "'none' for generic OpenFOAM. Defaults to built-in cardiacFoam. "
            "A colon always selects the import form. Either form executes "
            "the plugin's Python code."
        ),
    )
    parser.add_argument(
        "--entry",
        required=False,
        help=(
            "Entry name or relative workflow/case path to run "
            f"({', '.join(list_tutorials())}, genericCase)"
        ),
    )
    parser.add_argument(
        "--entry-kind",
        choices=list(ENTRY_KIND_VALUES),
        help="Optional entry classification override for --entry resolution.",
    )
    parser.add_argument(
        "--run-document",
        help=(
            "Path to an agent-authored RunDocument v2 JSON file. With "
            "action=run/step, executes the document's workflowDag/config "
            "instead of regenerating the plan from --entry. Mutually "
            "exclusive with --entry."
        ),
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Plan and print simulation cases without running OpenFOAM.",
    )
    parser.add_argument(
        "--strict",
        action="store_true",
        help="For action=plan/step/run, fail on incomplete machine-readable coverage.",
    )
    parser.add_argument(
        "--openfoam-bashrc",
        help=(
            "OpenFOAM bashrc to source for strict plan/step/run. Defaults to "
            "OPENFOAM_BASHRC, $WM_PROJECT_DIR/etc/bashrc, or a known local "
            "OpenFOAM install path."
        ),
    )
    parser.add_argument(
        "--step",
        help="Workflow step id to execute when action=step.",
    )
    parser.add_argument(
        "--tail-lines",
        type=int,
        default=200,
        help="For action=step/run --strict, number of log lines to include in failure_context (default 200).",
    )
    parser.add_argument(
        "--apply",
        metavar="OVERRIDES_JSON",
        help="action=step only: apply an override set (JSON list of "
             "{driver_path, value}) to the case, then rerun the step.",
    )
    parser.add_argument(
        "--continue-on-error",
        action="store_true",
        help="Continue executing remaining cases after a failure.",
    )
    parser.add_argument(
        "--config",
        help=(
            "Path to JSON file with make_spec overrides. Supports either a top-level "
            "entry map (keys: singleCell, niederer2012, manufacturedMonodomainPseudoECG, "
            "manufacturedBidomain, manufacturedBathBidomain, "
            "manufacturedEikonalECG, manufacturedMonodomainTotalLagrangianEM, "
            "manufacturedPurkinjeGraph, restitutionCurves, genericCase/randomCase) "
            "or a direct parameter object for the selected entry."
        ),
    )
    parser.add_argument(
        "--tutorials-root",
        help=(
            "Optional path to the tutorials folder. Defaults to '<repo>/tutorials' when present."
        ),
    )
    parser.add_argument(
        "--spec",
        help="Path to a sweep.json for action=sweep-plan/sweep-run.",
    )
    parser.add_argument(
        "--output-dir",
        help=(
            "Output directory for action=sweep-plan/sweep-run. Defaults to "
            "<repo>/.tmp/driverfoam/sweeps/<spec-name>; generated cases and "
            "logs never belong under tutorials/."
        ),
    )
    parser.add_argument(
        "--max-cases",
        type=int,
        default=200,
        help="Safety cap on expanded sweep case count (default 200).",
    )
    parser.add_argument(
        "--retry-failed",
        action="store_true",
        help="For action=sweep-run: rerun cases whose last recorded status was failed.",
    )
    parser.add_argument(
        "--fresh",
        action="store_true",
        help=(
            "For action=step/run/sweep-run: delete the resolved output "
            "directory before running, so the workflow executes as if no "
            "prior run existed. Use after a code/config change to guarantee "
            "a real rerun instead of silently resuming a stale "
            "workflow_state.json/sweep_manifest.json as 'completed'. "
            "Refuses to delete the filesystem root, your home directory, a "
            "too-shallow path, anything outside DRIVERFOAM_ALLOWED_RUNS_ROOT "
            "when set, or a directory with no recognizable driverFOAM "
            "artifact. No confirmation prompt -- treat --output-dir as fully "
            "disposable when passing this flag. Mutually exclusive with "
            "--retry-failed."
        ),
    )
    parser.add_argument(
        "--max-total-attempts",
        type=int,
        default=None,
        help=(
            "For action=run: whole-run ceiling on step executions (retry-storm "
            "guard). Default: unbounded (only per-step max_attempts applies)."
        ),
    )
    parser.add_argument(
        "--case-timeout-s",
        type=float,
        default=None,
        help=(
            "For action=sweep-run: wall-clock timeout per case subprocess; a case "
            "that exceeds it is marked failed and the sweep continues. Default: none."
        ),
    )

    return parser


def _load_spec_overrides(config_path: str, entry: str) -> dict:
    payload = json.loads(Path(config_path).read_text())
    if not isinstance(payload, dict):
        raise ValueError("Config file must contain a JSON object")

    normalized_requested = entry.strip().casefold()
    for key, value in payload.items():
        if key.casefold() == normalized_requested:
            if not isinstance(value, dict):
                raise ValueError(f"Config section '{key}' must be a JSON object")
            return _normalize_spec_overrides(value)

    known_tutorial_keys = {
        *(name.casefold() for name in list_tutorials()),
        "genericcase",
        "randomcase",
    }
    if any(key.casefold() in known_tutorial_keys for key in payload):
        raise KeyError(
            f"No config section found for entry '{entry}'. "
            f"Available config sections: {', '.join(payload.keys())}"
        )

    return _normalize_spec_overrides(payload)


def _normalize_spec_overrides(overrides: dict) -> dict:
    normalized = dict(overrides)

    case_dir_name = normalized.get("case_dir_name")
    setup_dir_name = normalized.get("setup_dir_name")
    if case_dir_name is not None and setup_dir_name is not None:
        if str(setup_dir_name) == default_setup_dir_name(str(case_dir_name)):
            normalized.pop("setup_dir_name")

    return normalized


_FLAG_ERRORS_BY_ACTION = {
    "describe": (
        ("dry_run", "--dry-run is not valid with action=describe"),
        ("continue_on_error", "--continue-on-error is not valid with action=describe"),
    ),
    "plan": (
        ("dry_run", "--dry-run is not valid with action=plan"),
        ("continue_on_error", "--continue-on-error is not valid with action=plan"),
    ),
    "step": (
        ("dry_run", "--dry-run is not valid with action=step"),
        ("continue_on_error", "--continue-on-error is not valid with action=step"),
    ),
    "run": (
        ("dry_run", "--dry-run is not valid with action=run"),
        ("continue_on_error", "--continue-on-error is not valid with action=run"),
    ),
}


def _validate_args(parser: argparse.ArgumentParser, args) -> None:
    for flag_name, message in _FLAG_ERRORS_BY_ACTION.get(args.action, ()):
        if getattr(args, flag_name):
            parser.error(message)
    if args.action not in {"plan", "step", "run"} and args.strict:
        parser.error("--strict is only valid with action=plan, action=step, or action=run")
    if args.openfoam_bashrc and args.action not in {"plan", "step", "run"}:
        parser.error("--openfoam-bashrc is only valid with action=plan, action=step, or action=run")
    if args.action != "step" and args.step:
        parser.error("--step is only valid with action=step")
    if args.apply is not None and args.action != "step":
        parser.error("--apply is only valid with action=step")
    if args.action not in {"step", "run"} and args.tail_lines != 200:
        parser.error("--tail-lines is only valid with action=step or action=run")
    if args.run_document and args.action not in {"run", "step"}:
        parser.error("--run-document is only valid with action=run or action=step")
    if args.run_document and args.entry:
        parser.error("--run-document and --entry are mutually exclusive")
    if args.run_document and (args.config or args.entry_kind or args.tutorials_root):
        parser.error("--config/--entry-kind/--tutorials-root are not valid with --run-document")
    if args.action in {"sweep-plan", "sweep-run"}:
        if args.entry:
            parser.error(f"--entry is not valid with action={args.action}; use --spec")
        if args.config or args.entry_kind or args.tutorials_root:
            parser.error(f"--config/--entry-kind/--tutorials-root are not valid with action={args.action}")
        if args.strict or args.run_document or args.step or args.apply is not None:
            parser.error(f"--strict/--run-document/--step/--apply are not valid with action={args.action}")
    if args.action in {"sweep-plan", "sweep-run"} and not args.spec:
        parser.error(f"action={args.action} requires --spec")
    # Sweep output defaults to the repository-local disposable workspace.
    if args.retry_failed and args.action != "sweep-run":
        parser.error("--retry-failed is only valid with action=sweep-run")
    if args.fresh and args.action not in {"step", "run", "sweep-run"}:
        parser.error("--fresh is only valid with action=step, action=run, or action=sweep-run")
    if args.fresh and args.retry_failed:
        parser.error("--fresh and --retry-failed are mutually exclusive")
    if args.max_cases != 200 and args.action not in {"sweep-plan", "sweep-run"}:
        parser.error("--max-cases is only valid with action=sweep-plan or action=sweep-run")
    if args.action not in {"sweep-plan", "sweep-run"} and (args.spec or args.output_dir):
        parser.error("--spec/--output-dir are only valid with action=sweep-plan or action=sweep-run")
    if not args.run_document and not args.entry and args.action not in {
        "sweep-plan", "sweep-run"
    }:
        parser.error("--entry is required (or use --run-document with action=run/step)")


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    _validate_args(parser, args)

    from .core.plugin_interface import (
        default_driver_context,
        generic_openfoam_context,
        load_plugin_context,
    )
    try:
        if args.plugin == "none":
            driver_context = generic_openfoam_context()
        elif args.plugin:
            driver_context = load_plugin_context(args.plugin)
        else:
            driver_context = default_driver_context()
    except Exception as exc:
        parser.error(f"Failed to load plugin {args.plugin!r}: {exc}")


    selected_entry = args.entry

    overrides = _load_spec_overrides(args.config, selected_entry) if args.config else None
    if args.tutorials_root:
        if overrides is None:
            overrides = {}
        overrides["tutorials_root"] = args.tutorials_root

    if args.action == "describe":
        print(
            json.dumps(
                describe_entry(
                    selected_entry,
                    entry_kind=args.entry_kind,
                    overrides=overrides,
                    config_path=args.config,
                    driver_context=driver_context,
                ),
                indent=2,
            )
        )
        return 0

    if args.action == "plan":
        if not args.strict:
            parser.error("action=plan currently requires --strict")
        report = strict_plan(
            selected_entry,
            entry_kind=args.entry_kind,
            overrides=overrides,
            config_path=args.config,
            openfoam_bashrc=args.openfoam_bashrc,
            driver_context=driver_context,
        )
        print(json.dumps(report.to_json(), indent=2))
        readiness = is_launchable(
            plan_status=report.status,
            environment_diagnostics=report.environment_diagnostics,
        )
        return 0 if readiness.structural_ok else 1

    if args.action == "step":
        if not (args.strict or args.run_document):
            parser.error("action=step requires --strict or --run-document")
        if not args.step:
            parser.error("action=step requires --step <id>")
        if args.run_document:
            return _run_document_dispatch(args, driver_context)
        context, failure_code = _context_from_entry(
            selected_entry=selected_entry,
            entry_kind=args.entry_kind,
            overrides=overrides,
            config_path=args.config,
            openfoam_bashrc=args.openfoam_bashrc,
            driver_context=driver_context,
            stage_for_execution=True,
            fresh=args.fresh,
        )
        if context is None:
            return failure_code
        return _dispatch_context(args, context)

    if args.action == "run":
        if not (args.strict or args.run_document):
            parser.error("action=run requires --strict or --run-document")
        if args.run_document:
            return _run_document_dispatch(args, driver_context)
        context, failure_code = _context_from_entry(
            selected_entry=selected_entry,
            entry_kind=args.entry_kind,
            overrides=overrides,
            config_path=args.config,
            openfoam_bashrc=args.openfoam_bashrc,
            driver_context=driver_context,
            stage_for_execution=True,
            fresh=args.fresh,
        )
        if context is None:
            return failure_code
        return _dispatch_context(args, context)

    if args.action == "sweep-plan":
        output_dir = args.output_dir or default_sweep_output_dir(args.spec)
        result = sweep_plan(
            args.spec,
            output_dir=output_dir,
            max_cases=args.max_cases,
            driver_context=driver_context,
        )
        print(json.dumps(result, indent=2))
        any_failed = any(case["status"] != "ok" for case in result["cases"])
        return 1 if any_failed else 0

    if args.action == "sweep-run":
        output_dir = args.output_dir or default_sweep_output_dir(args.spec)
        result = sweep_run(
            args.spec,
            output_dir=output_dir,
            max_cases=args.max_cases,
            retry_failed=args.retry_failed,
            case_timeout_s=args.case_timeout_s,
            fresh=args.fresh,
            driver_context=driver_context,
        )
        print(json.dumps(result, indent=2))
        return 1 if result["failed_count"] > 0 else 0

    raise AssertionError(f"unreachable: unhandled action {args.action!r}")
