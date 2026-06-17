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
from pathlib import Path

from .core.runtime.engine import DriverEngine
from .core.runtime.failure_context import build_failure_context
from .core.runtime.workflow_runner import run_workflow_step, _step_state_by_id
from .core.runtime.workflow_orchestrator import run_workflow
from .core.runtime.workflow_state import workflow_state_from_json
from .core.runtime.registry import ENTRY_KIND_VALUES, list_tutorials, load_entry_spec
from .introspection import describe_entry
from .specs.common import default_setup_dir_name
from .strict_planning import strict_plan, _utility_produces_by_command
from .core.runtime.run_document_exec import build_execution_inputs, load_run_document


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

    The single source of truth for the strict success contract, so step and run
    cannot drift. Decisions derive from status, never from a subprocess exit code.
    Any non-completed terminal state (pending, running, failed, skipped) is a
    failure at this boundary.
    """
    return "ok" if workflow_status == "completed" else "failed"


def _attach_failure_context(payload: dict, state, step_id: str | None, *, tail_lines: int) -> None:
    """Attach a failure_context bundle when the named step state is failed.

    Mutates ``payload`` in place. Shared by step and run so the two paths
    surface failures identically. Never persisted into workflow_state.json.
    """
    if step_id is None:
        return
    step_state = _step_state_by_id(state, step_id)
    if step_state.status == "failed":
        payload["failure_context"] = build_failure_context(step_state, max_lines=tail_lines)


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
    try:
        result = run_workflow_step(
            workflow_dag,
            workflow_state,
            step_id,
            case_root=case_root,
            log_dir=output_dir / "workflow_logs",
            state_path=state_path,
            expected_artifacts=expected_artifacts,
        )
    except Exception as exc:
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
    print(json.dumps(payload, indent=2))
    return 0 if status == "ok" else 1


def _execute_run(
    *,
    entry_label: str,
    workflow_dag: dict,
    planned_state,
    case_root: Path,
    output_dir: Path,
    expected_artifacts,
    tail_lines: int,
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
    _attach_failure_context(payload, workflow_state, workflow_state.failed_step_id, tail_lines=tail_lines)
    print(json.dumps(payload, indent=2))
    return 0 if status == "ok" else 1


def _run_document_dispatch(args) -> int:
    """Load + validate an agent-authored RunDocument and execute it.

    Shared by action=step and action=run; the two differ only in the final
    execution helper. Loading is strict (schema-validated via
    ``load_run_document``); a non-executable document returns its diagnostics.
    """
    try:
        run_doc = load_run_document(args.run_document)
    except Exception as exc:
        print(json.dumps({
            "status": "failed",
            "error": f"Could not load run document: {exc}",
            "run_document": args.run_document,
        }, indent=2))
        return 1
    inputs, diagnostics = build_execution_inputs(
        run_doc, utility_produces=_utility_produces_by_command(),
    )
    if inputs is None:
        print(json.dumps({
            "status": "failed",
            "run_document": args.run_document,
            "diagnostics": list(diagnostics),
        }, indent=2))
        return 1
    if args.action == "step":
        return _execute_step(
            entry_label=run_doc.name,
            step_id=args.step,
            workflow_dag=inputs.workflow_dag,
            planned_state=inputs.workflow_state,
            case_root=inputs.case_root,
            output_dir=inputs.output_dir,
            expected_artifacts=inputs.expected_artifacts,
            tail_lines=args.tail_lines,
        )
    return _execute_run(
        entry_label=run_doc.name,
        workflow_dag=inputs.workflow_dag,
        planned_state=inputs.workflow_state,
        case_root=inputs.case_root,
        output_dir=inputs.output_dir,
        expected_artifacts=inputs.expected_artifacts,
        tail_lines=args.tail_lines,
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Generic OpenFOAM tutorial automation driver")
    parser.add_argument(
        "action",
        choices=["sim", "post", "all", "describe", "plan", "step", "run"],
        help="Pipeline stage to execute",
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
        "--continue-on-error",
        action="store_true",
        help="Continue executing remaining cases after a failure.",
    )
    parser.add_argument(
        "--config",
        help=(
            "Path to JSON file with make_spec overrides. Supports either a top-level "
            "entry map (keys: singleCell, niederer2012, manufacturedFDA, "
            "manufacturedFDABidomain, manufacturedFDABathBidomain, "
            "manufacturedEikonalECG, manufacturedMonodomainTotalLagrangianEM, "
            "restitutionCurves, genericCase/randomCase) "
            "or a direct parameter object for the selected entry."
        ),
    )
    parser.add_argument(
        "--tutorials-root",
        help=(
            "Optional path to the tutorials folder. Defaults to '<repo>/tutorials' when present."
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


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    if args.action == "post" and args.dry_run:
        parser.error("--dry-run is not valid with action=post")
    if args.action == "describe" and args.dry_run:
        parser.error("--dry-run is not valid with action=describe")
    if args.action == "describe" and args.continue_on_error:
        parser.error("--continue-on-error is not valid with action=describe")
    if args.action == "plan" and args.dry_run:
        parser.error("--dry-run is not valid with action=plan")
    if args.action == "plan" and args.continue_on_error:
        parser.error("--continue-on-error is not valid with action=plan")
    if args.action == "step" and args.dry_run:
        parser.error("--dry-run is not valid with action=step")
    if args.action == "step" and args.continue_on_error:
        parser.error("--continue-on-error is not valid with action=step")
    if args.action == "run" and args.dry_run:
        parser.error("--dry-run is not valid with action=run")
    if args.action == "run" and args.continue_on_error:
        parser.error("--continue-on-error is not valid with action=run")
    if args.action not in {"plan", "step", "run"} and args.strict:
        parser.error("--strict is only valid with action=plan, action=step, or action=run")
    if args.action != "step" and args.step:
        parser.error("--step is only valid with action=step")
    if args.action not in {"step", "run"} and args.tail_lines != 200:
        parser.error("--tail-lines is only valid with action=step or action=run")
    if args.run_document and args.action not in {"run", "step"}:
        parser.error("--run-document is only valid with action=run or action=step")
    if args.run_document and args.entry:
        parser.error("--run-document and --entry are mutually exclusive")
    if args.run_document and (args.config or args.entry_kind or args.tutorials_root):
        parser.error("--config/--entry-kind/--tutorials-root are not valid with --run-document")
    if not args.run_document and not args.entry:
        parser.error("--entry is required (or use --run-document with action=run/step)")

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
        )
        print(json.dumps(report.to_json(), indent=2))
        return 0 if report.status == "ok" else 1

    if args.action == "step":
        if not (args.strict or args.run_document):
            parser.error("action=step requires --strict or --run-document")
        if not args.step:
            parser.error("action=step requires --step <id>")
        if args.run_document:
            return _run_document_dispatch(args)
        report = strict_plan(
            selected_entry,
            entry_kind=args.entry_kind,
            overrides=overrides,
            config_path=args.config,
        )
        if report.status != "ok":
            print(json.dumps(report.to_json(), indent=2))
            return 1
        if report.workflow_dag is None or report.workflow_state is None:
            print(json.dumps({
                "status": "failed",
                "error": "strict plan did not produce workflow_dag and workflow_state",
            }, indent=2))
            return 1
        return _execute_step(
            entry_label=selected_entry,
            step_id=args.step,
            workflow_dag=report.workflow_dag,
            planned_state=report.workflow_state,
            case_root=Path(report.launch["case_root"]),
            output_dir=Path(report.launch["output_dir"]),
            expected_artifacts=report.expected_artifacts,
            tail_lines=args.tail_lines,
        )

    if args.action == "run":
        if not (args.strict or args.run_document):
            parser.error("action=run requires --strict or --run-document")
        if args.run_document:
            return _run_document_dispatch(args)
        report = strict_plan(
            selected_entry,
            entry_kind=args.entry_kind,
            overrides=overrides,
            config_path=args.config,
        )
        if report.status != "ok":
            print(json.dumps(report.to_json(), indent=2))
            return 1
        if report.workflow_dag is None or report.workflow_state is None:
            print(json.dumps({
                "status": "failed",
                "error": "strict plan did not produce workflow_dag and workflow_state",
            }, indent=2))
            return 1
        return _execute_run(
            entry_label=selected_entry,
            workflow_dag=report.workflow_dag,
            planned_state=report.workflow_state,
            case_root=Path(report.launch["case_root"]),
            output_dir=Path(report.launch["output_dir"]),
            expected_artifacts=report.expected_artifacts,
            tail_lines=args.tail_lines,
        )

    spec = load_entry_spec(selected_entry, entry_kind=args.entry_kind, overrides=overrides)
    engine = DriverEngine(
        spec=spec,
        dry_run=args.dry_run,
        continue_on_error=args.continue_on_error,
        requested_action=args.action,
    )

    if args.action == "sim":
        engine.run_simulations()
    elif args.action == "post":
        engine.run_postprocess()
    else:
        engine.run_all()

    return 0
