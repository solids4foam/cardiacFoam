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
from .core.runtime.workflow_runner import run_workflow_step
from .core.runtime.workflow_orchestrator import run_workflow
from .core.runtime.workflow_state import workflow_state_from_json
from .core.runtime.registry import ENTRY_KIND_VALUES, list_tutorials, load_entry_spec
from .introspection import describe_entry
from .specs.common import default_setup_dir_name
from .strict_planning import strict_plan


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


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Generic OpenFOAM tutorial automation driver")
    parser.add_argument(
        "action",
        choices=["sim", "post", "all", "describe", "plan", "step", "run"],
        help="Pipeline stage to execute",
    )
    parser.add_argument(
        "--entry",
        required=True,
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
            "manufacturedEikonalECG, restitutionCurves, genericCase/randomCase) "
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
        if not args.strict:
            parser.error("action=step currently requires --strict")
        if not args.step:
            parser.error("action=step requires --step <id>")
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
        case_root = Path(report.launch["case_root"])
        output_dir = Path(report.launch["output_dir"])
        state_path = output_dir / "workflow_state.json"
        workflow_state = report.workflow_state
        if state_path.exists():
            try:
                workflow_state = workflow_state_from_json(json.loads(state_path.read_text()))
            except Exception as exc:
                print(json.dumps({
                    "status": "failed",
                    "entry": selected_entry,
                    "step": args.step,
                    "error": f"Could not read existing workflow state: {exc}",
                    "workflow_state_path": str(state_path),
                }, indent=2))
                return 1
        try:
            result = run_workflow_step(
                report.workflow_dag,
                workflow_state,
                args.step,
                case_root=case_root,
                log_dir=output_dir / "workflow_logs",
                state_path=state_path,
                expected_artifacts=report.expected_artifacts,
            )
        except Exception as exc:
            print(json.dumps({
                "status": "failed",
                "entry": selected_entry,
                "step": args.step,
                "error": str(exc),
                "workflow_state": report.workflow_state.to_json(),
            }, indent=2))
            return 1
        payload = {
            "status": "ok" if result.exit_code == 0 else "failed",
            "entry": selected_entry,
            "step": args.step,
            "exit_code": result.exit_code,
            "stdout_log": result.stdout_log,
            "stderr_log": result.stderr_log,
            "workflow_state_path": str(state_path),
            "workflow_state": result.state.to_json(),
        }
        print(json.dumps(payload, indent=2))
        return 0 if result.exit_code == 0 else 1

    if args.action == "run":
        if not args.strict:
            parser.error("action=run currently requires --strict")
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
        case_root = Path(report.launch["case_root"])
        output_dir = Path(report.launch["output_dir"])
        state_path = output_dir / "workflow_state.json"
        workflow_state = report.workflow_state
        if state_path.exists():
            try:
                workflow_state = workflow_state_from_json(json.loads(state_path.read_text()))
            except Exception as exc:
                print(json.dumps({
                    "status": "failed",
                    "entry": selected_entry,
                    "error": f"Could not read existing workflow state: {exc}",
                    "workflow_state_path": str(state_path),
                }, indent=2))
                return 1
        if workflow_state.status == "failed":
            print(json.dumps({
                "status": "failed",
                "entry": selected_entry,
                "error": "workflow_state is failed; use action=step to rerun a failed step explicitly",
                "workflow_state_path": str(state_path),
                "workflow_state": workflow_state.to_json(),
            }, indent=2))
            return 1
        try:
            outcome = run_workflow(
                report.workflow_dag,
                workflow_state,
                case_root=case_root,
                output_dir=output_dir,
                expected_artifacts=report.expected_artifacts,
                state_path=state_path,
            )
        except Exception as exc:
            try:
                error_state = workflow_state_from_json(json.loads(state_path.read_text()))
            except Exception:
                error_state = workflow_state
            print(json.dumps(_step_payload(
                status="failed",
                entry=selected_entry,
                step=error_state.current_step_id or workflow_state.current_step_id,
                workflow_state_path=state_path,
                workflow_state=error_state.to_json(),
                error=str(exc),
            ), indent=2))
            return 1
        workflow_state = outcome.state
        results = list(outcome.steps)
        status = "ok" if workflow_state.status == "completed" else "failed"
        if workflow_state.status == "pending" and workflow_state.current_step_id is None:
            status = "failed"
        payload = {
            "status": status,
            "entry": selected_entry,
            "steps": results,
            "workflow_state_path": str(state_path),
            "workflow_state": workflow_state.to_json(),
        }
        if workflow_state.status == "pending" and workflow_state.current_step_id is None:
            payload["error"] = "workflow_state is pending but has no current_step_id"
        print(json.dumps(payload, indent=2))
        return 0 if status == "ok" else 1

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
