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
#     sweep_runner
#
# Description
#     Orchestrates sweep-plan: expand, materialize, and audit each resolved case.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

from __future__ import annotations

import json
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from ...strict_planning import strict_plan
from ...sweep_derivation_catalog import get_derivation
from ...sweep_expansion import SweepValidationError, check_case_count_cap, expand_sweep
from ...sweep_materialize import materialize_case
from ...sweep_routing import route_case_values
from .sweep_manifest import (
    CaseManifestEntry,
    SweepManifest,
    compute_spec_hash,
    compute_override_hash,
    write_manifest,
    read_manifest,
)


def _load_spec(spec_path: str | Path) -> dict[str, Any]:
    return json.loads(Path(spec_path).read_text())


def sweep_plan(
    spec_path: str | Path,
    *,
    output_dir: str | Path,
    max_cases: int = 200,
) -> dict[str, Any]:
    sweep_spec = _load_spec(spec_path)
    check_case_count_cap(sweep_spec, max_cases=max_cases)

    output_dir = Path(output_dir)
    resolved_cases = expand_sweep(sweep_spec, get_derivation=get_derivation)
    base = sweep_spec.get("base", {})

    case_reports = []
    for case in resolved_cases:
        case_dir = output_dir / case.case_id
        try:
            routed = route_case_values(base=base, resolved_axis_values=case.resolved_axis_values)
            materialize_case(case_dir=case_dir, routed=routed)
        except (OSError, ValueError) as exc:
            case_reports.append(
                {
                    "case_id": case.case_id,
                    "resolved_axis_values": case.resolved_axis_values,
                    "status": "failed",
                    "materialization_error": str(exc),
                }
            )
            continue

        report = strict_plan(case.case_id, entry_kind="case_folder", overrides={"tutorials_root": str(output_dir)})
        report_payload = report.to_json()
        case_reports.append(
            {
                "case_id": case.case_id,
                "resolved_axis_values": case.resolved_axis_values,
                "status": report.status,
                "plan": report_payload,
            }
        )

    return {"case_count": len(resolved_cases), "cases": case_reports}


def _now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _workflow_state_path_from_run_document(run_document: dict[str, Any]) -> Path:
    try:
        output_dir = run_document["launch"]["outputDir"]
    except KeyError as exc:
        raise ValueError("strict_plan run_document is missing launch.outputDir") from exc
    return Path(output_dir) / "workflow_state.json"


def sweep_run(
    spec_path: str | Path,
    *,
    output_dir: str | Path,
    max_cases: int = 200,
    retry_failed: bool = False,
    case_timeout_s: float | None = None,
) -> dict[str, Any]:
    sweep_spec = _load_spec(spec_path)
    check_case_count_cap(sweep_spec, max_cases=max_cases)

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = output_dir / "sweep_manifest.json"

    spec_hash = compute_spec_hash(sweep_spec)
    existing_status_by_case: dict[str, str] = {}
    existing_entry_by_case = {}
    if manifest_path.exists():
        existing = read_manifest(manifest_path)
        if existing.sweep_spec_hash != spec_hash:
            raise SweepValidationError(
                "sweep.json has changed since this output directory was created "
                f"(hash mismatch: expected {existing.sweep_spec_hash}, got {spec_hash}); "
                "spec changed — use a fresh --output-dir or resolve the mismatch."
            )
        existing_status_by_case = {c.case_id: c.status for c in existing.cases}
        existing_entry_by_case = {c.case_id: c for c in existing.cases}

    resolved_cases = expand_sweep(sweep_spec, get_derivation=get_derivation)
    base = sweep_spec.get("base", {})

    manifest = SweepManifest(
        schema_version="1.0", sweep_spec_hash=spec_hash,
        created_at=_now(), updated_at=_now(), cases=[],
    )

    completed_count = 0
    failed_count = 0
    skipped_count = 0
    case_summaries: list[dict[str, Any]] = []

    for case in resolved_cases:
        case_dir = output_dir / case.case_id
        run_document_path = case_dir / "run_document.json"
        workflow_state_path = case_dir / "postProcessing" / "workflow_state.json"

        prior_status = existing_status_by_case.get(case.case_id)
        prior_entry = existing_entry_by_case.get(case.case_id)
        outcome = "fresh"
        materialization_error = None
        plan_error = None
        timeout_error = None

        routing_error: str | None = None
        try:
            routed = route_case_values(base=base, resolved_axis_values=case.resolved_axis_values)
        except (OSError, ValueError) as exc:
            # An unrecognized/unroutable axis (e.g. "dx") is a per-case
            # failure, not a crash of the whole sweep -- same treatment as a
            # materialize_case failure below.
            routed = {}
            routing_error = str(exc)

        if routing_error is not None:
            status = "failed"
            materialization_error = routing_error
            failed_count += 1
        elif prior_status == "completed":
            outcome = "skipped"
            skipped_count += 1
            completed_count += 1
            status = "completed"
            if prior_entry is not None:
                workflow_state_path = output_dir / prior_entry.workflow_state_path
                run_document_path = output_dir / prior_entry.run_document_path
        elif prior_status == "failed" and not retry_failed:
            status = "failed"
            failed_count += 1
            if prior_entry is not None:
                workflow_state_path = output_dir / prior_entry.workflow_state_path
                run_document_path = output_dir / prior_entry.run_document_path
        else:
            if prior_status == "failed" and retry_failed:
                outcome = "retried"
            status = "failed"
            try:
                materialize_case(case_dir=case_dir, routed=routed)
                report = strict_plan(case.case_id, entry_kind="case_folder", overrides={"tutorials_root": str(output_dir)})
                payload = report.to_json()
                if report.status != "ok":
                    plan_error = "strict_plan reported failed status"
                else:
                    run_document = payload["run_document"]
                    workflow_state_path = _workflow_state_path_from_run_document(run_document)
                    run_document_path.write_text(json.dumps(run_document, indent=2))
            except (OSError, ValueError) as exc:
                materialization_error = str(exc)
            except Exception as exc:
                plan_error = str(exc)
            else:
                if plan_error is None:
                    try:
                        result = subprocess.run(
                            [sys.executable, "-m", "openfoam_driver", "run", "--run-document", str(run_document_path)],
                            capture_output=True, text=True,
                            timeout=case_timeout_s,
                        )
                    except subprocess.TimeoutExpired as exc:
                        # A hung case must not block the whole serial sweep: mark
                        # it failed and continue. The manifest stays resumable.
                        status = "failed"
                        timeout_error = (
                            f"case exceeded timeout of {case_timeout_s}s "
                            f"and was terminated: {exc}"
                        )
                    else:
                        if workflow_state_path.exists():
                            state = json.loads(workflow_state_path.read_text())
                            status = state.get("status", "pending")
                        elif result.returncode != 0:
                            status = "failed"
                        else:
                            status = "pending"
            if status == "completed":
                completed_count += 1
            else:
                failed_count += 1

        case_summary = {
            "case_id": case.case_id,
            "status": status,
            "outcome": outcome,
            "run_document_path": str(run_document_path.relative_to(output_dir)),
            "workflow_state_path": str(workflow_state_path.relative_to(output_dir)),
        }
        if materialization_error is not None:
            case_summary["materialization_error"] = materialization_error
        if plan_error is not None:
            case_summary["plan_error"] = plan_error
        if timeout_error is not None:
            case_summary["timeout_error"] = timeout_error
        case_summaries.append(case_summary)

        manifest.cases.append(
            CaseManifestEntry(
                case_id=case.case_id,
                resolved_axis_values=case.resolved_axis_values,
                override_hash=compute_override_hash(routed),
                run_document_path=str(run_document_path.relative_to(output_dir)),
                workflow_state_path=str(workflow_state_path.relative_to(output_dir)),
                status=status,
                outcome=outcome,
                started_at=_now(),
                updated_at=_now(),
            )
        )
        manifest.updated_at = _now()
        write_manifest(manifest_path, manifest)

    return {
        "case_count": len(resolved_cases),
        "completed_count": completed_count,
        "failed_count": failed_count,
        "skipped_count": skipped_count,
        "cases": case_summaries,
    }
