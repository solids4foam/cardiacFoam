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
import os
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from openfoam_driver.core.strict_planning import strict_plan
from openfoam_driver.core.sweep.sweep_derivation_catalog import get_derivation
from openfoam_driver.core.sweep.sweep_expansion import SweepValidationError, check_case_count_cap, expand_sweep
from openfoam_driver.sweep_materialize import materialize_case
from openfoam_driver.sweep_routing import route_case_values, route_entry_case_values
from .fresh import ensure_fresh_output_dir
from .openfoam_environment import configure_plugin_environment
from .output_collection import collect_new_outputs, snapshot_postprocessing
from .postprocess_phase import build_sweep_context, run_postprocessing_module
from .registry import load_entry_spec
from .run_document_exec import _allowed_runs_root
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


def _entry_name(sweep_spec: dict[str, Any]) -> str | None:
    return sweep_spec.get("base", {}).get("entry")


def _relative_or_absolute(path: Path, base: Path) -> str:
    """Path relative to `base` when possible, else the absolute path.

    `run_document_path` is always written under the sweep's own
    `--output-dir` (safe to make relative). `workflow_state_path` is not:
    in entry mode it comes from `launch.outputDir`, which resolves to the
    target tutorial's own case_root/output_dir_name -- a directory tree
    entirely unrelated to the sweep's --output-dir (confirmed via a real,
    non-mocked sweep-run: Path.relative_to raised ValueError there). Record
    the absolute path in that case rather than crash the whole sweep over a
    manifest cosmetic.
    """
    try:
        return str(path.relative_to(base))
    except ValueError:
        return str(path)


def _is_stale_time_dir_name(name: str) -> bool:
    """True for an OpenFOAM time-directory name other than the literal '0'.

    '0' is a real, on-disk initial condition for some tutorials and must
    never be touched here. Every other numeric name (e.g. '0.0982143',
    written by a prior case's reconstructPar) is always solve *output*,
    never an input -- safe to clear before the next case reuses this
    case_root.
    """
    if name == "0":
        return False
    try:
        float(name)
    except ValueError:
        return False
    return True


def _clean_stale_time_directories(case_root: Path) -> None:
    """Remove a prior case's leftover reconstructed time directories.

    Entry-based sweeps reuse one shared case_root across cases (see
    _materialize_entry_case's docstring). decomposePar's -force flag already
    clears stale processor*/ dirs, but a case with no real 0/ (e.g. a
    manufactured-solution verifier, whose IC the solver computes rather than
    reads from disk) has nothing to anchor OpenFOAM's -time selector to --
    "-time <value>" matches the *nearest* existing time, not an exact one
    (confirmed via timeSelector.C), so a leftover time directory from the
    previous case would silently get decomposed onto this case's new mesh
    instead. Clearing them here, before the case that would otherwise read
    them, is the fix.
    """
    if not case_root.is_dir():
        return
    for child in case_root.iterdir():
        if child.is_dir() and _is_stale_time_dir_name(child.name):
            shutil.rmtree(child)


def _materialize_entry_case(
    entry: str,
    routed: dict[str, Any],
    *,
    staging_root: Path | None = None,
    driver_context=None,
) -> dict[str, Any]:
    """Materialize one entry-based sweep case via the tutorial's own spec.

    Entry-based sweeps target an existing registered tutorial whose
    apply_case()/build_cases() mutate a case root in place.  That root must be
    a disposable staging copy, never the checked-in tutorial directory.  The
    returned overrides point the entry at that staged case so strict planning
    and execution use exactly the same paths.

    ``staging_root`` is optional for compatibility with low-level callers and
    tests that provide an already-isolated fake spec.  Real sweep callers
    always pass it.

    Raises ValueError if the resolved overrides don't collapse to exactly one
    case -- the sweep model is one case per resolved axis combination.
    """
    spec = load_entry_spec(entry, overrides=routed, driver_context=driver_context)
    effective_routed = dict(routed)
    if staging_root is not None and spec.case_root.exists():
        source_case_root = Path(spec.case_root).resolve()
        staged_case_root = Path(staging_root).resolve()
        _stage_entry_case(source_case_root, staged_case_root)
        # ``make_spec`` resolves case_root as tutorials_root/case_dir_name.
        # Redirect both values together; changing only tutorials_root would
        # leave a nested original case_dir_name and recreate the source tree
        # below the scratch directory.
        effective_routed["tutorials_root"] = str(staged_case_root.parent)
        effective_routed["case_dir_name"] = staged_case_root.name
        staged_spec = load_entry_spec(
            entry, overrides=effective_routed, driver_context=driver_context,
        )
        # A real registered factory consumes tutorials_root/case_dir_name and
        # therefore returns the staged path.  Keep compatibility with test
        # doubles and third-party factories that intentionally return their
        # own fixed spec regardless of overrides.
        if Path(staged_spec.case_root).resolve() != staged_case_root:
            effective_routed = dict(routed)
        else:
            spec = staged_spec
    cases = spec.build_cases()
    if len(cases) != 1:
        raise ValueError(
            f"entry-based sweep axis combination resolved to {len(cases)} cases "
            f"for entry '{entry}'; expected exactly 1 -- add enough constraining "
            "overrides (e.g. 'solvers') to collapse this combination to a single case"
        )
    _clean_stale_time_directories(spec.case_root)
    spec.apply_case(spec.case_root, cases[0])
    return effective_routed


def _stage_entry_case(source_case_root: Path, staged_case_root: Path) -> None:
    """Copy a registered case into scratch storage without old run output.

    Registered tutorial folders contain source dictionaries and scripts next
    to OpenFOAM's generated mesh, time, processor, log, and post-processing
    trees.  Copying those generated trees would reintroduce the stale-state
    bug this staging boundary is meant to prevent, so the filter is explicit
    and conservative: keep authored inputs (including ``0/``) and omit only
    known derived artifacts.
    """
    if staged_case_root.exists():
        shutil.rmtree(staged_case_root)

    generated_dir_names = {
        "postProcessing", "logs", "workflow_logs", "cachedCasePostProcessing",
        "polyMesh", "archivedPostProcessing", "results", "data",
    }
    generated_file_names = {
        "workflow_state.json", "run_document.json", "sweep_manifest.json",
    }

    def ignore_generated(_directory: str, names: list[str]) -> set[str]:
        ignored: set[str] = set()
        for name in names:
            candidate = Path(_directory) / name
            # A previous driverFOAM case can have a descriptive directory name
            # (for example ``gauss_linear_40_*``) rather than a numeric
            # OpenFOAM time name. Its workflow markers are the reliable
            # boundary between authored tutorial content and generated case
            # content, so omit the whole directory when they are present.
            if candidate.is_dir() and any(
                (candidate / marker).exists()
                for marker in ("workflow_state.json", "workflow_logs", "run_document.json")
            ):
                ignored.add(name)
                continue
            if name in generated_dir_names or name in generated_file_names:
                ignored.add(name)
                continue
            if name.startswith("processor") or name.startswith("driverPostProcessingArchive"):
                ignored.add(name)
                continue
            if (
                name.startswith("log.")
                or name.endswith(".foam")
                or name.endswith(".msh")
                or (name.endswith(".geo") and not name.endswith(".geo.template"))
            ):
                ignored.add(name)
                continue
            path = Path(name)
            # Numeric OpenFOAM time directories are generated output.  Keep
            # the authored initial-condition directory ``0``.
            if path.name != "0":
                try:
                    float(path.name)
                except ValueError:
                    pass
                else:
                    ignored.add(name)
        return ignored

    if not source_case_root.is_dir():
        raise FileNotFoundError(f"Registered case root does not exist: {source_case_root}")
    staged_case_root.parent.mkdir(parents=True, exist_ok=True)
    shutil.copytree(
        source_case_root,
        staged_case_root,
        ignore=ignore_generated,
        symlinks=True,
    )


def sweep_plan(
    spec_path: str | Path,
    *,
    output_dir: str | Path,
    max_cases: int = 200,
    driver_context=None,
) -> dict[str, Any]:
    from ..compatibility import resolve_public_driver_context

    driver_context = resolve_public_driver_context(driver_context)
    sweep_spec = _load_spec(spec_path)
    check_case_count_cap(sweep_spec, max_cases=max_cases)

    output_dir = Path(output_dir)
    resolved_cases = expand_sweep(sweep_spec, get_derivation=get_derivation)
    base = sweep_spec.get("base", {})
    entry = _entry_name(sweep_spec)

    case_reports = []
    for case in resolved_cases:
        try:
            if entry is not None:
                routed = route_entry_case_values(base=base, resolved_axis_values=case.resolved_axis_values)
                routed = _materialize_entry_case(
                    entry,
                    routed,
                    staging_root=output_dir / "cases" / case.case_id,
                    driver_context=driver_context,
                )
            else:
                routed = route_case_values(
                    base=base,
                    resolved_axis_values=case.resolved_axis_values,
                    driver_context=driver_context,
                )
                materialize_case(case_dir=output_dir / case.case_id, routed=routed)
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

        if entry is not None:
            report = strict_plan(entry, overrides=routed, driver_context=driver_context)
        else:
            report = strict_plan(
                case.case_id,
                entry_kind="case_folder",
                overrides={"tutorials_root": str(output_dir)},
                driver_context=driver_context,
            )
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
    fresh: bool = False,
    task: str = "summarize",
    driver_context=None,
) -> dict[str, Any]:
    """`task` plays no part in the sweep loop itself -- expanding, routing,
    materializing, and running cases is fully deterministic and has no use
    for it. It is only consumed at the very end, handed to
    run_postprocessing_module: the sweep is task(sweep), no reasoning
    involved; the postprocess hand-off is where a task actually matters.
    """
    from ..compatibility import resolve_public_driver_context

    driver_context = resolve_public_driver_context(driver_context)
    execution_environment = configure_plugin_environment(
        os.environ,
        driver_context,
    ).env
    sweep_spec = _load_spec(spec_path)
    check_case_count_cap(sweep_spec, max_cases=max_cases)

    output_dir = Path(output_dir)
    fresh_error = ensure_fresh_output_dir(
        output_dir, fresh=fresh, allowed_root=_allowed_runs_root(),
    )
    if fresh_error is not None:
        raise SweepValidationError(fresh_error)
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
    entry = _entry_name(sweep_spec)

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
        case_record_path = case_dir / "case_record.json"

        prior_status = existing_status_by_case.get(case.case_id)
        prior_entry = existing_entry_by_case.get(case.case_id)
        outcome = "fresh"
        materialization_error = None
        plan_error = None
        timeout_error = None

        routing_error: str | None = None
        try:
            if entry is not None:
                routed = route_entry_case_values(base=base, resolved_axis_values=case.resolved_axis_values)
            else:
                routed = route_case_values(
                    base=base,
                    resolved_axis_values=case.resolved_axis_values,
                    driver_context=driver_context,
                )
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
                if entry is not None:
                    routed = _materialize_entry_case(
                        entry,
                        routed,
                        staging_root=output_dir / "cases" / case.case_id,
                        driver_context=driver_context,
                    )
                    report = strict_plan(entry, overrides=routed, driver_context=driver_context)
                else:
                    materialize_case(case_dir=case_dir, routed=routed)
                    report = strict_plan(
                        case.case_id,
                        entry_kind="case_folder",
                        overrides={"tutorials_root": str(output_dir)},
                        driver_context=driver_context,
                    )
                payload = report.to_json()
                if report.status != "ok":
                    plan_error = "strict_plan reported failed status"
                else:
                    run_document = payload["run_document"]
                    workflow_state_path = _workflow_state_path_from_run_document(run_document)
                    # In entry mode, case_dir (this sweep's own bookkeeping
                    # location for run_document.json) is unrelated to the
                    # tutorial's real case_root and is never created by
                    # _materialize_entry_case, unlike generic mode's
                    # materialize_case which creates it as a side effect.
                    run_document_path.parent.mkdir(parents=True, exist_ok=True)
                    run_document_path.write_text(json.dumps(run_document, indent=2))
            except (OSError, ValueError) as exc:
                materialization_error = str(exc)
            except Exception as exc:
                plan_error = str(exc)
            else:
                if plan_error is None:
                    # Entry-mode cases share one case_root, so OpenFOAM's own
                    # output (mesh, solved field time-dirs, postProcessing/)
                    # all lands in that one shared case_root/postProcessing/
                    # -- overwritten by each subsequent case, since OpenFOAM
                    # itself has no notion of driverFOAM's per-case
                    # output_dir_name. Snapshot it now so collect_new_outputs
                    # below (called once this case's own output_dir is known)
                    # can tell this case's own new/changed output apart from
                    # anything left over, and archive it into that same
                    # case's own output_dir_name folder -- the same directory
                    # workflow_state.json lives in -- so nothing needs a
                    # separate cache location to find it later.
                    archive_dir_name = (base.get("archive_dir_name") or "collectedOutput") if entry is not None else None
                    pp_before: dict[str, tuple[float, int]] = {}
                    if archive_dir_name:
                        case_root_for_archive = Path(run_document["launch"]["caseRoot"])
                        pp_before = snapshot_postprocessing(case_root_for_archive)
                    try:
                        if workflow_state_path.exists():
                            workflow_state_path.unlink()
                        result = subprocess.run(
                            [sys.executable, "-m", "openfoam_driver", "run", "--run-document", str(run_document_path)],
                            capture_output=True, text=True,
                            env=execution_environment,
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
                        if archive_dir_name and workflow_state_path.exists():
                            collect_new_outputs(
                                case_root_for_archive,
                                pp_before,
                                workflow_state_path.parent / archive_dir_name,
                                label=case.case_id,
                            )
            if status == "completed":
                completed_count += 1
            else:
                failed_count += 1

        case_summary = {
            "case_id": case.case_id,
            "status": status,
            "outcome": outcome,
            "run_document_path": str(run_document_path.relative_to(output_dir)),
            "workflow_state_path": _relative_or_absolute(workflow_state_path, output_dir),
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
                workflow_state_path=_relative_or_absolute(workflow_state_path, output_dir),
                status=status,
                outcome=outcome,
                started_at=_now(),
                updated_at=_now(),
                case_record_path=str(case_record_path.relative_to(output_dir)),
            )
        )
        manifest.updated_at = _now()
        write_manifest(manifest_path, manifest)

    context = build_sweep_context(output_dir)
    if failed_count == 0:
        postprocess = run_postprocessing_module(context, task=task).to_json()
    else:
        postprocess = {
            "status": "skipped",
            "message": f"sweep had {failed_count} failed case(s); postprocess not run",
        }

    return {
        "case_count": len(resolved_cases),
        "completed_count": completed_count,
        "failed_count": failed_count,
        "skipped_count": skipped_count,
        "cases": case_summaries,
        "postprocess": postprocess,
    }
