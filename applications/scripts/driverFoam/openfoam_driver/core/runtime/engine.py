"""DriverEngine: orchestrates tutorial sweeps and emits agent-readable artifacts.

Agent-facing contracts emitted under ``spec.output_dir`` (or ``setup_root``
during dry runs):

* ``run_manifest.json`` — single document, fully rewritten on every state
  change. Writes go through a sibling ``run_manifest.json.tmp`` file and an
  ``os.replace`` so polling readers never observe a torn document. Schema is
  versioned via the ``schema_version`` key; v2.x is additive-only.
* ``action_events.jsonl`` — append-only event log. Each event is a single
  ``json.dumps(event) + "\n"`` write, so a line is either fully present or
  not present at all. Agents may ``tail -F`` the file safely.
* ``run_report.md`` — human-readable companion to ``run_manifest.json``.
"""

from __future__ import annotations

import json
import os
import time
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
from uuid import uuid4

from .models import CaseConfig, TutorialSpec


@dataclass
class CaseResult:
    case_id: str
    status: str
    duration_s: float
    params: dict
    error: str | None = None
    index: int | None = None
    total_cases: int | None = None
    started_at_utc: str | None = None
    finished_at_utc: str | None = None


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


class DriverEngine:
    """Shared execution engine for all tutorial specs."""

    def __init__(
        self,
        spec: TutorialSpec,
        dry_run: bool = False,
        continue_on_error: bool = False,
        requested_action: str = "sim",
    ) -> None:
        self.spec = spec
        self.dry_run = dry_run
        self.continue_on_error = continue_on_error
        self.requested_action = requested_action
        self.run_id = f"{int(time.time())}-{uuid4().hex[:8]}"
        self.started_at_utc: str | None = None
        self.finished_at_utc: str | None = None

    def run_simulations(self) -> list[CaseResult]:
        cases = self.spec.build_cases()
        total = len(cases)
        results: list[CaseResult] = []
        self.started_at_utc = self.started_at_utc or _utc_now()

        print(f"Entry: {self._entry_name()}")
        print(f"Case root: {self.spec.case_root}")
        print(f"Planned simulations: {total}")
        self._write_action_event(
            "sim_started",
            {
                "entry": self._entry_name(),
                "entry_kind": self._entry_kind(),
                "total_cases": total,
                "dry_run": self.dry_run,
            },
        )
        self._write_manifest(
            results,
            status="running",
            current_case_id=None,
            total_cases=total,
            postprocess_status="not_started",
        )

        for idx, case in enumerate(cases, start=1):
            print("-" * 72)
            print(f"[{idx}/{total}] {case.case_id}")
            print(f"Params: {case.params}")
            self._write_action_event(
                "case_started",
                {
                    "entry": self._entry_name(),
                    "entry_kind": self._entry_kind(),
                    "case_id": case.case_id,
                    "index": idx,
                    "total": total,
                },
            )
            self._write_manifest(
                results,
                status="running",
                current_case_id=case.case_id,
                total_cases=total,
                postprocess_status="not_started",
            )

            if self.dry_run:
                results.append(
                    CaseResult(
                        case_id=case.case_id,
                        status="planned",
                        duration_s=0.0,
                        params=case.params,
                        index=idx,
                        total_cases=total,
                        started_at_utc=None,
                        finished_at_utc=None,
                    )
                )
                self._write_manifest(
                    results,
                    status="running",
                    current_case_id=None,
                    total_cases=total,
                    postprocess_status="not_started",
                )
                self._write_action_event(
                    "case_finished",
                    {
                        "entry": self._entry_name(),
                        "entry_kind": self._entry_kind(),
                        "case_id": case.case_id,
                        "status": "planned",
                    },
                )
                continue

            started_at_utc = _utc_now()
            start = time.time()
            try:
                self.spec.apply_case(self.spec.case_root, case)
                self.spec.run_case(self.spec.case_root, self.spec.setup_root, case)
                duration = time.time() - start
                results.append(
                    CaseResult(
                        case_id=case.case_id,
                        status="ok",
                        duration_s=duration,
                        params=case.params,
                        index=idx,
                        total_cases=total,
                        started_at_utc=started_at_utc,
                        finished_at_utc=_utc_now(),
                    )
                )
                self._write_manifest(
                    results,
                    status="running",
                    current_case_id=None,
                    total_cases=total,
                    postprocess_status="not_started",
                )
                self._write_action_event(
                    "case_finished",
                    {
                        "entry": self._entry_name(),
                        "entry_kind": self._entry_kind(),
                        "case_id": case.case_id,
                        "status": "ok",
                        "duration_s": duration,
                    },
                )
            except Exception as exc:
                duration = time.time() - start
                results.append(
                    CaseResult(
                        case_id=case.case_id,
                        status="failed",
                        duration_s=duration,
                        params=case.params,
                        error=str(exc),
                        index=idx,
                        total_cases=total,
                        started_at_utc=started_at_utc,
                        finished_at_utc=_utc_now(),
                    )
                )
                print(f"Simulation failed for {case.case_id}: {exc}")
                failure_status = "failed" if not self.continue_on_error else "running"
                self._write_manifest(
                    results,
                    status=failure_status,
                    current_case_id=None,
                    total_cases=total,
                    postprocess_status="not_started",
                    error=str(exc) if not self.continue_on_error else None,
                )
                self._write_action_event(
                    "case_finished",
                    {
                        "entry": self._entry_name(),
                        "entry_kind": self._entry_kind(),
                        "case_id": case.case_id,
                        "status": "failed",
                        "error": str(exc),
                    },
                )
                if not self.continue_on_error:
                    raise

        if not self.dry_run:
            self.spec.output_dir.mkdir(parents=True, exist_ok=True)
            if self.spec.collect_outputs is not None:
                self.spec.collect_outputs(self.spec.case_root, self.spec.output_dir)

        self.finished_at_utc = _utc_now()
        final_status = "planned" if self.dry_run else (
            "completed_with_failures" if any(item.status == "failed" for item in results) else "completed"
        )
        postprocess_status = "skipped" if self.dry_run else "not_started"
        self._write_manifest(
            results,
            status=final_status,
            current_case_id=None,
            total_cases=total,
            postprocess_status=postprocess_status,
            finished_at_utc=self.finished_at_utc,
        )
        self._write_action_event(
            "sim_finished",
            {
                "entry": self._entry_name(),
                "entry_kind": self._entry_kind(),
                "status": final_status,
            },
        )
        return results

    def run_postprocess(self) -> None:
        if self.spec.postprocess is None:
            print(f"No post-processing hook is defined for entry '{self._entry_name()}'.")
            return
        self.started_at_utc = self.started_at_utc or _utc_now()
        self.spec.output_dir.mkdir(parents=True, exist_ok=True)
        self._write_action_event(
            "postprocess_started",
            {"entry": self._entry_name(), "entry_kind": self._entry_kind()},
        )
        self._write_manifest(
            [],
            status="postprocessing",
            current_case_id=None,
            total_cases=0,
            postprocess_status="running",
        )
        try:
            if self.spec.collect_outputs is not None:
                self.spec.collect_outputs(self.spec.case_root, self.spec.output_dir)
            self.spec.postprocess(self.spec.setup_root, self.spec.output_dir)
        except Exception as exc:
            self.finished_at_utc = _utc_now()
            self._write_manifest(
                [],
                status="postprocess_failed",
                current_case_id=None,
                total_cases=0,
                postprocess_status="failed",
                error=str(exc),
                finished_at_utc=self.finished_at_utc,
            )
            self._write_action_event(
                "postprocess_finished",
                {
                    "entry": self._entry_name(),
                    "entry_kind": self._entry_kind(),
                    "status": "failed",
                    "error": str(exc),
                },
            )
            raise

        self.finished_at_utc = _utc_now()
        self._write_manifest(
            [],
            status="completed",
            current_case_id=None,
            total_cases=0,
            postprocess_status="completed",
            finished_at_utc=self.finished_at_utc,
        )
        self._write_action_event(
            "postprocess_finished",
            {"entry": self._entry_name(), "entry_kind": self._entry_kind(), "status": "ok"},
        )

    def run_all(self) -> list[CaseResult]:
        self._write_action_event(
            "all_started",
            {"entry": self._entry_name(), "entry_kind": self._entry_kind()},
        )
        results = self.run_simulations()
        if self.dry_run:
            print("Dry-run enabled. Post-processing was skipped.")
            self._write_action_event(
                "all_finished",
                {
                    "entry": self._entry_name(),
                    "entry_kind": self._entry_kind(),
                    "status": "planned",
                },
            )
            return results
        self._write_manifest(
            results,
            status="postprocessing",
            current_case_id=None,
            total_cases=len(results),
            postprocess_status="running",
        )
        try:
            self.spec.postprocess(self.spec.setup_root, self.spec.output_dir)
        except Exception as exc:
            self.finished_at_utc = _utc_now()
            self._write_manifest(
                results,
                status="postprocess_failed",
                current_case_id=None,
                total_cases=len(results),
                postprocess_status="failed",
                error=str(exc),
                finished_at_utc=self.finished_at_utc,
            )
            raise

        self.finished_at_utc = _utc_now()
        final_status = (
            "completed_with_failures" if any(item.status == "failed" for item in results) else "completed"
        )
        self._write_manifest(
            results,
            status=final_status,
            current_case_id=None,
            total_cases=len(results),
            postprocess_status="completed",
            finished_at_utc=self.finished_at_utc,
        )
        self._write_action_event(
            "all_finished",
            {
                "entry": self._entry_name(),
                "entry_kind": self._entry_kind(),
                "status": final_status,
            },
        )
        return results

    def _write_action_event(self, event_type: str, payload: dict | None = None) -> None:
        destination_root = self._manifest_destination_root()
        destination_root.mkdir(parents=True, exist_ok=True)
        events_path = destination_root / "action_events.jsonl"
        event = {
            "event": event_type,
            "timestamp_utc": _utc_now(),
            "run_id": self.run_id,
        }
        if payload is not None:
            event.update(payload)
        with events_path.open("a") as handle:
            handle.write(json.dumps(event) + "\n")

    def _manifest_destination_root(self) -> Path:
        return self.spec.setup_root if self.dry_run else self.spec.output_dir

    def _entry_name(self) -> str:
        return str(self.spec.metadata.get("entry_name", self.spec.name))

    def _entry_kind(self) -> str | None:
        value = self.spec.metadata.get("entry_kind")
        return str(value) if value is not None else None

    def _entry_path(self) -> str | None:
        value = self.spec.metadata.get("entry_path")
        return str(value) if value is not None else None

    def _source_type(self) -> str | None:
        value = self.spec.metadata.get("source_type")
        return str(value) if value is not None else None

    def _workflow_family(self) -> str | None:
        value = self.spec.metadata.get("workflow_family")
        return str(value) if value is not None else None

    def _write_human_report(self, manifest: dict, *, destination_root: Path) -> Path:
        destination_root.mkdir(parents=True, exist_ok=True)
        report_path = destination_root / "run_report.md"
        archived_logs_dir = Path(str(manifest["output_dir"])) / "logs"

        lines = [
            "# driverFoam Run Report",
            "",
            "## Summary",
            "",
            f"- Entry: `{manifest['entry']}`",
            f"- Entry kind: `{manifest.get('entry_kind')}`",
            f"- Entry path: `{manifest.get('entry_path')}`",
            f"- Requested action: `{manifest['requested_action']}`",
            f"- Run ID: `{manifest['run_id']}`",
            f"- Status: `{manifest['status']}`",
            f"- Post-process status: `{manifest['postprocess_status']}`",
            f"- Dry run: `{manifest['dry_run']}`",
            f"- Continue on error: `{manifest['continue_on_error']}`",
            f"- Case root: `{manifest['case_root']}`",
            f"- Setup root: `{manifest['setup_root']}`",
            f"- Output dir: `{manifest['output_dir']}`",
            f"- Started at (UTC): `{manifest['started_at_utc']}`",
            f"- Updated at (UTC): `{manifest['updated_at_utc']}`",
            f"- Finished at (UTC): `{manifest['finished_at_utc']}`",
            f"- Current case: `{manifest['current_case_id']}`",
            f"- Total cases: `{manifest['total_cases']}`",
            f"- Planned cases: `{manifest['planned_cases']}`",
            f"- Completed cases: `{manifest['completed_cases']}`",
            f"- Failed cases: `{manifest['failed_cases']}`",
        ]

        if manifest.get("error"):
            lines.extend(["", "## Run Error", "", "```text", str(manifest["error"]), "```"])

        results = manifest.get("results", [])
        lines.extend(["", "## Case Results", ""])
        if not results:
            lines.append("No case results recorded.")
        else:
            for result in results:
                lines.extend(
                    [
                        f"### {result['case_id']}",
                        f"- Status: `{result['status']}`",
                        f"- Index: `{result.get('index')}` / `{result.get('total_cases')}`",
                        f"- Duration (s): `{result['duration_s']:.6f}`",
                        f"- Started at (UTC): `{result.get('started_at_utc')}`",
                        f"- Finished at (UTC): `{result.get('finished_at_utc')}`",
                        f"- Parameters: `{json.dumps(result.get('params', {}), sort_keys=True)}`",
                    ]
                )
                if result.get("error"):
                    lines.extend(["- Error:", "", "```text", str(result["error"]), "```"])
                lines.append("")

        plots_manifest_path = manifest.get("plots_manifest_path")
        artifact_lines: list[str] = []
        if archived_logs_dir.exists():
            artifact_lines.append(f"- Archived logs: `{archived_logs_dir}`")
        if plots_manifest_path:
            artifact_lines.append(f"- Plots manifest: `{plots_manifest_path}`")
        if artifact_lines:
            lines.extend(["## Artifacts", "", *artifact_lines, ""])

        report_path.write_text("\n".join(lines).rstrip() + "\n")
        return report_path

    def _write_manifest(
        self,
        results: list[CaseResult],
        *,
        status: str,
        current_case_id: str | None,
        total_cases: int | None = None,
        postprocess_status: str = "not_started",
        error: str | None = None,
        finished_at_utc: str | None = None,
    ) -> None:
        destination_root = self._manifest_destination_root()
        plots_manifest = self.spec.output_dir / "plots.json"
        total = total_cases if total_cases is not None else len(results)
        manifest = {
            "schema_version": "2.2",
            "run_id": self.run_id,
            "requested_action": self.requested_action,
            "entry": self._entry_name(),
            "entry_kind": self._entry_kind(),
            "entry_path": self._entry_path(),
            "source_type": self._source_type(),
            "workflow_family": self._workflow_family(),
            "case_root": str(self.spec.case_root),
            "setup_root": str(self.spec.setup_root),
            "output_dir": str(self.spec.output_dir),
            "dry_run": self.dry_run,
            "continue_on_error": self.continue_on_error,
            "status": status,
            "postprocess_status": postprocess_status,
            "current_case_id": current_case_id,
            "started_at_utc": self.started_at_utc,
            "updated_at_utc": _utc_now(),
            "finished_at_utc": finished_at_utc,
            "total_cases": total,
            "planned_cases": sum(1 for item in results if item.status == "planned"),
            "completed_cases": sum(1 for item in results if item.status == "ok"),
            "failed_cases": sum(1 for item in results if item.status == "failed"),
            "error": error,
            "plots_manifest_path": str(plots_manifest) if plots_manifest.exists() else None,
            "results": [asdict(item) for item in results],
        }

        report_path = self._write_human_report(manifest, destination_root=destination_root)
        manifest["human_report_path"] = str(report_path)

        destination_root.mkdir(parents=True, exist_ok=True)
        manifest_path = destination_root / "run_manifest.json"
        tmp_path = manifest_path.with_name(manifest_path.name + ".tmp")
        tmp_path.write_text(json.dumps(manifest, indent=2))
        os.replace(tmp_path, manifest_path)
        print(f"Run manifest written to: {manifest_path}")
