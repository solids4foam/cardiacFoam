"""Machine-readable verification experiment contracts and legacy execution."""

from __future__ import annotations

import json
import re
import subprocess
from pathlib import Path


CONTRACT_PATH = Path(__file__).resolve().parent / "verification_experiments.json"
REPO_ROOT = Path(__file__).resolve().parents[3]
_ID = re.compile(r"^[a-z][a-z0-9_]*$")


def load_contracts(path: Path = CONTRACT_PATH) -> list[dict]:
    payload = json.loads(path.read_text())
    if payload.get("schema_version") != "1":
        raise ValueError("verification contract schema_version must be '1'")
    experiments = payload.get("experiments")
    if not isinstance(experiments, list) or not experiments:
        raise ValueError("verification contracts require a non-empty experiments list")
    seen: set[str] = set()
    for experiment in experiments:
        experiment_id = experiment.get("experiment_id", "")
        if not _ID.fullmatch(experiment_id):
            raise ValueError(f"invalid experiment_id: {experiment_id!r}")
        if experiment_id in seen:
            raise ValueError(f"duplicate experiment_id: {experiment_id}")
        seen.add(experiment_id)
        for field in ("description", "case_dir", "execution", "matrix", "observables", "aggregation"):
            if field not in experiment:
                raise ValueError(f"{experiment_id}: missing {field}")
        execution = experiment["execution"]
        if execution.get("kind") not in {"driver_sweep", "legacy_bash"}:
            raise ValueError(f"{experiment_id}: invalid execution kind")
        if execution["kind"] == "driver_sweep":
            specs = execution.get("driver_specs")
            if not isinstance(specs, list) or not specs:
                raise ValueError(f"{experiment_id}: driver_sweep requires driver_specs")
            # runner is optional here: reproduce_verification.sh runs
            # driver_sweep experiments straight from driver_specs +
            # aggregation.key, no wrapper script required.
        elif not execution.get("runner"):
            raise ValueError(f"{experiment_id}: legacy_bash requires execution runner")
        aggregation = experiment["aggregation"]
        expected_result = f"setup/results/{experiment_id}.csv"
        if aggregation.get("result") != expected_result:
            raise ValueError(
                f"{experiment_id}: result must be {expected_result}"
            )
    return experiments


def plan(experiment_id: str | None = None) -> dict:
    experiments = load_contracts()
    if experiment_id is not None:
        experiments = [item for item in experiments if item["experiment_id"] == experiment_id]
        if not experiments:
            raise KeyError(f"unknown verification experiment: {experiment_id}")
    planned = []
    for item in experiments:
        case_root = REPO_ROOT / item["case_dir"]
        execution = item["execution"]
        aggregation = item["aggregation"]
        runner = execution.get("runner")
        checks = {
            "case_dir": case_root.is_dir(),
            "runner": True if not runner else (case_root / runner).is_file(),
            "result": (case_root / aggregation["result"]).is_file(),
            "reference": (
                True if aggregation.get("reference") is None
                else (case_root / aggregation["reference"]).is_file()
            ),
        }
        if execution["kind"] == "driver_sweep":
            checks["driver_specs"] = all(
                (case_root / spec).is_file() for spec in execution["driver_specs"]
            )
        planned.append({**item, "checks": checks, "ready": all(checks.values())})
    return {"schema_version": "1", "experiments": planned}


def tsv_rows(experiment_id: str | None = None) -> str:
    rows = []
    experiments = load_contracts()
    if experiment_id is not None:
        experiments = [
            item for item in experiments if item["experiment_id"] == experiment_id
        ]
        if not experiments:
            raise KeyError(f"unknown verification experiment: {experiment_id}")
    for item in experiments:
        execution = item["execution"]
        aggregation = item["aggregation"]
        driver_specs = execution.get("driver_specs") or []
        rows.append("\t".join((
            item["experiment_id"], item["case_dir"], execution["kind"],
            execution.get("runner") or "-",
            driver_specs[0] if driver_specs else "-",
            aggregation.get("key") or "-", aggregation["result"],
            aggregation.get("reference") or "-",
        )))
    return "\n".join(rows)


def run(experiment_id: str) -> int:
    plan(experiment_id)  # validate before executing anything
    command = [str(REPO_ROOT / "reproduce_verification.sh"), experiment_id]
    return subprocess.run(command, cwd=REPO_ROOT).returncode
