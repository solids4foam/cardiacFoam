"""Agent-reproduces-regression check (gated on a built/sourced cardiacFoam).

The equivalence bar (per user, 2026-07-04): drive each regression case through
the **agent's own run path** and require the agent-produced outputs to match the
committed ``.reference`` within the case's own tolerances. The committed
reference is the ground truth — the hand-authored path is not re-run.

- Agent run (strict): ``foamctl run --strict --entry <name> --tutorials-root
  <staged>`` — the agent resolves the registered spec, plans it (non-mutating),
  and executes the case's workflow (solver + post). No dictionary overrides are
  applied; dict mutation lives only in the sweep path.
- Agent run (generic): the same via the case-folder path with
  ``--entry-kind case_folder`` for cases with no registered spec.

Then the agent's ``postProcessing`` outputs are checked against the committed
reference points.

Everything that touches a solver is gated by :func:`solver_available` and
returns a ``skipped`` result when cardiacFoam is not built/sourced, so this
module imports and its pure helpers unit-test anywhere.
"""
from __future__ import annotations

import os
import shutil
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path

import openfoam_driver
from openfoam_driver.specs.common import tutorials_root_default
from openfoam_driver.tests.regression_equivalence.registry import RegressionCase


# --------------------------------------------------------------------------- #
# Pure helpers (solver-free, unit-tested)
# --------------------------------------------------------------------------- #

@dataclass(frozen=True)
class ReferencePoint:
    data_file: str
    time: float
    variable: str
    expected: float
    tolerance: float


def parse_columnar_reference(text: str) -> list[ReferencePoint]:
    """Parse a `file time variable expected tolerance` reference file.

    Returns [] for reference files that don't follow this columnar layout
    (e.g. the bidomain `kind key metric ...` metric style), signalling the
    caller to fall back to the case's own regressionTest.sh as the gate.
    """
    points: list[ReferencePoint] = []
    for raw in text.splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        cols = line.split()
        if len(cols) < 5:
            return []
        data_file, time, variable, expected, tolerance = cols[:5]
        try:
            points.append(
                ReferencePoint(
                    data_file, float(time), variable,
                    float(expected), float(tolerance),
                )
            )
        except ValueError:
            # Non-numeric where numbers are expected -> not this layout.
            return []
    return points


def read_series_value(
    text: str, target_time: float, variable: str, *, time_atol: float = 1e-9
) -> float | None:
    """Return `variable` at the row whose first column is nearest `target_time`.

    Mirrors the awk extractor in the cases' regressionTest.sh: the first line is
    a header of column names; data rows follow with time in column 1. A match
    requires the nearest time to be within `time_atol` of the target.
    """
    lines = [ln for ln in text.splitlines() if ln.strip()]
    if not lines:
        return None
    header = lines[0].split()
    try:
        col = header.index(variable)
    except ValueError:
        return None
    best_diff = float("inf")
    best_val: float | None = None
    for row in lines[1:]:
        cells = row.split()
        if len(cells) <= col:
            continue
        try:
            t = float(cells[0])
            v = float(cells[col])
        except ValueError:
            continue
        diff = abs(t - target_time)
        if diff < best_diff:
            best_diff = diff
            best_val = v
    if best_val is None or best_diff > time_atol:
        return None
    return best_val


def values_agree(a: float, b: float, tolerance: float) -> bool:
    return abs(a - b) <= tolerance


@dataclass(frozen=True)
class ManufacturedReferencePoint:
    kind: str
    key: str
    metric: str
    expected: float
    tolerance: float


def parse_manufactured_reference(text: str) -> list[ManufacturedReferencePoint]:
    points = []
    for raw in text.splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        cols = line.split()
        if len(cols) < 5:
            return []
        kind, key, metric, expected, tolerance = cols[:5]
        if kind not in ("summary", "error", "pseudoECG"):
            return []
        try:
            points.append(
                ManufacturedReferencePoint(
                    kind, key, metric, float(expected), float(tolerance)
                )
            )
        except ValueError:
            return []
    return points


def find_manufactured_error_file(case_path: Path) -> Path | None:
    for f in case_path.glob("postProcessing/*.dat"):
        text = f.read_text(errors="ignore")
        if "manufactured-solution error summary" in text.lower() or "manufactured activation-time summary" in text.lower():
            return f
    for f in case_path.glob("processor*/postProcessing/*.dat"):
        text = f.read_text(errors="ignore")
        if "manufactured-solution error summary" in text.lower() or "manufactured activation-time summary" in text.lower():
            return f
    return None


def extract_summary_value(text: str, key: str) -> float | None:
    for line in text.splitlines():
        if key == "cells" and "Number of cells" in line:
            return float(line.split("=")[1].strip())
        if key == "finalTime" and "Final simulation time" in line:
            return float(line.split("=")[1].strip())
    return None


def extract_error_metric(text: str, key: str, metric: str) -> float | None:
    col = {"L1": 1, "L2": 2, "Linf": 3}.get(metric)
    if col is None: return None
    for line in text.splitlines():
        parts = line.split()
        if not parts: continue
        if parts[0] == key and len(parts) > col:
            try:
                return float(parts[col])
            except ValueError:
                pass
    return None


def find_pseudo_ecg_file(case_path: Path) -> Path | None:
    for p in [case_path / "postProcessing" / "eikonalECG.dat"] + list(case_path.glob("processor*/postProcessing/eikonalECG.dat")):
        if p.exists(): return p
    return None


def extract_pseudo_ecg_value(text: str, key: str) -> float | None:
    lines = text.splitlines()
    if not lines: return None
    header = lines[0].split()
    col = -1
    for i, h in enumerate(header):
        if h == key or h == f"numeric_{key}":
            col = i - 1 if header[0] == "#" else i
            break
    if col < 0: return None
    
    for line in reversed(lines):
        if line.startswith("#"): continue
        parts = line.split()
        if len(parts) > col:
            try:
                return float(parts[col])
            except ValueError:
                pass
    return None


def _check_manufactured_reference(case_path: Path, points: list[ManufacturedReferencePoint]) -> tuple[bool, str]:
    problems = []
    checks = 0
    error_file = find_manufactured_error_file(case_path)
    error_text = error_file.read_text(errors="ignore") if error_file else ""
    
    ecg_file = find_pseudo_ecg_file(case_path)
    ecg_text = ecg_file.read_text(errors="ignore") if ecg_file else ""
    
    for p in points:
        checks += 1
        val = None
        if p.kind == "summary":
            if error_text: val = extract_summary_value(error_text, p.key)
        elif p.kind == "error":
            if error_text: val = extract_error_metric(error_text, p.key, p.metric)
        elif p.kind == "pseudoECG":
            if ecg_text: val = extract_pseudo_ecg_value(ecg_text, p.key)
            
        if val is None:
            problems.append(f"{p.kind} {p.key} {p.metric}: missing in agent output")
        elif not values_agree(val, p.expected, p.tolerance):
            problems.append(f"{p.kind} {p.key} {p.metric}: agent={val} expected={p.expected} tol={p.tolerance}")
            
    if problems:
        return False, "\n".join(problems)
    return True, f"{checks} reference points reproduced within tolerance"


# --------------------------------------------------------------------------- #
# Gated agent-driven orchestration
# --------------------------------------------------------------------------- #

# Output write times carry a small offset from the requested grid (e.g.
# 1.5000010 rather than 1.5000000), so we sample the row nearest each reference
# time; this tolerance exceeds that offset while staying below the write
# interval so the nearest row is unambiguous.
TIME_MATCH_ATOL = 1e-2


@dataclass(frozen=True)
class ReproResult:
    case_dir: str
    driver: str
    # "reproduced": agent output matches committed reference within tolerance.
    # "mismatch":   agent ran but some reference point is off/missing.
    # "run_failed": the agent run itself returned non-zero.
    # "unsupported_ref": reference format not parsed by this harness.
    # "skipped":    no solver, or case not agent-addressable.
    status: str
    detail: str


def solver_available() -> bool:
    return bool(shutil.which("cardiacFoam")) and bool(os.environ.get("WM_PROJECT_DIR"))


def _stage_tutorials_root(case: RegressionCase) -> tuple[Path, Path]:
    """Copy the case into a throwaway tutorials root; return (root, case_path)."""
    src = tutorials_root_default() / case.case_dir
    root = Path(tempfile.mkdtemp(prefix="regressioneq_")) / "tutorials"
    case_path = root / case.case_dir
    case_path.parent.mkdir(parents=True, exist_ok=True)
    shutil.copytree(src, case_path)
    # Remove any committed outputs so we prove the agent produced fresh ones.
    for stale in ("postProcessing", "workflow_state.json", "workflow_logs"):
        p = case_path / stale
        shutil.rmtree(p, ignore_errors=True) if p.is_dir() else p.unlink(missing_ok=True)
    for log_file in case_path.glob("log.*"):
        log_file.unlink(missing_ok=True)
    return root, case_path


def _drive_agent(case: RegressionCase, driver: str, tutorials_root: Path) -> subprocess.CompletedProcess:
    """Invoke `foamctl run --strict` through the agent CLI on the staged case."""
    pkg_parent = str(Path(openfoam_driver.__file__).resolve().parent.parent)
    env = os.environ.copy()
    env["PYTHONPATH"] = pkg_parent + os.pathsep + env.get("PYTHONPATH", "")
    if driver == "strict":
        entry_args = ["--entry", case.entry_name]
    else:
        entry_args = ["--entry", case.case_dir, "--entry-kind", "case_folder"]
    argv = [
        sys.executable, "-m", "openfoam_driver", "run", "--strict",
        *entry_args, "--tutorials-root", str(tutorials_root),
    ]
    return subprocess.run(argv, env=env, capture_output=True, text=True)


def check_reference(case_path: Path, reference_text: str) -> tuple[bool, str]:
    """Check agent outputs under `case_path` against a reference file.

    Returns (ok, detail). Supports columnar layout (`file time variable
    expected tolerance`) or manufactured layout (`kind key metric expected
    tolerance`). Non-supported references raise NotImplementedError.
    """
    points = parse_columnar_reference(reference_text)
    if points:
        problems: list[str] = []
        checks = 0
        for p in points:
            checks += 1
            f = case_path / p.data_file
            val = (
                read_series_value(f.read_text(), p.time, p.variable, time_atol=TIME_MATCH_ATOL)
                if f.exists() else None
            )
            if val is None:
                problems.append(f"{p.data_file} {p.variable}@{p.time}: missing in agent output")
            elif not values_agree(val, p.expected, p.tolerance):
                problems.append(
                    f"{p.data_file} {p.variable}@{p.time}: agent={val} "
                    f"expected={p.expected} tol={p.tolerance}"
                )
        if problems:
            return False, "\n".join(problems)
        return True, f"{checks} reference points reproduced within tolerance"

    man_points = parse_manufactured_reference(reference_text)
    if man_points:
        return _check_manufactured_reference(case_path, man_points)

    raise NotImplementedError("reference format is not supported")


def verify_reproduction(case: RegressionCase, *, driver: str) -> ReproResult:
    """Drive `case` through the agent and check outputs vs committed reference."""
    if not solver_available():
        return ReproResult(
            case.case_dir, driver, "skipped",
            "cardiacFoam not built/sourced (WM_PROJECT_DIR unset or binary absent)",
        )
    if driver == "generic" and not case.generic_addressable:
        return ReproResult(
            case.case_dir, driver, "skipped",
            "case layout not addressable by agent discovery",
        )

    reference_text = (
        tutorials_root_default() / case.case_dir / case.reference_file
    ).read_text()
    root, case_path = _stage_tutorials_root(case)
    try:
        proc = _drive_agent(case, driver, root)
        if proc.returncode != 0:
            return ReproResult(
                case.case_dir, driver, "run_failed",
                f"agent run rc={proc.returncode}; stderr tail:\n{proc.stderr[-1500:]}",
            )
        try:
            ok, detail = check_reference(case_path, reference_text)
        except NotImplementedError as exc:
            return ReproResult(
                case.case_dir, driver, "unsupported_ref",
                f"agent run ok, but {exc} (reference {case.reference_file})",
            )
        return ReproResult(
            case.case_dir, driver, "reproduced" if ok else "mismatch", detail
        )
    finally:
        shutil.rmtree(root.parent, ignore_errors=True)
