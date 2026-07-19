#!/usr/bin/env python3
"""Summarise monodomain corrector-study runs and compare final Vm fields."""

from __future__ import annotations

import csv
import hashlib
import math
from pathlib import Path
import re
import statistics
import sys


def read_metadata(path: Path) -> dict[str, str]:
    result: dict[str, str] = {}
    for line in path.read_text().splitlines():
        if "=" in line:
            key, value = line.split("=", 1)
            result[key] = value
    return result


def parse_vm(path: Path) -> list[float]:
    text = path.read_text()
    match = re.search(
        r"internalField\s+nonuniform\s+List<scalar>\s+(\d+)\s*\((.*?)\)\s*;",
        text,
        flags=re.S,
    )
    if match:
        expected = int(match.group(1))
        values = [float(token) for token in match.group(2).split()]
        if len(values) != expected:
            raise ValueError(f"{path}: expected {expected} values, found {len(values)}")
        return values

    uniform = re.search(r"internalField\s+uniform\s+([^;]+);", text)
    if uniform:
        return [float(uniform.group(1))]
    raise ValueError(f"Could not parse internalField from {path}")


def parse_mms(path: Path) -> dict[str, float | int]:
    text = path.read_text()
    vm = re.search(
        r"^Vm\s+([0-9.eE+-]+)\s+([0-9.eE+-]+)\s+([0-9.eE+-]+)",
        text,
        flags=re.M,
    )
    dx = re.search(r"Grid spacing \(dx\)\s+=\s+([0-9.eE+-]+)", text)
    steps = re.search(r"Number of steps\s+=\s+(\d+)", text)
    if not (vm and dx and steps):
        raise ValueError(f"Could not parse MMS summary {path}")
    return {
        "vm_l1": float(vm.group(1)),
        "vm_l2": float(vm.group(2)),
        "vm_linf": float(vm.group(3)),
        "dx": float(dx.group(1)),
        "n_steps": int(steps.group(1)),
    }


def parse_log(path: Path) -> dict[str, float | int]:
    text = path.read_text(errors="replace")
    timings = re.findall(
        r"ExecutionTime\s*=\s*([0-9.eE+-]+)\s+s\s+ClockTime\s*=\s*([0-9.eE+-]+)",
        text,
    )
    execution, clock = (float(value) for value in timings[-1]) if timings else (math.nan, math.nan)
    return {
        "vm_solves": len(re.findall(r"Solving for Vm(?:Final)?\b", text)),
        "execution_time_s": execution,
        "clock_time_s": clock,
    }


def field_delta(values: list[float], reference: list[float]) -> tuple[float, float]:
    if len(values) != len(reference):
        raise ValueError(f"Field sizes differ: {len(values)} != {len(reference)}")
    differences = [a - b for a, b in zip(values, reference)]
    l2 = math.sqrt(sum(value * value for value in differences) / len(differences))
    linf = max(abs(value) for value in differences)
    return l2, linf


def main() -> None:
    if len(sys.argv) != 2:
        raise SystemExit("usage: summarize_corrector_study.py RESULTS_DIR")
    results_dir = Path(sys.argv[1]).resolve()

    rows: list[dict[str, object]] = []
    fields: dict[str, list[float]] = {}
    for metadata_path in sorted(results_dir.glob("*/metadata.env")):
        run_dir = metadata_path.parent
        metadata = read_metadata(metadata_path)
        mms = parse_mms(run_dir / "summary.dat")
        log = parse_log(run_dir / "log.cardiacFoam")
        vm_path = run_dir / "Vm"
        fields[metadata["run_id"]] = parse_vm(vm_path)
        rows.append(
            {
                **metadata,
                **mms,
                **log,
                "vm_sha256": hashlib.sha256(vm_path.read_bytes()).hexdigest(),
            }
        )

    if not rows:
        raise SystemExit(f"No runs found below {results_dir}")

    references: dict[tuple[int, str], dict[str, object]] = {}
    resolutions = sorted({int(row["resolution"]) for row in rows})
    for resolution in resolutions:
        outer_candidates = [
            row
            for row in rows
            if int(row["resolution"]) == resolution
            and row["mode"] == "outer"
            and int(row["repeat"]) == 1
        ]
        references[(resolution, "outer")] = max(
            outer_candidates, key=lambda row: int(row["n_outer"])
        )
        nonorth_candidates = [
            row
            for row in rows
            if int(row["resolution"]) == resolution
            and row["mode"] == "nonorth"
            and int(row["repeat"]) == 1
        ]
        references[(resolution, "nonorth")] = min(
            nonorth_candidates, key=lambda row: int(row["n_nonorth"])
        )

    for row in rows:
        key = (int(row["resolution"]), str(row["mode"]))
        reference = references[key]
        delta_l2, delta_linf = field_delta(
            fields[str(row["run_id"])], fields[str(reference["run_id"])]
        )
        row["reference_run"] = reference["run_id"]
        row["delta_l2_to_reference"] = delta_l2
        row["delta_linf_to_reference"] = delta_linf
        row["delta_l2_over_mms_l2"] = delta_l2 / float(reference["vm_l2"])

    fieldnames = [
        "run_id",
        "mode",
        "resolution",
        "n_outer",
        "n_nonorth",
        "repeat",
        "delta_t",
        "end_time",
        "latest_time",
        "dx",
        "n_steps",
        "vm_solves",
        "vm_l1",
        "vm_l2",
        "vm_linf",
        "execution_time_s",
        "clock_time_s",
        "harness_wall_time_s",
        "vm_sha256",
        "reference_run",
        "delta_l2_to_reference",
        "delta_linf_to_reference",
        "delta_l2_over_mms_l2",
    ]
    with (results_dir / "raw_results.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    groups: dict[tuple[int, str, int, int], list[dict[str, object]]] = {}
    for row in rows:
        key = (
            int(row["resolution"]),
            str(row["mode"]),
            int(row["n_outer"]),
            int(row["n_nonorth"]),
        )
        groups.setdefault(key, []).append(row)

    lines = [
        "# Corrector-study summary",
        "",
        "| N | mode | outer | non-orth | repeats | Vm solves/step | Vm L2 MMS | delta L2/reference | delta/MMS | median clock (s) | field hashes identical across repeats |",
        "|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|:---:|",
    ]
    for key in sorted(groups):
        group = groups[key]
        first = min(group, key=lambda row: int(row["repeat"]))
        clocks = [float(row["harness_wall_time_s"]) for row in group]
        hashes_identical = len({str(row["vm_sha256"]) for row in group}) == 1
        solves_per_step = float(first["vm_solves"]) / int(first["n_steps"])
        lines.append(
            "| {n} | {mode} | {outer} | {nonorth} | {repeats} | {solves:.2f} | "
            "{mms:.6e} | {delta:.6e} | {ratio:.3e} | {clock:.3f} | {same} |".format(
                n=key[0],
                mode=key[1],
                outer=key[2],
                nonorth=key[3],
                repeats=len(group),
                solves=solves_per_step,
                mms=float(first["vm_l2"]),
                delta=float(first["delta_l2_to_reference"]),
                ratio=float(first["delta_l2_over_mms_l2"]),
                clock=statistics.median(clocks),
                same="yes" if hashes_identical else "no",
            )
        )

    lines.extend(
        [
            "",
            "The `outer` reference is the highest outer-corrector count at each N. ",
            "The `nonorth` reference is `nNonOrthogonalCorrectors=0` at each N. ",
            "Use `delta/MMS <= 1e-2` as the predeclared corrector-convergence threshold.",
            "",
            "## Automated checks",
            "",
        ]
    )

    for resolution in resolutions:
        outer_groups = []
        for key, group in groups.items():
            if key[0] == resolution and key[1] == "outer":
                first = min(group, key=lambda row: int(row["repeat"]))
                outer_groups.append(first)
        outer_groups.sort(key=lambda row: int(row["n_outer"]))

        passing = None
        for index, row in enumerate(outer_groups[:-1]):
            next_row = outer_groups[index + 1]
            correction_ratio = float(row["delta_l2_over_mms_l2"])
            next_error_change = abs(float(next_row["vm_l2"]) - float(row["vm_l2"])) / float(
                next_row["vm_l2"]
            )
            if correction_ratio <= 0.01 and next_error_change <= 0.01:
                passing = int(row["n_outer"])
                break

        if passing is None:
            lines.append(
                f"- N={resolution}: no non-reference outer count passes both 1% criteria."
            )
        else:
            lines.append(
                f"- N={resolution}: smallest outer count passing both 1% criteria is {passing}."
            )

        nonorth_rows = [
            row
            for row in rows
            if int(row["resolution"]) == resolution
            and row["mode"] == "nonorth"
            and int(row["repeat"]) == 1
        ]
        nonorth_reference = references[(resolution, "nonorth")]
        nonorth_active = any(
            int(row["n_nonorth"]) > 0
            and
            (
                float(row["delta_linf_to_reference"]) > 0.0
                or row["vm_sha256"] != nonorth_reference["vm_sha256"]
            )
            for row in nonorth_rows
        )
        solve_counts_match = all(
            int(row["vm_solves"])
            == int(row["n_steps"])
            * int(row["n_outer"])
            * (int(row["n_nonorth"]) + 1)
            for row in nonorth_rows
        )
        lines.append(
            f"- N={resolution}: nNonOrthogonalCorrectors active check: "
            f"{'PASS' if nonorth_active and solve_counts_match else 'FAIL'} "
            f"(field changed={nonorth_active}, solve counts match={solve_counts_match})."
        )

    lines.append("")
    (results_dir / "summary.md").write_text("\n".join(lines))
    print("\n".join(lines))


if __name__ == "__main__":
    main()
