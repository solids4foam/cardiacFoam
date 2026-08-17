#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import json
import re
from pathlib import Path


def parse_probe_file(path: Path) -> tuple[list[tuple[float, float, float]], list[float], list[list[float]]]:
    if not path.is_file():
        raise FileNotFoundError(f"Probe file not found: {path}")

    positions: list[tuple[float, float, float]] = []
    numeric_rows: list[list[float]] = []
    probe_pattern = re.compile(r"#\s*Probe\s+\d+\s+\(([^)]+)\)")

    with path.open("r", encoding="ascii") as handle:
        for raw_line in handle:
            line = raw_line.strip()

            match = probe_pattern.match(line)
            if match:
                coords = tuple(float(value) for value in match.group(1).split())
                if len(coords) != 3:
                    raise ValueError(f"Unexpected probe coordinates: {line}")
                positions.append(coords)
                continue

            if not line or line.startswith("#"):
                continue

            numeric_rows.append([float(token) for token in line.split()])

    if not numeric_rows:
        raise ValueError(f"No numeric probe data found in {path}")

    times = [row[0] for row in numeric_rows]
    vm_data = []
    for i in range(1, len(numeric_rows[0])):
        vm_data.append([row[i] for row in numeric_rows])

    if positions and len(positions) != len(vm_data):
        raise ValueError(
            "Probe position count does not match sampled value count: "
            f"{len(positions)} vs {len(vm_data)}"
        )

    return positions, times, vm_data


def detect_s2_crossings(times: list[float], vm_data: list[list[float]]) -> list[float]:
    s2_crossings = []
    for p in range(len(vm_data)):
        crossings = []
        for t in range(1, len(times)):
            if vm_data[p][t-1] < 0 and vm_data[p][t] >= 0:
                crossings.append(times[t])
        
        # Debounce crossings
        clean = [crossings[0]] if crossings else []
        for c in crossings[1:]:
            if c - clean[-1] > 0.05:
                clean.append(c)
                
        if len(clean) < 2:
            raise ValueError(f"Probe {p} did not activate a second time for S2 beat! Crossings detected: {clean}")
            
        s2_crossings.append(clean[-1])
        
    return s2_crossings


def compute_segment_cv(
    positions: list[tuple[float, float, float]],
    activation_times: list[float],
) -> list[tuple[int, int, float, float, float]]:
    segments: list[tuple[int, int, float, float, float]] = []

    for index in range(len(activation_times) - 1):
        x0 = positions[index][0]
        x1 = positions[index + 1][0]
        t0 = activation_times[index]
        t1 = activation_times[index + 1]
        dx = x1 - x0
        dt = t1 - t0

        if t0 < 0 or t1 < 0:
            raise ValueError("At least one probe never activated.")
        if dt <= 0:
            raise ValueError(
                f"Non-positive activation delay between probes {index} and {index + 1}: {dt}"
            )

        segments.append((index, index + 1, dx, dt, dx / dt))

    return segments


def build_summary(probe_path: Path) -> dict:
    positions, times, vm_data = parse_probe_file(probe_path)
    s2_activation_times = detect_s2_crossings(times, vm_data)
    segments = compute_segment_cv(positions, s2_activation_times)

    central_dx = positions[3][0] - positions[1][0]
    central_dt = s2_activation_times[3] - s2_activation_times[1]
    if central_dt <= 0:
        raise ValueError("Central-cable activation delay is non-positive.")
    central_cv = central_dx / central_dt

    return {
        "probe_positions_m": [
            {"probe_index": index, "x": xyz[0], "y": xyz[1], "z": xyz[2]}
            for index, xyz in enumerate(positions)
        ],
        "activation_times_s": [
            {"probe_index": index, "activation_time_s": activation_time}
            for index, activation_time in enumerate(s2_activation_times)
        ],
        "segments": [
            {
                "start_probe": start,
                "end_probe": end,
                "dx_m": dx,
                "dt_s": dt,
                "cv_m_per_s": cv,
            }
            for start, end, dx, dt, cv in segments
        ],
        "central_cv": {
            "start_probe": 1,
            "end_probe": 3,
            "dx_m": central_dx,
            "dt_s": central_dt,
            "cv_m_per_s": central_cv,
        },
    }


def write_summary_files(
    *,
    case_dir: Path,
    output_dir: Path | None = None,
    case_id: str | None = None,
) -> dict:
    probe_files = list((case_dir / "postProcessing" / "cableProbes").glob("*/Vm"))
    if not probe_files:
        raise FileNotFoundError(f"No Vm probe files found in {case_dir}/postProcessing/cableProbes/")
        
    # Sort by numeric time directory, take the latest
    probe_files.sort(key=lambda p: float(p.parent.name))
    target_probe = probe_files[-1]
    
    summary = build_summary(target_probe)
    
    # Write summary payload
    output_dir.mkdir(parents=True, exist_ok=True)
    json_path = output_dir / f"{case_id}_cv_summary.json"
    activation_path = output_dir / f"{case_id}_activation_times.csv"
    segment_path = output_dir / f"{case_id}_segments.csv"

    json_payload = {"case_id": case_id, **summary}
    json_path.write_text(json.dumps(json_payload, indent=2), encoding="ascii")

    with activation_path.open("w", newline="", encoding="ascii") as handle:
        writer = csv.writer(handle)
        writer.writerow(["probe_index", "x_mm", "s2_activation_time_ms"])
        for position, activation in zip(summary["probe_positions_m"], summary["activation_times_s"]):
            writer.writerow(
                [
                    position["probe_index"],
                    f"{1e3 * position['x']:.6f}",
                    f"{1e3 * activation['activation_time_s']:.6f}",
                ]
            )

    with segment_path.open("w", newline="", encoding="ascii") as handle:
        writer = csv.writer(handle)
        writer.writerow(["start_probe", "end_probe", "dx_mm", "dt_ms", "s2_cv_m_per_s"])
        for segment in summary["segments"]:
            writer.writerow(
                [
                    segment["start_probe"],
                    segment["end_probe"],
                    f"{1e3 * segment['dx_m']:.6f}",
                    f"{1e3 * segment['dt_s']:.6f}",
                    f"{segment['cv_m_per_s']:.8f}",
                ]
            )
    return summary

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--case-dir", default=".")
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--case-id", default=None)
    args = parser.parse_args()

    case_dir = Path(args.case_dir).resolve()
    case_id = args.case_id
    if case_id is None:
        sentinel = case_dir / ".driverfoam_case_id"
        if not sentinel.exists():
            raise SystemExit(
                f"--case-id not supplied and sentinel {sentinel} not found."
            )
        case_id = sentinel.read_text().strip()

    write_summary_files(
        case_dir=case_dir,
        output_dir=Path(args.output_dir).resolve(),
        case_id=case_id,
    )
