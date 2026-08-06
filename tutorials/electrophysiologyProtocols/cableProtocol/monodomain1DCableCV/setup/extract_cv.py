#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import json
import re
from pathlib import Path


CASE_DIR = Path(__file__).resolve().parents[1]


def probe_file(case_dir: Path) -> Path:
    return case_dir / "postProcessing" / "cableProbes" / "0" / "activationTime"


def summary_file(case_dir: Path) -> Path:
    return case_dir / "postProcessing" / "cv_summary.txt"


def parse_probe_file(path: Path) -> tuple[list[tuple[float, float, float]], list[float]]:
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

    last_row = numeric_rows[-1]
    if len(last_row) < 2:
        raise ValueError(f"Unexpected final probe row: {last_row}")

    sampled_values = last_row[1:]
    if positions and len(positions) != len(sampled_values):
        raise ValueError(
            "Probe position count does not match sampled value count: "
            f"{len(positions)} vs {len(sampled_values)}"
        )

    return positions, sampled_values


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
            raise ValueError("At least one probe never activated; increase endTime.")
        if dt <= 0:
            raise ValueError(
                f"Non-positive activation delay between probes {index} and {index + 1}: {dt}"
            )

        segments.append((index, index + 1, dx, dt, dx / dt))

    return segments


def format_mm(value_m: float) -> str:
    return f"{1e3 * value_m:.3f}"


def format_ms(value_s: float) -> str:
    return f"{1e3 * value_s:.3f}"


def build_summary(case_dir: Path) -> dict:
    positions, activation_times = parse_probe_file(probe_file(case_dir))
    segments = compute_segment_cv(positions, activation_times)

    central_dx = positions[3][0] - positions[1][0]
    central_dt = activation_times[3] - activation_times[1]
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
            for index, activation_time in enumerate(activation_times)
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


def render_summary_text(summary: dict) -> str:
    lines = [
        "monodomain1DCableCV summary",
        "",
        "Probe activation times:",
    ]

    positions = summary["probe_positions_m"]
    activation_times = summary["activation_times_s"]
    for position, activation in zip(positions, activation_times):
        lines.append(
            f"  probe {position['probe_index']}: x = {format_mm(position['x'])} mm, "
            f"activationTime = {format_ms(activation['activation_time_s'])} ms"
        )

    lines.extend(["", "Segment conduction velocities:"])
    for segment in summary["segments"]:
        lines.append(
            f"  probes {segment['start_probe']}->{segment['end_probe']}: "
            f"dx = {format_mm(segment['dx_m'])} mm, "
            f"dt = {format_ms(segment['dt_s'])} ms, "
            f"CV = {segment['cv_m_per_s']:.4f} m/s"
        )

    central = summary["central_cv"]
    lines.extend(
        [
            "",
            "Central calibration CV:",
            (
                f"  probes {central['start_probe']}->{central['end_probe']} "
                "(5 mm to 15 mm): "
                f"dx = {format_mm(central['dx_m'])} mm, "
                f"dt = {format_ms(central['dt_s'])} ms, "
                f"CV = {central['cv_m_per_s']:.4f} m/s"
            ),
        ]
    )
    return "\n".join(lines) + "\n"


def write_summary_files(
    *,
    case_dir: Path,
    output_dir: Path | None = None,
    case_id: str | None = None,
) -> dict:
    summary = build_summary(case_dir)
    text = render_summary_text(summary)

    case_summary_file = summary_file(case_dir)
    case_summary_file.parent.mkdir(parents=True, exist_ok=True)
    case_summary_file.write_text(text, encoding="ascii")

    if output_dir is not None and case_id is not None:
        output_dir.mkdir(parents=True, exist_ok=True)
        json_path = output_dir / f"{case_id}_cv_summary.json"
        activation_path = output_dir / f"{case_id}_activation_times.csv"
        segment_path = output_dir / f"{case_id}_segments.csv"

        json_payload = {"case_id": case_id, **summary}
        json_path.write_text(json.dumps(json_payload, indent=2), encoding="ascii")

        with activation_path.open("w", newline="", encoding="ascii") as handle:
            writer = csv.writer(handle)
            writer.writerow(["probe_index", "x_mm", "activation_time_ms"])
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
            writer.writerow(["start_probe", "end_probe", "dx_mm", "dt_ms", "cv_m_per_s"])
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

def main() -> None:
    parser = argparse.ArgumentParser(description="Compute 1D cable CV summary.")
    parser.add_argument("--case-dir", default=str(CASE_DIR))
    parser.add_argument("--output-dir")
    parser.add_argument("--case-id")
    args = parser.parse_args()

    case_dir = Path(args.case_dir).resolve()
    output_dir = Path(args.output_dir).resolve() if args.output_dir else None
    summary = write_summary_files(
        case_dir=case_dir,
        output_dir=output_dir,
        case_id=args.case_id,
    )
    print(render_summary_text(summary), end="")


if __name__ == "__main__":
    main()
