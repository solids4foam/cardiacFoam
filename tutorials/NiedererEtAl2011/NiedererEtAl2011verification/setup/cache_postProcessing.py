from __future__ import annotations

import csv
import math
import re
from pathlib import Path

_PROBE_LINE = re.compile(r"^# Probe (\d+) \(([^)]+)\)")


def _parse_probes_file(path: Path) -> tuple[list[tuple[str, float, float, float]], list[float]]:
    """Parse an OpenFOAM `probes` functionObject output file.

    Format: a `# Probe N (x y z)` header line per sample location, then a
    `# Time ...` header and one data row per sampled time. Only the file's
    last data row is used -- niederer_2012.py's DAG always samples with
    `postProcess -func <name> -latestTime`, so there is never more than one.
    """
    probes: list[tuple[str, float, float, float]] = []
    last_row: list[float] | None = None
    for line in path.read_text().splitlines():
        match = _PROBE_LINE.match(line)
        if match:
            index, coords = match.groups()
            x, y, z = (float(v) for v in coords.split())
            probes.append((f"P{index}", x, y, z))
            continue
        if not line.strip() or line.startswith("#"):
            continue
        last_row = [float(token) for token in line.split()]
    if last_row is None:
        raise ValueError(f"No data rows found in probes file: {path}")
    return probes, last_row[1:]  # drop the leading Time column


def _latest_probes_file(function_object_dir: Path) -> Path | None:
    """Find the most recent time directory's `activationTime` sample file."""
    if not function_object_dir.is_dir():
        return None
    time_dirs: list[tuple[float, Path]] = []
    for child in function_object_dir.iterdir():
        if not child.is_dir():
            continue
        try:
            time_dirs.append((float(child.name), child))
        except ValueError:
            continue
    if not time_dirs:
        return None
    _, latest = max(time_dirs, key=lambda item: item[0])
    candidate = latest / "activationTime"
    return candidate if candidate.is_file() else None


def _write_points_csv(
    probes: list[tuple[str, float, float, float]], values: list[float], destination: Path,
) -> None:
    with destination.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["Label", "Points:0", "Points:1", "Points:2", "activationTime"])
        for (label, x, y, z), value in zip(probes, values):
            writer.writerow([label, x, y, z, value])


def _write_line_csv(
    probes: list[tuple[str, float, float, float]], values: list[float], destination: Path,
) -> None:
    origin = probes[0][1:] if probes else (0.0, 0.0, 0.0)
    with destination.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["activationTime", "arc_length", "Points_0", "Points_1", "Points_2"])
        for (label, x, y, z), value in zip(probes, values):
            arc_length = math.dist(origin, (x, y, z))
            writer.writerow([value, arc_length, x, y, z])


def run_postprocessing(
    *,
    output_dir: str,
    setup_root: str | None = None,
    archive_subdir: str = "sweepCases",
    **_: object,
):
    """Convert this case's archived raw probe samples into labeled CSVs.

    driverFOAM's sweep loop (base.archive_dir_name in the sweep spec)
    archives each case's raw postProcessing/ output -- OpenFOAM's own
    `probes` functionObject format, coordinates in a comment header, one
    data row per sampled time -- directly into that case's own output_dir,
    under output_dir/<archive_subdir>/, before the next case's run
    overwrites the live (shared) postProcessing/ folder. That raw format
    isn't what line_postProcessing.py/points_postProcessing.py/
    table_summary.py read, though -- they need a labeled CSV. This reads
    this one case's archived raw samples and writes
    output_dir/points.csv and output_dir/line.csv.

    Operates on exactly one case's own output_dir per call (matching the
    PostprocessingProtocol contract) -- no cross-case scanning, no separate
    cache location. Each case already has its own directory; that's the
    only disambiguation needed.
    """
    output_path = Path(output_dir)
    archive_dir = output_path / archive_subdir

    if not archive_dir.is_dir():
        print(f"No archived raw postProcessing found at: {archive_dir}")
        return []

    artifacts: list[dict] = []

    points_raw = _latest_probes_file(archive_dir / "Niedererpoints")
    if points_raw is not None:
        probes, values = _parse_probes_file(points_raw)
        destination = output_path / "points.csv"
        _write_points_csv(probes, values, destination)
        artifacts.append({
            "path": str(destination), "label": "points", "kind": "table", "format": "csv",
        })

    lines_raw = _latest_probes_file(archive_dir / "Niedererlines")
    if lines_raw is not None:
        probes, values = _parse_probes_file(lines_raw)
        destination = output_path / "line.csv"
        _write_line_csv(probes, values, destination)
        artifacts.append({
            "path": str(destination), "label": "line", "kind": "table", "format": "csv",
        })

    print(f"Converted {len(artifacts)} raw probe sample file(s) for {output_path}.")
    return artifacts
