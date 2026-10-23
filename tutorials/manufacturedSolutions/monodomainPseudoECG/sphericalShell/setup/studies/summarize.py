#----------------------------------------------------------------------------#
# Module
#     summarize
#
# Description
#     Convergence tables of the monodomainPseudoECG sphericalShell studies from the
#     verifier summaries postProcessing/3D_<mesh>_<N>_cells.dat of the run
#     copies in a results folder. Run copies are named
#     <sealed|stock>_<hex|tet>_N<n>_s<seed>[_T2|_dthalf|_serial] (wall
#     gradient, mesh, N, gmsh seed, long window, halved time step or serial
#     run). For every wall gradient and window it prints the error table
#     (means over seeds), the order fitted by least squares to log(error)
#     against log(h) over every run of the three finest levels with its
#     standard error, the change under the halved time step, the largest
#     relative difference between a serial and the parallel run over all rows,
#     and the spread over seeds.
#
# Usage
#     python3 summarize.py hexConvergence/results
#----------------------------------------------------------------------------#

from __future__ import annotations

import argparse
import math
import re
import statistics
from pathlib import Path

COLUMNS = (
    ("Vm", "L2"),
    ("Vm", "Linf"),
    ("Vm_boundary", "L2"),
    ("gradVm_boundary", "L2"),
    ("gradReconVm_boundary", "L2"),
)
SPREAD_ROWS = ("Vm",)
NORMS = {"L1": 0, "L2": 1, "Linf": 2}
RUN_RE = re.compile(
    r"^(sealed|stock)_(hex|tet)_N(\d+)_s(\d+)(?:_(T2|dthalf|serial))?$"
)


def read_summary(path: Path) -> tuple[float, dict[str, list[float]]]:
    h = None
    rows: dict[str, list[float]] = {}
    for line in path.read_text().splitlines():
        words = line.split()
        if line.startswith("# h "):
            h = float(words[2])
        elif not line.startswith("#") and len(words) == 4:
            rows[words[0]] = [float(w) for w in words[1:]]
    if h is None:
        raise SystemExit(f"{path}: no '# h' line")
    return h, rows


def load(results: Path) -> dict[tuple, list]:
    runs: dict[tuple, list] = {}
    for run in sorted(results.iterdir()):
        match = RUN_RE.match(run.name)
        summaries = sorted(run.glob("postProcessing/3D_*_cells.dat"))
        if not match or not summaries:
            continue
        gradient, mesh, n, seed, window = match.groups()
        key = (gradient, mesh, window or "")
        runs.setdefault(key, []).append(
            (int(n), int(seed), *read_summary(summaries[0]))
        )
    return runs


def value(rows: dict, column: tuple[str, str]) -> float:
    return rows[column[0]][NORMS[column[1]]]


def fitted_order(points: list[tuple[float, float]]) -> tuple[float, float]:
    x = [math.log(h) for h, _ in points]
    y = [math.log(e) for _, e in points]
    n = len(x)
    xm, ym = sum(x)/n, sum(y)/n
    sxx = sum((xi - xm)**2 for xi in x)
    slope = sum((xi - xm)*(yi - ym) for xi, yi in zip(x, y))/sxx
    if n <= 2:
        return slope, float("nan")
    residual = sum((yi - ym - slope*(xi - xm))**2 for xi, yi in zip(x, y))
    return slope, math.sqrt(residual/(n - 2)/sxx)


def table(runs: list) -> str:
    header = " | ".join(f"{r} {norm}" for r, norm in COLUMNS)
    levels = sorted({n for n, *_ in runs})
    lines = [
        f"| N | h | {header} |",
        "|---" * (len(COLUMNS) + 2) + "|",
    ]
    for level in levels:
        at = [r for r in runs if r[0] == level]
        h = statistics.mean(r[2] for r in at)
        cells = " | ".join(
            f"{statistics.mean(value(r[3], c) for r in at):.3e}"
            for c in COLUMNS
        )
        lines.append(f"| {level} | {h:.4f} | {cells} |")
    finest = levels[-3:]
    tail = [r for r in runs if r[0] in finest]
    orders = " | ".join(
        "%.2f ± %.2f" % fitted_order([(r[2], value(r[3], c)) for r in tail])
        for c in COLUMNS
    )
    lines += [
        "",
        f"| N | {header} |",
        "|---" * (len(COLUMNS) + 1) + "|",
        f"| {', '.join(map(str, finest))} | {orders} |",
    ]
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("results", type=Path, help="folder of run copies")
    runs = load(parser.parse_args().results)
    if not runs:
        raise SystemExit("no run copies with a verifier summary found")

    for gradient in ("sealed", "stock"):
        for window, label in (("", "T = 0.2"), ("T2", "T = 2")):
            for mesh in ("hex", "tet"):
                key = (gradient, mesh, window)
                if len({r[0] for r in runs.get(key, [])}) >= 3:
                    print(f"\n### {gradient}, {mesh}, {label}\n")
                    print(table(runs[key]))

    for (gradient, mesh, window), halved in sorted(runs.items()):
        if window != "dthalf":
            continue
        for n, seed, _, rows in halved:
            full = [
                r for r in runs.get((gradient, mesh, ""), [])
                if r[0] == n and r[1] == seed
            ]
            if full:
                change = ", ".join(
                    f"{row} L2 {100*(rows[row][1]/full[0][3][row][1] - 1):+.2f} %"
                    for row in SPREAD_ROWS
                )
                print(f"\nHalved dt, {gradient} {mesh} N {n} seed {seed}: {change}")

    for (gradient, mesh, window), serial in sorted(runs.items()):
        if window != "serial":
            continue
        for n, seed, _, rows in serial:
            full = [
                r for r in runs.get((gradient, mesh, ""), [])
                if r[0] == n and r[1] == seed
            ]
            if full:
                difference = max(
                    abs(a/b - 1)
                    for row in rows if row in full[0][3]
                    for a, b in zip(rows[row], full[0][3][row]) if b != 0
                )
                print(
                    f"\nSerial against parallel, {gradient} {mesh} N {n} seed {seed}: "
                    f"largest relative difference {difference:.1e}"
                )

    for (gradient, mesh, window), at in sorted(runs.items()):
        seeds = {r[1] for r in at}
        if window or len(seeds) < 2:
            continue
        for level in sorted({r[0] for r in at}):
            values = [r for r in at if r[0] == level]
            if len(values) < 2:
                continue
            spread = ", ".join(
                f"{row} L2 {100*(max(v[3][row][1] for v in values)/min(v[3][row][1] for v in values) - 1):.2f} %"
                for row in SPREAD_ROWS
            )
            print(f"Seed spread, {gradient} {mesh} N {level}: {spread}")


if __name__ == "__main__":
    main()
