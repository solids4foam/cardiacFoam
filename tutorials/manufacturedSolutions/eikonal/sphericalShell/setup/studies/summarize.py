#----------------------------------------------------------------------------#
# Module
#     summarize
#
# Description
#     Convergence tables of the eikonalECG sphericalShell studies from the
#     verifier summaries postProcessing/3D_<mesh>_<N>_cells.dat and the
#     solver logs of the run copies in a results folder. Run copies are named
#     <sealed|stock>_<hex|tet>_N<n>_s<seed>[_ex][_exact|_best|_cap|_parallel]
#     (wall gradient, mesh, N, gmsh seed, explicit-metric formulation, one
#     iteration from the exact field, stopped at the iteration of the smallest
#     outer residual, stopped at the iteration cap of a run that does not
#     converge, or a parallel run). For every wall gradient, formulation and
#     window it prints the error table with the outer iterations, the smallest
#     outer residual and the iteration where it occurs, and, except for the
#     capped runs, the pairwise order between every level N and the level 2N
#     with the ratio r of their mean h, and the order fitted by least squares
#     to log(error) against log(h) over every run of the three finest levels
#     with its standard error; for a parallel run, the largest relative
#     difference from the serial run over all rows.
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
    ("activationTime", "L2"),
    ("activationTime", "Linf"),
    ("activationTime_boundary", "L2"),
    ("gradActivationTime_boundary", "L2"),
    ("gradReconActivationTime_boundary", "L2"),
)
NORMS = {"L1": 0, "L2": 1, "Linf": 2}
RUN_RE = re.compile(
    r"^(sealed|stock)_(hex|tet)_N(\d+)_s(\d+)(_ex)?(?:_(exact|best|cap|parallel))?$"
)
RESIDUAL_RE = re.compile(r"Initial residual = ([0-9.eE+-]+),")


def read_summary(path: Path, advection: bool) -> tuple[float, dict[str, list[float]]]:
    h = None
    rows: dict[str, list[float]] = {}
    for line in path.read_text().splitlines():
        words = line.split()
        if line.startswith("# h "):
            h = float(words[2])
        elif line.startswith("# eikonalAdvectionDiffusionApproach "):
            if (words[2] in ("true", "on", "yes", "1")) != advection:
                raise SystemExit(f"{path}: the run name and the solved formulation differ")
        elif not line.startswith("#") and len(words) == 4:
            rows[words[0]] = [float(w) for w in words[1:]]
    if h is None:
        raise SystemExit(f"{path}: no '# h' line")
    return h, rows


def read_residuals(run: Path) -> list[float]:
    log = run / "log.cardiacFoam"
    if not log.exists():
        return []
    return [float(r) for r in RESIDUAL_RE.findall(log.read_text())]


def load(results: Path) -> dict[tuple, list]:
    runs: dict[tuple, list] = {}
    for run in sorted(results.iterdir()):
        match = RUN_RE.match(run.name)
        summaries = sorted(run.glob("postProcessing/3D_*_cells.dat"))
        if not match or not summaries:
            continue
        gradient, mesh, n, seed, explicit, window = match.groups()
        key = (gradient, mesh, "explicit" if explicit else "advectionDiffusion",
               window or "converged")
        runs.setdefault(key, []).append(
            (int(n), int(seed), *read_summary(summaries[0], not explicit), read_residuals(run))
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


def history(residuals: list[float]) -> str:
    if not residuals:
        return "- | - | -"
    low = min(residuals)
    return f"{len(residuals)} | {low:.1e} @ {residuals.index(low) + 1} | {residuals[-1]:.1e}"


def table(runs: list, capped: bool = False) -> str:
    header = " | ".join(f"{r} {norm}" for r, norm in COLUMNS)
    levels = sorted({n for n, *_ in runs})
    lines = [
        f"| N | h | {header} | outer its | min residual | last residual |",
        "|---" * (len(COLUMNS) + 5) + "|",
    ]
    for level in levels:
        at = [r for r in runs if r[0] == level]
        h = statistics.mean(r[2] for r in at)
        cells = " | ".join(
            f"{statistics.mean(value(r[3], c) for r in at):.3e}"
            for c in COLUMNS
        )
        lines.append(f"| {level} | {h:.4f} | {cells} | {history(at[0][4])} |")
    pairs = [(n, 2*n) for n in levels if 2*n in levels]
    if pairs and not capped:
        mean_h = {n: statistics.mean(r[2] for r in runs if r[0] == n) for n in levels}
        lines += [
            "",
            f"| N | r | {header} |",
            "|---" * (len(COLUMNS) + 2) + "|",
        ]
        for coarse, fine in pairs:
            ratio = mean_h[coarse]/mean_h[fine]
            orders = " | ".join(
                "%.2f" % (
                    math.log(
                        statistics.mean(value(r[3], c) for r in runs if r[0] == coarse)
                        / statistics.mean(value(r[3], c) for r in runs if r[0] == fine)
                    )/math.log(ratio)
                )
                for c in COLUMNS
            )
            lines.append(f"| {coarse} to {fine} | {ratio:.3f} | {orders} |")
    if len(levels) >= 3 and not capped:
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

    for key in sorted(runs):
        gradient, mesh, formulation, window = key
        if window == "parallel":
            continue
        print(f"\n### {gradient}, {mesh}, {formulation}, {window}\n")
        print(table(runs[key], window == "cap"))

    for (gradient, mesh, formulation, window), parallel in sorted(runs.items()):
        if window != "parallel":
            continue
        for n, seed, _, rows, _ in parallel:
            serial = [
                r for r in runs.get((gradient, mesh, formulation, "converged"), [])
                if r[0] == n and r[1] == seed
            ]
            if serial:
                difference = max(
                    abs(a/b - 1)
                    for row in rows if row in serial[0][3]
                    for a, b in zip(rows[row], serial[0][3][row]) if b != 0
                )
                print(
                    f"\nParallel against serial, {gradient} {mesh} {formulation} N {n} "
                    f"seed {seed}: largest relative difference {difference:.1e}"
                )

    for (gradient, mesh, formulation, window), at in sorted(runs.items()):
        seeds = {r[1] for r in at}
        if len(seeds) < 2:
            continue
        for level in sorted({r[0] for r in at}):
            values = [r for r in at if r[0] == level]
            if len(values) < 2:
                continue
            spread = ", ".join(
                f"{c[0]} {c[1]} {100*(max(value(v[3], c) for v in values)/min(value(v[3], c) for v in values) - 1):.2f} %"
                for c in COLUMNS[:1]
            )
            print(f"Seed spread, {gradient} {mesh} {formulation} {window} N {level}: {spread}")


if __name__ == "__main__":
    main()
