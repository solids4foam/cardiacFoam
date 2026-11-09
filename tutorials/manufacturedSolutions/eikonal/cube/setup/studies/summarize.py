#----------------------------------------------------------------------------#
# Module
#     summarize
#
# Description
#     Convergence tables of the eikonalECG studies from the verifier summaries
#     of the run copies in a results folder: the activation-time summary
#     postProcessing/manufacturedEikonalActivationTime.dat and the eikonal ECG
#     summary postProcessing/manufacturedEikonalECGSummary_ECG.dat. Run copies
#     are grouped by their name with the token N<n> (the level) removed; every
#     other token names the configuration. For every configuration it prints
#     the activation-time L2 and Linf, the largest and mean over E1-E5 of the
#     ECG L2 and Linf errors, the outer iterations, the bulk/boundary split of
#     the activation-time L2 and of the reconstruction of the exact gradient
#     (gradReconActivationTime), and the order between consecutive levels of
#     h = 1/n, with n the cells per direction on hexahedra and
#     round(cbrt(cells)) on tetrahedra.
#
# Usage
#     python3 summarize.py tetConvergence/results
#----------------------------------------------------------------------------#

from __future__ import annotations

import argparse
import math
import re
from pathlib import Path

ELECTRODES = ("E1", "E2", "E3", "E4", "E5")
N_RE = re.compile(r"(?:^|_)N(\d+)(?=_|$)")


def read_rows(path: Path) -> dict[str, list[float]]:
    rows: dict[str, list[float]] = {}
    for line in path.read_text().splitlines():
        words = line.split()
        if len(words) < 4 or words[0].startswith("#"):
            continue
        try:
            rows.setdefault(words[0], [float(w) for w in words[1:4]])
        except ValueError:
            continue
    return rows


def read_header(path: Path) -> tuple[int, int]:
    text = path.read_text()
    cells = int(re.search(r"Number of cells = (\d+)", text).group(1))
    dim = int(re.search(r"# dimension (\d)D", text).group(1))
    return cells, dim


def spacing(cells: int, dim: int) -> float:
    n = round(cells ** (1.0/dim))
    if n**dim != cells:
        n = round(cells ** (1.0/3.0))
    return 1.0/n


def outer_iterations(run: Path) -> int | None:
    log = run / "log.cardiacFoam"
    if not log.exists():
        return None
    return sum(1 for line in log.open() if "Solving for activationTime" in line)


def load(results: Path) -> dict[str, list[dict]]:
    groups: dict[str, list[dict]] = {}
    for run in sorted(results.iterdir()):
        summary = run / "postProcessing" / "manufacturedEikonalActivationTime.dat"
        n = N_RE.search(run.name)
        if not summary.exists() or not n:
            continue
        cells, dim = read_header(summary)
        ecg = run / "postProcessing" / "manufacturedEikonalECGSummary_ECG.dat"
        ecgRows = read_rows(ecg) if ecg.exists() else {}
        electrodes = [ecgRows[e] for e in ELECTRODES if e in ecgRows]
        config = N_RE.sub("", run.name).strip("_") or "default"
        groups.setdefault(config, []).append({
            "N": int(n.group(1)),
            "cells": cells,
            "h": spacing(cells, dim),
            "rows": read_rows(summary),
            "ecgL2": [e[1] for e in electrodes],
            "ecgLinf": [e[2] for e in electrodes],
            "outer": outer_iterations(run),
        })
    return groups


def order(coarse: float | None, fine: float | None, hCoarse: float, hFine: float) -> str:
    if not coarse or not fine or coarse <= 0 or fine <= 0 or hCoarse == hFine:
        return ""
    return f"{math.log(coarse/fine)/math.log(hCoarse/hFine):.2f}"


def markdown(header: list[str], runs: list[dict], columns) -> str:
    lines = ["| " + " | ".join(header) + " |", "|" + "---|"*len(header)]
    previous = None
    for run in runs:
        lines.append("| " + " | ".join(columns(run, previous)) + " |")
        previous = run
    return "\n".join(lines)


def value(run: dict, row: str, norm: int) -> float | None:
    return run["rows"][row][norm] if row in run["rows"] else None


def fmt(x: float | None) -> str:
    return f"{x:.4e}" if x is not None else ""


def field_table(runs: list[dict]) -> str:
    def columns(run, previous):
        L2, Linf = value(run, "activationTime", 1), value(run, "activationTime", 2)
        cells = [str(run["N"]), str(run["cells"]), f"{run['h']:.4g}", fmt(L2), fmt(Linf)]
        if previous:
            cells += [order(value(previous, "activationTime", 1), L2, previous["h"], run["h"]),
                      order(value(previous, "activationTime", 2), Linf, previous["h"], run["h"])]
        else:
            cells += ["", ""]
        cells.append(str(run["outer"]) if run["outer"] is not None else "")
        return cells
    return markdown(["N", "cells", "h", "psi L2", "psi Linf", "order L2", "order Linf",
                     "outer iterations"], runs, columns)


def ecg_table(runs: list[dict]) -> str:
    def stats(run):
        L2, Linf = run["ecgL2"], run["ecgLinf"]
        if not L2:
            return [None]*4
        return [max(L2), sum(L2)/len(L2), max(Linf), sum(Linf)/len(Linf)]

    def columns(run, previous):
        now = stats(run)
        cells = [str(run["N"])] + [fmt(x) for x in now]
        before = stats(previous) if previous else [None]*4
        cells += [order(before[i], now[i], previous["h"], run["h"]) if previous else ""
                  for i in (0, 2)]
        return cells
    return markdown(["N", "ECG max L2", "ECG mean L2", "ECG max Linf", "ECG mean Linf",
                     "order max L2", "order max Linf"], runs, columns)


def split_table(runs: list[dict]) -> str:
    names = ("activationTime_bulk", "activationTime_boundary", "gradReconActivationTime_bulk",
             "gradReconActivationTime_boundary", "gradReconActivationTime")

    def columns(run, previous):
        cells = [str(run["N"])]
        for name in names:
            cells.append(fmt(value(run, name, 1)))
        for name in names:
            cells.append(order(value(previous, name, 1), value(run, name, 1),
                               previous["h"], run["h"]) if previous else "")
        return cells
    short = ("psi bulk", "psi boundary", "gradRecon bulk", "gradRecon boundary", "gradRecon")
    return markdown(["N"] + [f"{s} L2" for s in short] + [f"order {s}" for s in short],
                    runs, columns)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("results", type=Path, help="folder of run copies")
    groups = load(parser.parse_args().results)
    if not groups:
        raise SystemExit("no run copies with a verifier summary found")
    for config in sorted(groups):
        runs = sorted(groups[config], key=lambda r: -r["h"])
        print(f"\n### {config}\n")
        print(field_table(runs))
        if any(run["ecgL2"] for run in runs):
            print()
            print(ecg_table(runs))
        print()
        print(split_table(runs))


if __name__ == "__main__":
    main()
