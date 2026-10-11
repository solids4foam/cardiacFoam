#----------------------------------------------------------------------------#
# Module
#     summarize
#
# Description
#     Convergence tables of the monodomain1D3D studies from the verifier
#     summaries of the run copies in a results folder: the 3D field summary
#     postProcessing/3D_<n>_cells.dat, the graph summary
#     postProcessing/graph_1D_<k>_nodes.dat and the PVJ source diagnostics
#     verification/coupled1D3DMonodomain_diagnostics.csv. Run copies are
#     grouped by their name with the level removed: the token N<n> is the level
#     of a coupled ladder (h = 1/n), else the token nodes<k> the level of a
#     graph ladder (h = 1/(k - 1), two branches of length 1/2); every other
#     token names the configuration. For every configuration it prints the L2
#     and Linf errors of the 3D Vm (coupled runs), of the graph Vm1D and u1, the
#     L1 PVJ source error (coupled runs), and the order between consecutive
#     levels.
#
# Usage
#     python3 summarize.py coupledConvergence/results
#----------------------------------------------------------------------------#

from __future__ import annotations

import argparse
import csv
import math
import re
from pathlib import Path

N_RE = re.compile(r"(?:^|_)N(\d+)(?=_|$)")
NODES_RE = re.compile(r"(?:^|_)nodes(\d+)(?=_|$)")


def read_rows(path: Path) -> dict[str, list[float]]:
    rows: dict[str, list[float]] = {}
    for line in path.read_text().splitlines():
        words = line.split()
        if len(words) != 4:
            continue
        try:
            rows.setdefault(words[0], [float(w) for w in words[1:4]])
        except ValueError:
            continue
    return rows


def read_source_error(path: Path) -> float | None:
    if not path.exists():
        return None
    rows = list(csv.DictReader(path.open()))
    return float(rows[-1]["sourceErrorL1"]) if rows else None


def load(results: Path) -> dict[str, list[tuple]]:
    groups: dict[str, list[tuple]] = {}
    for run in sorted(results.iterdir()):
        graph = sorted(run.glob("postProcessing/graph_1D_*_nodes.dat"))
        if not graph:
            continue
        tissue = sorted(run.glob("postProcessing/3D_*_cells.dat"))
        n, nodes = N_RE.search(run.name), NODES_RE.search(run.name)
        if n:
            level, h = int(n.group(1)), 1/int(n.group(1))
            config = N_RE.sub("", run.name)
        elif nodes:
            level, h = int(nodes.group(1)), 1/(int(nodes.group(1)) - 1)
            config = NODES_RE.sub("", run.name)
        else:
            continue
        config = config.strip("_") or "default"
        groups.setdefault(config, []).append((
            level,
            h,
            read_rows(tissue[0]) if tissue else {},
            read_rows(graph[0]),
            read_source_error(run/"verification"/"coupled1D3DMonodomain_diagnostics.csv"),
        ))
    return groups


def order(coarse: float | None, fine: float | None, hCoarse: float, hFine: float) -> str:
    if not coarse or not fine or coarse <= 0 or fine <= 0 or hCoarse == hFine:
        return ""
    return f"{math.log(coarse/fine)/math.log(hCoarse/hFine):.2f}"


def table(runs: list[tuple]) -> str:
    runs = sorted(runs, key=lambda r: -r[1])
    coupled = any(r[2] for r in runs)
    columns = []
    if coupled:
        columns += [("3D Vm L2", lambda r: r[2]["Vm"][1]), ("3D Vm Linf", lambda r: r[2]["Vm"][2])]
    columns += [("1D Vm L2", lambda r: r[3]["Vm1D"][1]), ("1D Vm Linf", lambda r: r[3]["Vm1D"][2]),
                ("1D u1 L2", lambda r: r[3]["u1"][1])]
    if coupled:
        columns += [("PVJ source error L1", lambda r: r[4])]
    header = ["level", "h"] + [c[0] for c in columns] + [f"order {c[0]}" for c in columns]
    lines = ["| " + " | ".join(header) + " |", "|" + "---|"*len(header)]
    previous = None
    for run in runs:
        values = [c[1](run) for c in columns]
        cells = [str(run[0]), f"{run[1]:.6g}"]
        cells += [f"{v:.4e}" if v is not None else "" for v in values]
        if previous:
            cells += [order(v0, v, previous[0], run[1]) for v0, v in zip(previous[1], values)]
        else:
            cells += [""]*len(columns)
        lines.append("| " + " | ".join(cells) + " |")
        previous = (run[1], values)
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("results", type=Path, help="folder of run copies")
    groups = load(parser.parse_args().results)
    if not groups:
        raise SystemExit("no run copies with a verifier summary found")
    for config in sorted(groups):
        print(f"\n### {config}\n")
        print(table(groups[config]))


if __name__ == "__main__":
    main()
