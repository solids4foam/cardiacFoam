#----------------------------------------------------------------------------#
# Module
#     summarize
#
# Description
#     Tables of the bathBidomain studies from the outputs of the run copies in
#     a results folder: the verifier summary postProcessing/<d>D_<n>_cells.dat
#     and the interface diagnostics postProcessing/bathBidomainInterfaceMetrics
#     .csv (solved field) and bathBidomainInterfaceMetricsExactField.csv
#     (exact phiE). Run copies are grouped by their name with the level
#     removed: the token dt<value> is the level of a temporal ladder, else the
#     token N<n> the level of a spatial ladder; every other token names the
#     configuration. For every configuration it prints the field errors (Vm,
#     phiE, phiI L2 and Linf) and the interface rows (L2), with the order
#     between consecutive levels of h = (cells/3)^(-1/d) (cells and solution
#     directions d of the whole mesh from log.checkMesh, so 1/N on
#     hexahedra) or of dt. With
#     --relative-to <configuration> it prints instead the relative change of
#     every other configuration against that one at the same level.
#
# Usage
#     python3 summarize.py cartesianConvergence/results
#     python3 summarize.py coupling/results --relative-to predictor_tet
#----------------------------------------------------------------------------#

from __future__ import annotations

import argparse
import csv
import math
import re
from pathlib import Path

FIELDS = ("Vm", "phiE", "phiI")
INTERFACE = (
    "heartPhiE", "bathPhiE", "x0FluxJump", "x1FluxJump",
    "x0IntracellularLeak", "x1IntracellularLeak",
    "x0AssembledFlux", "x1AssembledFlux",
)
N_RE = re.compile(r"(?:^|_)N(\d+)(?=_|$)")
DT_RE = re.compile(r"(?:^|_)dt([0-9.eE+-]+?)(?=_|$)")


def read_summary(path: Path) -> dict[str, list[float]]:
    rows = {}
    for line in path.read_text().splitlines():
        words = line.split()
        if len(words) == 4 and not line.startswith("#"):
            try:
                rows[words[0]] = [float(w) for w in words[1:]]
            except ValueError:
                continue
    return rows


def read_metrics(path: Path) -> dict[str, float]:
    if not path.exists():
        return {}
    row = next(csv.DictReader(path.open()))
    return {k: float(row[k + "_L2"]) for k in INTERFACE if k + "_L2" in row}


def read_mesh(run: Path) -> tuple[int | None, int | None]:
    log = run / "log.checkMesh"
    text = log.read_text() if log.exists() else ""
    cells = re.search(r"^\s*cells:\s*(\d+)", text, re.M)
    dirs = re.search(r"Mesh has (\d) solution", text)
    return (int(cells.group(1)) if cells else None,
            int(dirs.group(1)) if dirs else None)


def load(results: Path) -> dict[str, list[dict]]:
    groups: dict[str, list[dict]] = {}
    for run in sorted(results.iterdir()):
        summaries = sorted(run.glob("postProcessing/*D_*_cells.dat"))
        if not summaries:
            continue
        dt, n = DT_RE.search(run.name), N_RE.search(run.name)
        fields = read_summary(summaries[0])
        if dt:
            level = spacing = float(dt.group(1))
            config = DT_RE.sub("", run.name)
        elif n:
            cells, dim = read_mesh(run)
            level = int(n.group(1))
            spacing = (cells/3)**(-1/dim) if cells and dim else 1/level
            config = N_RE.sub("", run.name)
        else:
            continue
        groups.setdefault(config.strip("_") or "default", []).append({
            "level": level,
            "h": spacing,
            "fields": fields,
            "solved": read_metrics(run/"postProcessing"/"bathBidomainInterfaceMetrics.csv"),
            "exact": read_metrics(run/"postProcessing"/"bathBidomainInterfaceMetricsExactField.csv"),
        })
    for runs in groups.values():
        runs.sort(key=lambda r: -r["h"])
    return groups


def order(coarse: float, fine: float, hCoarse: float, hFine: float) -> str:
    if coarse <= 0 or fine <= 0 or hCoarse == hFine:
        return ""
    return f"{math.log(coarse/fine)/math.log(hCoarse/hFine):.2f}"


def table(header: list[str], rows: list[list[str]]) -> str:
    lines = ["| " + " | ".join(header) + " |", "|" + "---|"*len(header)]
    lines += ["| " + " | ".join(r) + " |" for r in rows]
    return "\n".join(lines)


def field_table(runs: list[dict]) -> str:
    header = ["level", "h or dt"]
    for f in FIELDS:
        header += [f"{f} L2", f"{f} Linf", "order L2"]
    rows, previous = [], None
    for r in runs:
        cells = [str(r["level"]), f"{r['h']:.4g}"]
        for f in FIELDS:
            e = r["fields"].get(f)
            p = previous["fields"].get(f) if previous else None
            cells += [f"{e[1]:.4e}", f"{e[2]:.4e}"] if e else ["", ""]
            cells.append(order(p[1], e[1], previous["h"], r["h"]) if e and p else "")
        rows.append(cells)
        previous = r
    return table(header, rows)


def interface_table(runs: list[dict], source: str) -> str:
    keys = [k for k in INTERFACE if any(k in r[source] for r in runs)]
    if not keys:
        return ""
    header = ["level"] + [f"{k} L2" for k in keys] + [f"order {k}" for k in keys]
    rows, previous = [], None
    for r in runs:
        values = [f"{r[source][k]:.4e}" if k in r[source] else "" for k in keys]
        orders = [
            order(previous[source][k], r[source][k], previous["h"], r["h"])
            if previous and k in previous[source] and k in r[source] else ""
            for k in keys
        ]
        rows.append([str(r["level"])] + values + orders)
        previous = r
    return table(header, rows)


def relative_table(base: list[dict], other: list[dict]) -> str:
    reference = {r["level"]: r for r in base}
    keys = [f"{f} L2" for f in FIELDS] + [f"{k} L2" for k in INTERFACE]
    rows = []
    for r in other:
        b = reference.get(r["level"])
        if not b:
            continue
        cells = [str(r["level"])]
        for f in FIELDS:
            x, y = r["fields"].get(f), b["fields"].get(f)
            cells.append(f"{100*(x[1]/y[1] - 1):+.2f}%" if x and y and y[1] else "")
        for k in INTERFACE:
            x, y = r["solved"].get(k), b["solved"].get(k)
            cells.append(f"{100*(x/y - 1):+.2f}%" if x is not None and y else "")
        rows.append(cells)
    return table(["level"] + keys, rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("results", type=Path, help="folder of run copies")
    parser.add_argument("--relative-to", metavar="CONFIGURATION",
                        help="print relative changes against this configuration")
    args = parser.parse_args()
    groups = load(args.results)
    if not groups:
        raise SystemExit("no run copies with a verifier summary found")
    if args.relative_to:
        if args.relative_to not in groups:
            raise SystemExit(f"no configuration {args.relative_to} in {sorted(groups)}")
        for config in sorted(groups):
            if config != args.relative_to:
                print(f"\n### {config} against {args.relative_to}\n")
                print(relative_table(groups[args.relative_to], groups[config]))
        return
    for config in sorted(groups):
        print(f"\n### {config}\n")
        print(field_table(groups[config]))
        for source, title in (("solved", "interface, solved field"),
                              ("exact", "interface, exact phiE")):
            text = interface_table(groups[config], source)
            if text:
                print(f"\n{title}:\n")
                print(text)


if __name__ == "__main__":
    main()
