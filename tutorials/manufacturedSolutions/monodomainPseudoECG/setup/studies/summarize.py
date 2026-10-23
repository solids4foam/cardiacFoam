#----------------------------------------------------------------------------#
# Module
#     summarize
#
# Description
#     Convergence tables of the monodomainPseudoECG studies from the verifier
#     summaries of the run copies in a results folder: the field summary
#     postProcessing/<d>D_<n>_cells.dat and, in 3D, the manufactured pseudo-ECG
#     summary postProcessing/manufacturedPseudoECGSummary_ECG.dat. Run copies
#     are grouped by their name with the level removed: the token dt<value> is
#     the level of a temporal ladder, else the token N<n> the level of a spatial
#     ladder; every other token names the configuration. For every
#     configuration it prints Vm, u1 and u2 L2 and Linf, the largest pseudo-ECG
#     Linf error over E1-E5, and the order between consecutive levels (of
#     h = 1/n for hexahedra, the verifier's h on tetrahedra, or dt).
#
# Usage
#     python3 summarize.py cartesianConvergence/results
#----------------------------------------------------------------------------#

from __future__ import annotations

import argparse
import math
import re
from pathlib import Path

FIELDS = ("Vm", "u1", "u2")
ELECTRODES = ("E1", "E2", "E3", "E4", "E5")
N_RE = re.compile(r"(?:^|_)N(\d+)(?=_|$)")
DT_RE = re.compile(r"(?:^|_)dt([0-9.eE+-]+?)(?=_|$)")


def read_rows(path: Path) -> dict[str, list[float]]:
    rows: dict[str, list[float]] = {}
    for line in path.read_text().splitlines():
        words = line.split()
        if len(words) < 4:
            continue
        try:
            rows.setdefault(words[0], [float(w) for w in words[1:4]])
        except ValueError:
            continue
    return rows


def read_spacing(path: Path) -> float | None:
    match = re.search(r"Grid spacing \(dx\)\s*=\s*([0-9.eE+-]+)", path.read_text())
    return float(match.group(1)) if match else None


def load(results: Path) -> dict[str, list[tuple]]:
    groups: dict[str, list[tuple]] = {}
    for run in sorted(results.iterdir()):
        summaries = sorted(run.glob("postProcessing/*D_*_cells.dat"))
        if not summaries:
            continue
        dt, n = DT_RE.search(run.name), N_RE.search(run.name)
        if dt:
            level = spacing = float(dt.group(1))
            config = DT_RE.sub("", run.name)
        elif n:
            level, spacing = int(n.group(1)), read_spacing(summaries[0])
            config = N_RE.sub("", run.name)
        else:
            continue
        config = config.strip("_") or "default"
        ecg = run / "postProcessing" / "manufacturedPseudoECGSummary_ECG.dat"
        ecgRows = read_rows(ecg) if ecg.exists() else {}
        ecgLinf = max((ecgRows[e][2] for e in ELECTRODES if e in ecgRows), default=None)
        groups.setdefault(config, []).append((level, spacing, read_rows(summaries[0]), ecgLinf))
    return groups


def order(coarse: float, fine: float, hCoarse: float, hFine: float) -> str:
    if coarse <= 0 or fine <= 0 or hCoarse == hFine:
        return ""
    return f"{math.log(coarse/fine)/math.log(hCoarse/hFine):.2f}"


def table(runs: list[tuple]) -> str:
    runs = sorted(runs, key=lambda r: -r[1])
    header = ["level", "h or dt"]
    for f in FIELDS:
        header += [f"{f} L2", f"{f} Linf"]
    header += ["pseudo-ECG max Linf", "order Vm L2", "order Vm Linf", "order ECG"]
    lines = ["| " + " | ".join(header) + " |", "|" + "---|"*len(header)]
    previous = None
    for level, spacing, rows, ecg in runs:
        cells = [str(level), f"{spacing:.6g}"]
        for f in FIELDS:
            cells += [f"{rows[f][1]:.4e}", f"{rows[f][2]:.4e}"] if f in rows else ["", ""]
        cells.append(f"{ecg:.4e}" if ecg is not None else "")
        if previous:
            _, h0, rows0, ecg0 = previous
            cells += [order(rows0["Vm"][1], rows["Vm"][1], h0, spacing),
                      order(rows0["Vm"][2], rows["Vm"][2], h0, spacing),
                      order(ecg0, ecg, h0, spacing) if ecg and ecg0 else ""]
        else:
            cells += ["", "", ""]
        lines.append("| " + " | ".join(cells) + " |")
        previous = (level, spacing, rows, ecg)
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
