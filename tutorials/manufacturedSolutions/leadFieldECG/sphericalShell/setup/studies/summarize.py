#----------------------------------------------------------------------------#
# Module
#     summarize
#
# Description
#     Tables of the leadFieldECG/sphericalShell studies from the verifier
#     outputs of the run copies in a results folder, named <hex|tet>_N<n>.
#     For each mesh family it prints, per level: the static lead-field error
#     against the exact potential at each electrode cell centroid; the
#     largest gap between the static lead field and the static direct solve;
#     the RMS and maximum error of the direct static potential field over
#     the whole mesh, the heart, the torso and the interface band; and the
#     orders between consecutive levels of h = cbrt(mean cell volume).
#
# Usage
#     python3 summarize.py curvedConvergence/results
#----------------------------------------------------------------------------#

from __future__ import annotations

import argparse
import math
import re
from pathlib import Path

RUN_RE = re.compile(r"^(hex|tet)_N(\d+)$")
REGIONS = ("phiStatic", "phiStatic_heart", "phiStatic_torso",
           "phiStatic_interfaceBand")


def read_run(run: Path):
    post = run / "postProcessing"
    series = (post / "manufacturedLeadFieldECG_ECG.dat").read_text().splitlines()
    header = [l for l in series if l.startswith("#")][-1].split()[1:]
    values = dict(zip(header, map(float, [l for l in series
                                          if l.strip() and not l.startswith("#")][-1].split())))
    summary = (post / "manufacturedLeadFieldECGSummary_ECG.dat").read_text()
    rows = {l.split()[0]: list(map(float, l.split()[1:]))
            for l in summary.splitlines() if l.startswith("phiStatic")}
    cells = int(re.search(r"^# cells (\d+)", summary, re.M).group(1))
    electrodes = [k[len("leadFieldStatic_"):] for k in header
                  if k.startswith("leadFieldStatic_")]
    return values, rows, cells, electrodes


def order(a: float, b: float, ha: float, hb: float) -> str:
    if a <= 0 or b <= 0 or ha == hb:
        return ""
    return f"{math.log(a/b)/math.log(ha/hb):.2f}"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("results", type=Path, help="folder of run copies")
    results = parser.parse_args().results
    runs: dict[str, list] = {}
    for run in sorted(results.iterdir()):
        match = RUN_RE.match(run.name)
        if not match or not (run / "postProcessing").is_dir():
            continue
        values, rows, cells, electrodes = read_run(run)
        volume = (4.0/3.0)*math.pi*1.5**3
        h = (volume/cells)**(1.0/3.0)
        runs.setdefault(match.group(1), []).append(
            (int(match.group(2)), h, values, rows, electrodes))

    for mesh in ("hex", "tet"):
        if mesh not in runs:
            continue
        levels = sorted(runs[mesh], key=lambda r: r[0])
        electrodes = levels[0][4]
        print(f"\n### {mesh}: static lead-field error at the electrode cell centroid\n")
        print("| N | h | " + " | ".join(electrodes)
              + " | lead field - direct |")
        print("|---" * (len(electrodes) + 3) + "|")
        prev = None
        for n, h, v, rows, _ in levels:
            errs = [abs(v[f"leadFieldStatic_{e}"] - v[f"exactCentroid_{e}"])
                    for e in electrodes]
            gap = max(abs(v[f"leadFieldStatic_{e}"] - v[f"directStatic_{e}"])
                      for e in electrodes)
            cells = " | ".join(
                f"{x:.3e}" + (f" ({order(prev[2][i], x, prev[1], h)})" if prev else "")
                for i, x in enumerate(errs))
            print(f"| {n} | {h:.4f} | {cells} | {gap:.1e} |")
            prev = (n, h, errs)

        print(f"\n### {mesh}: error of the direct static potential field\n")
        print("| N | " + " | ".join(f"{r} RMS | {r} max" for r in REGIONS) + " |")
        print("|---" * (2*len(REGIONS) + 1) + "|")
        prev = None
        for n, h, v, rows, _ in levels:
            cells = []
            for r in REGIONS:
                for k in (1, 2):
                    x = rows[r][k]
                    o = f" ({order(prev[1][r][k], x, prev[0], h)})" if prev else ""
                    cells.append(f"{x:.3e}{o}")
            print(f"| {n} | " + " | ".join(cells) + " |")
            prev = (h, rows)


if __name__ == "__main__":
    main()
