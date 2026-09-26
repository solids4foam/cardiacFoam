#!/usr/bin/env python3
#----------------------------------------------------------------------------#
# License
#     This file is part of cardiacFoam.
#
#     cardiacFoam is free software: you can redistribute it and/or modify it
#     under the terms of the GNU General Public License as published by the
#     Free Software Foundation, either version 3 of the License, or (at your
#     option) any later version.
#
#     cardiacFoam is distributed in the hope that it will be useful, but
#     WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#     General Public License for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with cardiacFoam.  If not, see <http://www.gnu.org/licenses/>.
#
# Script
#     postprocess_springStiffness
#----------------------------------------------------------------------------#
"""Summarise finished springSupportedSlab cases.

Every subdirectory of --input-dir that holds a finished case is read. The
script writes springStiffness_summary.csv and springStiffness.png to
--output-dir.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
import re

import numpy as np

try:
    import matplotlib

    matplotlib.use("Agg", force=True)
    import matplotlib.pyplot as plt
except ModuleNotFoundError:
    plt = None


def readParameter(case: Path, name: str) -> float:
    text = (case / "system" / "caseParameters").read_text()
    return float(re.search(rf"^\s*{name}\s+([^;]+);", text, re.M).group(1))


def slabExtent(case: Path) -> np.ndarray:
    """Slab extent (Lx, Ly, Lz) [m] from blockMeshDict."""
    text = (case / "system" / "blockMeshDict").read_text()
    scale = float(re.search(r"^\s*scale\s+([^;]+);", text, re.M).group(1))
    block = re.search(r"vertices\s*\((.*?)\n\);", text, re.S).group(1)
    pts = np.array([[float(v) for v in p.split()]
                    for p in re.findall(r"\(([^()]*)\)", block)])
    return (pts.max(axis=0) - pts.min(axis=0))*scale


def youngsModulus(case: Path) -> float:
    text = (case / "constant" / "solid" / "mechanicalProperties").read_text()
    return float(re.search(r"\bE\s+E\s+\[[^]]*\]\s+([^;]+);", text).group(1))


def readVectorSeries(path: Path) -> tuple[np.ndarray, np.ndarray]:
    t, v = [], []
    for line in path.read_text().splitlines():
        m = re.match(r"\s*([0-9.eE+-]+)\s+\(([^)]*)\)", line)
        if m:
            t.append(float(m.group(1)))
            v.append([float(a) for a in m.group(2).split()])
    return np.array(t), np.array(v)


def readColumns(path: Path) -> np.ndarray:
    rows = [[float(a) for a in line.split()]
            for line in path.read_text().splitlines()
            if line.strip() and not line.startswith("#")]
    return np.array(rows)


def readCase(case: Path) -> dict:
    pp = case / "postProcessing"
    tD, dMin = readVectorSeries(next(pp.glob("solid/D_xMin/*/surfaceFieldValue.dat")))
    _, dMax = readVectorSeries(next(pp.glob("solid/D_xMax/*/surfaceFieldValue.dat")))
    n = min(len(dMin), len(dMax))
    force = readColumns(next(pp.glob("*/solidForcesxMin.dat")))
    ta = readColumns(next(pp.glob("Taprobes/solid/*/Ta")))
    vm = readColumns(next(pp.glob("Vmprobes/electro/*/Vm")))

    k = readParameter(case, "kEnds")
    ext = slabExtent(case)
    A0, length = ext[1]*ext[2], ext[0]
    t = tD[:n]
    dxMin, dxMax = dMin[:n, 0], dMax[:n, 0]
    fStress = np.interp(t, force[:, 0], force[:, 1])
    fSpring = -k*A0*dxMin
    forceScale = ta[:, 1:].max()*A0

    return {
        "name": case.name,
        "kEnds": k,
        "deltaT": readParameter(case, "deltaT"),
        "t": t,
        "shortening": (dxMin - dxMax)/length,
        "fStress": fStress,
        "fSpring": fSpring,
        "ta": ta,
        "vm": vm,
        "kTransition": 2*youngsModulus(case)/length,
        "peakShortening": float(np.max((dxMin - dxMax)/length)),
        "peakForce": float(np.max(np.abs(fStress))),
        "springLawError": float(np.max(np.abs(fStress - fSpring))/forceScale),
    }


def plot(cases: list[dict], outFile: Path) -> None:
    fig, ax = plt.subplots(2, 2, figsize=(12, 8.5))
    (aTa, aS), (aF, aK) = ax
    ref = min(cases, key=lambda c: (abs(np.log10(c["kEnds"]) - 7), -c["deltaT"]))
    for i, label in ((1, "Ta, x = 0.75 mm"), (4, "Ta, x = 10 mm"), (5, "Ta, x = 19.25 mm")):
        aTa.plot(ref["ta"][:, 0]*1e3, ref["ta"][:, i]/1e3, label=label)
    aVm = aTa.twinx()
    for i, ls in ((1, ":"), (3, "--")):
        aVm.plot(ref["vm"][:, 0]*1e3, ref["vm"][:, i], ls, color="gray", lw=0.8)
    aVm.set_ylabel("Vm at x = 0.75 (:) and 19.25 mm (--)")
    aTa.set(xlabel="time [ms]", ylabel="Ta [kPa]", title=f"Activation ({ref['name']})")

    for c in cases:
        label = f"{c['name']} (k = {c['kEnds']:.0e}, dt = {c['deltaT']:.0e})"
        style = "--" if c["deltaT"] < 1e-5 else "-"
        aS.plot(c["t"]*1e3, c["shortening"]*100, style, label=label)
        line, = aF.plot(c["t"]*1e3, c["fStress"]*1e3, style, label=f"{c['name']}: from stress")
        step = max(1, len(c["t"])//30)
        aF.plot(c["t"][::step]*1e3, c["fSpring"][::step]*1e3, "o", ms=3,
                color=line.get_color(), label=f"{c['name']}: -k A0 <Dx>")

    main = sorted((c for c in cases if c["deltaT"] >= 1e-5), key=lambda c: c["kEnds"])
    aK.semilogx([c["kEnds"] for c in main], [c["peakShortening"]*100 for c in main],
                "o-", label="dt = 1e-5 s")
    for c in cases:
        if c["deltaT"] < 1e-5:
            aK.semilogx(c["kEnds"], c["peakShortening"]*100, "s", mfc="none",
                        label=f"dt = {c['deltaT']:.0e} s")
    aK.axvline(ref["kTransition"], ls=":", color="gray", label="k = 2E/L")

    aS.set(xlabel="time [ms]", ylabel="shortening [% of L]", title="Slab shortening")
    aF.set(xlabel="time [ms]", ylabel="axial force on xMin [mN]",
           title="Spring law: force from stress vs -k A0 <Dx>")
    aK.set(xlabel="spring stiffness kEnds [Pa/m]", ylabel="peak shortening [% of L]",
           title="Peak shortening vs spring stiffness")
    for a in ax.flat:
        a.grid(alpha=0.3)
        a.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(outFile, dpi=130)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, required=True,
                        help="directory whose subdirectories are finished cases")
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    cases = [readCase(c) for c in sorted(args.input_dir.iterdir())
             if (c / "postProcessing" / "solid" / "D_xMin").is_dir()]
    if not cases:
        raise SystemExit(f"no finished cases found under {args.input_dir}")
    cases.sort(key=lambda c: (-c["kEnds"], -c["deltaT"]))
    args.output_dir.mkdir(parents=True, exist_ok=True)

    summary = args.output_dir / "springStiffness_summary.csv"
    with summary.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["case", "kEnds_Pa_per_m", "deltaT_s", "peak_shortening_percent",
                    "peak_end_force_mN", "spring_law_error_over_TaMax_A0"])
        for c in cases:
            w.writerow([c["name"], f"{c['kEnds']:.3e}", f"{c['deltaT']:.1e}",
                        f"{c['peakShortening']*100:.3f}", f"{c['peakForce']*1e3:.3f}",
                        f"{c['springLawError']:.2e}"])
    print(summary.read_text())

    if plt is not None:
        plot(cases, args.output_dir / "springStiffness.png")


if __name__ == "__main__":
    main()
