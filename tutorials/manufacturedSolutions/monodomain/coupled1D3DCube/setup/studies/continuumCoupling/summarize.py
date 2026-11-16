#!/usr/bin/env python3
"""Tables and observed orders for the continuum-reference 1D-3D coupling runs.

Usage: summarize.py <results dir> <prefix> [<prefix> ...]
       summarize.py <results dir> --scan <prefix>   (runs <prefix>_c<c>)
       summarize.py <results dir> --dt <prefix>     (runs <prefix>_k<k>, dt halving)

Reads <results dir>/<prefix>_N<n>/verification/coupled1D3DMonodomain_local.csv
(final-time errors local to each support), ..._map.csv (discrete map against
the continuum map, every diagnostic kernel) and ..._steps.csv (per-step
junction charge on both sides), and the whole-domain error summaries in
log.cardiacFoam. Orders are log(e_prev/e)/log(N/N_prev) between successive
levels, with h proportional to 1/N.
"""

import csv
import math
import re
import sys
from pathlib import Path


def rows(path):
    with open(path) as f:
        return list(csv.DictReader(f))


def level(name):
    m = re.search(r"_N(\d+)$", name)
    return int(m.group(1)) if m else None


def log_errors(log):
    """Whole-domain L2 of Vm (tissue) and Vm1D (graph) from log.cardiacFoam."""
    out = {}
    if not log.exists():
        return out
    for line in log.read_text().splitlines():
        parts = line.split()
        if len(parts) == 4 and parts[0] in ("Vm", "Vm1D"):
            out[parts[0]] = float(parts[2])
    return out


DIVERGED = 1e3  # an error above this (or not finite) marks an unstable run


def diverged(r):
    for k in ("supportL2", "currentError", "tissueL2"):
        v = abs(float(r.get(k) or "nan"))
        if not math.isfinite(v) or v > DIVERGED:
            return True
    return False


def order(a, b, na, nb):
    try:
        return math.log(abs(a)/abs(b))/math.log(nb/na)
    except (ValueError, ZeroDivisionError):
        return float("nan")


def fit(levels, values, last=4):
    """Least-squares order and its standard error over the finest levels."""
    pts = [(math.log(1.0/n), math.log(abs(v))) for n, v in zip(levels, values) if v][-last:]
    if len(pts) < 3:
        return float("nan"), float("nan")
    k = len(pts)
    mx = sum(x for x, _ in pts)/k
    my = sum(y for _, y in pts)/k
    sxx = sum((x - mx)**2 for x, _ in pts)
    p = sum((x - mx)*(y - my) for x, y in pts)/sxx
    res = sum((y - my - p*(x - mx))**2 for x, y in pts)
    return p, math.sqrt(res/(k - 2)/sxx)


def table(title, header, data):
    print(f"\n### {title}\n")
    print("| " + " | ".join(header) + " |")
    print("|" + "---|"*len(header))
    for r in data:
        print("| " + " | ".join(r) + " |")


def fmt(x):
    return f"{x:.3e}" if isinstance(x, float) else str(x)


def summarize(root, prefix):
    cases = sorted(
        (c for c in root.glob(f"{prefix}_N*") if level(c.name)),
        key=lambda c: level(c.name),
    )
    local = []
    for c in cases:
        f = c/"verification"/"coupled1D3DMonodomain_local.csv"
        if not f.exists():
            continue
        r = [x for x in rows(f) if x["terminal"] == "0"][0]
        if diverged(r):
            print(f"{c.name}: diverged (unstable), left out of the table")
            continue
        steps = rows(c/"verification"/"coupled1D3DMonodomain_steps.csv")
        imb = max(abs(float(s["relativeImbalance"])) for s in steps)
        logs = log_errors(c/"log.cardiacFoam")
        local.append((level(c.name), r, imb, logs))

    keys = [
        ("supportL2", "support L2"),
        ("supportLinf", "support Linf"),
        ("shellL2", "2R L2"),
        ("tissueGatheredError", "gathered V3"),
        ("graphTerminalError", "graph V1(t)"),
        ("currentError", "current"),
        ("farL2", "far-field L2"),
    ]
    data = []
    prev = None
    for n, r, imb, logs in local:
        vals = [abs(float(r.get(k) or "nan")) for k, _ in keys]
        row = [str(n), r["supportCells"]]
        for i, v in enumerate(vals):
            row.append(fmt(v) + (f" ({order(prev[1][i], v, prev[0], n):.2f})" if prev else ""))
        row.append(fmt(logs.get("Vm", float("nan"))))
        row.append(fmt(logs.get("Vm1D", float("nan"))))
        row.append(f"{imb:.1e}")
        data.append(row)
        prev = (n, vals)
    if len(local) >= 3:
        row = ["LS fit", ""]
        for i, (k, _) in enumerate(keys):
            o, e = fit([n for n, *_ in local], [abs(float(r.get(k) or "nan")) for _, r, *_ in local])
            row.append(f"{o:.2f} +- {e:.2f}")
        data.append(row + ["", "", ""])
    if local:
        table(
            f"{prefix}: kernel {local[0][1]['kernel']}, R {local[0][1]['R']}, "
            f"rPvj {local[0][1]['rPvj']}, terminal 0 (order in brackets)",
            ["N", "support cells"] + [k[1] for k in keys]
            + ["tissue L2 (all)", "graph L2 (all)", "max charge imbalance"],
            data,
        )

    # Map convergence, every diagnostic kernel, terminal 0
    maps = {}
    for c in cases:
        f = c/"verification"/"coupled1D3DMonodomain_map.csv"
        if f.exists():
            for r in rows(f):
                if r["terminal"] == "0":
                    maps.setdefault(r["kernel"], []).append((level(c.name), r))
    for kernel, lv in maps.items():
        data = []
        prev = None
        for n, r in lv:
            vals = [
                abs(float(r["normalisationError"])),
                abs(float(r["gatherError"])),
                float(r["depositL1"]),
                float(r["depositL2rel"]),
            ]
            row = [str(n), r["supportCells"]]
            for i, v in enumerate(vals):
                row.append(fmt(v) + (f" ({order(prev[1][i], v, prev[0], n):.2f})" if prev else ""))
            row.append(f"{float(r['meshIntegratedNormalisationError']):.1e}")
            data.append(row)
            prev = (n, vals)
        if len(lv) >= 3:
            row = ["LS fit", ""]
            for key in ("normalisationError", "gatherError", "depositL1", "depositL2rel"):
                o, e = fit([n for n, _ in lv], [abs(float(r[key])) for _, r in lv])
                row.append(f"{o:.2f} +- {e:.2f}")
            data.append(row + [""])
        table(
            f"{prefix}: discrete map vs continuum, kernel {kernel}, terminal 0",
            ["N", "support cells", "Wh/Wc - 1", "<F>h - <F>c",
             "deposit L1", "deposit L2 rel", "sum omega/Wc - 1"],
            data,
        )


def scan(root, prefix):
    """Errors against R/h at a fixed mesh (runs <prefix>_c<c>)."""
    data = []
    for c in sorted(root.glob(f"{prefix}_c*"), key=lambda c: float(c.name.split("_c")[-1])):
        f = c/"verification"/"coupled1D3DMonodomain_local.csv"
        if not f.exists():
            data.append([c.name.split("_c")[-1], "no output (failed)"] + [""]*6)
            continue
        r = [x for x in rows(f) if x["terminal"] == "0"][0]
        if diverged(r):
            data.append([c.name.split("_c")[-1], r["R"], r["supportCells"], "diverged (unstable)"] + [""]*4)
            continue
        data.append([
            c.name.split("_c")[-1], r["R"], r["supportCells"],
            fmt(abs(float(r["supportL2"]))), fmt(abs(float(r["shellL2"]))),
            fmt(abs(float(r["currentError"]))), fmt(abs(float(r["graphTerminalError"]))),
            fmt(abs(float(r["farL2"]))),
        ])
    table(
        f"{prefix}: errors against R/h at fixed mesh, terminal 0",
        ["R/h", "R", "support cells", "support L2", "2R L2", "current",
         "graph V1(t)", "far-field L2"],
        data,
    )


def dtstudy(root, prefix):
    """Fixed-mesh time-step refinement: signed final-time errors and the
    self-convergence order log2(|q_k - q_k+1| / |q_k+1 - q_k+2|)."""
    keys = ["graphTerminalError", "tissueGatheredError", "currentError", "supportL2", "farL2"]
    runs = []
    for c in sorted(root.glob(f"{prefix}_k*"), key=lambda c: int(c.name.split("_k")[-1])):
        f = c/"verification"/"coupled1D3DMonodomain_local.csv"
        if f.exists():
            r = [x for x in rows(f) if x["terminal"] == "0"][0]
            steps = rows(c/"verification"/"coupled1D3DMonodomain_steps.csv")
            runs.append((float(steps[0]["dt"]), [float(r[k]) for k in keys]))
    data = []
    for j, (dt, v) in enumerate(runs):
        row = [f"{dt:.4e}"] + [f"{x:+.4e}" for x in v]
        if j >= 2:
            row += [
                f"{math.log2(abs(runs[j-2][1][i] - runs[j-1][1][i])/abs(runs[j-1][1][i] - v[i])):.2f}"
                if runs[j-1][1][i] != v[i] else "nan"
                for i in range(len(keys))
            ]
        else:
            row += [""]*len(keys)
        data.append(row)
    table(
        f"{prefix}: fixed mesh, dt halving, terminal 0 (signed values; self-convergence orders)",
        ["dt"] + keys + [f"order {k}" for k in keys],
        data,
    )


def main():
    if len(sys.argv) < 3:
        print(__doc__)
        sys.exit(2)
    root = Path(sys.argv[1])
    if sys.argv[2] == "--dt":
        for prefix in sys.argv[3:]:
            dtstudy(root, prefix)
        return
    if sys.argv[2] == "--scan":
        for prefix in sys.argv[3:]:
            scan(root, prefix)
        return
    for prefix in sys.argv[2:]:
        summarize(root, prefix)


if __name__ == "__main__":
    main()
