#!/usr/bin/env python3
"""Interior-terminal junction runs on the sealed shell (tetrahedra): tables and observed orders.

Usage: summarize_interior.py <results dir>

Reads the runs written by the "Terminals inside the wall" loop of README.md:
  in<kernel>_s<seed>_N<n>   unidirectional, terminals at r = 0.75 (lin, cos, uni, gau)
  inoff_s<seed>_N<n>        the same, junction suppressed (R_pvj = 1e6)
  inbi_s<seed>_N<n>         bidirectional
  enduni_s1_N<n>            unidirectional, terminals on the endocardium
  indton_k<k>, indtoff_k<k> N = 14, dt = 0.1/(22 2^k), active and suppressed
Values are means over the seeds present at each level and over the two
terminals; h is the mean cell size written by the shell verifier
(postProcessing/3D_tet_<N>_cells.dat). Orders are least-squares fits of
log(error) on log(h) over all levels, with their standard error, and the
pairwise orders between successive levels.
"""

import csv
import math
import re
import sys
from collections import defaultdict
from pathlib import Path


def local_rows(run):
    with open(run / "verification/coupled1D3DMonodomain_local.csv") as f:
        rows = list(csv.DictReader(f))
    # the file holds one row per terminal per write; keep the last time
    t = max(float(r["t"]) for r in rows)
    last = {}
    for r in rows:
        if float(r["t"]) == t:
            last[int(r["terminal"])] = r
    return [last[k] for k in sorted(last)]


def map_rows(run, kernel):
    with open(run / "verification/coupled1D3DMonodomain_map.csv") as f:
        return [r for r in csv.DictReader(f) if r["kernel"] == kernel]


def mean_h(run):
    for p in (run / "postProcessing").glob("3D_tet_*_cells.dat"):
        m = re.search(r"^# h ([0-9.eE+-]+)", p.read_text(), re.M)
        if m:
            return float(m.group(1))
    return float("nan")


def fit(hs, es):
    pts = [(math.log(h), math.log(abs(e))) for h, e in zip(hs, es) if e and math.isfinite(e)]
    k = len(pts)
    if k < 3:
        return float("nan"), float("nan")
    mx = sum(x for x, _ in pts)/k
    my = sum(y for _, y in pts)/k
    sxx = sum((x - mx)**2 for x, _ in pts)
    p = sum((x - mx)*(y - my) for x, y in pts)/sxx
    res = sum((y - my - p*(x - mx))**2 for x, y in pts)
    return p, math.sqrt(res/(k - 2)/sxx)


QUANT = [("supportL2", "support L2"), ("currentError", "junction current"),
         ("graphTerminalError", "graph V1D at the terminal"), ("farL2", "far field L2"),
         ("tissueL2", "whole-domain L2")]


def ladder(root, prefix, kernel):
    """{N: {quantity: mean, 'h': mean h, 'seeds': n, 'cells': mean support cells, 'deposit': mean L1}}"""
    by = defaultdict(list)
    for run in root.glob(f"{prefix}_s*_N*"):
        m = re.search(r"_s(\d+)_N(\d+)$", run.name)
        if m and (run / "verification/coupled1D3DMonodomain_local.csv").exists():
            by[int(m.group(2))].append(run)
    out = {}
    for N in sorted(by):
        vals = defaultdict(list)
        for run in by[N]:
            rows = local_rows(run)
            for q, _ in QUANT:
                vals[q] += [abs(float(r[q])) for r in rows]
            vals["Ic"] += [abs(float(r["Ic"])) for r in rows]
            vals["cells"] += [float(r["supportCells"]) for r in rows]
            vals["charge"] += [float(r["maxChargeImbalance"]) for r in rows]
            vals["h"].append(mean_h(run))
            vals["deposit"] += [float(r["depositL1"]) for r in map_rows(run, kernel)]
            vals["seedspread"].append(sum(abs(float(r["currentError"])) for r in rows)/len(rows))
        out[N] = {k: sum(v)/len(v) for k, v in vals.items() if v}
        out[N]["seeds"] = len(by[N])
        s = vals["seedspread"]
        out[N]["spread"] = (max(s) - min(s))/(sum(s)/len(s)) if len(s) > 1 else float("nan")
    return out


def print_ladder(title, lad):
    Ns = list(lad)
    hs = [lad[N]["h"] for N in Ns]
    print(f"\n### {title}\n")
    print("| N | meshes | h | support cells | " + " | ".join(n for _, n in QUANT) + " | deposit L1 |")
    print("|" + "---|"*(5 + len(QUANT)))
    for N in Ns:
        L = lad[N]
        print(f"| {N} | {L['seeds']} | {L['h']:.4f} | {L['cells']:.0f} | "
              + " | ".join(f"{L.get(q, float('nan')):.3e}" for q, _ in QUANT)
              + f" | {L.get('deposit', float('nan')):.3e} |")
    row = []
    for q, _ in QUANT + [("deposit", "")]:
        p, s = fit(hs, [lad[N].get(q, float("nan")) for N in Ns])
        pw = [math.log(lad[a][q]/lad[b][q])/math.log(lad[a]["h"]/lad[b]["h"])
              if lad[a].get(q) and lad[b].get(q) else float("nan") for a, b in zip(Ns, Ns[1:])]
        row.append(f"{p:.2f} ± {s:.2f} ({', '.join(f'{x:.2f}' for x in pw)})")
    print("| fit (pairwise) | | | | " + " | ".join(row) + " |")
    print(f"\nexact current |Ic| ~ {lad[Ns[-1]]['Ic']:.4f}; seed spread of the current error at "
          f"N = {Ns[-1]}: {100*lad[Ns[-1]]['spread']:.1f} %; max charge imbalance "
          f"{max(lad[N]['charge'] for N in Ns):.1e}")


def main():
    root = Path(sys.argv[1])
    lads = {}
    for prefix, kernel, title in [("inlin", "linear", "Interior, unidirectional, linear"),
                                  ("incos", "cosine", "Interior, unidirectional, cosine (seed 1)"),
                                  ("inuni", "uniform", "Interior, unidirectional, top-hat (seed 1)"),
                                  ("ingau", "gaussian", "Interior, unidirectional, truncated Gaussian (seed 1)"),
                                  ("inbi", "linear", "Interior, bidirectional, linear"),
                                  ("inoff", "linear", "Interior, unidirectional, junction suppressed"),
                                  ("enduni", "linear", "Endocardial, unidirectional, linear (seed 1)")]:
        lad = ladder(root, prefix, kernel)
        if lad:
            lads[prefix] = lad
            print_ladder(title, lad)

    if "inlin" in lads and "inoff" in lads:
        print("\n### Suppression: error with the junction / error with it suppressed (linear, unidirectional)\n")
        print("| N | within R | far field | whole domain |")
        print("|---|---|---|---|")
        for N in lads["inlin"]:
            if N in lads["inoff"]:
                a, b = lads["inlin"][N], lads["inoff"][N]
                print(f"| {N} | {a['supportL2']/b['supportL2']:.1f} | {a['farL2']/b['farL2']:.2f} | "
                      f"{a['tissueL2']/b['tissueL2']:.2f} |")

    dt = {}
    for tag in ("indton", "indtoff"):
        for run in sorted(root.glob(f"{tag}_k*")):
            k = int(run.name.rsplit("_k", 1)[1])
            rows = local_rows(run)
            dt.setdefault(tag, {})[k] = {q: [float(r[q]) for r in rows] for q in
                                         ("currentError", "tissueGatheredError", "graphTerminalError", "supportL2")}
    for tag, d in dt.items():
        ks = sorted(d)
        print(f"\n### Time step, N = 14, {tag} (self-convergence orders, terminal 0)\n")
        for q in ("currentError", "tissueGatheredError", "supportL2"):
            v = [d[k][q][0] for k in ks]
            # self-convergence: differences of successive solutions are not available from errors alone,
            # so report the orders of the error itself and of the successive error differences
            diffs = [v[i] - v[i + 1] for i in range(len(v) - 1)]
            od = [math.log(abs(diffs[i]/diffs[i + 1]))/math.log(2) for i in range(len(diffs) - 1)
                  if diffs[i + 1]]
            print(f"- {q}: errors {', '.join(f'{x:.4e}' for x in v)}; orders of successive differences "
                  f"{', '.join(f'{x:.2f}' for x in od)}")


if __name__ == "__main__":
    main()
