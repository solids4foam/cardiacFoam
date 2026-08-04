#!/usr/bin/env python3
"""Per-electrode convergence of the manufactured pseudo-ECG.

Reads the manufacturedPseudoECGSummary.dat files preserved by
run_per_electrode_tet.sh and reports, for each electrode separately, the
Linf error against the high-order reference kernel and its observed order.

The point is to test whether the scattered tetrahedral pseudo-ECG rate
reported for the max-over-electrodes diagnostic is driven by one electrode.
E3 sits 0.200 from the domain against 0.350-0.550 for the others, so its
1/|x - r_e| kernel is the sharpest and its midpoint quadrature error the
largest.  If the scatter is an E3 effect the remaining four should be clean.

Usage
    ./analyse_per_electrode.py <dir-with-summaries> [N ...]
"""

import sys
import os
import math

# Distance from each electrode to the unit cube, for reporting alongside the
# rates.  Kept here rather than recomputed so the script stays standalone.
GEOMETRY = {
    "E1": ((-0.50, 0.50, 0.50), 0.500),
    "E2": ((1.50, 0.50, 0.50), 0.500),
    "E3": ((1.20, 0.23, 0.61), 0.200),
    "E4": ((1.35, 0.74, 0.28), 0.350),
    "E5": ((1.55, 0.41, 0.83), 0.550),
}


# Effective spacing per nominal level, h = nCells^(-1/3).  Orders must be
# computed against the actual spacing ratio, not against the nominal N ratio:
# a gmsh Delaunay mesh at lc=1/N does not contain exactly N^3 cells, and using
# log(N2/N1) instead of log(h1/h2) shifts the reported order by ~0.2 at the
# coarse end.  Override with --h if a different ladder is used.
DEFAULT_H = {10: 0.0588235, 20: 0.030303, 40: 0.0151515, 80: 0.00757576}


def parse_summary(path):
    """Return {electrode: Linf_err_ref} from a verifier summary file."""
    out = {}
    header_seen = False
    with open(path) as fh:
        for line in fh:
            parts = line.split()
            if not parts:
                continue
            if parts[0] == "Electrode":
                header_seen = True
                continue
            if not header_seen:
                continue
            # Electrode  L1_err_ref  L2_err_ref  Linf_err_ref  ...
            if len(parts) >= 4 and parts[0].startswith("E"):
                try:
                    out[parts[0]] = float(parts[3])
                except ValueError:
                    pass
    if not out:
        raise SystemExit(f"no per-electrode rows parsed from {path}")
    return out


def main():
    if len(sys.argv) < 2:
        raise SystemExit(__doc__)
    d = sys.argv[1]
    levels = [int(a) for a in sys.argv[2:]] or [10, 20, 40, 80]

    data = {}
    for N in levels:
        p = os.path.join(d, f"manufacturedPseudoECGSummary_N{N}.dat")
        if os.path.exists(p):
            data[N] = parse_summary(p)
        else:
            print(f"  (missing N={N}: {p})")
    if len(data) < 2:
        raise SystemExit("need at least two levels to compute an order")

    have = sorted(data)
    names = sorted(set().union(*(set(v) for v in data.values())))

    print(f"\nLinf error against the reference kernel, per electrode")
    hdr = f"{'elec':>5} {'dist':>6}" + "".join(f"{('N=%d' % N):>13}" for N in have)
    print(hdr)
    print("-" * len(hdr))
    for e in names:
        dist = GEOMETRY.get(e, (None, float("nan")))[1]
        row = f"{e:>5} {dist:>6.3f}"
        for N in have:
            row += f"{data[N].get(e, float('nan')):>13.4e}"
        print(row)

    print(f"\nObserved order per refinement pair (h halves each level)")
    hdr = f"{'elec':>5} {'dist':>6}" + "".join(
        f"{('%d->%d' % (have[i], have[i+1])):>10}" for i in range(len(have) - 1)
    )
    print(hdr)
    print("-" * len(hdr))
    orders = {}
    for e in names:
        row = f"{e:>5} {GEOMETRY.get(e, (None, float('nan')))[1]:>6.3f}"
        os_ = []
        for i in range(len(have) - 1):
            a, b = data[have[i]].get(e), data[have[i + 1]].get(e)
            if a and b and b > 0:
                p = math.log(a / b) / math.log(DEFAULT_H[have[i]] / DEFAULT_H[have[i + 1]])
                os_.append(p)
                row += f"{p:>10.2f}"
            else:
                os_.append(float("nan"))
                row += f"{'--':>10}"
        orders[e] = os_
        print(row)

    # max-over-electrodes, i.e. what the paper currently reports
    print()
    row = f"{'MAX':>5} {'--':>6}"
    finest = []
    for i in range(len(have) - 1):
        a = max(data[have[i]].values())
        b = max(data[have[i + 1]].values())
        p = math.log(a / b) / math.log(DEFAULT_H[have[i]] / DEFAULT_H[have[i + 1]])
        finest.append(p)
        row += f"{p:>10.2f}"
    print(row + "   <-- the reported diagnostic")

    print("\nWhich electrode sets the maximum at each level")
    for N in have:
        e = max(data[N], key=lambda k: data[N][k])
        print(f"  N={N:<4} {e}  ({data[N][e]:.4e})")

    print("\nInterpretation")
    last = {e: orders[e][-1] for e in names if orders[e] and not math.isnan(orders[e][-1])}
    if last:
        worst = min(last, key=last.get)
        others = [v for k, v in last.items() if k != worst]
        print(f"  finest-pair orders: " + ", ".join(f"{k}={v:.2f}" for k, v in sorted(last.items())))
        if others and last[worst] < min(others) - 0.3:
            print(f"  {worst} is a clear outlier; the other electrodes are tighter.")
            print(f"  -> the reported max-over-electrodes rate is an {worst} effect.")
        else:
            print("  no single electrode is a clear outlier; the scatter is shared.")


if __name__ == "__main__":
    main()
