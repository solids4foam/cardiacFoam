#!/usr/bin/env python3
"""Spatial localisation of the manufactured eikonal activation-time error.

Discriminates two explanations for the observed ~1.5 convergence order:

  H1  boundary-sourced.  The reconstruction stencil degrades at the domain
      edge; the elliptic solve and the outward characteristic sweep then
      spread that error through the interior.  Signature: |e| grows with
      distance travelled from the seed planes, i.e. an accumulation ramp.

  H2  interior truncation.  The bulk gradient is only first order, and
      elliptic regularity lifts the solved field to ~1.5 everywhere.  The
      boundary is irrelevant to the rate.  Signature: no systematic trend
      with distance from the seeds.

Two independent coordinates separate the hypotheses:

  d_seed = min(x, y, z)              small only near the three seed planes
  d_wall = min(x, y, z, L-x, L-y, L-z)   small near ALL six faces

A ramp in d_seed with a flat profile in d_wall is propagation.  Elevation at
small d_wall with no d_seed trend is a local boundary effect.  Flat in both
is interior truncation.

Usage
    postProcess -func writeCellCentres          # writes Cx, Cy, Cz
    ./analyse_error_localisation.py [timeDir] [nbins]

With no timeDir given, resolves the latest numeric time directory under the
current directory itself (same convention as `foamListTimes -latestTime`),
so this can run as a driverFOAM workflow_dag step -- a static command with
no data flow from a preceding step's output -- immediately after a
`postProcess -func writeCellCentres -latestTime` step.
"""

import sys
import math
import os


def _latest_time_dir(case_root="."):
    candidates = []
    for name in os.listdir(case_root):
        try:
            value = float(name)
        except ValueError:
            continue
        candidates.append((value, name))
    if not candidates:
        raise SystemExit(f"no numeric time directories found under {case_root}")
    return max(candidates)[1]


def read_internal_field(path):
    """Parse an OpenFOAM ASCII volScalarField internalField into a list."""
    with open(path) as fh:
        text = fh.read()

    i = text.find("internalField")
    if i < 0:
        raise SystemExit(f"no internalField in {path}")
    seg = text[i:]

    if seg.lstrip()[13:].lstrip().startswith("uniform"):
        val = float(seg.split("uniform", 1)[1].split(";")[0])
        return None, val

    open_paren = seg.find("(")
    close_paren = seg.find(")", open_paren)
    body = seg[open_paren + 1:close_paren]
    return [float(t) for t in body.split()], None


def load(time_dir, name):
    path = os.path.join(time_dir, name)
    if not os.path.exists(path):
        raise SystemExit(
            f"missing {path}\n"
            "Run 'postProcess -func writeCellCentres' and make sure the case "
            "was run with writeErrorField yes."
        )
    vals, uniform = read_internal_field(path)
    if vals is None:
        raise SystemExit(f"{name} is uniform; expected a cellwise field")
    return vals


def stats(pairs, nbins, label):
    """Bin |e| by coordinate and report mean |e| per bin."""
    lo = min(p[0] for p in pairs)
    hi = max(p[0] for p in pairs)
    width = (hi - lo) / nbins if hi > lo else 1.0

    sums = [0.0] * nbins
    counts = [0] * nbins
    for coord, err in pairs:
        b = min(int((coord - lo) / width), nbins - 1)
        sums[b] += abs(err)
        counts[b] += 1

    print(f"\n  mean |e| binned by {label}")
    print(f"  {'bin range':>22} {'n cells':>10} {'mean |e|':>14} {'ratio':>8}")
    first = None
    for b in range(nbins):
        if counts[b] == 0:
            continue
        m = sums[b] / counts[b]
        if first is None:
            first = m
        r = m / first if first else float("nan")
        rng = f"[{lo + b*width:.4f}, {lo + (b+1)*width:.4f})"
        print(f"  {rng:>22} {counts[b]:>10,} {m:>14.6e} {r:>8.2f}")

    # Pearson and Spearman against |e|
    n = len(pairs)
    mx = sum(p[0] for p in pairs) / n
    my = sum(abs(p[1]) for p in pairs) / n
    num = sum((p[0] - mx) * (abs(p[1]) - my) for p in pairs)
    dx = math.sqrt(sum((p[0] - mx) ** 2 for p in pairs))
    dy = math.sqrt(sum((abs(p[1]) - my) ** 2 for p in pairs))
    pear = num / (dx * dy) if dx > 0 and dy > 0 else float("nan")

    order_c = sorted(range(n), key=lambda i: pairs[i][0])
    order_e = sorted(range(n), key=lambda i: abs(pairs[i][1]))
    rc = [0] * n
    re_ = [0] * n
    for rank, i in enumerate(order_c):
        rc[i] = rank
    for rank, i in enumerate(order_e):
        re_[i] = rank
    mr = (n - 1) / 2.0
    num = sum((rc[i] - mr) * (re_[i] - mr) for i in range(n))
    den = sum((rc[i] - mr) ** 2 for i in range(n))
    spear = num / den if den > 0 else float("nan")

    print(f"  Pearson  r(|e|, {label}) = {pear:+.4f}")
    print(f"  Spearman r(|e|, {label}) = {spear:+.4f}")
    return pear, spear


def main():
    time_dir = sys.argv[1] if len(sys.argv) > 1 else _latest_time_dir()
    nbins = int(sys.argv[2]) if len(sys.argv) > 2 else 10

    err = load(time_dir, "activationTimeError")
    cx = load(time_dir, "Cx")
    cy = load(time_dir, "Cy")
    cz = load(time_dir, "Cz")

    if not (len(err) == len(cx) == len(cy) == len(cz)):
        raise SystemExit("field length mismatch")

    lx, ly, lz = max(cx), max(cy), max(cz)
    print(f"cells            : {len(err):,}")
    print(f"domain extent    : x<={lx:.4f} y<={ly:.4f} z<={lz:.4f}")
    print(f"mean |e|         : {sum(abs(e) for e in err)/len(err):.6e}")
    print(f"max  |e|         : {max(abs(e) for e in err):.6e}")

    seed_pairs = []
    wall_pairs = []
    for i, e in enumerate(err):
        x, y, z = cx[i], cy[i], cz[i]
        seed_pairs.append((min(x, y, z), e))
        wall_pairs.append((min(x, y, z, lx - x, ly - y, lz - z), e))

    print("\n" + "=" * 68)
    print("H1 test: does error accumulate along the propagation path?")
    print("=" * 68)
    ps, ss = stats(seed_pairs, nbins, "d_seed = min(x,y,z)")

    print("\n" + "=" * 68)
    print("Boundary-proximity test: is error elevated near ANY wall?")
    print("=" * 68)
    pw, sw = stats(wall_pairs, nbins, "d_wall")

    # d_seed and d_wall are correlated by construction (d_seed >= d_wall is
    # false in general, but both are small near the three min faces), so
    # neither marginal trend is decisive on its own.  The controlled test is
    # whether a d_seed ramp survives among cells held at comparable distance
    # from every wall.  Restrict to interior cells and re-test.
    interior_cut = 0.15 * min(lx, ly, lz)
    interior = [
        (seed_pairs[i][0], err[i])
        for i in range(len(err))
        if wall_pairs[i][0] > interior_cut
    ]

    print("\n" + "=" * 68)
    print("Controlled test: d_seed ramp among interior cells only")
    print(f"(cells with d_wall > {interior_cut:.4f}; "
          f"{len(interior):,} of {len(err):,})")
    print("=" * 68)
    if len(interior) < 100:
        print("  too few interior cells to test")
        si = 0.0
    else:
        _, si = stats(interior, max(4, nbins // 2), "d_seed | interior")

    print("\n" + "=" * 68)
    print("Interpretation")
    print("=" * 68)
    print(f"  Spearman d_seed (all cells)      : {ss:+.4f}")
    print(f"  Spearman d_wall (all cells)      : {sw:+.4f}")
    print(f"  Spearman d_seed (interior only)  : {si:+.4f}   <-- controlled")
    print()

    # Optional machine-readable row, so a ladder over N can be archived rather
    # than living only in this script's stdout. Appends, writing the header
    # only when the file does not yet exist.
    csv_path = os.environ.get("LOCALISATION_CSV")
    if csv_path:
        new = not os.path.exists(csv_path)
        with open(csv_path, "a") as fh:
            if new:
                fh.write(
                    "N,scheme,n_cells,interior_cut,"
                    "spearman_d_seed_all,spearman_d_wall_all,"
                    "spearman_d_seed_interior,mean_abs_err,max_abs_err\n"
                )
            fh.write(
                "{},{},{},{:.6g},{:+.4f},{:+.4f},{:+.4f},{:.6e},{:.6e}\n".format(
                    os.environ.get("LOCALISATION_N", ""),
                    os.environ.get("LOCALISATION_SCHEME", ""),
                    len(err),
                    interior_cut,
                    ss,
                    sw,
                    si,
                    sum(abs(e) for e in err) / len(err),
                    max(abs(e) for e in err),
                )
            )
        print(f"  appended a row to {csv_path}")
    if abs(si) > 0.2:
        print("  A propagation ramp survives with wall proximity held out.")
        print("  -> supports H1: error accumulates along the sweep.")
    elif abs(sw) > 0.2 and abs(si) < 0.1:
        print("  Trend disappears once wall proximity is held out.")
        print("  -> local boundary effect, NOT propagation.")
    elif abs(ss) < 0.1 and abs(sw) < 0.1:
        print("  No systematic spatial trend in any coordinate.")
        print("  -> consistent with H2: interior truncation error.")
    else:
        print("  Mixed / weak signal; neither hypothesis clearly favoured.")

    # The two effects are not exclusive: a strong near-wall elevation and an
    # interior propagation ramp can coexist, and here they do.  Report the
    # near-wall contrast separately rather than forcing a single verdict.
    if abs(sw) > 0.2:
        print()
        print("  A near-wall elevation is also present "
              f"(Spearman {sw:+.4f} against d_wall).")
        print("  Local boundary error and downstream accumulation are")
        print("  independent effects; both can operate at once.")
    print()
    print("  (Rank correlation is the robust statistic here; Pearson assumes")
    print("   a linear relationship the error field need not satisfy.)")


if __name__ == "__main__":
    main()
