#----------------------------------------------------------------------------#
# Module
#     volume_weighted_norm
#
# Description
#     Recomputes the manufactured-monodomain Vm error in both the
#     cell-count-averaged L2 (what the verifier .dat reports) and the
#     volume-weighted L2 (sqrt(sum(e^2 V)/sum(V))), the proper discrete
#     L2(Omega) norm on a non-uniform tet mesh.
#
#     Self-validating: the exact field is reconstructed from the analytic
#     manufactured solution V_ex = sqrt(1+t)*cos(pi x)cos(2 pi y)cos(3 pi z)
#     (Appendix C, 3D), and the recomputed *unweighted* L2 is printed next to
#     the verifier's own .dat L2. If they agree, the exact-field
#     reconstruction is correct and the volume-weighted value is trustworthy.
#
# Usage
#     volume_weighted_norm.py <caseDir> <timeName> [--dat <verifier.dat>]
#----------------------------------------------------------------------------#

from __future__ import annotations

import argparse
import math
import re
from pathlib import Path


def _internal_block(text: str) -> str:
    i = text.index("internalField")
    j = text.index("(", i)
    depth = 0
    for k in range(j, len(text)):
        if text[k] == "(":
            depth += 1
        elif text[k] == ")":
            depth -= 1
            if depth == 0:
                return text[j + 1 : k]
    raise ValueError("unterminated internalField list")


def read_scalar_field(path: Path) -> list[float]:
    body = _internal_block(path.read_text())
    return [float(x) for x in body.split()]


def read_vector_field(path: Path) -> list[tuple[float, float, float]]:
    body = _internal_block(path.read_text())
    trip = re.findall(
        r"\(\s*([-\d.eE+]+)\s+([-\d.eE+]+)\s+([-\d.eE+]+)\s*\)", body
    )
    return [(float(a), float(b), float(c)) for a, b, c in trip]


def dat_l2(path: Path) -> float | None:
    m = re.search(r"^Vm\s+[\d.eE+-]+\s+([\d.eE+-]+)", path.read_text(), re.MULTILINE)
    return float(m.group(1)) if m else None


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("case", type=Path)
    ap.add_argument("time", type=str)
    ap.add_argument("--dat", type=Path, default=None)
    args = ap.parse_args()

    t = float(args.time)
    tdir = args.case / args.time
    Vm = read_scalar_field(tdir / "Vm")
    V = read_scalar_field(tdir / "V")
    C = read_vector_field(tdir / "C")
    n = len(Vm)
    assert len(V) == n and len(C) == n, (n, len(V), len(C))

    amp = math.sqrt(1.0 + t)
    se2 = 0.0        # sum e^2
    se2v = 0.0       # sum e^2 * V
    sv = 0.0         # sum V
    linf = 0.0
    for i in range(n):
        x, y, z = C[i]
        vex = amp * math.cos(math.pi * x) * math.cos(2 * math.pi * y) * math.cos(3 * math.pi * z)
        e = Vm[i] - vex
        e2 = e * e
        se2 += e2
        se2v += e2 * V[i]
        sv += V[i]
        linf = max(linf, abs(e))

    l2_unw = math.sqrt(se2 / n)
    l2_vw = math.sqrt(se2v / sv)

    ref = dat_l2(args.dat) if args.dat else None
    check = ""
    if ref is not None:
        rel = abs(l2_unw - ref) / ref
        check = f"  (verifier .dat L2={ref:.6e}, rel.diff {rel:.1%} {'OK' if rel < 0.02 else 'MISMATCH!'})"

    print(f"  t={t:.6g}  nCells={n}")
    print(f"  L2 unweighted   = {l2_unw:.6e}{check}")
    print(f"  L2 volume-weight= {l2_vw:.6e}")
    print(f"  Linf            = {linf:.6e}")
    # machine-readable last line: t nCells l2_unw l2_vw linf
    print(f"CSV {t:.6g} {n} {l2_unw:.6e} {l2_vw:.6e} {linf:.6e}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
