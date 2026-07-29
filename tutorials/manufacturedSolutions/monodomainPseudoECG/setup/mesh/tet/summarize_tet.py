#----------------------------------------------------------------------------#
# Module
#     summarize_tet
#
# Description
#     Parses manufacturedFDAMonodomainVerifier .dat summaries plus checkMesh
#     logs from the tetrahedral-mesh monodomain MMS sweep into a single
#     order-of-accuracy table. Observed order is computed from the effective
#     mesh spacing dx = 1/cbrt(nCells) reported by the verifier (not an
#     assumed factor-of-two refinement), which is the correct abscissa for an
#     unstructured mesh. Reuses checkmesh_parse for the non-orthogonality
#     metric and adds max skewness, both reported alongside the rate so the
#     mesh quality that the correction has to cope with is explicit.
#----------------------------------------------------------------------------#

from __future__ import annotations

import argparse
import csv
import re
import sys
from dataclasses import dataclass
from pathlib import Path

sys.path.insert(
    0, str(Path(__file__).resolve().parents[4] / "applications/scripts/paperI_results")
)
from checkmesh_parse import parse_checkmesh_log  # noqa: E402
from schema import observed_order  # noqa: E402

_VM_ROW_RE = re.compile(
    r"^Vm\s+([\d.eE+-]+)\s+([\d.eE+-]+)\s+([\d.eE+-]+)", re.MULTILINE
)
_DX_RE = re.compile(r"Grid spacing \(dx\)\s*=\s*([\d.eE+-]+)")
_SKEW_RE = re.compile(r"Max skewness\s*=\s*([\d.eE+-]+)", re.IGNORECASE)


@dataclass
class Row:
    n: int
    dx: float
    l2: float
    linf: float
    max_non_ortho: float | None
    max_skewness: float | None


def parse_dat(path: Path) -> tuple[float, float, float]:
    text = path.read_text()
    vm = _VM_ROW_RE.search(text)
    dx = _DX_RE.search(text)
    if not (vm and dx):
        raise ValueError(f"could not parse dx/Vm norms from {path}")
    # Vm row is L1, L2, Linf
    return float(dx.group(1)), float(vm.group(2)), float(vm.group(3))


def parse_skewness(text: str) -> float | None:
    m = _SKEW_RE.search(text)
    return float(m.group(1)) if m else None


def collect(results_dir: Path, resolutions: list[int]) -> list[Row]:
    rows: list[Row] = []
    for n in resolutions:
        d = results_dir / str(n)
        dx, l2, linf = parse_dat(d / "summary.dat")
        log = (d / "log.checkMesh").read_text()
        cm = parse_checkmesh_log(log)
        rows.append(
            Row(n, dx, l2, linf, cm.max_non_orthogonality, parse_skewness(log))
        )
    # coarse -> fine (largest dx first) so orders read top-to-bottom
    rows.sort(key=lambda r: r.dx, reverse=True)
    return rows


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("results_dir", type=Path)
    ap.add_argument("--resolutions", type=int, nargs="+", required=True)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()

    rows = collect(args.results_dir, args.resolutions)

    header = [
        "N", "dx", "L2_Vm", "p_L2", "Linf_Vm", "p_Linf",
        "maxNonOrtho_deg", "maxSkewness",
    ]
    print("  ".join(f"{h:>15}" for h in header))
    with args.out.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(header)
        for i, r in enumerate(rows):
            if i == 0:
                p2 = pinf = None
            else:
                prev = rows[i - 1]
                p2 = observed_order(prev.l2, r.l2, prev.dx, r.dx)
                pinf = observed_order(prev.linf, r.linf, prev.dx, r.dx)
            cells = [
                r.n,
                f"{r.dx:.6g}",
                f"{r.l2:.6e}",
                "--" if p2 is None else f"{p2:.2f}",
                f"{r.linf:.6e}",
                "--" if pinf is None else f"{pinf:.2f}",
                "--" if r.max_non_ortho is None else f"{r.max_non_ortho:.1f}",
                "--" if r.max_skewness is None else f"{r.max_skewness:.2f}",
            ]
            print("  ".join(f"{str(c):>15}" for c in cells))
            w.writerow(cells)
    print(f"\nwrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
