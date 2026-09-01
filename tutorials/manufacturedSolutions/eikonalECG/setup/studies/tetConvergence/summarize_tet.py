#----------------------------------------------------------------------------#
# Module
#     summarize_tet
#
# Description
#     Parses manufacturedEikonalVerifier .dat summaries plus checkMesh logs
#     from the tetrahedral-mesh eikonal MMS sweep into a single
#     order-of-accuracy table. Unlike the monodomain verifier, the eikonal
#     verifier's .dat has no "Grid spacing (dx)" line -- it only reports
#     "Number of cells = N" -- so the unit-cube effective spacing is computed
#     as h = nCells^(-1/3), without rounding the cube root. It reports the
#     realised pairwise h ratio alongside every observed order, rather than
#     assuming the nominal Gmsh lc ladder is exactly a factor of two. It also
#     retains the requested N separately from the realised cell count.
#----------------------------------------------------------------------------#

from __future__ import annotations

import argparse
import csv
import math
import re
from dataclasses import dataclass
from pathlib import Path

_NCELLS_RE = re.compile(r"Number of cells\s*=\s*(\d+)")
_ACTIVATIONTIME_ROW_RE = re.compile(
    r"^activationTime\s+([\d.eE+-]+)\s+([\d.eE+-]+)\s+([\d.eE+-]+)", re.MULTILINE
)
_NON_ORTHO_RE = re.compile(
    r"Mesh non-orthogonality\s+Max:\s*([\d.eE+-]+)", re.IGNORECASE
)
_SKEW_RE = re.compile(r"Max skewness\s*=\s*([\d.eE+-]+)", re.IGNORECASE)


@dataclass
class Row:
    requested_n: int
    n_cells: int
    h_eff: float
    l2: float
    linf: float
    ecg_l2: float | None
    ecg_linf: float | None
    max_non_ortho: float | None
    max_skewness: float | None


def observed_order(
    error_coarse: float, error_fine: float, h_coarse: float, h_fine: float
) -> float:
    return math.log(error_coarse / error_fine) / math.log(h_coarse / h_fine)


def parse_dat(path: Path) -> tuple[int, float, float, float]:
    text = path.read_text()
    ncells = _NCELLS_RE.search(text)
    row = _ACTIVATIONTIME_ROW_RE.search(text)
    if not (ncells and row):
        raise ValueError(f"could not parse cell count/activationTime norms from {path}")
    n_cells = int(ncells.group(1))
    h_eff = n_cells ** (-1.0 / 3.0)
    # activationTime row is L1, L2, Linf
    return n_cells, h_eff, float(row.group(2)), float(row.group(3))


def parse_ecg_dat(path: Path) -> tuple[float, float] | None:
    """Manufactured pseudo-ECG summary: one header block (title, samples,
    dimension, qChecks, qReference, k, column-name row -- 7 lines) then one
    row per electrode with columns Electrode L1_err_ref L2_err_ref
    Linf_err_ref ... Reports the worst (max) L2/Linf error over electrodes,
    same convention as run_eikonal_tet.sh's ecg_metrics().
    """
    if not path.exists():
        return None
    lines = path.read_text().splitlines()
    header_idx = next(
        (i for i, ln in enumerate(lines) if ln.startswith("Electrode")), None
    )
    if header_idx is None:
        return None
    l2 = linf = 0.0
    for ln in lines[header_idx + 1 :]:
        parts = ln.split()
        if len(parts) < 4:
            continue
        l2 = max(l2, float(parts[2]))
        linf = max(linf, float(parts[3]))
    return l2, linf


def parse_skewness(text: str) -> float | None:
    m = _SKEW_RE.search(text)
    return float(m.group(1)) if m else None


def parse_max_non_orthogonality(text: str) -> float | None:
    m = _NON_ORTHO_RE.search(text)
    return float(m.group(1)) if m else None


def collect(results_dir: Path, resolutions: list[int]) -> list[Row]:
    rows: list[Row] = []
    for requested_n in resolutions:
        d = results_dir / str(requested_n)
        n_cells, h_eff, l2, linf = parse_dat(d / "summary.dat")
        ecg = parse_ecg_dat(d / "pseudoECG_summary.dat")
        log = (d / "log.checkMesh").read_text()
        rows.append(
            Row(
                requested_n, n_cells, h_eff, l2, linf,
                ecg[0] if ecg else None, ecg[1] if ecg else None,
                parse_max_non_orthogonality(log), parse_skewness(log),
            )
        )
    # Coarse -> fine (largest realised h first) so orders read top-to-bottom.
    rows.sort(key=lambda r: r.h_eff, reverse=True)
    return rows




def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("results_dir", type=Path)
    ap.add_argument("--resolutions", type=int, nargs="+", required=True)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()

    rows = collect(args.results_dir, args.resolutions)

    header = [
        "N_requested", "nCells", "h_eff", "h_ratio", "L2_activationTime", "p_L2", "Linf_activationTime", "p_Linf",
        "ECG_L2", "p_ECG_L2", "ECG_Linf", "p_ECG_Linf",
        "maxNonOrtho_deg", "maxSkewness",
    ]
    print("  ".join(f"{h:>19}" for h in header))
    with args.out.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(header)
        for i, r in enumerate(rows):
            if i == 0:
                h_ratio = p2 = pinf = pecg2 = pecginf = None
            else:
                prev = rows[i - 1]
                h_ratio = prev.h_eff / r.h_eff
                p2 = observed_order(prev.l2, r.l2, prev.h_eff, r.h_eff)
                pinf = observed_order(prev.linf, r.linf, prev.h_eff, r.h_eff)
                pecg2 = (
                    observed_order(prev.ecg_l2, r.ecg_l2, prev.h_eff, r.h_eff)
                    if prev.ecg_l2 is not None and r.ecg_l2 is not None
                    else None
                )
                pecginf = (
                    observed_order(prev.ecg_linf, r.ecg_linf, prev.h_eff, r.h_eff)
                    if prev.ecg_linf is not None and r.ecg_linf is not None
                    else None
                )
            cells = [
                r.requested_n,
                r.n_cells,
                f"{r.h_eff:.9g}",
                "--" if h_ratio is None else f"{h_ratio:.6g}",
                f"{r.l2:.6e}",
                "--" if p2 is None else f"{p2:.2f}",
                f"{r.linf:.6e}",
                "--" if pinf is None else f"{pinf:.2f}",
                "--" if r.ecg_l2 is None else f"{r.ecg_l2:.6e}",
                "--" if pecg2 is None else f"{pecg2:.2f}",
                "--" if r.ecg_linf is None else f"{r.ecg_linf:.6e}",
                "--" if pecginf is None else f"{pecginf:.2f}",
                "--" if r.max_non_ortho is None else f"{r.max_non_ortho:.1f}",
                "--" if r.max_skewness is None else f"{r.max_skewness:.2f}",
            ]
            print("  ".join(f"{str(c):>19}" for c in cells))
            w.writerow(cells)
    print(f"\nwrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
