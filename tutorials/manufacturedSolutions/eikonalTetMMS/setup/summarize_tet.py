#----------------------------------------------------------------------------#
# Module
#     summarize_tet
#
# Description
#     Parses manufacturedEikonalVerifier .dat summaries plus checkMesh logs
#     from the tetrahedral-mesh eikonal MMS sweep into a single
#     order-of-accuracy table. Unlike the monodomain verifier, the eikonal
#     verifier's .dat has no "Grid spacing (dx)" line -- it only reports
#     "Number of cells = N" -- so the effective spacing
#     dx = 1/round(cbrt(nCells)) is computed here directly, same formula the
#     monodomain summarizer uses for its unstructured mesh. Reuses
#     checkmesh_parse for the non-orthogonality metric and adds max
#     skewness, both reported alongside the rate so the mesh quality that
#     the gradient scheme has to cope with is explicit.
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

_NCELLS_RE = re.compile(r"Number of cells\s*=\s*(\d+)")
_ACTIVATIONTIME_ROW_RE = re.compile(
    r"^activationTime\s+([\d.eE+-]+)\s+([\d.eE+-]+)\s+([\d.eE+-]+)", re.MULTILINE
)
_SKEW_RE = re.compile(r"Max skewness\s*=\s*([\d.eE+-]+)", re.IGNORECASE)


@dataclass
class Row:
    n: int
    dx: float
    l2: float
    linf: float
    ecg_l2: float | None
    ecg_linf: float | None
    max_non_ortho: float | None
    max_skewness: float | None


def parse_dat(path: Path) -> tuple[float, float, float]:
    text = path.read_text()
    ncells = _NCELLS_RE.search(text)
    row = _ACTIVATIONTIME_ROW_RE.search(text)
    if not (ncells and row):
        raise ValueError(f"could not parse cell count/activationTime norms from {path}")
    n_cells = int(ncells.group(1))
    dx = 1.0 / round(n_cells ** (1.0 / 3.0))
    # activationTime row is L1, L2, Linf
    return dx, float(row.group(2)), float(row.group(3))


def parse_ecg_dat(path: Path) -> tuple[float, float] | None:
    """Manufactured pseudo-ECG summary: one header block (title, samples,
    dimension, qChecks, qReference, k, column-name row -- 7 lines) then one
    row per electrode with columns Electrode L1_err_ref L2_err_ref
    Linf_err_ref ... Reports the worst (max) L2/Linf error over electrodes,
    same convention as run_scheme_study.sh's ecg_metrics().
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


def collect(results_dir: Path, resolutions: list[int]) -> list[Row]:
    rows: list[Row] = []
    for n in resolutions:
        d = results_dir / str(n)
        dx, l2, linf = parse_dat(d / "summary.dat")
        ecg = parse_ecg_dat(d / "pseudoECG_summary.dat")
        log = (d / "log.checkMesh").read_text()
        cm = parse_checkmesh_log(log)
        rows.append(
            Row(
                n, dx, l2, linf,
                ecg[0] if ecg else None, ecg[1] if ecg else None,
                cm.max_non_orthogonality, parse_skewness(log),
            )
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
        "N", "dx", "L2_activationTime", "p_L2", "Linf_activationTime", "p_Linf",
        "ECG_L2", "p_ECG_L2", "ECG_Linf", "p_ECG_Linf",
        "maxNonOrtho_deg", "maxSkewness",
    ]
    print("  ".join(f"{h:>19}" for h in header))
    with args.out.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(header)
        for i, r in enumerate(rows):
            if i == 0:
                p2 = pinf = pecg2 = pecginf = None
            else:
                prev = rows[i - 1]
                p2 = observed_order(prev.l2, r.l2, prev.dx, r.dx)
                pinf = observed_order(prev.linf, r.linf, prev.dx, r.dx)
                pecg2 = (
                    observed_order(prev.ecg_l2, r.ecg_l2, prev.dx, r.dx)
                    if prev.ecg_l2 is not None and r.ecg_l2 is not None
                    else None
                )
                pecginf = (
                    observed_order(prev.ecg_linf, r.ecg_linf, prev.dx, r.dx)
                    if prev.ecg_linf is not None and r.ecg_linf is not None
                    else None
                )
            cells = [
                r.n,
                f"{r.dx:.6g}",
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
