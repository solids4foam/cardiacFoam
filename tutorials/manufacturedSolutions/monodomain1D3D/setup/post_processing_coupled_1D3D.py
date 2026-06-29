#!/usr/bin/env python3
"""Post-process the joint 1D-3D manufactured solution convergence sweep.

Reads per-resolution output directories produced by run_coupled_1D3D_sweep.sh
and writes unified convergence tables (errors + rates) for:
  - 3D myocardium: Vm, u1, u2
  - 1D graph:      Vm1D, u1, u2
  - Coupling:      sourceL1, sourceL2, totalAbsCurrent
"""

from __future__ import annotations

import argparse
import csv
import math
import re
from pathlib import Path


# ---------------------------------------------------------------------------
# Regex patterns
# ---------------------------------------------------------------------------

# postProcessing/3D_10_cells_implicit.dat
_3D_FILE_RE = re.compile(r"3D_(?P<n>\d+)_cells_(?P<alg>\w+)\.dat$")

# postProcessing/graph_1D_11_nodes.dat
_GRAPH_FILE_RE = re.compile(r"graph_1D_(?P<nodes>\d+)_nodes\.dat$")

# Error lines common to both dat files
_ERROR_LINE_RE = re.compile(
    r"^(?P<field>\S+)\s+"
    r"(?P<l1>[-+0-9.eE]+)\s+"
    r"(?P<l2>[-+0-9.eE]+)\s+"
    r"(?P<linf>[-+0-9.eE]+)\s*$"
)

# ---------------------------------------------------------------------------
# Parsing helpers
# ---------------------------------------------------------------------------

def _parse_dat(path: Path) -> dict[str, tuple[float, float, float]]:
    """Return {field: (L1, L2, Linf)} from a verifier .dat file."""
    result: dict[str, tuple[float, float, float]] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        m = _ERROR_LINE_RE.match(line.strip())
        if m:
            result[m.group("field")] = (
                float(m.group("l1")),
                float(m.group("l2")),
                float(m.group("linf")),
            )
    return result


def _parse_coupling_csv(path: Path) -> dict[str, float]:
    """Return scalar fields from the coupling diagnostics CSV."""
    with path.open(newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            return {k: float(v) for k, v in row.items() if k != "t"}
    return {}


def _safe_rate(e_coarse: float, e_fine: float, h_coarse: float, h_fine: float) -> str:
    if e_coarse <= 0.0 or e_fine <= 0.0:
        return ""
    return f"{math.log(e_coarse / e_fine) / math.log(h_coarse / h_fine):.4f}"


# ---------------------------------------------------------------------------
# Collection
# ---------------------------------------------------------------------------

def collect(output_dir: Path) -> list[dict]:
    rows = []

    for n_dir in sorted(output_dir.iterdir(), key=lambda p: int(p.name) if p.name.isdigit() else 0):
        if not n_dir.is_dir() or not n_dir.name.isdigit():
            continue

        N = int(n_dir.name)
        h = 1.0 / N

        # 3D errors
        mesh_errors: dict[str, tuple[float, float, float]] = {}
        for f in n_dir.glob("3D_*_cells_*.dat"):
            mesh_errors = _parse_dat(f)
            break

        # 1D graph errors
        graph_errors: dict[str, tuple[float, float, float]] = {}
        for f in n_dir.glob("graph_1D_*_nodes.dat"):
            m = _GRAPH_FILE_RE.match(f.name)
            graph_errors = _parse_dat(f)
            n_nodes = int(m.group("nodes")) if m else 0
            break
        else:
            n_nodes = 0

        # Coupling diagnostics
        coupling: dict[str, float] = {}
        coupling_csv = n_dir / "coupling_diagnostics.csv"
        if coupling_csv.exists():
            coupling = _parse_coupling_csv(coupling_csv)

        if not mesh_errors and not graph_errors:
            continue

        row: dict = {"N": N, "h": h, "nodes_1D": n_nodes}

        for field, tag in (("Vm", "3D_Vm"), ("u1", "3D_u1"), ("u2", "3D_u2")):
            if field in mesh_errors:
                l1, l2, linf = mesh_errors[field]
                row[f"L1_{tag}"]   = l1
                row[f"L2_{tag}"]   = l2
                row[f"Linf_{tag}"] = linf

        for field, tag in (("Vm1D", "1D_Vm"), ("u1", "1D_u1"), ("u2", "1D_u2")):
            if field in graph_errors:
                l1, l2, linf = graph_errors[field]
                row[f"L1_{tag}"]   = l1
                row[f"L2_{tag}"]   = l2
                row[f"Linf_{tag}"] = linf

        for key in ("sourceL1", "sourceL2", "totalAbsCurrent"):
            if key in coupling:
                row[f"coupling_{key}"] = coupling[key]

        rows.append(row)

    return sorted(rows, key=lambda r: r["N"])


# ---------------------------------------------------------------------------
# Convergence rates
# ---------------------------------------------------------------------------

_RATE_TARGETS = [
    "Linf_3D_Vm", "Linf_3D_u1", "Linf_3D_u2",
    "Linf_1D_Vm", "Linf_1D_u1", "Linf_1D_u2",
    "coupling_sourceL1", "coupling_totalAbsCurrent",
]


def compute_rates(rows: list[dict]) -> list[dict]:
    rates = []
    for coarse, fine in zip(rows, rows[1:]):
        rate_row: dict = {
            "N_coarse": coarse["N"],
            "N_fine":   fine["N"],
            "h_coarse": coarse["h"],
            "h_fine":   fine["h"],
        }
        for key in _RATE_TARGETS:
            if key in coarse and key in fine:
                rate_row[f"rate_{key}"] = _safe_rate(
                    coarse[key], fine[key], coarse["h"], fine["h"]
                )
        rates.append(rate_row)
    return rates


# ---------------------------------------------------------------------------
# CSV output
# ---------------------------------------------------------------------------

def _write_csv(rows: list[dict], path: Path) -> None:
    if not rows:
        return
    fieldnames = list(rows[0].keys())
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {path}")


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def run_postprocessing(*, output_dir: str, **_) -> list[dict]:
    out = Path(output_dir).resolve()
    rows = collect(out)

    if not rows:
        print(f"No coupled sweep output found under {out}")
        return []

    rates = compute_rates(rows)

    summary_csv = out / "coupled_convergence_summary.csv"
    rates_csv   = out / "coupled_convergence_rates.csv"

    _write_csv(rows,  summary_csv)
    _write_csv(rates, rates_csv)

    # Human-readable table to stdout
    print("\nCoupled 1D-3D convergence summary:")
    print(f"{'N':>4}  {'h':>7}  "
          f"{'Linf_3D_Vm':>12}  {'Linf_1D_Vm':>12}  {'coupling_srcL1':>14}")
    for r in rows:
        print(
            f"{r['N']:>4}  {r['h']:>7.4f}  "
            f"{r.get('Linf_3D_Vm', float('nan')):>12.4e}  "
            f"{r.get('Linf_1D_Vm', float('nan')):>12.4e}  "
            f"{r.get('coupling_sourceL1', float('nan')):>14.4e}"
        )

    if rates:
        print("\nConvergence rates (Linf, consecutive pairs):")
        print(f"{'N_c->N_f':>12}  "
              f"{'rate_3D_Vm':>10}  {'rate_1D_Vm':>10}  {'rate_coupling':>13}")
        for r in rates:
            print(
                f"{r['N_coarse']:>4}->{r['N_fine']:<4}  "
                f"{r.get('rate_Linf_3D_Vm', ''):>10}  "
                f"{r.get('rate_Linf_1D_Vm', ''):>10}  "
                f"{r.get('rate_coupling_sourceL1', ''):>13}"
            )

    return [
        {"path": str(summary_csv), "label": "Coupled 1D-3D convergence summary", "kind": "table", "format": "csv"},
        {"path": str(rates_csv),   "label": "Coupled 1D-3D convergence rates",   "kind": "table", "format": "csv"},
    ]


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path(__file__).resolve().parents[1] / "outputs" / "coupled1D3DConvergence",
    )
    args = parser.parse_args()
    run_postprocessing(output_dir=str(args.output_dir))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
