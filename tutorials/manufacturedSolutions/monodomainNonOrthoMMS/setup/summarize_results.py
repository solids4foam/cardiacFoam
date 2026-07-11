#----------------------------------------------------------------------------#
# Module
#     summarize_results
#
# Description
#     Parses manufacturedFDAMonodomainVerifier .dat summaries plus
#     checkMesh logs from the non-orthogonal distortion sweep into a
#     single order-of-accuracy table.
#----------------------------------------------------------------------------#

from __future__ import annotations

import argparse
import csv
import math
import re
import sys
from dataclasses import dataclass
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from checkmesh_parse import parse_checkmesh_log  # noqa: E402

_VM_ROW_RE = re.compile(
    r"^Vm\s+([\d.eE+-]+)\s+([\d.eE+-]+)\s+([\d.eE+-]+)", re.MULTILINE
)
_DX_RE = re.compile(r"Grid spacing \(dx\)\s*=\s*([\d.eE+-]+)")
_DT_RE = re.compile(r"Time step \(dt\)\s*=\s*([\d.eE+-]+)")


@dataclass
class DatSummary:
    dx: float
    dt: float
    l1_vm: float
    l2_vm: float
    linf_vm: float


def parse_dat_file(path: Path) -> DatSummary:
    text = path.read_text()
    vm_match = _VM_ROW_RE.search(text)
    dx_match = _DX_RE.search(text)
    dt_match = _DT_RE.search(text)
    if not (vm_match and dx_match and dt_match):
        raise ValueError(f"could not parse expected fields from {path}")
    return DatSummary(
        dx=float(dx_match.group(1)),
        dt=float(dt_match.group(1)),
        l1_vm=float(vm_match.group(1)),
        l2_vm=float(vm_match.group(2)),
        linf_vm=float(vm_match.group(3)),
    )


def convergence_order(errors: list[float], refinement_ratio: float = 2.0) -> list[float | None]:
    orders: list[float | None] = [None]
    for coarse, fine in zip(errors[:-1], errors[1:]):
        if coarse <= 0.0 or fine <= 0.0:
            orders.append(None)
            continue
        orders.append(math.log(coarse / fine) / math.log(refinement_ratio))
    return orders


def build_summary_rows(
    results_dir: Path, amplitudes: list[str], resolutions: list[int]
) -> list[dict]:
    rows: list[dict] = []
    for amplitude in amplitudes:
        errors_by_n = []
        for n in resolutions:
            dat_path = results_dir / amplitude / f"3D_{n}_cells_implicit.dat"
            checkmesh_path = results_dir / amplitude / f"log.checkMesh.{n}"
            summary = parse_dat_file(dat_path)
            mesh = parse_checkmesh_log(checkmesh_path.read_text())
            errors_by_n.append(summary.l2_vm)
            rows.append(
                {
                    "A": amplitude,
                    "N": n,
                    "dx": summary.dx,
                    "dt": summary.dt,
                    "L2_Vm": summary.l2_vm,
                    "max_non_orthogonality": mesh.max_non_orthogonality,
                }
            )
        orders = convergence_order(errors_by_n)
        for row, order in zip(rows[-len(resolutions):], orders):
            row["order_Vm"] = order
    return rows


def write_summary_csv(rows: list[dict], out_path: Path) -> None:
    fieldnames = ["A", "N", "dx", "dt", "L2_Vm", "order_Vm", "max_non_orthogonality"]
    with out_path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def build_arg_parser() -> argparse.ArgumentParser:
    """Construct the CLI argument parser.

    NOTE: --amplitudes must stay type=str. The sweep orchestrator
    (run_nonortho_sweep.sh) names result directories with literal shell
    tokens like "0.10"; parsing amplitudes as float would silently mangle
    that into 0.1 (str(float("0.10")) == "0.1"), breaking the directory
    lookup in build_summary_rows.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("results_dir", type=Path)
    parser.add_argument("--amplitudes", type=str, nargs="+", required=True)
    parser.add_argument("--resolutions", type=int, nargs="+", default=[10, 20, 40, 80])
    parser.add_argument("--out", type=Path, default=Path("results/summary.csv"))
    return parser


if __name__ == "__main__":
    args = build_arg_parser().parse_args()

    rows = build_summary_rows(args.results_dir, args.amplitudes, args.resolutions)
    write_summary_csv(rows, args.out)
    print(f"wrote {len(rows)} rows to {args.out}")
