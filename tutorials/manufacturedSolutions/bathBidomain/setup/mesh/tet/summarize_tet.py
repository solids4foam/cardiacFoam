"""Summarize the bath-bidomain tetrahedral convergence sweep."""

from __future__ import annotations

import argparse
import csv
import math
import re
from dataclasses import dataclass
from pathlib import Path

FIELDS = ("Vm", "phiE", "phiI")
NORMS = ("L1", "L2", "Linf")


@dataclass
class Row:
    nominal_n: int
    total_cells: int
    heart_cells: int
    bath_cells: int
    interface_faces: int
    h_heart: float
    h_global: float
    max_non_ortho: float
    average_non_ortho: float
    max_skewness: float
    min_volume: float
    strict_small_determinant_cells: int
    phi_e_iterations_mean: float
    phi_e_iterations_max: int
    vm_iterations_mean: float
    vm_iterations_max: int
    execution_time_s: float
    errors: dict[str, dict[str, float]]


def key_values(path: Path) -> dict[str, str]:
    values = {}
    for line in path.read_text().splitlines():
        if "=" in line:
            key, value = line.split("=", 1)
            values[key.strip()] = value.strip()
    return values


def required_float(pattern: str, text: str, description: str) -> float:
    match = re.search(pattern, text, flags=re.IGNORECASE)
    if not match:
        raise ValueError(f"could not parse {description}")
    return float(match.group(1).rstrip("."))


def parse_errors(path: Path) -> dict[str, dict[str, float]]:
    errors = {}
    for line in path.read_text().splitlines():
        parts = line.split()
        if parts and parts[0] in FIELDS and len(parts) >= 4:
            errors[parts[0]] = dict(zip(NORMS, map(float, parts[1:4])))
    if set(errors) != set(FIELDS):
        raise ValueError(f"missing field errors in {path}: {set(FIELDS) - set(errors)}")
    return errors


def solver_iterations(text: str, field: str) -> tuple[float, int]:
    values = [
        int(value)
        for value in re.findall(
            rf"Solving for {field},.*?No Iterations\s+(\d+)", text
        )
    ]
    if not values:
        raise ValueError(f"no {field} linear-solver iterations found")
    return sum(values) / len(values), max(values)


def collect(results_dir: Path, resolutions: list[int]) -> list[Row]:
    rows = []
    for nominal_n in resolutions:
        case = results_dir / f"N{nominal_n}"
        manifest = key_values(case / "mesh_manifest.txt")
        checkmesh = (case / "log.checkMesh").read_text()
        solver_log_path = case / "log.cardiacFoam"
        summary_path = case / "summary.dat"
        if not solver_log_path.exists():
            solver_log_path = case / "smoke" / "log.cardiacFoam"
        if not summary_path.exists():
            summaries = sorted(
                (case / "smoke" / "postProcessing").glob(
                    "bathBidomain_3D_*_cells_implicit.dat"
                )
            )
            if not summaries:
                raise FileNotFoundError(f"no manufactured summary under {case}")
            summary_path = summaries[0]
        solver_log = solver_log_path.read_text()

        heart_cells = int(manifest["cells_myocardium"])
        bath_cells = int(manifest["cells_bath"])
        total_cells = heart_cells + bath_cells
        phi_mean, phi_max = solver_iterations(solver_log, "phiE")
        vm_mean, vm_max = solver_iterations(solver_log, "Vm")
        execution_times = [
            float(value)
            for value in re.findall(r"ExecutionTime\s*=\s*([\d.eE+-]+)\s+s", solver_log)
        ]

        rows.append(
            Row(
                nominal_n=nominal_n,
                total_cells=total_cells,
                heart_cells=heart_cells,
                bath_cells=bath_cells,
                interface_faces=int(manifest["internal_interface_faces"]),
                h_heart=(1.0 / heart_cells) ** (1.0 / 3.0),
                h_global=(3.0 / total_cells) ** (1.0 / 3.0),
                max_non_ortho=required_float(
                    r"Mesh non-orthogonality Max:\s*([\d.eE+-]+)",
                    checkmesh,
                    "maximum non-orthogonality",
                ),
                average_non_ortho=required_float(
                    r"Mesh non-orthogonality Max:\s*[\d.eE+-]+\s+average:\s*([\d.eE+-]+)",
                    checkmesh,
                    "average non-orthogonality",
                ),
                max_skewness=required_float(
                    r"Max skewness\s*=\s*([\d.eE+-]+)", checkmesh, "maximum skewness"
                ),
                min_volume=required_float(
                    r"Min volume\s*=\s*([\d.eE+-]+)", checkmesh, "minimum volume"
                ),
                strict_small_determinant_cells=int(
                    manifest.get("strict_small_determinant_cells") or 0
                ),
                phi_e_iterations_mean=phi_mean,
                phi_e_iterations_max=phi_max,
                vm_iterations_mean=vm_mean,
                vm_iterations_max=vm_max,
                execution_time_s=max(execution_times),
                errors=parse_errors(summary_path),
            )
        )
    return sorted(rows, key=lambda row: row.h_heart, reverse=True)


def observed_order(coarse_error: float, fine_error: float, coarse_h: float, fine_h: float) -> float:
    return math.log(coarse_error / fine_error) / math.log(coarse_h / fine_h)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("results_dir", type=Path)
    parser.add_argument("--resolutions", type=int, nargs="+", required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()

    rows = collect(args.results_dir, args.resolutions)
    header = [
        "N_nominal", "cells_total", "cells_myocardium", "cells_bath",
        "interface_faces", "h_heart", "h_global", "maxNonOrtho_deg",
        "avgNonOrtho_deg", "maxSkewness", "minVolume",
        "strictSmallDeterminantCells", "phiE_iters_mean", "phiE_iters_max",
        "Vm_iters_mean", "Vm_iters_max", "executionTime_s",
    ]
    for field in FIELDS:
        for norm in NORMS:
            header.extend((f"{norm}_{field}", f"p_{norm}_{field}"))

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=header)
        writer.writeheader()
        for index, row in enumerate(rows):
            output = {
                "N_nominal": row.nominal_n,
                "cells_total": row.total_cells,
                "cells_myocardium": row.heart_cells,
                "cells_bath": row.bath_cells,
                "interface_faces": row.interface_faces,
                "h_heart": f"{row.h_heart:.9g}",
                "h_global": f"{row.h_global:.9g}",
                "maxNonOrtho_deg": f"{row.max_non_ortho:.6g}",
                "avgNonOrtho_deg": f"{row.average_non_ortho:.6g}",
                "maxSkewness": f"{row.max_skewness:.6g}",
                "minVolume": f"{row.min_volume:.9g}",
                "strictSmallDeterminantCells": row.strict_small_determinant_cells,
                "phiE_iters_mean": f"{row.phi_e_iterations_mean:.3f}",
                "phiE_iters_max": row.phi_e_iterations_max,
                "Vm_iters_mean": f"{row.vm_iterations_mean:.3f}",
                "Vm_iters_max": row.vm_iterations_max,
                "executionTime_s": f"{row.execution_time_s:.6g}",
            }
            for field in FIELDS:
                for norm in NORMS:
                    error = row.errors[field][norm]
                    output[f"{norm}_{field}"] = f"{error:.9e}"
                    if index:
                        coarse = rows[index - 1]
                        rate = observed_order(
                            coarse.errors[field][norm], error, coarse.h_heart, row.h_heart
                        )
                        output[f"p_{norm}_{field}"] = f"{rate:.4f}"
                    else:
                        output[f"p_{norm}_{field}"] = ""
            writer.writerow(output)

    print(f"wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
