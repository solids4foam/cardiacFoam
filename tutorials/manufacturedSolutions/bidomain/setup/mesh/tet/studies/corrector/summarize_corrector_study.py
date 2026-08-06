#!/usr/bin/env python3
"""Summarise the standalone tetrahedral bidomain corrector study."""

from __future__ import annotations

import csv
import math
from pathlib import Path
import re
import sys


FIELDS = ("Vm", "phiE_gauge", "phiI_gauge")


def metadata(path: Path) -> dict[str, str]:
    return dict(line.split("=", 1) for line in path.read_text().splitlines())


def parse_errors(path: Path) -> dict[str, float]:
    text = path.read_text()
    result: dict[str, float] = {}
    for field in FIELDS:
        match = re.search(
            rf"^{field}\s+[0-9.eE+-]+\s+([0-9.eE+-]+)\s+[0-9.eE+-]+",
            text,
            flags=re.M,
        )
        if not match:
            raise ValueError(f"Could not parse {field} from {path}")
        result[f"{field}_l2"] = float(match.group(1))
    dx = re.search(r"Grid spacing \(dx\)\s+=\s+([0-9.eE+-]+)", text)
    if not dx:
        raise ValueError(f"Could not parse grid spacing from {path}")
    result["dx"] = float(dx.group(1))
    return result


def parse_field(path: Path) -> list[float]:
    text = path.read_text()
    match = re.search(
        r"internalField\s+nonuniform\s+List<scalar>\s+(\d+)\s*\((.*?)\)\s*;",
        text,
        flags=re.S,
    )
    if match:
        expected = int(match.group(1))
        values = [float(token) for token in match.group(2).split()]
        if len(values) != expected:
            raise ValueError(f"{path}: expected {expected} values, found {len(values)}")
        return values

    uniform = re.search(r"internalField\s+uniform\s+([^;]+);", text)
    if uniform:
        return [float(uniform.group(1))]
    raise ValueError(f"Could not parse internalField from {path}")


def field_delta(values: list[float], reference: list[float]) -> float:
    if len(values) != len(reference):
        raise ValueError(f"Field sizes differ: {len(values)} != {len(reference)}")
    differences = [a - b for a, b in zip(values, reference)]
    return math.sqrt(sum(value * value for value in differences)/len(differences))


def solve_counts(path: Path) -> tuple[int, int]:
    text = path.read_text(errors="replace")
    return (
        len(re.findall(r"Solving for phiE\b", text)),
        len(re.findall(r"Solving for Vm\b", text)),
    )


def main() -> None:
    if len(sys.argv) != 2:
        raise SystemExit("usage: summarize_corrector_study.py RESULTS_DIR")
    root = Path(sys.argv[1]).resolve()
    rows: list[dict[str, object]] = []
    fields: dict[tuple[str, str], list[float]] = {}

    for meta_path in sorted(root.glob("*/metadata.env")):
        run_dir = meta_path.parent
        row: dict[str, object] = metadata(meta_path)
        row.update(parse_errors(run_dir / "summary.dat"))
        phi_solves, vm_solves = solve_counts(run_dir / "log.cardiacFoam")
        row["phiE_solves"] = phi_solves
        row["Vm_solves"] = vm_solves
        for name in ("Vm", "phiE"):
            path = run_dir / name
            fields[(str(row["run_id"]), name)] = parse_field(path) if path.exists() else []
        rows.append(row)

    if not rows:
        raise SystemExit(f"No results found below {root}")

    baseline = {
        int(row["resolution"]): row
        for row in rows
        if row["variant"] == "baseline"
    }
    for row in rows:
        ref = baseline[int(row["resolution"])]
        for field in FIELDS:
            key = f"{field}_l2"
            row[f"{field}_change_percent"] = 100.0 * (
                float(row[key])/float(ref[key]) - 1.0
            )

    # Outer-corrector convergence. The MMS norm is not a valid selection
    # criterion: an under-iterated field can score better through accidental
    # cancellation against spatial truncation error. Measure instead how far
    # each count sits from the most-iterated field on the same mesh.
    outer_reference: dict[int, dict[str, object]] = {}
    for row in rows:
        if int(row["n_nonorth"]) != 0:
            continue
        resolution = int(row["resolution"])
        current = outer_reference.get(resolution)
        if current is None or int(row["n_outer"]) > int(current["n_outer"]):
            outer_reference[resolution] = row

    for row in rows:
        ref = outer_reference.get(int(row["resolution"]))
        if ref is None or int(row["n_nonorth"]) != 0:
            continue
        row["outer_reference_run"] = ref["run_id"]
        for name, mms_key in (("Vm", "Vm_l2"), ("phiE", "phiE_gauge_l2")):
            values = fields[(str(row["run_id"]), name)]
            reference = fields[(str(ref["run_id"]), name)]
            if not values or not reference:
                continue
            delta = field_delta(values, reference)
            row[f"{name}_delta_l2_to_outer_reference"] = delta
            row[f"{name}_delta_over_mms_l2"] = delta/float(ref[mms_key])

    fieldnames = sorted({key for row in rows for key in row})
    with (root / "summary.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, restval="")
        writer.writeheader()
        writer.writerows(rows)

    lines = [
        "# Standalone bidomain corrector study",
        "",
        "Changes are relative to one outer sweep and no additional "
        "non-orthogonal assembly on the same mesh.",
        "",
        "| N | variant | equation solves/step | Vm L2 | change | phiE L2 | change | phiI L2 | change |",
        "|---:|---|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in sorted(rows, key=lambda item: (int(item["resolution"]), str(item["variant"]))):
        steps = int(row["n_steps"])
        solves_per_step = (int(row["phiE_solves"]) + int(row["Vm_solves"]))/steps
        lines.append(
            "| {resolution} | {variant} | {solves:.0f} | {vm:.6e} | {vmc:+.2f}% "
            "| {pe:.6e} | {pec:+.2f}% | {pi:.6e} | {pic:+.2f}% |".format(
                resolution=row["resolution"],
                variant=row["variant"],
                solves=solves_per_step,
                vm=float(row["Vm_l2"]),
                vmc=float(row["Vm_change_percent"]),
                pe=float(row["phiE_gauge_l2"]),
                pec=float(row["phiE_gauge_change_percent"]),
                pi=float(row["phiI_gauge_l2"]),
                pic=float(row["phiI_gauge_change_percent"]),
            )
        )

    ladder = [row for row in rows if "Vm_delta_over_mms_l2" in row]
    if ladder:
        lines.extend([
            "",
            "## Outer-corrector convergence",
            "",
            "Field difference from the most-iterated same-mesh field, divided by "
            "that field's MMS L2 error. Predeclared acceptance: <= 1%.",
            "",
            "| N | outer | reference | Vm delta/MMS | phiE delta/MMS |",
            "|---:|---:|---|---:|---:|",
        ])
        for row in sorted(ladder, key=lambda item: (int(item["resolution"]), int(item["n_outer"]))):
            lines.append(
                "| {resolution} | {outer} | {ref} | {vm:.3e} | {pe:.3e} |".format(
                    resolution=row["resolution"],
                    outer=row["n_outer"],
                    ref=row["outer_reference_run"],
                    vm=float(row["Vm_delta_over_mms_l2"]),
                    pe=float(row.get("phiE_delta_over_mms_l2", math.nan)),
                )
            )

        lines.extend(["", "### Smallest passing outer count", ""])
        by_resolution: dict[int, list[dict[str, object]]] = {}
        for row in ladder:
            by_resolution.setdefault(int(row["resolution"]), []).append(row)
        for resolution, group in sorted(by_resolution.items()):
            group.sort(key=lambda item: int(item["n_outer"]))
            passing = [
                row for row in group
                if int(row["n_outer"]) < int(outer_reference[resolution]["n_outer"])
                and float(row["Vm_delta_over_mms_l2"]) <= 0.01
                and float(row.get("phiE_delta_over_mms_l2", math.inf)) <= 0.01
            ]
            verdict = (
                f"outer = {passing[0]['n_outer']}" if passing
                else "none below the reference count"
            )
            lines.append(f"- N={resolution}: smallest outer count passing both 1% criteria is {verdict}.")

    lines.extend(["", "## Two-level observed orders", ""])
    by_variant: dict[str, list[dict[str, object]]] = {}
    for row in rows:
        by_variant.setdefault(str(row["variant"]), []).append(row)
    for variant, group in sorted(by_variant.items()):
        group.sort(key=lambda item: int(item["resolution"]))
        rates = []
        for coarse, fine in zip(group, group[1:]):
            h_ratio = float(coarse["dx"])/float(fine["dx"])
            rate = math.log(float(coarse["Vm_l2"])/float(fine["Vm_l2"]))/math.log(h_ratio)
            rates.append(f"{coarse['resolution']}->{fine['resolution']}: {rate:.2f}")
        lines.append(f"- `{variant}` Vm L2: " + ", ".join(rates))

    (root / "summary.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
