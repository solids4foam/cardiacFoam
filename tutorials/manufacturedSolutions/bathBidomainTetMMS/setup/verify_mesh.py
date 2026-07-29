"""Verify the imported two-zone conformal bath-bidomain mesh."""

from __future__ import annotations

import argparse
import re
from pathlib import Path


def strip_comments(text: str) -> str:
    text = re.sub(r"/\*.*?\*/", "", text, flags=re.DOTALL)
    return re.sub(r"//.*", "", text)


def read_label_list(path: Path) -> list[int]:
    text = strip_comments(path.read_text())
    match = re.search(r"\b\d+\s*\((.*?)\)\s*$", text, flags=re.DOTALL)
    if not match:
        raise ValueError(f"cannot parse label list from {path}")
    return [int(value) for value in re.findall(r"\b\d+\b", match.group(1))]


def read_cell_zones(path: Path) -> dict[str, set[int]]:
    text = strip_comments(path.read_text())
    zones: dict[str, set[int]] = {}
    pattern = re.compile(
        r"([A-Za-z_][A-Za-z0-9_]*)\s*\{[^{}]*?cellLabels\s+List<label>\s+"
        r"\d+\s*\((.*?)\)\s*;?\s*\}",
        flags=re.DOTALL,
    )
    for name, body in pattern.findall(text):
        zones[name] = {int(value) for value in re.findall(r"\b\d+\b", body)}
    return zones


def read_boundary_names(path: Path) -> set[str]:
    text = strip_comments(path.read_text())
    return set(
        re.findall(
            r"\b([A-Za-z_][A-Za-z0-9_]*)\s*\{\s*type\s+[^;]+;",
            text,
            flags=re.DOTALL,
        )
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("case", type=Path)
    args = parser.parse_args()

    mesh = args.case / "constant" / "polyMesh"
    owner = read_label_list(mesh / "owner")
    neighbour = read_label_list(mesh / "neighbour")
    zones = read_cell_zones(mesh / "cellZones")
    patches = read_boundary_names(mesh / "boundary")

    expected_zones = {"myocardium", "bath"}
    if set(zones) != expected_zones:
        raise SystemExit(f"expected zones {expected_zones}, found {set(zones)}")

    myocardium = zones["myocardium"]
    bath = zones["bath"]
    overlap = myocardium & bath
    if overlap:
        raise SystemExit(f"cellZones overlap in {len(overlap)} cells")

    all_zone_cells = myocardium | bath
    n_cells = len(all_zone_cells)
    expected_cells = set(range(max(all_zone_cells) + 1))
    missing = expected_cells - all_zone_cells
    if missing:
        raise SystemExit(f"{len(missing)} cells are outside myocardium/bath zones")
    if len(expected_cells) != n_cells:
        raise SystemExit("cellZone labels are not a contiguous zero-based partition")

    expected_patches = {"xMin", "xMax", "sides"}
    if patches != expected_patches:
        raise SystemExit(f"expected patches {expected_patches}, found {patches}")

    interface_faces = 0
    for face, neighbour_cell in enumerate(neighbour):
        owner_cell = owner[face]
        if (owner_cell in myocardium) != (neighbour_cell in myocardium):
            interface_faces += 1

    if interface_faces == 0:
        raise SystemExit("no internal myocardium--bath interface faces found")

    print(f"cells_total={n_cells}")
    print(f"cells_myocardium={len(myocardium)}")
    print(f"cells_bath={len(bath)}")
    print(f"internal_interface_faces={interface_faces}")
    print("patches=" + ",".join(sorted(patches)))
    print("zone_partition=PASS")
    print("internal_interface=PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
