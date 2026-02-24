from __future__ import annotations

import argparse
import json
from pathlib import Path

from .mesh_context import build_segment_mesh_context, save_mesh_ascii
from .vtk_bullseye_workbench import map_percentages_to_segments


def _parse_segment_values(raw: object, label: str) -> dict[int, float]:
    if raw is None:
        return {}
    if not isinstance(raw, dict):
        raise ValueError(f"'{label}' must be a JSON object mapping segment IDs to values.")

    out: dict[int, float] = {}
    for key, value in raw.items():
        try:
            seg_id = int(key)
        except Exception as exc:  # pragma: no cover - defensive
            raise ValueError(f"Invalid segment id in '{label}': {key!r}") from exc
        try:
            out[seg_id] = float(value)
        except Exception as exc:  # pragma: no cover - defensive
            raise ValueError(f"Invalid value for segment {key!r} in '{label}': {value!r}") from exc
    return out


def load_solution_json(path: Path) -> tuple[dict[int, float], dict[int, float], dict[int, float]]:
    with open(path, "r", encoding="utf-8") as f:
        payload = json.load(f)
    if not isinstance(payload, dict):
        raise ValueError("Solution JSON must be a top-level object.")
    purkinje = _parse_segment_values(payload.get("purkinje"), "purkinje")
    pmj = _parse_segment_values(payload.get("pmj"), "pmj")
    thickness_bundles = _parse_segment_values(
        payload.get("thickness_bundles"),
        "thickness_bundles",
    )
    return purkinje, pmj, thickness_bundles


def build_arg_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        description="Apply a saved unified-workbench solution JSON to a tagged VTK mesh."
    )
    ap.add_argument("--vtk", type=Path, required=True, help="Input tagged VTK/VTU/VTP mesh.")
    ap.add_argument(
        "--solution-json",
        type=Path,
        required=True,
        help="Path to JSON produced by unified_workbench_gui (typically *.solution.json).",
    )
    ap.add_argument("--output", type=Path, required=True, help="Output mesh path.")
    ap.add_argument(
        "--segment-field",
        default=None,
        help="Anatomical segment field. If omitted, auto-detected.",
    )
    ap.add_argument(
        "--segment-association",
        choices=["auto", "point", "cell"],
        default="auto",
        help="Whether segment IDs are in point_data or cell_data.",
    )
    ap.add_argument(
        "--purkinje-field",
        default="purkinje_density",
        help="Output scalar field name for Purkinje density.",
    )
    ap.add_argument(
        "--pmj-field",
        default="pmj_per_mm2",
        help="Output scalar field name for PMJ /mm2.",
    )
    ap.add_argument(
        "--thickness-bundles-field",
        default="thickness_bundles",
        help="Output scalar field name for Thickness Bundles.",
    )
    ap.add_argument(
        "--endocardial-area-field",
        default="endocardial_surface_area_by_tag",
        help="Output scalar field name for automatic endocardial area-by-tag values.",
    )
    ap.add_argument(
        "--endocardial-selector-field",
        default=None,
        help=(
            "Optional point/cell field to select endocardial cells for area "
            "computation. If omitted, auto-detects endo_surface tags and "
            "falls back to all cells."
        ),
    )
    ap.add_argument(
        "--endocardial-selector-value",
        default="1",
        help="Selector value used with --endocardial-selector-field (default: 1).",
    )
    return ap


def main() -> None:
    args = build_arg_parser().parse_args()

    if not args.vtk.exists():
        raise FileNotFoundError(f"Input VTK not found: {args.vtk}")
    if not args.solution_json.exists():
        raise FileNotFoundError(f"Solution JSON not found: {args.solution_json}")

    try:
        import numpy as np
        import pyvista as _pv  # noqa: F401
    except Exception as exc:
        raise RuntimeError(
            "Applying solution JSON requires numpy + pyvista."
        ) from exc

    purkinje_map, pmj_map, thickness_bundles_map = load_solution_json(args.solution_json)
    ctx = build_segment_mesh_context(
        args.vtk,
        segment_field=args.segment_field,
        segment_association=args.segment_association,
        endocardial_selector_field=args.endocardial_selector_field,
        endocardial_selector_value=args.endocardial_selector_value,
    )
    mesh = ctx.mesh
    selection = ctx.selection
    segment_ids = np.asarray(ctx.segment_ids).reshape(-1).astype(int)

    purkinje_values = map_percentages_to_segments(segment_ids, purkinje_map, default_value=0.0)
    pmj_values = map_percentages_to_segments(segment_ids, pmj_map, default_value=0.0)
    thickness_bundles_values = map_percentages_to_segments(
        segment_ids,
        thickness_bundles_map,
        default_value=0.0,
    )
    area_result = ctx.area_result

    out = mesh.copy(deep=True)
    if selection.association == "point":
        out.point_data[str(args.purkinje_field)] = purkinje_values
        out.point_data[str(args.pmj_field)] = pmj_values
        out.point_data[str(args.thickness_bundles_field)] = thickness_bundles_values
        out.point_data[str(args.endocardial_area_field)] = area_result.mapped_values
    else:
        out.cell_data[str(args.purkinje_field)] = purkinje_values
        out.cell_data[str(args.pmj_field)] = pmj_values
        out.cell_data[str(args.thickness_bundles_field)] = thickness_bundles_values
        out.cell_data[str(args.endocardial_area_field)] = area_result.mapped_values

    save_mesh_ascii(out, args.output)

    print(f"Read input mesh: {args.vtk}")
    print(f"Used segment field: {selection.field_name} ({selection.association})")
    print(f"Loaded solution JSON: {args.solution_json}")
    print(f"Wrote: {args.output}")
    print(
        "Added fields: "
        f"{args.purkinje_field} (Purkinje), "
        f"{args.pmj_field} (PMJ /mm2), "
        f"{args.thickness_bundles_field} (Thickness Bundles), "
        f"{args.endocardial_area_field} (Endocardial area by tag)"
    )
    print("Endocardial surface area by tag:")
    for seg_id in sorted(area_result.area_by_segment):
        print(f"  {seg_id}: {area_result.area_by_segment[seg_id]:.6f}")
    print(f"Endocardial selector source: {area_result.selector_source}")


if __name__ == "__main__":
    main()
