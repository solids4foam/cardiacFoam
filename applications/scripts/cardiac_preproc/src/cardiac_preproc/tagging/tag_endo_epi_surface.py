"""Tag endocardial/epicardial surfaces from UVC fields.

Inputs:
- Legacy ASCII VTK (volume or surface) with POINT_DATA:
  - uvc_transmural
  - uvc_intraventricular

Outputs:
- outputs/<input>_endo_epi_surface.vtk (surface tags)
- outputs/<input>_endo_epi_volume.vtk (volume-mapped tags)
- outputs/<input>_complementary_endocardium_epicardium_label.vtk (extra labels)

Main:
- Tags LV/RV endo, epi, shared boundary, and overlap diagnostics.
"""
import argparse
import sys
from pathlib import Path

import numpy as np
import pyvista as pv

PROJECT_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(PROJECT_ROOT / "src"))
sys.path.insert(0, str(PROJECT_ROOT / "configs"))

from cardiac_preproc.utils.vtk_utils import inspect_fields, remove_blank_lines  # noqa: E402
from cardiac_preproc.utils.vtk_convert_arrays_to_fields import convert_vtk_file  # noqa: E402

try:
    from tag_endo_epi_surface_config import (  # type: ignore
        COMPLEMENTARY_SURFACE_OUTPUT,
        EXTRACT_VOLUME_OUTPUT,
        SHARED_BOUNDARY_SEED,
    )
except ImportError:
    COMPLEMENTARY_SURFACE_OUTPUT = True
    EXTRACT_VOLUME_OUTPUT = True
    SHARED_BOUNDARY_SEED = None


def surface_point_ids(surface: pv.DataSet, mesh: pv.DataSet) -> np.ndarray:
    if "vtkOriginalPointIds" in surface.point_data:
        return surface.point_data["vtkOriginalPointIds"].astype(int)

    ids = np.zeros(surface.n_points, dtype=int)
    for i, p in enumerate(surface.points):
        ids[i] = mesh.find_closest_point(p)
    return ids


def map_surface_tags_to_volume(
    mesh: pv.DataSet,
    surface: pv.DataSet,
    tag_names: list[str],
) -> pv.DataSet:
    point_ids = surface_point_ids(surface, mesh)

    for tag in tag_names:
        if tag not in surface.cell_data:
            continue
        surface_point_mask = np.zeros(mesh.n_points, dtype=bool)
        tagged_cells = surface.cell_data[tag] > 0
        for cell_id in np.where(tagged_cells)[0]:
            pts = surface.get_cell(cell_id).point_ids
            surface_point_mask[point_ids[pts]] = True

        cell_tags = np.zeros(mesh.n_cells, dtype=int)
        for i in range(mesh.n_cells):
            pts = mesh.get_cell(i).point_ids
            if np.any(surface_point_mask[pts]):
                cell_tags[i] = 1
        mesh.cell_data[tag] = cell_tags

    return mesh


def find_unused_point_ids(mesh: pv.DataSet) -> np.ndarray:
    used = np.zeros(mesh.n_points, dtype=bool)
    for cell_id in range(mesh.n_cells):
        used[mesh.get_cell(cell_id).point_ids] = True
    return np.where(~used)[0]


def report_unused_points(
    mesh: pv.DataSet,
    label: str,
    output_path: Path | None,
) -> None:
    unused_ids = find_unused_point_ids(mesh)
    if unused_ids.size == 0:
        print(f"{label}: no unused points detected.")
        return

    print(f"{label}: found {unused_ids.size} unused points.")
    print(f"{label}: unused point coordinates:\n{mesh.points[unused_ids]}")
    if output_path is not None:
        points = mesh.points[unused_ids]
        pv.PolyData(points).save(str(output_path), binary=False)
        print(f"{label}: unused points written to {output_path.resolve()}")
def tag_surface(
    mesh: pv.DataSet,
    transmural_field: str,
    intraventricular_field: str,
    transmural_min: float,
    transmural_max: float,
    transmural_eq: float,
    transmural_eps: float,
    ventricle: str,
    lv_value: int,
    rv_value: int,
    mode: str,
    tag_name: str,
    extract_surface: bool = True,
) -> pv.DataSet:
    surface = mesh.extract_surface() if extract_surface else mesh

    if transmural_field not in surface.point_data:
        raise RuntimeError(
            f"Missing point_data '{transmural_field}' on surface. "
            "Check field names with inspect_fields."
        )
    if intraventricular_field not in surface.point_data:
        raise RuntimeError(
            f"Missing point_data '{intraventricular_field}' on surface. "
            "Check field names with inspect_fields."
        )

    uvc_trans = surface.point_data[transmural_field]
    uvc_intra = surface.point_data[intraventricular_field]

    if transmural_eq is not None:
        is_endo = np.abs(uvc_trans - transmural_eq) <= transmural_eps
    else:
        is_endo = (uvc_trans >= transmural_min) & (uvc_trans <= transmural_max)
    if ventricle == "lv":
        is_endo &= uvc_intra == lv_value
    elif ventricle == "rv":
        is_endo &= uvc_intra == rv_value
    else:
        is_endo &= (uvc_intra == lv_value) | (uvc_intra == rv_value)

    cell_tags = np.zeros(surface.n_cells, dtype=int)
    for i in range(surface.n_cells):
        point_ids = surface.get_cell(i).point_ids
        if mode == "all":
            tagged = np.all(is_endo[point_ids])
        else:
            tagged = np.any(is_endo[point_ids])
        cell_tags[i] = 1 if tagged else 0

    surface.cell_data[tag_name] = cell_tags
    return surface


def tag_epi_surface(
    surface: pv.DataSet,
    transmural_field: str,
    transmural_eq: float,
    transmural_eps: float,
    mode: str,
    tag_name: str,
) -> pv.DataSet:
    if transmural_field not in surface.point_data:
        raise RuntimeError(
            f"Missing point_data '{transmural_field}' on surface. "
            "Check field names with inspect_fields."
        )

    uvc_trans = surface.point_data[transmural_field]
    is_epi = np.abs(uvc_trans - transmural_eq) <= transmural_eps

    cell_tags = np.zeros(surface.n_cells, dtype=int)
    for i in range(surface.n_cells):
        point_ids = surface.get_cell(i).point_ids
        if mode == "all":
            tagged = np.all(is_epi[point_ids])
        else:
            tagged = np.any(is_epi[point_ids])
        cell_tags[i] = 1 if tagged else 0

    surface.cell_data[tag_name] = cell_tags
    return surface


def tag_shared_boundary(
    surface: pv.DataSet,
    transmural_field: str,
    intraventricular_field: str,
    rv_transmural_eq: float,
    rv_transmural_eps: float,
    lv_transmural_eq: float,
    lv_transmural_eps: float,
    lv_value: int,
    rv_value: int,
    tag_name: str,
) -> pv.DataSet:
    if transmural_field not in surface.point_data:
        raise RuntimeError(
            f"Missing point_data '{transmural_field}' on surface. "
            "Check field names with inspect_fields."
        )
    if intraventricular_field not in surface.point_data:
        raise RuntimeError(
            f"Missing point_data '{intraventricular_field}' on surface. "
            "Check field names with inspect_fields."
        )

    uvc_trans = surface.point_data[transmural_field]
    uvc_intra = surface.point_data[intraventricular_field]

    is_rv_endo = (uvc_intra == rv_value) & (
        np.abs(uvc_trans - rv_transmural_eq) <= rv_transmural_eps
    )
    is_lv_epi = (uvc_intra == lv_value) & (
        np.abs(uvc_trans - lv_transmural_eq) <= lv_transmural_eps
    )

    cell_tags = np.zeros(surface.n_cells, dtype=int)
    for i in range(surface.n_cells):
        point_ids = surface.get_cell(i).point_ids
        has_rv_endo = np.any(is_rv_endo[point_ids])
        has_lv_epi = np.any(is_lv_epi[point_ids])
        cell_tags[i] = 1 if (has_rv_endo and has_lv_epi) else 0

    surface.cell_data[tag_name] = cell_tags
    return surface


def parse_vector(value: str) -> np.ndarray:
    parts = [p.strip() for p in value.split(",") if p.strip()]
    if len(parts) != 3:
        raise ValueError(f"Expected 3 comma-separated values, got: '{value}'")
    return np.array([float(p) for p in parts])


def cell_adjacency(surface: pv.DataSet):
    edge_to_cells = {}
    for i in range(surface.n_cells):
        pts = surface.get_cell(i).point_ids
        n = len(pts)
        for a in range(n):
            e = tuple(sorted((pts[a], pts[(a + 1) % n])))
            edge_to_cells.setdefault(e, []).append(i)

    adj = [[] for _ in range(surface.n_cells)]
    for cells in edge_to_cells.values():
        if len(cells) < 2:
            continue
        for i in cells:
            for j in cells:
                if i != j:
                    adj[i].append(j)
    return adj


def tag_inside_shared_boundary(
    surface: pv.DataSet,
    seed_point: np.ndarray,
    intraventricular_field: str,
    transmural_field: str,
    lv_value: int,
    lv_transmural_eq: float,
    lv_transmural_eps: float,
    shared_tag_name: str,
    inside_tag_name: str,
) -> pv.DataSet:
    if shared_tag_name not in surface.cell_data:
        raise RuntimeError(f"Missing cell_data '{shared_tag_name}' on surface.")
    if intraventricular_field not in surface.point_data:
        raise RuntimeError(
            f"Missing point_data '{intraventricular_field}' on surface."
        )
    if transmural_field not in surface.point_data:
        raise RuntimeError(f"Missing point_data '{transmural_field}' on surface.")

    uvc_trans = surface.point_data[transmural_field]
    uvc_intra = surface.point_data[intraventricular_field]
    is_lv_epi_pt = (uvc_intra == lv_value) & (
        np.abs(uvc_trans - lv_transmural_eq) <= lv_transmural_eps
    )

    lv_epi_cell = np.zeros(surface.n_cells, dtype=bool)
    for i in range(surface.n_cells):
        point_ids = surface.get_cell(i).point_ids
        lv_epi_cell[i] = np.any(is_lv_epi_pt[point_ids])

    centers = surface.cell_centers().points
    start = np.linalg.norm(centers - seed_point, axis=1).argmin()

    shared = surface.cell_data[shared_tag_name].astype(bool)
    inside = np.zeros(surface.n_cells, dtype=int)
    adj = cell_adjacency(surface)
    stack = [start]
    visited = np.zeros(surface.n_cells, dtype=bool)
    visited[start] = True

    while stack:
        c = stack.pop()
        if shared[c]:
            continue
        if not lv_epi_cell[c]:
            continue
        inside[c] = 1
        for nb in adj[c]:
            if not visited[nb]:
                visited[nb] = True
                stack.append(nb)

    surface.cell_data[inside_tag_name] = inside
    return surface


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Tag endocardial surface faces using UVC fields."
    )
    parser.add_argument("--input", required=True, metavar="INPUT_VTK", help="Input VTK mesh.")
    parser.add_argument(
        "--output",
        default=None,
        metavar="OUTPUT_VTK",
        help="Output surface VTK with cell tags.",
    )
    parser.add_argument(
        "--output-volume",
        default=None,
        metavar="OUTPUT_VOLUME_VTK",
        help="Output VTK containing volume-mapped tags.",
    )
    parser.add_argument(
        "--transmural-field",
        default="uvc_transmural",
        help="Point data field for transmural coordinate.",
    )
    parser.add_argument(
        "--intraventricular-field",
        default="uvc_intraventricular",
        help="Point data field for intraventricular tag.",
    )
    parser.add_argument(
        "--transmural-min",
        type=float,
        default=0.0,
        help="Minimum transmural value (endo threshold).",
    )
    parser.add_argument(
        "--transmural-max",
        type=float,
        default=0.1,
        help="Maximum transmural value (endo threshold).",
    )
    parser.add_argument(
        "--transmural-eq",
        type=float,
        default=None,
        help="Exact transmural target (overrides min/max).",
    )
    parser.add_argument(
        "--transmural-eps",
        type=float,
        default=1e-6,
        help="Tolerance for transmural-eq comparison.",
    )
    parser.add_argument(
        "--epi-eq",
        type=float,
        default=1.0,
        help="Exact transmural target for epicardium.",
    )
    parser.add_argument(
        "--epi-eps",
        type=float,
        default=1e-6,
        help="Tolerance for epicardium transmural-eq comparison.",
    )
    parser.add_argument(
        "--lv-value",
        type=int,
        default=-1,
        help="LV tag value in intraventricular field.",
    )
    parser.add_argument(
        "--rv-value",
        type=int,
        default=1,
        help="RV tag value in intraventricular field.",
    )
    parser.add_argument(
        "--mode",
        choices=["any", "all"],
        default="any",
        help="Tag a face if any or all points meet criteria.",
    )
    parser.add_argument(
        "--tag-name",
        default="endo_surface",
        help="Cell data name for the tag.",
    )
    parser.add_argument(
        "--epi-tag-name",
        default="epi_surface",
        help="Cell data name for epicardial tag.",
    )
    parser.add_argument(
        "--shared-tag-name",
        default="shared_boundary",
        help="Cell data name for shared boundary tag.",
    )
    parser.add_argument(
        "--seed-point",
        default=SHARED_BOUNDARY_SEED,
        help="Seed point as 'x,y,z' to grow inside-shared region on LV epi.",
    )
    parser.add_argument(
        "--inside-tag-name",
        default="inside_shared_boundary",
        help="Cell data name for inside-boundary tag.",
    )
    parser.add_argument(
        "--no-surface-output",
        action="store_true",
        help="Skip writing surface outputs.",
    )
    parser.add_argument(
        "--no-volume-output",
        action="store_true",
        help="Skip writing volume-mapped output.",
    )
    args = parser.parse_args()

    print("\n-------------------------\nInspecting Input Fields\n-------------------------")
    inspect_fields(args.input)
    if args.seed_point is None:
        print(
            "Seed point: not set. "
            "Inside-shared boundary growth is disabled. "
            "If boundary artifacts appear, set --seed-point (see tag_endo_epi_surface_config.py)."
        )

    print("\n-------------------------\nTagging Endocardium (LV/RV)\n-------------------------")
    mesh = pv.read(args.input)
    def _postprocess_vtk_inplace(path: Path) -> None:
        convert_vtk_file(str(path), str(path))
        remove_blank_lines(str(path), str(path))
        print(f"Converted file written to {path} (blank lines removed)")
        print(f"Post-processed output: {path}")

    def _write_surface_output(surface: pv.DataSet, output_path: Path):
        output_path.parent.mkdir(parents=True, exist_ok=True)

        surface.save(str(output_path), binary=False)
        print(f"Tagged surface written to {output_path.resolve()}")

        # Keep the original surface output untouched.

    surface = tag_surface(
        mesh=mesh,
        transmural_field=args.transmural_field,
        intraventricular_field=args.intraventricular_field,
        transmural_min=args.transmural_min,
        transmural_max=args.transmural_max,
        transmural_eq=args.transmural_eq,
        transmural_eps=args.transmural_eps,
        ventricle="lv",
        lv_value=args.lv_value,
        rv_value=args.rv_value,
        mode=args.mode,
        tag_name=f"{args.tag_name}_lv",
        extract_surface=True,
    )
    surface = tag_surface(
        mesh=surface,
        transmural_field=args.transmural_field,
        intraventricular_field=args.intraventricular_field,
        transmural_min=args.transmural_min,
        transmural_max=args.transmural_max,
        transmural_eq=args.transmural_eq,
        transmural_eps=args.transmural_eps,
        ventricle="rv",
        lv_value=args.lv_value,
        rv_value=args.rv_value,
        mode=args.mode,
        tag_name=f"{args.tag_name}_rv",
        extract_surface=False,
    )
    lv_count = int(np.count_nonzero(surface.cell_data[f"{args.tag_name}_lv"]))
    rv_count = int(np.count_nonzero(surface.cell_data[f"{args.tag_name}_rv"]))
    print(f"Tagged LV endo cells: {lv_count}")
    print(f"Tagged RV endo cells: {rv_count}")
    surface.cell_data[args.tag_name] = np.maximum(
        surface.cell_data[f"{args.tag_name}_lv"],
        surface.cell_data[f"{args.tag_name}_rv"],
    )
    print(f"Tagged endo total cells: {int(np.count_nonzero(surface.cell_data[args.tag_name]))}")

    print("\n-------------------------\nTagging Epicardium\n-------------------------")
    surface = tag_epi_surface(
        surface=surface,
        transmural_field=args.transmural_field,
        transmural_eq=args.epi_eq,
        transmural_eps=args.epi_eps,
        mode=args.mode,
        tag_name=args.epi_tag_name,
    )
    print(f"Tagged epi cells: {int(np.count_nonzero(surface.cell_data[args.epi_tag_name]))}")

    print("\n-------------------------\nTagging Shared Boundary\n-------------------------")
    surface = tag_shared_boundary(
        surface=surface,
        transmural_field=args.transmural_field,
        intraventricular_field=args.intraventricular_field,
        rv_transmural_eq=0.0,
        rv_transmural_eps=args.transmural_eps,
        lv_transmural_eq=1.0,
        lv_transmural_eps=args.epi_eps,
        lv_value=args.lv_value,
        rv_value=args.rv_value,
        tag_name=args.shared_tag_name,
    )
    print(f"Tagged shared-boundary cells: {int(np.count_nonzero(surface.cell_data[args.shared_tag_name]))}")
    if args.seed_point is not None:
        seed = parse_vector(args.seed_point)
        rv_before = int(np.count_nonzero(surface.cell_data[args.tag_name + "_rv"]))
        surface = tag_inside_shared_boundary(
            surface=surface,
            seed_point=seed,
            intraventricular_field=args.intraventricular_field,
            transmural_field=args.transmural_field,
            lv_value=args.lv_value,
            lv_transmural_eq=1.0,
            lv_transmural_eps=args.epi_eps,
            shared_tag_name=args.shared_tag_name,
            inside_tag_name=args.inside_tag_name,
        )
        surface.cell_data[args.tag_name + "_rv"] = np.maximum(
            surface.cell_data[args.tag_name + "_rv"],
            surface.cell_data[args.inside_tag_name],
        )
        rv_after = int(np.count_nonzero(surface.cell_data[args.tag_name + "_rv"]))
        surface.cell_data[args.tag_name] = np.maximum(
            surface.cell_data[f"{args.tag_name}_lv"],
            surface.cell_data[f"{args.tag_name}_rv"],
        )
        print(
            "Seed point applied: RV endo cells changed from "
            f"{rv_before} to {rv_after}."
        )
        if args.epi_tag_name in surface.cell_data:
            surface.cell_data[args.epi_tag_name] = np.where(
                surface.cell_data[args.inside_tag_name] > 0,
                0,
                surface.cell_data[args.epi_tag_name],
            )
    if args.epi_tag_name in surface.cell_data and args.shared_tag_name in surface.cell_data:
        surface.cell_data[args.epi_tag_name] = np.where(
            surface.cell_data[args.shared_tag_name] > 0,
            0,
            surface.cell_data[args.epi_tag_name],
        )

    if args.epi_tag_name in surface.cell_data and args.tag_name in surface.cell_data:
        overlap = np.logical_and(
            surface.cell_data[args.epi_tag_name] > 0,
            surface.cell_data[args.tag_name] > 0,
        )
        surface.cell_data["epi_endo_overlap"] = overlap.astype(int)
        overlap_count = int(np.count_nonzero(overlap))
        if overlap_count > 0:
            print(
                f"WARNING: Epi/endo overlap detected: {overlap_count} cells. "
                "Tagged as 'epi_endo_overlap'."
            )

    # Preserve a full surface copy for complementary output (keeps all tags).
    complementary_surface = surface.copy(deep=True)

    if args.inside_tag_name in surface.cell_data:
        del surface.cell_data[args.inside_tag_name]
    if args.shared_tag_name in surface.cell_data:
        del surface.cell_data[args.shared_tag_name]
    if "epi_endo_overlap" in surface.cell_data:
        del surface.cell_data["epi_endo_overlap"]

    surface = surface.triangulate()
    complementary_surface = complementary_surface.triangulate()

    input_stem = Path(args.input).stem
    extract_volume_output = EXTRACT_VOLUME_OUTPUT and not args.no_volume_output

    output_path = (
        Path(args.output)
        if args.output
        else (PROJECT_ROOT / "outputs" / f"{input_stem}_endo_epi_surface.vtk")
    )
    surface_unstructured = surface.cast_to_unstructured_grid()
    unused_path = output_path.with_name(f"{input_stem}_endo_epi_surface_unused_points.vtk")
    report_unused_points(
        surface_unstructured,
        "Main surface",
        unused_path,
    )
    surface_unstructured = surface_unstructured.remove_unused_points()
    _write_surface_output(surface_unstructured, output_path)
    convert_vtk_file(str(output_path), str(output_path))
    remove_blank_lines(str(output_path), str(output_path))
    print(f"Converted and cleaned {output_path}")

    extract_complementary_output = COMPLEMENTARY_SURFACE_OUTPUT and not args.no_surface_output
    if extract_complementary_output:
        complementary_path = output_path.with_name(
            f"{input_stem}_endo_epi_surface_complementary.vtk"
        )
        complementary_path.parent.mkdir(parents=True, exist_ok=True)
        complementary_unused = complementary_path.with_name(
            f"{input_stem}_endo_epi_surface_complementary_unused_points.vtk"
        )
        report_unused_points(
            complementary_surface,
            "Complementary surface",
            complementary_unused,
        )
        complementary_surface = complementary_surface.remove_unused_points()
        complementary_surface.save(str(complementary_path), binary=False)
        print(f"Complementary labels written to {complementary_path.resolve()}")
        convert_vtk_file(str(complementary_path), str(complementary_path))
        remove_blank_lines(str(complementary_path), str(complementary_path))
        print(f"Converted and cleaned {complementary_path}")
    else:
        print("Complementary output disabled (no complementary file written).")

    if extract_volume_output:
        volume_path = (
            Path(args.output_volume)
            if args.output_volume
            else (PROJECT_ROOT / "outputs" / f"{input_stem}_endo_epi_volume.vtk")
        )
        volume_path.parent.mkdir(parents=True, exist_ok=True)
        volume_mesh = mesh.copy(deep=True)
        volume_mesh = map_surface_tags_to_volume(
            volume_mesh,
            surface,
            [
                f"{args.tag_name}_lv",
                f"{args.tag_name}_rv",
                args.tag_name,
                args.epi_tag_name,
                "epi_endo_overlap",
            ],
        )
        volume_mesh.save(str(volume_path), binary=False)
        print(f"Volume-mapped labels written to {volume_path.resolve()}")
    else:
        print("Volume output disabled (no volume file written).")

    print("\n-------------------------\nInspecting Final Outputs\n-------------------------")
    inspect_fields(str(output_path))
    if extract_complementary_output:
        inspect_fields(str(complementary_path))
    if extract_volume_output:
        inspect_fields(str(volume_path))


if __name__ == "__main__":
    main()
