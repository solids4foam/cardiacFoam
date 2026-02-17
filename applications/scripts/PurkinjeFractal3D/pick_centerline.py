"""Pick a seed on the endocardium and trace a centerline by following normals."""

import argparse
from pathlib import Path

import numpy as np
import pyvista as pv


def build_endo_surface(mesh: pv.DataSet, transmural_field: str, eps: float) -> pv.PolyData:
    if transmural_field not in mesh.point_data:
        raise KeyError(f"Missing point data field '{transmural_field}'.")
    transmural = np.asarray(mesh.point_data[transmural_field]).reshape(-1)
    endo_mask = np.isclose(transmural, 0.0, atol=eps)
    if not np.any(endo_mask):
        raise RuntimeError("No endocardial points found (transmural ~ 0).")
    surface = mesh.extract_surface().triangulate()
    endo_surface = surface.extract_points(endo_mask, adjacent_cells=True)
    endo_surface = endo_surface.extract_surface().triangulate()
    if endo_surface.n_cells == 0:
        raise RuntimeError("No endocardial surface cells found.")
    endo_surface = endo_surface.compute_normals(point_normals=True, cell_normals=False)
    return endo_surface


def average_normal(endo_surface: pv.PolyData, point: np.ndarray, radius: float) -> np.ndarray:
    distances = np.linalg.norm(endo_surface.points - point, axis=1)
    idx = np.where(distances <= radius)[0]
    if idx.size == 0:
        idx = np.array([int(distances.argmin())])
    normals = endo_surface.point_data["Normals"][idx]
    n = normals.mean(axis=0)
    n /= np.linalg.norm(n)
    return n


def trace_line(
    endo_surface: pv.PolyData,
    start: np.ndarray,
    step_size: float,
    max_steps: int,
    radius: float,
    normal_dot_min: float,
) -> np.ndarray:
    points = [start]
    direction = average_normal(endo_surface, start, radius)
    for _ in range(max_steps):
        candidate = points[-1] + direction * step_size
        closest_id = endo_surface.find_closest_point(candidate)
        next_point = endo_surface.points[closest_id]
        next_dir = average_normal(endo_surface, next_point, radius)
        if np.dot(next_dir, direction) < 0:
            next_dir = -next_dir
        if np.dot(next_dir, direction) < normal_dot_min:
            break
        direction = next_dir
        points.append(next_point)
    return np.array(points)


def to_polyline(points: np.ndarray) -> pv.PolyData:
    poly = pv.PolyData()
    poly.points = points
    lines = np.hstack([[points.shape[0]], np.arange(points.shape[0])])
    poly.lines = lines
    return poly


def pick_point(mesh: pv.PolyData) -> np.ndarray:
    picked = {"point": None}

    def _callback(point):
        picked["point"] = np.array(point)

    plotter = pv.Plotter()
    plotter.add_mesh(mesh, color="white", opacity=0.3, show_edges=False)
    plotter.enable_point_picking(callback=_callback, show_message=True, use_mesh=True)
    plotter.show()
    if picked["point"] is None:
        raise RuntimeError("No point picked.")
    return picked["point"]


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Pick a seed and trace a centerline along endocardial normals."
    )
    parser.add_argument("--input", required=True, help="Input VTK/VTU mesh.")
    parser.add_argument(
        "--transmural-field",
        default="uvc_transmural",
        help="Point data field for transmural coordinate.",
    )
    parser.add_argument(
        "--longitudinal-field",
        default="uvc_longitudinal",
        help="Point data field for longitudinal coordinate.",
    )
    parser.add_argument(
        "--transmural-eps",
        type=float,
        default=1e-3,
        help="Tolerance for transmural ~ 0 endocardium.",
    )
    parser.add_argument(
        "--step-size",
        type=float,
        default=2.0,
        help="Step size along the normal direction.",
    )
    parser.add_argument(
        "--max-steps",
        type=int,
        default=200,
        help="Maximum steps in each direction.",
    )
    parser.add_argument(
        "--radius",
        type=float,
        default=5.0,
        help="Neighborhood radius for averaging normals.",
    )
    parser.add_argument(
        "--normal-dot-min",
        type=float,
        default=0.9,
        help="Minimum dot product to keep normal direction consistent.",
    )
    parser.add_argument(
        "--output",
        default=None,
        help="Output polyline VTK (defaults to <input>_centerline.vtk).",
    )
    parser.add_argument(
        "--plot",
        action="store_true",
        help="Plot the mesh, sampled normals, and centerline.",
    )
    args = parser.parse_args()

    mesh = pv.read(args.input)
    if args.longitudinal_field not in mesh.point_data:
        raise KeyError(f"Missing point data field '{args.longitudinal_field}'.")
    endo = build_endo_surface(mesh, args.transmural_field, args.transmural_eps)

    seed = pick_point(endo)
    forward = trace_line(
        endo,
        seed,
        args.step_size,
        args.max_steps,
        args.radius,
        args.normal_dot_min,
    )
    backward = trace_line(
        endo,
        seed,
        -args.step_size,
        args.max_steps,
        args.radius,
        args.normal_dot_min,
    )
    centerline = np.vstack([backward[::-1][:-1], forward])
    polyline = to_polyline(centerline)

    longitudinal = np.asarray(mesh.point_data[args.longitudinal_field]).reshape(-1)
    ids = [mesh.find_closest_point(p) for p in centerline]
    line_long = longitudinal[ids]
    polyline.point_data["uvc_longitudinal"] = line_long

    base_idx = int(np.argmax(line_long))
    apex_idx = int(np.argmin(line_long))
    print(f"Base point: {centerline[base_idx]} (long={line_long[base_idx]})")
    print(f"Apex point: {centerline[apex_idx]} (long={line_long[apex_idx]})")

    output_path = (
        Path(args.output)
        if args.output
        else Path(args.input).with_name(f"{Path(args.input).stem}_centerline.vtk")
    )
    polyline.save(str(output_path), binary=False)
    print(f"Centerline written to: {output_path}")

    if args.plot:
        plotter = pv.Plotter()
        plotter.add_mesh(endo, color="white", opacity=0.25)
        plotter.add_mesh(polyline, color="blue", line_width=4)
        plotter.add_points(centerline[base_idx], color="green", point_size=15)
        plotter.add_points(centerline[apex_idx], color="red", point_size=15)
        plotter.show()

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
