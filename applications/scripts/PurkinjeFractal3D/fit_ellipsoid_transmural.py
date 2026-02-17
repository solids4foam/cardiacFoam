"""Fit a canonical open ellipsoid volume and write a VTK with uvc_transmural."""

import argparse
from pathlib import Path

import numpy as np
import pyvista as pv


def fit_ellipsoid(points: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return (center, axes, rotation) for a PCA-aligned ellipsoid fit."""
    center = points.mean(axis=0)
    centered = points - center
    _, _, vh = np.linalg.svd(centered, full_matrices=False)
    rotation = vh.T
    local = centered @ rotation
    axes = np.max(np.abs(local), axis=0)
    axes = np.where(axes == 0, 1.0, axes)
    return center, axes, rotation


def build_ellipsoid_volume(
    center: np.ndarray,
    axes: np.ndarray,
    rotation: np.ndarray,
    theta_res: int,
    phi_res: int,
    radial_res: int,
    expansion: float,
) -> pv.PolyData:
    theta = np.linspace(0.0, 2.0 * np.pi, theta_res)
    phi = np.linspace(0.0, np.pi, phi_res)
    radial = np.linspace(0.0, 1.0, radial_res)

    theta_grid, phi_grid, radial_grid = np.meshgrid(
        theta, phi, radial, indexing="xy"
    )

    axes_x = axes[0] + expansion * radial_grid
    axes_y = axes[1] + expansion * radial_grid
    axes_z = axes[2] + expansion * radial_grid

    x = axes_x * np.sin(phi_grid) * np.cos(theta_grid)
    y = axes_y * np.sin(phi_grid) * np.sin(theta_grid)
    z = axes_z * np.cos(phi_grid)

    local = np.column_stack((x.ravel(), y.ravel(), z.ravel()))
    points = local @ rotation.T + center

    mesh = pv.StructuredGrid(
        x.reshape(phi_res, theta_res, radial_res),
        y.reshape(phi_res, theta_res, radial_res),
        z.reshape(phi_res, theta_res, radial_res),
    )
    mesh.points = points

    # Transmural scalar for the volume (0 = inner, 1 = outer).
    mesh.point_data["uvc_transmural"] = radial_grid.ravel()
    return mesh


def parse_vector(value: str) -> np.ndarray:
    parts = [p.strip() for p in value.split(",") if p.strip()]
    if len(parts) != 3:
        raise ValueError("Expected 3 comma-separated values for a vector.")
    return np.array([float(p) for p in parts])


def compute_axis_from_endo_normals(
    mesh: pv.DataSet,
    endo_mask: np.ndarray,
    sample_count: int,
) -> tuple[np.ndarray, np.ndarray]:
    surface = mesh.extract_surface().triangulate()
    endo_surface = surface.extract_points(endo_mask, adjacent_cells=True)
    if endo_surface.n_cells == 0:
        raise RuntimeError("No endocardial surface cells found for normals.")
    endo_surface = endo_surface.extract_surface().triangulate()
    if endo_surface.n_cells == 0:
        raise RuntimeError("No endocardial surface cells found for normals.")
    endo_surface = endo_surface.compute_normals(
        point_normals=True,
        cell_normals=False,
    )
    normals = endo_surface.point_data.get("Normals")
    if normals is None or normals.size == 0:
        raise RuntimeError("No normals computed for endocardial surface.")

    mean_normal = normals.mean(axis=0)
    mean_normal /= np.linalg.norm(mean_normal)

    rng = np.random.default_rng(1234)
    count = min(sample_count, endo_surface.n_points)
    sample_idx = rng.choice(endo_surface.n_points, size=count, replace=False)
    sample_points = endo_surface.points[sample_idx]
    sample_normals = normals[sample_idx]
    return mean_normal, np.column_stack((sample_points, sample_normals))


def find_apex(points: np.ndarray, longitudinal: np.ndarray, tol: float) -> np.ndarray:
    mask = np.isclose(longitudinal, 0.0, atol=tol)
    if not np.any(mask):
        idx = int(np.argmin(longitudinal))
        return points[idx]
    return points[mask].mean(axis=0)


def select_axis(rotation: np.ndarray, axes: np.ndarray, mode: str) -> np.ndarray:
    if mode == "auto":
        idx = int(np.argmax(axes))
    else:
        idx = int(mode)
    axis = rotation[:, idx]
    axis /= np.linalg.norm(axis)
    return axis


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Fit an ellipsoid and write a VTK surface with uvc_transmural."
    )
    parser.add_argument("--input", required=True, help="Input VTK/VTU mesh file.")
    parser.add_argument(
        "--output",
        default=None,
        help="Output VTK file (defaults to <input>_ellipsoid_volume.vtk).",
    )
    parser.add_argument(
        "--expansion",
        type=float,
        default=9.0,
        help="Outward expansion from the inner surface (same units as mesh).",
    )
    parser.add_argument(
        "--theta-res",
        type=int,
        default=120,
        help="Azimuthal resolution for ellipsoid volume.",
    )
    parser.add_argument(
        "--phi-res",
        type=int,
        default=60,
        help="Polar resolution for ellipsoid volume.",
    )
    parser.add_argument(
        "--radial-res",
        type=int,
        default=16,
        help="Number of layers through the wall (transmural).",
    )
    parser.add_argument(
        "--cut-distance",
        type=float,
        required=True,
        help="Distance along the axis from the origin to place the opening cut plane.",
    )
    parser.add_argument(
        "--cut-origin",
        default=None,
        help="Origin for the cut plane distance as 'x,y,z' (defaults to ellipsoid center).",
    )
    parser.add_argument(
        "--cut-side",
        choices=["base", "apex"],
        default="base",
        help="Keep base-side or apex-side of the cut plane.",
    )
    parser.add_argument(
        "--longitudinal-field",
        default="uvc_longitudinal",
        help="Point data field used to locate the apex (near 0).",
    )
    parser.add_argument(
        "--transmural-field",
        default="uvc_transmural",
        help="Point data field used to select endocardium (near 0).",
    )
    parser.add_argument(
        "--longitudinal-eps",
        type=float,
        default=1e-3,
        help="Tolerance for longitudinal ~ 0 when locating the apex.",
    )
    parser.add_argument(
        "--transmural-eps",
        type=float,
        default=1e-3,
        help="Tolerance for transmural ~ 0 when selecting endocardium points.",
    )
    parser.add_argument(
        "--plot-axis",
        action="store_true",
        help="Plot the mesh, apex, and fitted axis in PyVista.",
    )
    parser.add_argument(
        "--normal-samples",
        type=int,
        default=300,
        help="Number of endocardial normals to sample for plotting.",
    )
    args = parser.parse_args()

    input_path = Path(args.input)
    mesh = pv.read(str(input_path))
    points = mesh.points
    if points.size == 0:
        raise RuntimeError(f"No points found in {input_path}")

    if args.longitudinal_field not in mesh.point_data:
        raise KeyError(
            f"Missing point data field '{args.longitudinal_field}' in {input_path}"
        )
    if args.transmural_field not in mesh.point_data:
        raise KeyError(
            f"Missing point data field '{args.transmural_field}' in {input_path}"
        )

    longitudinal = np.asarray(mesh.point_data[args.longitudinal_field]).reshape(-1)
    transmural = np.asarray(mesh.point_data[args.transmural_field]).reshape(-1)
    apex = find_apex(points, longitudinal, args.longitudinal_eps)

    endo_mask = np.isclose(transmural, 0.0, atol=args.transmural_eps)
    if not np.any(endo_mask):
        raise RuntimeError(
            "No endocardial points found (transmural ~ 0). "
            "Increase --transmural-eps."
        )
    axis, normal_samples = compute_axis_from_endo_normals(
        mesh,
        endo_mask,
        args.normal_samples,
    )
    apex_center = apex

    if args.plot_axis:
        plotter = pv.Plotter()
        plotter.add_mesh(mesh, color="white", opacity=0.2, show_edges=False)
        plotter.add_points(apex_center, color="red", point_size=15, label="apex")
        plotter.add_lines(
            np.vstack([apex_center - axis * 50.0, apex_center + axis * 50.0]),
            color="blue",
            width=4,
            label="axis",
        )
        if normal_samples.size > 0:
            pts = normal_samples[:, :3]
            vecs = normal_samples[:, 3:]
            plotter.add_arrows(pts, vecs, mag=5.0, color="orange")
        plotter.show()

    center, axes, rotation = fit_ellipsoid(points)
    ellipsoid = build_ellipsoid_volume(
        center,
        axes,
        rotation,
        args.theta_res,
        args.phi_res,
        args.radial_res,
        args.expansion,
    )

    origin = center if args.cut_origin is None else parse_vector(args.cut_origin)
    projections = (mesh.points - origin) @ axis
    base_proj = projections.min()
    apex_proj = projections.max()
    print(f"Axis direction: {axis}")
    print(f"Cut origin: {origin}")
    print(f"Base (min proj): {base_proj}")
    print(f"Apex (max proj): {apex_proj}")

    plane_center = origin + axis * args.cut_distance
    signed = (ellipsoid.points - plane_center) @ axis
    if args.cut_side == "base":
        keep = signed <= 0
    else:
        keep = signed >= 0
    ellipsoid = ellipsoid.extract_points(keep, adjacent_cells=False)

    output_path = (
        Path(args.output)
        if args.output
        else input_path.with_name(f"{input_path.stem}_ellipsoid_volume.vtk")
    )
    ellipsoid.save(str(output_path), binary=False)
    print(f"Ellipsoid written to: {output_path}")
    print(f"Center: {center}")
    print(f"Axes: {axes}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
