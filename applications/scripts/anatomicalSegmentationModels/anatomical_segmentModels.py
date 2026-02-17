
#!/usr/bin/env python3
import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pyvista as pv
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401
from src.bullseye_imageCreator import (
    draw_all_bullseyes,
    draw_bullseye_reference,
)

LONGITUDINAL_FIELD = "uvc_longitudinal"
LONGITUDINAL_EPS = 0.01
ANGULAR_DIV_EPS = 0.08  # radians; boundary-line thickness in angle-division plot
GROOVE_PROJ_EPS = 0.06  # radians; thickness for groove-angle projection on surface
# Hardcoded plotting caps (planes mode) to reduce render time.
PLOT_CAP_PLANE_POINTS = 800
PLOT_CAP_Z_DIVISION_POINTS = 1200
PLOT_CAP_THETA_DIVISION_POINTS = 600
PLOT_CAP_SEGMENT_POINTS = 500
PLOT_CAP_INTERFACE_COMPONENT_POINTS = 1200
PLOT_CAP_GROOVE_INTERFACE_POINTS = 1000

# Anatomy-based circumferential segmentation utilities
from src.segment_layout import (  # noqa: E402
    GrooveAngles,
    build_ring_boundaries,
    cw_is_between,
    groove_angles_debug_from_intraventricular_interface,
    groove_angles_from_point_mask,
    interface_components_from_intraventricular,
    segment_ids_from_theta_layout,
)

try:
    import configs.lv_pca_axis_config as cfg  # type: ignore
except Exception:
    cfg = None
try:
    import configs.lv_division_config as dcfg  # type: ignore
except Exception:
    dcfg = None

z_apical_apex = float(getattr(cfg, "z_apical_apex", 0.08))
z_apical_mid = float(getattr(cfg, "z_apical_mid", 1.0 / 3.0))
z_mid_basal = float(getattr(cfg, "z_mid_basal", 2.0 / 3.0))
z_cap_valves = float(getattr(cfg, "z_cap_valves", 0.9))
LV_INTRAVENTRICULAR_TAG = int(getattr(cfg, "LV_INTRAVENTRICULAR_TAG", 1))
RV_INTRAVENTRICULAR_TAG = int(getattr(cfg, "RV_INTRAVENTRICULAR_TAG", -1))
ROT_START = float(getattr(cfg, "ROT_START", -np.pi))

# Anatomy-angle defaults. CLI values still override these when provided.
ANT_GROOVE_THETA = getattr(cfg, "ant_groove_theta", None)
POST_GROOVE_THETA = getattr(cfg, "post_groove_theta", None)
ANT_GROOVE_FIELD = getattr(cfg, "ant_groove_field", None)
POST_GROOVE_FIELD = getattr(cfg, "post_groove_field", None)
ANGLE_REFERENCE_IMAGE = getattr(cfg, "angle_reference_image", "aha_lv17_reference.png")
SEPTAL_CLICK_XYZ = getattr(cfg, "septal_click_xyz", None)
BULLSEYE_PLOTS_DIR = Path(__file__).resolve().parent / "inputs" / "bullseyePlots"
OUTPUT_DIR = Path(__file__).resolve().parent / "outputs"
# Fixed 2D bullseye orientation for visualization (independent from groove angles).
BULLSEYE_START_ANGLE_DEG = 120.0

# Division/segmentation scheme configuration (separate file).
DIVISION_REFERENCE = getattr(dcfg, "division_reference", "aha17")
DIVISION_REFERENCES = getattr(dcfg, "division_references", {})
_DEFAULT_DIVISION_REF = {
    "ring_layout": {"basal": [2, 4], "mid": [2, 4], "apical": [1, 3]},
    "ring_id_map": {"basal": [1, 2, 3, 4, 5, 6], "mid": [7, 8, 9, 10, 11, 12], "apical": [13, 14, 15, 16]},
    "apex_id": 17,
    "use_apex_cap": True,
    "colors": [
        "purple", "cyan", "yellow", "magenta", "lime", "orange",
        "red", "green", "blue", "pink", "brown", "gray",
        "gold", "navy", "teal", "olive", "black",
    ],
    "z_division_color": "white",
    "theta_division_color": "black",
}
DIVISION_CFG = DIVISION_REFERENCES.get(DIVISION_REFERENCE, _DEFAULT_DIVISION_REF)
RING_LAYOUT = DIVISION_CFG.get("ring_layout", _DEFAULT_DIVISION_REF["ring_layout"])
ANGLE_SEG_PALETTE = list(DIVISION_CFG.get("colors", _DEFAULT_DIVISION_REF["colors"]))
Z_DIVISION_COLOR = str(DIVISION_CFG.get("z_division_color", _DEFAULT_DIVISION_REF["z_division_color"]))
THETA_DIVISION_COLOR = str(
    DIVISION_CFG.get("theta_division_color", _DEFAULT_DIVISION_REF["theta_division_color"])
)
SEGMENT_ARTERY_MAP = dict(DIVISION_CFG.get("segment_artery_map", {}))
ARTERY_COLORS = dict(DIVISION_CFG.get("artery_colors", {}))
COLOR_SCHEME_DEFAULT = str(
    DIVISION_CFG.get("color_scheme_default", "segment")
).lower()


def flow(msg: str) -> None:
    """Minimal terminal flow message."""
    print(f"----------\n{msg}\n----------")


def show_reference_figures(
    division_reference: str,
    division_cfg: dict,
) -> None:
    """Show and save bullseye reference figures with fixed standard orientation."""
    BULLSEYE_PLOTS_DIR.mkdir(parents=True, exist_ok=True)
    out_seg = BULLSEYE_PLOTS_DIR / f"{division_reference}_bullseye.png"
    try:
        draw_bullseye_reference(
            division_reference,
            division_cfg,
            out_path=str(out_seg),
            start_angle_deg=BULLSEYE_START_ANGLE_DEG,
            color_mode="segment",
            show_labels=True,
        )
        if division_reference == "aha17" and division_cfg.get("artery_colors") and division_cfg.get("segment_artery_map"):
            out_art = BULLSEYE_PLOTS_DIR / f"{division_reference}_bullseye_artery.png"
            draw_bullseye_reference(
                division_reference,
                division_cfg,
                out_path=str(out_art),
                start_angle_deg=BULLSEYE_START_ANGLE_DEG,
                color_mode="artery",
                show_labels=True,
            )
        plt.show(block=False)
    except Exception:
        return

def _clear_active_renderer(plotter: pv.Plotter):
    """Clear only current subplot renderer actors to avoid background renderer bugs."""
    try:
        plotter.clear_actors()
        return
    except Exception:
        pass
    try:
        ren = plotter.renderer
        for name in list(getattr(ren, "actors", {}).keys()):
            try:
                plotter.remove_actor(name, reset_camera=False, render=False)
            except Exception:
                pass
    except Exception:
        pass


def _build_ring_id_map(layout: dict, explicit_map: dict | None) -> dict:
    """Create ring->segment-id map; use explicit map when valid, else auto-generate."""
    out = {}
    next_id = 1
    for ring in ("basal", "mid", "apical"):
        if ring not in layout:
            raise ValueError(f"Division ring_layout missing required ring '{ring}'.")
        septal_parts, freewall_parts = layout[ring]
        n = int(septal_parts) + int(freewall_parts)
        ids = None if explicit_map is None else explicit_map.get(ring)
        if ids is not None and len(ids) == n:
            out[ring] = [int(x) for x in ids]
            next_id = max(next_id, max(out[ring]) + 1)
        else:
            out[ring] = list(range(next_id, next_id + n))
            next_id += n
    return out


RING_ID_MAP = _build_ring_id_map(RING_LAYOUT, DIVISION_CFG.get("ring_id_map", None))
_max_ring_id = max(i for vals in RING_ID_MAP.values() for i in vals)
USE_APEX_CAP = bool(DIVISION_CFG.get("use_apex_cap", _DEFAULT_DIVISION_REF["use_apex_cap"]))
_apex_id_raw = DIVISION_CFG.get("apex_id", _max_ring_id + 1)
APEX_SEGMENT_ID = None if _apex_id_raw is None else int(_apex_id_raw)

def token_stream(f, bufsize=1 << 20):
    """Yield whitespace-separated tokens from a large text file efficiently."""
    buf = ""
    while True:
        chunk = f.read(bufsize)
        if not chunk:
            break
        buf += chunk
        parts = buf.split()
        if buf and not buf[-1].isspace():
            buf = parts.pop()  # keep last partial token
        else:
            buf = ""
        for p in parts:
            yield p
    if buf:
        for p in buf.split():
            yield p


def read_vtk_legacy_offsets_connectivity(path: Path):
    """
    Reads VTK legacy ASCII unstructured grid with OFFSETS/CONNECTIVITY layout.

    Returns:
      pts  (N,3) float64
      conn (Ncells,4) int32   (tet4 connectivity)
      tags (Ncells,) int32    (CELL_DATA SCALARS tags)
    """
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        ts = token_stream(f)

        # ---- POINTS ----
        for tok in ts:
            if tok == "POINTS":
                break
        else:
            raise ValueError("POINTS not found")

        n_points = int(next(ts))
        _pts_dtype = next(ts)  # float/double (ignored)
        pts = np.fromiter(
            (float(next(ts)) for _ in range(3 * n_points)),
            dtype=np.float64,
            count=3 * n_points,
        ).reshape((n_points, 3))

        # ---- CELLS ----
        for tok in ts:
            if tok == "CELLS":
                break
        else:
            raise ValueError("CELLS not found")

        # In your file: CELLS N_offsets N_connectivity
        n_offsets = int(next(ts))
        n_conn = int(next(ts))
        n_cells = n_offsets - 1

        # ---- OFFSETS ----
        for tok in ts:
            if tok == "OFFSETS":
                break
        else:
            raise ValueError("OFFSETS not found")
        _off_dtype = next(ts)  # ignored
        _offsets = np.fromiter(
            (int(next(ts)) for _ in range(n_offsets)),
            dtype=np.int64,
            count=n_offsets,
        )

        # ---- CONNECTIVITY ----
        for tok in ts:
            if tok == "CONNECTIVITY":
                break
        else:
            raise ValueError("CONNECTIVITY not found")
        _conn_dtype = next(ts)  # ignored
        conn = np.fromiter(
            (int(next(ts)) for _ in range(n_conn)),
            dtype=np.int32,
            count=n_conn,
        ).reshape((n_cells, 4))

        # ---- CELL_TYPES (ignore) ----
        for tok in ts:
            if tok == "CELL_TYPES":
                break
        else:
            raise ValueError("CELL_TYPES not found")
        n_ct = int(next(ts))
        _ = np.fromiter((int(next(ts)) for _ in range(n_ct)), dtype=np.int32, count=n_ct)

        # ---- CELL_DATA ----
        for tok in ts:
            if tok == "CELL_DATA":
                break
        else:
            raise ValueError("CELL_DATA not found")
        n_cd = int(next(ts))
        if n_cd != n_cells:
            raise ValueError(f"CELL_DATA size {n_cd} != n_cells {n_cells}")

        # ---- Find SCALARS tags ----
        found = False
        for tok in ts:
            if tok == "SCALARS":
                name = next(ts)
                _dtype = next(ts)
                peek = next(ts)

                # legacy may include optional ncomp
                if peek.isdigit():
                    _ncomp = int(peek)
                    tok = next(ts)  # should be LOOKUP_TABLE
                else:
                    tok = peek

                if name == "tags":
                    found = True
                    if tok != "LOOKUP_TABLE":
                        for tok in ts:
                            if tok == "LOOKUP_TABLE":
                                break
                    _ = next(ts)  # "default"
                    tags_f = np.fromiter(
                        (float(next(ts)) for _ in range(n_cells)),
                        dtype=np.float32,
                        count=n_cells,
                    )
                    tags = tags_f.astype(np.int32)
                    break

        if not found:
            raise ValueError("SCALARS tags not found in CELL_DATA")

    return pts, conn, tags


def tet_volumes_and_centroids(pts: np.ndarray, tet_conn: np.ndarray):
    """Vectorized volume and centroid for tet4 cells."""
    p0 = pts[tet_conn[:, 0]]
    p1 = pts[tet_conn[:, 1]]
    p2 = pts[tet_conn[:, 2]]
    p3 = pts[tet_conn[:, 3]]

    vol = np.abs(np.einsum("ij,ij->i", np.cross(p1 - p0, p2 - p0), (p3 - p0))) / 6.0
    cen = (p0 + p1 + p2 + p3) / 4.0
    return vol, cen


def volume_weighted_pca_axis(
    pts: np.ndarray,
    conn: np.ndarray,
    tags: np.ndarray,
    tag: int = 1,
    chunk_cells: int = 200000,
):
    """
    Volume-weighted PCA using tet centroids as mass-lumps.

    Returns:
      mu    (3,) centroid
      v1    (3,) principal eigenvector (unit) for largest eigenvalue
      evals (3,) eigenvalues sorted desc
      conn_tag (Ntets,4) connectivity for filtered cells
    """
    mask = tags == tag
    conn_tag = conn[mask]
    if conn_tag.size == 0:
        raise ValueError(f"No cells found with tag={tag}")

    # centroid
    sumV = 0.0
    sumVc = np.zeros(3, dtype=np.float64)
    for i in range(0, conn_tag.shape[0], chunk_cells):
        c = conn_tag[i : i + chunk_cells]
        V, C = tet_volumes_and_centroids(pts, c)
        sumV += float(V.sum())
        sumVc += (V[:, None] * C).sum(axis=0)
    mu = sumVc / sumV

    # covariance
    Cov = np.zeros((3, 3), dtype=np.float64)
    for i in range(0, conn_tag.shape[0], chunk_cells):
        c = conn_tag[i : i + chunk_cells]
        V, C = tet_volumes_and_centroids(pts, c)
        r = C - mu
        Cov += (r.T * V) @ r
    Cov /= sumV

    # eigen
    evals, evecs = np.linalg.eigh(Cov)
    order = np.argsort(evals)[::-1]
    evals = evals[order]
    evecs = evecs[:, order]
    v1 = evecs[:, 0]
    v1 = v1 / np.linalg.norm(v1)

    return mu, v1, evals, conn_tag


def plot_interface_components(
    plotter: pv.Plotter,
    mesh: pv.DataSet,
    intraventricular: np.ndarray,
    lv_tag: float,
    rv_tag: float,
    point_size: float = 10.0,
) -> int:
    """Overlay detected LV/RV interface components on an existing plotter."""
    comps = interface_components_from_intraventricular(
        mesh=mesh,
        intraventricular=intraventricular,
        lv_tag=lv_tag,
        rv_tag=rv_tag,
        min_component_points=10,
    )
    if not comps:
        return 0

    colors = ["deepskyblue", "limegreen", "gold", "tomato", "fuchsia", "white"]
    rng = np.random.default_rng(0)
    for i, comp in enumerate(comps[:6]):
        pts = mesh.points[comp]
        if pts.shape[0] > PLOT_CAP_INTERFACE_COMPONENT_POINTS:
            idx = rng.choice(
                pts.shape[0], size=PLOT_CAP_INTERFACE_COMPONENT_POINTS, replace=False
            )
            pts = pts[idx]
        color = colors[i % len(colors)]
        plotter.add_points(pts, color=color, point_size=point_size)
        center = pts.mean(axis=0)
        plotter.add_points(np.asarray([center]), color="black", point_size=point_size + 4.0)
    return len(comps)


def circular_distance(theta: np.ndarray, center: float) -> np.ndarray:
    """Shortest signed angular distance magnitude in radians."""
    return np.abs((theta - center + np.pi) % (2.0 * np.pi) - np.pi)


def seg_id_from_area_name(area_name: str) -> int:
    """Extract segment id from '..._seg<id>' area name."""
    try:
        return int(area_name.rsplit("_seg", 1)[1])
    except Exception:
        return -1


def color_for_segment(seg_id: int, mode: str) -> str:
    """Resolve display color for one segment id and color mode."""
    if mode == "artery":
        art = SEGMENT_ARTERY_MAP.get(int(seg_id))
        if art is not None and art in ARTERY_COLORS:
            return str(ARTERY_COLORS[art])
    return ANGLE_SEG_PALETTE[(int(seg_id) - 1) % len(ANGLE_SEG_PALETTE)]


def theta_from_direction_reference(
    lv_points: np.ndarray,
    lv_theta: np.ndarray,
    lv_center: np.ndarray,
    direction_vec: np.ndarray,
) -> float | None:
    """Map a spatial direction to nearest LV rotational theta."""
    v = np.asarray(direction_vec, dtype=float)
    n = np.linalg.norm(v)
    if n <= 1.0e-12:
        return None
    v = v / n
    rel = lv_points - lv_center
    rel_n = np.linalg.norm(rel, axis=1)
    valid = rel_n > 1.0e-12
    if not np.any(valid):
        return None
    rel_u = rel[valid] / rel_n[valid][:, None]
    dots = rel_u @ v
    idx = int(np.argmax(dots))
    return float(lv_theta[valid][idx])


def orient_grooves_with_septal_theta(grooves: GrooveAngles, septal_theta: float) -> GrooveAngles:
    """Ensure septal arc (anterior->posterior CW) contains septal_theta."""
    test = np.array([septal_theta], dtype=float)
    if bool(cw_is_between(test, grooves.anterior, grooves.posterior)[0]):
        return grooves
    return GrooveAngles(grooves.posterior, grooves.anterior)


def plot_each_division(mesh: pv.DataSet, bins) -> None:
    """Open one plot window per longitudinal z-band, showing rotational subdivisions."""
    if not bins:
        return
    grouped = {}
    for pts, _bin_center, seg_id, area_name in bins:
        # area_name pattern: lv_<z_band>_seg<id>
        z_name = area_name.rsplit("_seg", 1)[0]
        grouped.setdefault(z_name, []).append((pts, int(seg_id)))
    for z_name in sorted(grouped):
        p = pv.Plotter()
        p.add_mesh(mesh, color="white", opacity=0.08, show_edges=False)
        rng = np.random.default_rng(0)
        seg_labels = []
        for pts, seg_id in grouped[z_name]:
            color = color_for_segment(int(seg_id), COLOR_SCHEME_DEFAULT)
            pts_plot = pts
            if pts_plot.shape[0] > PLOT_CAP_SEGMENT_POINTS:
                idx = rng.choice(pts_plot.shape[0], size=PLOT_CAP_SEGMENT_POINTS, replace=False)
                pts_plot = pts_plot[idx]
            p.add_points(pts_plot, color=color, point_size=6)
            seg_labels.append(f"seg{int(seg_id)}")
        p.add_text(f"Z band: {z_name}", position="upper_left", font_size=12)
        p.add_text(
            "Rotation divisions: " + ", ".join(sorted(seg_labels)),
            position="lower_left",
            font_size=10,
        )
        p.show()


def plot_groove_diagnostics(
    mesh: pv.DataSet,
    rotational: np.ndarray,
    intrav: np.ndarray,
    grooves: GrooveAngles,
    groove_debug: dict,
    z_values: np.ndarray,
    z_min_exclusive: float | None,
    z_cap_max: float,
    lv_tag: float,
    rv_tag: float,
) -> None:
    """Plot points used for groove estimation and groove-angle surface projections."""
    p = pv.Plotter()
    p.add_mesh(mesh, color="white", opacity=0.08, show_edges=False)

    comps = interface_components_from_intraventricular(
        mesh=mesh,
        intraventricular=intrav,
        longitudinal_z=z_values,
        z_min_exclusive=z_min_exclusive,
        z_cap_max=z_cap_max,
        lv_tag=lv_tag,
        rv_tag=rv_tag,
        min_component_points=10,
    )
    for comp in comps:
        pts = mesh.points[comp]
        if pts.shape[0] > PLOT_CAP_GROOVE_INTERFACE_POINTS:
            idx = np.random.default_rng(0).choice(
                pts.shape[0], size=PLOT_CAP_GROOVE_INTERFACE_POINTS, replace=False
            )
            pts = pts[idx]
        p.add_points(pts, color="lightgray", point_size=4)

    ant_ids = np.asarray(groove_debug.get("anterior_point_ids", []), dtype=np.int64)
    post_ids = np.asarray(groove_debug.get("posterior_point_ids", []), dtype=np.int64)
    if ant_ids.size > 0:
        p.add_points(mesh.points[ant_ids], color="red", point_size=9)
    if post_ids.size > 0:
        p.add_points(mesh.points[post_ids], color="blue", point_size=9)

    lv_mask = np.isclose(intrav, lv_tag, atol=1.0e-8) & (z_values <= z_cap_max)
    rv_mask = np.isclose(intrav, rv_tag, atol=1.0e-8) & (z_values <= z_cap_max)
    if z_min_exclusive is not None:
        lv_mask = lv_mask & (z_values > z_min_exclusive)
        rv_mask = rv_mask & (z_values > z_min_exclusive)
    lv_pts = mesh.points[lv_mask]
    lv_theta = rotational[lv_mask]
    if lv_pts.shape[0] > 0:
        ant_proj = circular_distance(lv_theta, grooves.anterior) <= GROOVE_PROJ_EPS
        post_proj = circular_distance(lv_theta, grooves.posterior) <= GROOVE_PROJ_EPS
        if np.any(ant_proj):
            p.add_points(lv_pts[ant_proj], color="darkred", point_size=5)
        if np.any(post_proj):
            p.add_points(lv_pts[post_proj], color="darkblue", point_size=5)
        lv_center = lv_pts.mean(axis=0)
        p.add_points(np.asarray([lv_center]), color="black", point_size=14)
        if np.any(ant_proj):
            ant_center = lv_pts[ant_proj].mean(axis=0)
            p.add_mesh(pv.Line(lv_center, ant_center), color="darkred", line_width=4)
        if np.any(post_proj):
            post_center = lv_pts[post_proj].mean(axis=0)
            p.add_mesh(pv.Line(lv_center, post_center), color="darkblue", line_width=4)
        if np.any(rv_mask):
            rv_center = mesh.points[rv_mask].mean(axis=0)
            p.add_points(np.asarray([rv_center]), color="magenta", point_size=14)
            p.add_mesh(pv.Line(lv_center, rv_center), color="magenta", line_width=3)

    p.add_text(
        (
            "Groove Diagnostics\n"
            "gray=interface points, red/blue=points used for groove estimate,\n"
            "dark red/blue=projection of each groove angle on LV surface\n"
            f"method={groove_debug.get('method', 'unknown')}, "
            f"z_min_exclusive={z_min_exclusive if z_min_exclusive is not None else 'none'}, "
            f"z_cap_max={z_cap_max:.3f}, "
            f"ant={grooves.anterior:.3f}, post={grooves.posterior:.3f}"
        ),
        position="upper_left",
        font_size=10,
    )
    p.show()


def process_ventricle_centerline_and_slices(
    plotter: pv.Plotter,
    mesh: pv.DataSet,
    values: np.ndarray,
    intrav: np.ndarray,
    vent_label: str,
    vent_value: float,
    targets: list[float],
    longitudinal_eps: float,
    apex_active: bool,
    z_apical_apex: float,
    plane_color: str,
    point_color: str,
    seed: int,
) -> tuple[np.ndarray, np.ndarray, list[str]]:
    """Shared LV/RV extraction for z-min and longitudinal target slice points."""
    vent_centers = []
    legend_lines = []
    vent_mask = np.isclose(intrav, vent_value, atol=longitudinal_eps)
    if np.any(vent_mask):
        vent_values = values[vent_mask]
        z_min = float(vent_values.min())
        min_mask = np.isclose(values, z_min, atol=longitudinal_eps) & vent_mask
        min_pts = mesh.points[min_mask]
        if min_pts.shape[0] == 0:
            vent_pts = mesh.points[vent_mask]
            z_idx = int(np.argmin(np.abs(vent_values - z_min)))
            min_center = vent_pts[z_idx]
        else:
            min_center = min_pts.mean(axis=0)
        vent_centers.append(min_center)
        legend_lines.append(f"{vent_label} z_min={z_min:.3f}: included in centerline")
        plotter.add_points(np.asarray([min_center]), color=plane_color, point_size=14)

    for target in targets:
        mask = np.isclose(values, target, atol=longitudinal_eps) & vent_mask
        pts = mesh.points[mask]
        if pts.shape[0] == 0:
            if apex_active and np.isclose(target, z_apical_apex):
                vent_values = values[vent_mask]
                if vent_values.size > 0:
                    nearest = float(vent_values[np.argmin(np.abs(vent_values - target))])
                    mask = np.isclose(values, nearest, atol=longitudinal_eps) & vent_mask
                    pts = mesh.points[mask]
                if pts.shape[0] == 0:
                    continue
            else:
                continue

        center = pts.mean(axis=0)
        if pts.shape[0] > PLOT_CAP_PLANE_POINTS:
            rng = np.random.default_rng(seed)
            idx = rng.choice(pts.shape[0], size=PLOT_CAP_PLANE_POINTS, replace=False)
            pts_plot = pts[idx]
        else:
            pts_plot = pts
        plotter.add_points(pts_plot, color=point_color, point_size=6)
        vent_centers.append(center)
        legend_lines.append(f"{vent_label} z={target:.3f}: points={pts.shape[0]}")

    return np.asarray(vent_centers), vent_mask, legend_lines


def point_tags_to_cell_tags(mesh: pv.DataSet, point_tags: np.ndarray) -> np.ndarray:
    """Project point tags to cell tags by majority vote over nonzero point tags."""
    n_cells = int(mesh.n_cells)
    out = np.zeros(n_cells, dtype=np.int32)

    conn = getattr(mesh, "cell_connectivity", None)
    offs = getattr(mesh, "offset", None)
    if conn is not None and offs is not None:
        conn = np.asarray(conn, dtype=np.int64).reshape(-1)
        offs = np.asarray(offs, dtype=np.int64).reshape(-1)
        if offs.size == n_cells + 1:
            for cid in range(n_cells):
                ids = conn[offs[cid]:offs[cid + 1]]
                vals = point_tags[ids]
                nz = vals[vals > 0]
                if nz.size > 0:
                    out[cid] = int(np.bincount(nz).argmax())
            return out
        if offs.size == n_cells:
            for cid in range(n_cells):
                start = offs[cid]
                end = offs[cid + 1] if cid + 1 < n_cells else conn.size
                ids = conn[start:end]
                vals = point_tags[ids]
                nz = vals[vals > 0]
                if nz.size > 0:
                    out[cid] = int(np.bincount(nz).argmax())
            return out

    # Fallback for dataset types without connectivity arrays.
    for cid in range(n_cells):
        cell = mesh.get_cell(cid)
        ids = np.asarray(cell.point_ids, dtype=np.int64)
        if ids.size == 0:
            continue
        vals = point_tags[ids]
        nz = vals[vals > 0]
        if nz.size > 0:
            out[cid] = int(np.bincount(nz).argmax())
    return out


def main():
    ap = argparse.ArgumentParser(
        description="Anatomical segmentation and tagging for VTK meshes."
    )
    ap.add_argument("vtk", nargs="?", type=Path, help="Path to VTK legacy ASCII file")
    ap.add_argument(
        "--input",
        type=Path,
        default=None,
        help="Input mesh file path (same as positional 'vtk').",
    )
    ap.add_argument("--seed", type=int, default=0, help="RNG seed for subsampling (default: 0)")
    ap.add_argument(
        "--longitudinal-field",
        default=LONGITUDINAL_FIELD,
        help="Point data field for longitudinal coordinate.",
    )
    ap.add_argument(
        "--intraventricular-field",
        default="uvc_intraventricular",
        help="Point data field for LV/RV selection.",
    )
    ap.add_argument(
        "--rotational-field",
        default="uvc_rotational",
        help="Point data field for rotational coordinate (LV bands).",
    )

    # Anatomy-based circumferential segmentation (bullseye)
    ap.add_argument(
        "--ant-groove-theta",
        type=float,
        default=ANT_GROOVE_THETA,
        help=(
            "Anterior interventricular groove angle (radians) in the same convention as the rotational field. "
            "If not given, can be estimated from --ant-groove-field."
        ),
    )
    ap.add_argument(
        "--post-groove-theta",
        type=float,
        default=POST_GROOVE_THETA,
        help=(
            "Posterior interventricular groove angle (radians) in the same convention as the rotational field. "
            "If not given, can be estimated from --post-groove-field."
        ),
    )
    ap.add_argument(
        "--ant-groove-field",
        default=ANT_GROOVE_FIELD,
        help=(
            "Optional point-data mask/label field identifying the anterior interventricular groove points. "
            "Nonzero values are treated as True."
        ),
    )
    ap.add_argument(
        "--post-groove-field",
        default=POST_GROOVE_FIELD,
        help=(
            "Optional point-data mask/label field identifying the posterior interventricular groove points. "
            "Nonzero values are treated as True."
        ),
    )
    ap.add_argument(
        "--z-apical-apex",
        type=float,
        default=z_apical_apex,
        help="Longitudinal value for apex plane (used if apex cap is enabled in division config).",
    )
    ap.add_argument(
        "--z-apical-mid",
        type=float,
        default=z_apical_mid,
        help="Longitudinal value for apical-mid plane.",
    )
    ap.add_argument(
        "--z-mid-basal",
        type=float,
        default=z_mid_basal,
        help="Longitudinal value for mid-basal plane.",
    )
    ap.add_argument(
        "--z-cap-valves",
        type=float,
        default=z_cap_valves,
        help="Longitudinal value for cap-valves plane.",
    )
    ap.add_argument(
        "--longitudinal-eps",
        type=float,
        default=LONGITUDINAL_EPS,
        help="Tolerance for longitudinal value selection.",
    )
    ap.add_argument(
        "--plane-count",
        type=int,
        default=None,
        choices=[3, 4],
        help="Number of longitudinal planes to use (3 or 4).",
    )
    ap.add_argument(
        "--debug-steps",
        action="store_true",
        help="Enable all debug-step plots (groove diagnostics + each z-band divisions).",
    )
    ap.add_argument(
        "--ventricle",
        choices=["both", "lv", "rv"],
        default="both",
        help="Select which ventricle(s) to process/plot for centerlines and z divisions.",
    )
    ap.add_argument(
        "--draw-bullseye",
        action="store_true",
        help="Draw bullseye for the active division reference before mesh plots.",
    )
    ap.add_argument(
        "--draw-all-bullseyes",
        action="store_true",
        help="Draw all available bullseyes from division config before mesh plots.",
    )
    ap.add_argument(
        "--bullseye-outdir",
        default=None,
        help="Optional output directory to save bullseye images.",
    )
    ap.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Output mesh path with exported point-data tags (default: <input>_tagged.vtk).",
    )
    ap.add_argument(
        "--septal-click",
        type=float,
        nargs=3,
        default=SEPTAL_CLICK_XYZ,
        metavar=("X", "Y", "Z"),
        help=(
            "Fallback septal wall point (x y z). Used to orient groove ordering when RV "
            "reference is unavailable."
        ),
    )
    args = ap.parse_args()
    input_path = args.input if args.input is not None else args.vtk
    if input_path is None:
        ap.error("Provide an input mesh via positional 'vtk' or --input <file>.")
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    vent_suffix = "" if args.ventricle == "both" else f"_{args.ventricle}"
    out_name = (
        args.output.name
        if args.output is not None
        else f"{input_path.stem}_{DIVISION_REFERENCE}{vent_suffix}_tagged{input_path.suffix}"
    )
    output_path = OUTPUT_DIR / out_name

    if args.draw_all_bullseyes:
        flow("drawing all bullseye references")
        bullseye_outdir = args.bullseye_outdir or BULLSEYE_PLOTS_DIR
        draw_all_bullseyes(
            DIVISION_REFERENCES,
            output_dir=bullseye_outdir,
            show=True,
        )
    elif args.draw_bullseye and DIVISION_REFERENCE in DIVISION_REFERENCES:
        flow("drawing active bullseye reference")
        out_path = None
        out_dir = Path(args.bullseye_outdir) if args.bullseye_outdir else BULLSEYE_PLOTS_DIR
        out_dir.mkdir(parents=True, exist_ok=True)
        out_path = out_dir / f"{DIVISION_REFERENCE}_bullseye.png"
        draw_bullseye_reference(
            DIVISION_REFERENCE,
            DIVISION_REFERENCES[DIVISION_REFERENCE],
            out_path=str(out_path) if out_path else None,
        )
        plt.show()

    if True:
        flow("loading mesh and coordinate fields")
        flow(f"bullseye model: {DIVISION_REFERENCE}")
        flow(f"ventricle selection: {args.ventricle}")
        mesh = pv.read(str(input_path))
        apex_active = USE_APEX_CAP and (APEX_SEGMENT_ID is not None)
        z_min_exclusive = args.z_apical_apex if apex_active else None
        field = args.longitudinal_field
        if field not in mesh.point_data:
            raise KeyError(f"Missing point data field '{field}' in {input_path}")
        if args.intraventricular_field not in mesh.point_data:
            raise KeyError(
                f"Missing point data field '{args.intraventricular_field}' in {input_path}"
            )
        values = np.asarray(mesh.point_data[field]).reshape(-1)
        intrav = np.asarray(mesh.point_data[args.intraventricular_field]).reshape(-1)
        if args.rotational_field not in mesh.point_data:
            raise KeyError(
                f"Missing point data field '{args.rotational_field}' in {input_path}"
            )
        rotational = np.asarray(mesh.point_data[args.rotational_field]).reshape(-1)
        n_points = mesh.points.shape[0]
        # Exported point tag: 0 means "not assigned".
        tagged_segment_id = np.zeros(n_points, dtype=np.int32)

        flow("trying to read intraventricular coordinates for groove definitions")
        grooves = None
        groove_debug = {"method": "unknown", "source": "unknown"}
        try:
            grooves, groove_debug = groove_angles_debug_from_intraventricular_interface(
                mesh=mesh,
                rotational_theta=rotational,
                intraventricular=intrav,
                longitudinal_z=values,
                z_min_exclusive=z_min_exclusive,
                z_cap_max=args.z_cap_valves,
                rot_start=ROT_START,
                lv_tag=float(LV_INTRAVENTRICULAR_TAG),
                rv_tag=float(RV_INTRAVENTRICULAR_TAG),
            )
        except Exception:
            grooves = None

        if grooves is None:
            flow("did not find intraventricular groove interface, reading from config file")
            if args.ant_groove_theta is not None and args.post_groove_theta is not None:
                grooves = GrooveAngles(args.ant_groove_theta, args.post_groove_theta)
                groove_debug = {"method": "explicit_thetas", "source": "config_file"}
            else:
                if (
                    args.ant_groove_field
                    and args.post_groove_field
                    and args.ant_groove_field in mesh.point_data
                    and args.post_groove_field in mesh.point_data
                ):
                    ant_mask = np.asarray(mesh.point_data[args.ant_groove_field]).reshape(-1) != 0
                    post_mask = np.asarray(mesh.point_data[args.post_groove_field]).reshape(-1) != 0
                    cap_mask = values <= args.z_cap_valves
                    if z_min_exclusive is not None:
                        cap_mask = cap_mask & (values > z_min_exclusive)
                    ant_mask = ant_mask & cap_mask
                    post_mask = post_mask & cap_mask
                    ant_theta = groove_angles_from_point_mask(mesh.points, rotational, ant_mask)
                    post_theta = groove_angles_from_point_mask(mesh.points, rotational, post_mask)
                    grooves = GrooveAngles(ant_theta, post_theta)
                    groove_debug = {
                        "method": "mask_fields",
                        "source": "config_file",
                        "anterior_point_ids": np.flatnonzero(ant_mask),
                        "posterior_point_ids": np.flatnonzero(post_mask),
                    }

        if grooves is None:
            raise ValueError(
                "Anatomy-based bullseye needs groove angles. Provide either: "
                "--ant-groove-theta and --post-groove-theta (radians), or "
                "--ant-groove-field and --post-groove-field (point-data masks), "
                "or an intraventricular field that defines LV/RV interface."
            )
        groove_method = str(groove_debug.get("method", "unknown"))
        groove_source = str(groove_debug.get("source", "unknown"))
        uses_intraventricular = groove_source == "intraventricular_interface"
        flow(f"groove source uses intraventricular coordinates: {'yes' if uses_intraventricular else 'no'}")
        flow(f"groove source: {groove_source}")
        flow(f"groove source method: {groove_method}")

        flow("orienting grooves with rv or septal reference")
        # Re-orient groove ordering so septal arc is RV-facing if RV is available.
        cap_mask = values <= args.z_cap_valves
        if z_min_exclusive is not None:
            cap_mask = cap_mask & (values > z_min_exclusive)
        lv_mask = np.isclose(intrav, float(LV_INTRAVENTRICULAR_TAG), atol=1.0e-8) & cap_mask
        rv_mask = np.isclose(intrav, float(RV_INTRAVENTRICULAR_TAG), atol=1.0e-8) & cap_mask
        septal_theta_ref = None
        if np.any(lv_mask) and np.any(rv_mask):
            lv_center = mesh.points[lv_mask].mean(axis=0)
            rv_center = mesh.points[rv_mask].mean(axis=0)
            septal_theta_ref = theta_from_direction_reference(
                lv_points=mesh.points[lv_mask],
                lv_theta=rotational[lv_mask],
                lv_center=lv_center,
                direction_vec=(rv_center - lv_center),
            )
            if septal_theta_ref is not None:
                grooves = orient_grooves_with_septal_theta(grooves, septal_theta_ref)
        elif args.septal_click is not None and np.any(lv_mask):
            lv_center = mesh.points[lv_mask].mean(axis=0)
            click = np.asarray(args.septal_click, dtype=float)
            septal_theta_ref = theta_from_direction_reference(
                lv_points=mesh.points[lv_mask],
                lv_theta=rotational[lv_mask],
                lv_center=lv_center,
                direction_vec=(click - lv_center),
            )
            if septal_theta_ref is not None:
                grooves = orient_grooves_with_septal_theta(grooves, septal_theta_ref)

        if apex_active and args.z_apical_apex > args.z_apical_mid:
            raise ValueError("z_apical_apex must be <= z_apical_mid.")
        if args.z_apical_mid > args.z_mid_basal:
            raise ValueError("z_apical_mid must be <= z_mid_basal.")
        if args.z_mid_basal >= args.z_cap_valves:
            raise ValueError("z_cap_valves must be > z_mid_basal.")

        flow("extracting centerlines and longitudinal planes")
        targets = [args.z_apical_mid, args.z_mid_basal, args.z_cap_valves]
        if apex_active:
            apex_value = args.z_apical_apex
            targets = [apex_value] + [t for t in targets if t != apex_value]
        plane_count = args.plane_count
        if plane_count is None:
            plane_count = 4 if apex_active else 3
        targets = targets[:plane_count]

        plotter = pv.Plotter()
        plotter.add_mesh(mesh, color="white", opacity=0.2, show_edges=False)
        plotter_angles = pv.Plotter()
        plotter_angles.add_mesh(mesh, color="white", opacity=0.1, show_edges=False)
        rng = np.random.default_rng(args.seed)

        plane_centers = []
        legend_lines = []
        n_interface_components = plot_interface_components(
            plotter=plotter,
            mesh=mesh,
            intraventricular=intrav,
            lv_tag=float(LV_INTRAVENTRICULAR_TAG),
            rv_tag=float(RV_INTRAVENTRICULAR_TAG),
            point_size=8.0,
        )
        if n_interface_components > 0:
            legend_lines.append(f"Interface components: {n_interface_components}")
        vent_full_centers = {}
        seg_band_centers = None
        seg_bins = []
        area_counts = {}
        angle_division_points = 0
        angle_z_lines = []
        angle_theta_lines = []
        vent_specs_all = (
            ("LV", float(LV_INTRAVENTRICULAR_TAG), "blue", "red"),
            ("RV", float(RV_INTRAVENTRICULAR_TAG), "green", "orange"),
        )
        if args.ventricle == "lv":
            vent_specs = [vent_specs_all[0]]
        elif args.ventricle == "rv":
            vent_specs = [vent_specs_all[1]]
        else:
            vent_specs = list(vent_specs_all)
        seg_target_labels = {spec[0] for spec in vent_specs}

        for vent_label, vent_value, plane_color, point_color in vent_specs:
            vent_band_centers = []
            vent_centroids = []
            vent_targets = targets
            if vent_label == "RV":
                vent_targets = [t for t in targets if t != args.z_apical_apex]

            vent_centers, vent_mask, vent_legend = process_ventricle_centerline_and_slices(
                plotter=plotter,
                mesh=mesh,
                values=values,
                intrav=intrav,
                vent_label=vent_label,
                vent_value=vent_value,
                targets=vent_targets,
                longitudinal_eps=args.longitudinal_eps,
                apex_active=apex_active,
                z_apical_apex=args.z_apical_apex,
                plane_color=plane_color,
                point_color=point_color,
                seed=args.seed,
            )
            legend_lines.extend(vent_legend)
            if np.any(vent_mask):
                vent_full_centers[vent_label] = mesh.points[vent_mask].mean(axis=0)

            if len(vent_centers) >= 2:
                centers = np.array(vent_centers)
                centers = centers[np.argsort(centers[:, 2])]
                line = pv.PolyData()
                line.points = centers
                line.lines = np.hstack([[centers.shape[0]], np.arange(centers.shape[0])])
                plotter.add_mesh(line, color=plane_color, line_width=4)
                plotter.add_points(centers, color=plane_color, point_size=10)

            if vent_label in seg_target_labels and len(targets) >= 2:
                flow(f"computing {vent_label.lower()} ring bands and theta segment bins")
                sorted_targets = np.array(targets)
                sorted_targets = np.sort(sorted_targets)
                vent_mask = np.isclose(intrav, vent_value, atol=args.longitudinal_eps)
                vent_values = values[vent_mask]
                if vent_values.size == 0:
                    continue
                # Include the full apical region below the first target (e.g. below z_apical_mid).
                band_edges = np.concatenate(([float(vent_values.min())], sorted_targets))
                # Draw longitudinal division lines used by the segmentation.
                for z_cut in sorted_targets:
                    z_mask = np.isclose(values, z_cut, atol=args.longitudinal_eps) & vent_mask
                    if not np.any(z_mask):
                        continue
                    z_pts = mesh.points[z_mask]
                    if z_pts.shape[0] > PLOT_CAP_Z_DIVISION_POINTS:
                        idx = rng.choice(
                            z_pts.shape[0], size=PLOT_CAP_Z_DIVISION_POINTS, replace=False
                        )
                        z_pts = z_pts[idx]
                    angle_z_lines.append(z_pts.copy())
                    plotter_angles.add_points(z_pts, color=Z_DIVISION_COLOR, point_size=7)
                # Segment colors are resolved per segment id and active color scheme.

                for i in range(len(band_edges) - 1):
                    z_lo = band_edges[i]
                    z_hi = band_edges[i + 1]
                    if apex_active and np.isclose(z_hi, args.z_apical_apex):
                        long_name = "apex_to_apical_apex"
                    elif apex_active and np.isclose(z_lo, args.z_apical_apex) and np.isclose(z_hi, args.z_apical_mid):
                        long_name = "apical_apex_to_apical_mid"
                    elif np.isclose(z_hi, args.z_apical_mid):
                        long_name = "apex_to_apical_mid"
                    elif np.isclose(z_hi, args.z_mid_basal):
                        long_name = "apical_mid_to_mid_basal"
                    elif np.isclose(z_hi, args.z_cap_valves):
                        long_name = "mid_basal_to_cap_valves"
                    else:
                        long_name = f"long_{i+1}"
                    band_mask = (
                        (values >= z_lo - args.longitudinal_eps)
                        & (values <= z_hi + args.longitudinal_eps)
                        & vent_mask
                    )
                    if not np.any(band_mask):
                        continue
                    band_center = mesh.points[band_mask].mean(axis=0)
                    vent_band_centers.append(band_center)

                    # Determine ring for this longitudinal band
                    if z_hi <= args.z_apical_mid + args.longitudinal_eps:
                        ring = "apical"
                    elif z_hi <= args.z_mid_basal + args.longitudinal_eps:
                        ring = "mid"
                    else:
                        ring = "basal"

                    band_pts = mesh.points[band_mask]
                    band_theta = rotational[band_mask]
                    is_apex_cap_band = apex_active and (
                        z_hi <= (args.z_apical_apex + args.longitudinal_eps)
                    )
                    if is_apex_cap_band:
                        # Apex cap has no angular subdivision.
                        seg_ids = np.full(band_theta.shape, APEX_SEGMENT_ID, dtype=np.int32)
                    else:
                        seg_ids = segment_ids_from_theta_layout(
                            band_theta,
                            ring=ring,
                            grooves=grooves,
                            ring_layout=RING_LAYOUT,
                            ring_id_map=RING_ID_MAP,
                        )
                        # Draw circumferential boundary lines used by angle segmentation.
                        septal_parts, freewall_parts = RING_LAYOUT[ring]
                        boundaries = build_ring_boundaries(
                            grooves,
                            septal_parts=int(septal_parts),
                            freewall_parts=int(freewall_parts),
                        )
                        boundary_angles = sorted({float(v[0]) for v in boundaries})
                        for theta_b in boundary_angles:
                            bmask = circular_distance(band_theta, theta_b) <= ANGULAR_DIV_EPS
                            if not np.any(bmask):
                                continue
                            bpts = band_pts[bmask]
                            if bpts.shape[0] > PLOT_CAP_THETA_DIVISION_POINTS:
                                idx = rng.choice(
                                    bpts.shape[0],
                                    size=PLOT_CAP_THETA_DIVISION_POINTS,
                                    replace=False,
                                )
                                bpts = bpts[idx]
                            angle_theta_lines.append(bpts.copy())
                            plotter_angles.add_points(
                                bpts, color=THETA_DIVISION_COLOR, point_size=8
                            )

                    for seg_id in np.unique(seg_ids):
                        if seg_id < 0:
                            continue
                        seg_mask = seg_ids == seg_id
                        if not np.any(seg_mask):
                            continue
                        pts_seg = band_pts[seg_mask]
                        seg_center = pts_seg.mean(axis=0)
                        area_name = f"{vent_label.lower()}_{long_name}_seg{int(seg_id)}"
                        area_counts[area_name] = int(pts_seg.shape[0])
                        vent_centroids.append(seg_center)

                        color = color_for_segment(int(seg_id), COLOR_SCHEME_DEFAULT)
                        pts_plot = pts_seg
                        if pts_plot.shape[0] > PLOT_CAP_SEGMENT_POINTS:
                            idx = rng.choice(
                                pts_plot.shape[0], size=PLOT_CAP_SEGMENT_POINTS, replace=False
                            )
                            pts_plot = pts_plot[idx]
                        plotter_angles.add_points(pts_plot, color=color, point_size=5)
                        angle_division_points += int(pts_plot.shape[0])
                        seg_bins.append((pts_plot.copy(), seg_center, int(seg_id), area_name))
                        band_ids = np.flatnonzero(band_mask)
                        point_ids = band_ids[seg_mask]
                        if point_ids.size > 0:
                            tagged_segment_id[point_ids] = int(seg_id)

            if vent_band_centers:
                band_centers = np.array(vent_band_centers)
                plotter.add_points(band_centers, color="white", point_size=14)
                seg_band_centers = band_centers
            if vent_centroids:
                centroids = np.array(vent_centroids)
                plotter.add_points(centroids, color="white", point_size=12)
                seg_band_centers = centroids

        if legend_lines:
            plotter.add_text(
                "\n".join(legend_lines),
                position="upper_left",
                font_size=12,
            )
        # 1) First: centerlines + plane divisions.
        flow("showing centerlines and planes")
        plotter.show()
        # 2) Then: groove diagnostics (if requested).
        if args.debug_steps:
            flow("plotting groove diagnostics")
            if args.ventricle == "rv":
                dbg_primary = float(RV_INTRAVENTRICULAR_TAG)
                dbg_secondary = float(LV_INTRAVENTRICULAR_TAG)
            else:
                dbg_primary = float(LV_INTRAVENTRICULAR_TAG)
                dbg_secondary = float(RV_INTRAVENTRICULAR_TAG)
            plot_groove_diagnostics(
                mesh=mesh,
                rotational=rotational,
                intrav=intrav,
                grooves=grooves,
                groove_debug=groove_debug,
                z_values=values,
                z_min_exclusive=z_min_exclusive,
                z_cap_max=args.z_cap_valves,
                lv_tag=dbg_primary,
                rv_tag=dbg_secondary,
            )
        # 3) Then: angle-based divisions (if available/requested).
        if angle_division_points > 0:
            flow("plotting angle-based divisions")
            available_modes = ["segment"]
            if ARTERY_COLORS and SEGMENT_ARTERY_MAP:
                available_modes.append("artery")
            mode0 = COLOR_SCHEME_DEFAULT if COLOR_SCHEME_DEFAULT in available_modes else "segment"
            state = {"mode_idx": available_modes.index(mode0), "bg_idx": 0}
            bg_colors = ["white", "black", "dimgray", "midnightblue"]

            def _render_angle_window():
                mode = available_modes[state["mode_idx"]]
                bg = bg_colors[state["bg_idx"] % len(bg_colors)]
                _clear_active_renderer(plotter_angles)
                plotter_angles.set_background(bg)
                plotter_angles.add_mesh(mesh, color="white", opacity=0.1, show_edges=False)
                for pts in angle_z_lines:
                    plotter_angles.add_points(pts, color=Z_DIVISION_COLOR, point_size=7)
                for pts in angle_theta_lines:
                    plotter_angles.add_points(pts, color=THETA_DIVISION_COLOR, point_size=8)
                for pts, _center, seg_id, _name in seg_bins:
                    color = color_for_segment(int(seg_id), mode)
                    plotter_angles.add_points(pts, color=color, point_size=5)
                plotter_angles.add_text(
                    "Angle-Based Ventricle Division",
                    position="upper_left",
                    font_size=11,
                )

            def _toggle_mode():
                state["mode_idx"] = (state["mode_idx"] + 1) % len(available_modes)
                _render_angle_window()

            def _toggle_bg():
                state["bg_idx"] = (state["bg_idx"] + 1) % len(bg_colors)
                _render_angle_window()

            if DIVISION_REFERENCE == "aha17":
                plotter_angles.add_key_event("t", _toggle_mode)
                flow("angle plot toggle: press 't' to change color scheme")
            else:
                plotter_angles.add_key_event("c", _toggle_mode)
                plotter_angles.add_key_event("b", _toggle_bg)
            _render_angle_window()
            if DIVISION_REFERENCE in DIVISION_REFERENCES:
                if DIVISION_REFERENCE == "aha17":
                    flow("opening bullseye reference images")
                else:
                    flow("opening bullseye reference image")
                show_reference_figures(
                    DIVISION_REFERENCE,
                    DIVISION_REFERENCES[DIVISION_REFERENCE],
                )
            plotter_angles.show()
        if args.debug_steps and seg_bins:
            flow("plotting each z-band division")
            plot_each_division(mesh, seg_bins)

        flow("exporting vtk cell-data tags")
        out_mesh = mesh.copy(deep=True)
        tagged_cell_id = point_tags_to_cell_tags(out_mesh, tagged_segment_id)
        out_mesh.cell_data["anatomical_tag"] = tagged_cell_id
        out_mesh.save(str(output_path))
        flow(f"saved tagged vtk: {output_path}")
        return


if __name__ == "__main__":
    main()
