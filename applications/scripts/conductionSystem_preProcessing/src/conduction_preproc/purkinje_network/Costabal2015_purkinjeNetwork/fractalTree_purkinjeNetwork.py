import argparse
import importlib.util
import sys
from pathlib import Path

# Force Python to prefer the local folder (the one that contains ./fractal_tree/)
HERE = Path(__file__).resolve().parent
INPUT_DIR = HERE / "input"
OUTPUT_DIR = HERE / "outputs"
INPUT_DIR.mkdir(exist_ok=True)
OUTPUT_DIR.mkdir(exist_ok=True)
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE / "src" / "src"))

import inspect

import numpy as np
import logging
import meshio
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

from fractal_tree import generate_fractal_tree, FractalTreeParameters, Mesh
from fractal_tree.branch import Branch


logging.basicConfig(level=logging.INFO)
np.random.seed(1234)




# ------------------------- helpers -------------------------
def _parse_vector(value):
    if isinstance(value, (list, tuple, np.ndarray)):
        if len(value) != 3:
            raise ValueError(f"Expected 3 values, got {len(value)}")
        return np.array([float(v) for v in value])
    parts = [p.strip() for p in str(value).split(",") if p.strip()]
    if len(parts) != 3:
        raise ValueError(f"Expected 3 values, got {len(parts)} in '{value}'")
    return np.array([float(p) for p in parts])


def load_config(path: Path):
    if not path.exists():
        raise FileNotFoundError(path)
    spec = importlib.util.spec_from_file_location("fractalTree_config", path)
    module = importlib.util.module_from_spec(spec)
    assert spec and spec.loader
    spec.loader.exec_module(module)
    return module


def _line_to_dir(start, end):
    v0 = _parse_vector(start)
    v1 = _parse_vector(end)
    d = v1 - v0
    n = np.linalg.norm(d)
    if n == 0:
        raise ValueError("Line start/end are identical; cannot build direction.")
    return d / n


def _resolve_seed(line_cfg, name, fallback):
    if line_cfg is None:
        return fallback
    if hasattr(line_cfg, name):
        return _parse_vector(getattr(line_cfg, name))
    if hasattr(line_cfg, "base_seed"):
        return _parse_vector(getattr(line_cfg, "base_seed"))
    return fallback

def get_surface_connectivity(msh, physical_tag, cell_type="triangle"):
    """Robustly extract ONLY the triangles whose gmsh:physical tag matches physical_tag."""
    conns = []
    for block, phys in zip(msh.cells, msh.cell_data["gmsh:physical"]):
        if block.type != cell_type:
            continue
        phys = np.asarray(phys).reshape(-1)
        mask = phys == physical_tag
        if np.any(mask):
            conns.append(block.data[mask])
    if not conns:
        raise RuntimeError(f"No {cell_type} cells found for physical tag {physical_tag}")
    return np.vstack(conns)


def get_surface_connectivity_from_cell_data(msh, tag_name, cell_type="triangle"):
    conns = []
    if tag_name not in msh.cell_data:
        raise RuntimeError(f"Missing cell_data tag '{tag_name}'.")
    for block, tags in zip(msh.cells, msh.cell_data[tag_name]):
        if block.type != cell_type:
            continue
        tag_vals = np.asarray(tags).reshape(-1)
        mask = tag_vals != 0
        if np.any(mask):
            conns.append(block.data[mask])
    if not conns:
        raise RuntimeError(f"No {cell_type} cells found for tag {tag_name}")
    return np.vstack(conns)


def get_all_triangles(msh):
    conns = []
    for block in msh.cells:
        if block.type != "triangle":
            continue
        conns.append(block.data)
    if not conns:
        raise RuntimeError("No triangle cells found in base mesh.")
    return np.vstack(conns)


def meshio_from_pyvista(surface):
    if surface.n_cells == 0:
        raise RuntimeError("Input surface has no cells.")
    if surface.faces.size == 0:
        surface = surface.extract_surface()
    if surface.faces.size == 0:
        raise RuntimeError("Unable to extract surface faces from input mesh.")

    surface = surface.triangulate()
    faces = surface.faces.reshape(-1, 4)
    tri = faces[:, 1:4]
    cell_data = {name: [np.asarray(arr)] for name, arr in surface.cell_data.items()}
    return meshio.Mesh(points=surface.points, cells=[("triangle", tri)], cell_data=cell_data)


def surface_points(verts, conn):
    ids = np.unique(conn.ravel())
    return verts[ids]


def plot_debug(surfaces, seed_lv, seed_rv, init_dir_lv, init_dir_rv, max_points):
    fig = plt.figure(figsize=(10, 9))
    ax = fig.add_subplot(111, projection="3d")

    # Surfaces
    rng = np.random.default_rng(1234)
    for name, (s_verts, conn) in surfaces.items():
        pts = surface_points(s_verts, conn)
        if max_points > 0 and pts.shape[0] > max_points:
            idx = rng.choice(pts.shape[0], size=max_points, replace=False)
            pts = pts[idx]
        ax.scatter(pts[:, 0], pts[:, 1], pts[:, 2], s=1, alpha=0.15, label=name)

    # Seeds
    ax.scatter(*seed_lv, s=150, label="seed LV")
    ax.scatter(*seed_rv, s=150, label="seed RV")

    # Initial direction vectors at seeds (these are what go into FractalTreeParameters)
    ax.quiver(seed_lv[0], seed_lv[1], seed_lv[2],
              init_dir_lv[0], init_dir_lv[1], init_dir_lv[2],
              length=0.9, normalize=True)
    ax.text(seed_lv[0], seed_lv[1], seed_lv[2], "LV init_dir")

    ax.quiver(seed_rv[0], seed_rv[1], seed_rv[2],
              init_dir_rv[0], init_dir_rv[1], init_dir_rv[2],
              length=0.9, normalize=True)
    ax.text(seed_rv[0], seed_rv[1], seed_rv[2], "RV init_dir")

    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")
    ax.legend(loc="upper right")
    ax.set_title("Tagged surfaces + demo-style seeds + init direction vectors")
    ax.set_box_aspect([1, 1, 1])
    plt.tight_layout()
    plt.show()


def view_final_trees(lv_path: Path, rv_path: Path) -> None:
    import shutil
    import subprocess

    paraview_bin = shutil.which("paraview")
    if paraview_bin:
        subprocess.Popen([paraview_bin, str(lv_path), str(rv_path)])
        return

    try:
        subprocess.Popen(["open", "-a", "ParaView", "--args", str(lv_path), str(rv_path)])
    except Exception as exc:
        print(f"Viewer skipped (ParaView not found): {exc}")




# ------------------------- main -------------------------

parser = argparse.ArgumentParser(
    description="Generate a Costabal2015-style Purkinje fractal tree network."
)
parser.add_argument(
    "--input",
    default=str(Path(__file__).resolve().parents[4] / "inputs" / "meshes" / "biv_ellipsoid.msh"),
    help="Path to the input mesh file (defaults to inputs/meshes/biv_ellipsoid.msh).",
)
parser.add_argument(
    "--config",
    default=str(
        Path(__file__).resolve().parents[4]
        / "configs"
        / "purkinjeFractalTree_config.py"
    ),
    help="Path to combined config file (defaults to configs/purkinjeFractalTree_config.py).",
)
parser.add_argument(
    "--debug-plot",
    action="store_true",
    help="Plot tagged surfaces and seed points before generation.",
)
parser.add_argument(
    "--no-view-final",
    action="store_true",
    help="Skip opening the final LV/RV VTUs in an interactive viewer.",
)
args = parser.parse_args()

# Load config
cfg = load_config(Path(args.config))
Branch.set_debug_grow(bool(getattr(cfg, "debug_grow", False)))
Branch.set_debug_dir(str(getattr(cfg, "debug_dir", HERE / "debug_attempts")))
Branch.set_debug_every(int(getattr(cfg, "debug_every", 50)))

# 1) Get/create file (demo mesh)
path = Path(args.input)
if not path.is_file():
    import cardiac_geometries
    cardiac_geometries.gmsh.biv_ellipsoid(path, char_length=0.1)

# 2) Read mesh
try:
    msh = meshio.read(path)
except Exception as exc:
    import pyvista as pv

    surface = pv.read(path)
    msh = meshio_from_pyvista(surface)
    print(f"Loaded mesh via PyVista fallback ({path.name}): {exc}")
verts = msh.points

if msh.field_data:
    print("\nfield_data (name -> [tag_id, dim]):")
    for k, v in msh.field_data.items():
        print(f"  {k}: {v}")
else:
    print("\ncell_data fields:")
    for k in msh.cell_data.keys():
        print(f"  {k}")

# 3) Extract surfaces
default_tag_names = {
    "ENDO_LV": "ENDO_LV",
    "ENDO_RV": "ENDO_RV",
    "EPI": "EPI",
}
tag_names = default_tag_names
if hasattr(cfg, "SURFACE_TAGS"):
    tag_names = {**default_tag_names, **getattr(cfg, "SURFACE_TAGS")}

print("Note: Only ENDO_LV and ENDO_RV are required; EPI is used for plotting only.")

def prune_surface_vertices(points, conn, label):
    used = np.unique(conn)
    if used.size == points.shape[0]:
        return points, conn
    remap = -np.ones(points.shape[0], dtype=int)
    remap[used] = np.arange(used.size)
    pruned_points = points[used]
    pruned_conn = remap[conn]
    print(
        f"{label}: pruned {points.shape[0] - used.size} unused points "
        f"(kept {used.size})."
    )
    return pruned_points, pruned_conn


surfaces = {}
for name, required in (("ENDO_LV", True), ("ENDO_RV", True), ("EPI", False)):
    tag_key = tag_names.get(name, name)
    has_field = msh.field_data and tag_key in msh.field_data
    has_cell = tag_key in msh.cell_data
    if not (has_field or has_cell):
        if required:
            raise RuntimeError(f"Missing required surface tag '{tag_key}' for {name}.")
        print(f"Skipping optional surface {name} ({tag_key}): not found.")
        continue
    if has_field:
        tag_id = msh.field_data[tag_key][0]
        conn = get_surface_connectivity(msh, tag_id, cell_type="triangle")
    else:
        conn = get_surface_connectivity_from_cell_data(msh, tag_key, cell_type="triangle")
    surf_points, surf_conn = prune_surface_vertices(verts, conn, name)
    surfaces[name] = (surf_points, surf_conn)
    print(f"{name} ({tag_key}): triangles = {surf_conn.shape[0]}")

lv_verts, conn_lv = surfaces["ENDO_LV"]
rv_verts, conn_rv = surfaces["ENDO_RV"]

# -------------------------
# Config-driven seeds + initial directions
# -------------------------
init_node_lv = _parse_vector(getattr(cfg, "lv_seed", "0.0,1.0,0.0"))
idx_lv = np.linalg.norm(lv_verts - init_node_lv, axis=1).argmin()
seed_lv = lv_verts[idx_lv]

init_node_rv = _parse_vector(getattr(cfg, "rv_seed", "0.0,1.5,0.15"))
idx_rv = np.linalg.norm(rv_verts - init_node_rv, axis=1).argmin()
seed_rv = rv_verts[idx_rv]

# Initial direction vectors (can be overridden by line-config)
init_node_lv = _resolve_seed(cfg, "lv_seed", init_node_lv)
init_node_rv = _resolve_seed(cfg, "rv_seed", init_node_rv)

idx_lv = np.linalg.norm(lv_verts - init_node_lv, axis=1).argmin()
seed_lv = lv_verts[idx_lv]
idx_rv = np.linalg.norm(rv_verts - init_node_rv, axis=1).argmin()
seed_rv = rv_verts[idx_rv]

init_dir_lv = _parse_vector(getattr(cfg, "lv_init_dir", "1.0,0.0,0.0"))
init_dir_rv = _parse_vector(getattr(cfg, "rv_init_dir", "1.0,0.0,0.0"))

if hasattr(cfg, "lv_line_start") and hasattr(cfg, "lv_line_end"):
    init_dir_lv = _line_to_dir(cfg.lv_line_start, cfg.lv_line_end)
elif hasattr(cfg, "lv_line_end"):
    init_dir_lv = _line_to_dir(init_node_lv, cfg.lv_line_end)

if hasattr(cfg, "rv_line_start") and hasattr(cfg, "rv_line_end"):
    init_dir_rv = _line_to_dir(cfg.rv_line_start, cfg.rv_line_end)
elif hasattr(cfg, "rv_line_end"):
    init_dir_rv = _line_to_dir(init_node_rv, cfg.rv_line_end)

print("\nSeeds + init dirs (demo-style):")
print("seed_lv:", seed_lv)
print("seed_rv:", seed_rv)
print("init_dir_lv:", init_dir_lv)
print("init_dir_rv:", init_dir_rv)
print("lv_seed_source:", init_node_lv)
print("rv_seed_source:", init_node_rv)
print("lv_seed_mesh_idx:", idx_lv)
print("rv_seed_mesh_idx:", idx_rv)
print("lv_triangles:", conn_lv.shape[0])
print("rv_triangles:", conn_rv.shape[0])

# Plot everything so you can verify visually before generation
if args.debug_plot:
    max_plot_points = int(getattr(cfg, "plot_max_points", 2000))
    plot_debug(surfaces, seed_lv, seed_rv, init_dir_lv, init_dir_rv, max_plot_points)
    raise SystemExit(0)

# 4) Build fractal-tree meshes
mesh_lv = Mesh(verts=lv_verts, connectivity=conn_lv, init_node=seed_lv)
mesh_rv = Mesh(verts=rv_verts, connectivity=conn_rv, init_node=seed_rv)

# 5) Parameters (match the demo look: straight line + branching)
param_lv = FractalTreeParameters(
    filename=str(OUTPUT_DIR / getattr(cfg, "lv_filename", "biv-line-lv")),
    init_length=float(getattr(cfg, "lv_init_length", 3.0)),
    N_it=int(getattr(cfg, "lv_N_it", 15)),
    length=float(getattr(cfg, "lv_length", 0.25)),
    initial_direction=init_dir_lv,

)
param_lv.l_segment = float(getattr(cfg, "lv_l_segment", 0.01))
param_lv.branch_angle = float(getattr(cfg, "lv_branch_angle", 0.15))
param_lv.repulsitivity = float(getattr(cfg, "lv_repulsitivity", 0.1))
param_lv.save_frames = bool(getattr(cfg, "lv_save_frames", True))
print(
    "LV params: length=%g init_length=%g l_segment=%g num_segments=%d"
    % (
        param_lv.length,
        param_lv.init_length,
        param_lv.l_segment,
        int(param_lv.length / param_lv.l_segment),
    )
)


param_rv = FractalTreeParameters(
    filename=str(OUTPUT_DIR / getattr(cfg, "rv_filename", "biv-line-rv")),
    init_length=float(getattr(cfg, "rv_init_length", 4.0)),
    N_it=int(getattr(cfg, "rv_N_it", 25)),
    length=float(getattr(cfg, "rv_length", 0.15)),
    initial_direction=init_dir_rv,
)
param_rv.l_segment = float(getattr(cfg, "rv_l_segment", 0.01))
param_rv.branch_angle = float(getattr(cfg, "rv_branch_angle", 0.15))
param_rv.repulsitivity = float(getattr(cfg, "rv_repulsitivity", 0.1))
param_rv.save_frames = bool(getattr(cfg, "rv_save_frames", True))
print(
    "RV params: length=%g init_length=%g l_segment=%g num_segments=%d"
    % (
        param_rv.length,
        param_rv.init_length,
        param_rv.l_segment,
        int(param_rv.length / param_rv.l_segment),
    )
)

# 6) Generate trees
generate_fractal_tree(mesh_lv, param_lv)
generate_fractal_tree(mesh_rv, param_rv)

input_stem = Path(args.input).stem
geom_out = OUTPUT_DIR / f"{input_stem}_geometry.vtu"   # ParaView-friendly
meshio.write(geom_out, msh)
print(f"Wrote geometry mesh: {geom_out.resolve()}")

frames_root = Path(__file__).resolve().parents[4] / "outputs" / "videoEditorPurkinje"
frames_lv_dir = frames_root / f"frames_{Path(param_lv.filename).name}"
frames_rv_dir = frames_root / f"frames_{Path(param_rv.filename).name}"

print("\nOutputs:")
print(f"  Tree files (LV/RV): {OUTPUT_DIR.resolve()}")
if param_lv.save_frames:
    print(f"  LV frames: {frames_lv_dir.resolve()}")
if param_rv.save_frames:
    print(f"  RV frames: {frames_rv_dir.resolve()}")
print(f"  Geometry mesh: {geom_out.resolve()}")

print("\nDone. Check biv-line-lv* and biv-line-rv* outputs.")

if not args.no_view_final:
    lv_vtu = OUTPUT_DIR / f"{getattr(cfg, 'lv_filename', 'biv-line-lv')}.vtu"
    rv_vtu = OUTPUT_DIR / f"{getattr(cfg, 'rv_filename', 'biv-line-rv')}.vtu"
    view_final_trees(lv_vtu, rv_vtu)
