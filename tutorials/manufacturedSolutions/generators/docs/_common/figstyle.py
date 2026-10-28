"""Shared style and helpers of the figure scripts.

The conventions follow the paper's figures: one sans-serif family, 12 pt labels,
10 pt ticks, a fixed colour per quantity, light-grey geometry with a cut-away,
dark-red rods for fibre directions, blue dots for cube electrodes and red
triangles for named electrodes.  Scalar fields use one perceptually uniform
colormap (viridis) with a colourbar.  Fields are rendered with pyvista
off-screen into an image and placed in a matplotlib panel, so colourbars and
symbols are matplotlib mathtext.  Figures show the problem only: geometry,
exact fields, conductivity frames and sources.
"""

import argparse
import importlib.util
import os
import re
import shutil
import sys
import tempfile

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pyvista as pv
from matplotlib import cm
from matplotlib.colors import Normalize
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator
from mpl_toolkits.mplot3d import proj3d

pv.OFF_SCREEN = True

HERE = os.path.dirname(os.path.abspath(__file__))
GENERATORS = os.path.abspath(os.path.join(HERE, "..", ".."))
REFERENCE_DIR = os.path.join(GENERATORS, "reference")

CMAP = "viridis"

# colour per quantity (line plots of the quantity)
VM = "#1f77b4"
VM_LIGHT = "#9ecae1"
PSI = "#8c564b"
PSI_LIGHT = "#c49c94"
PSI_DARK = "#5e3a33"
ECG = "#9467bd"
ECG_LIGHT = "#c5b0d5"
ECG_DARK = "#5b3f7c"
NEUTRAL = "#7f7f7f"

# geometry
MESH_FACE = "#e3ebf3"
MESH_EDGE = "#5d6d7e"
FIBRE = "#c0392b"
FIBRE_DEPTH = ((0.96, 0.72, 0.69), (0.57, 0.17, 0.13))   # endocardium light, epicardium dark
SHEET = "#2c3e50"
NORMAL = "#95a5a6"
GRAPH = "#d35400"
PIN = "#d35400"
KERNEL = "#16a085"
DOMAIN = "#555555"
COHORT = "#1f77b4"
NAMED = "#b62c35"
GREY = "#808080"
DPI = 300
DRAFT_DPI = 80

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
    "mathtext.fontset": "dejavusans",
    "font.size": 10,
    "axes.titlesize": 12,
    "axes.labelsize": 12,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "legend.fontsize": 10,
    "axes.linewidth": 0.8,
    "lines.linewidth": 1.5,
    "figure.dpi": 100,
    "savefig.dpi": DPI,
})


class Settings:
    """Resolution of one run: final or draft."""

    def __init__(self, draft=False):
        self.draft = draft
        self.dpi = DRAFT_DPI if draft else DPI
        self.window = (360, 320) if draft else (1100, 960)
        self.n_cube = 24 if draft else 72
        self.n_shell = (6, 24, 36) if draft else (14, 90, 140)
        self.ecg_order = 12 if draft else 48
        self.n_time = 12 if draft else 90


def load_reference(name):
    """The emitted NumPy module reference/<name>.py (numpy only)."""
    path = os.path.join(REFERENCE_DIR, name + ".py")
    spec = importlib.util.spec_from_file_location("ref_" + name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def read_electrodes(path):
    """{name: position} of the electrodePositions block of an electroProperties file."""
    with open(path) as fh:
        text = fh.read()
    block = re.search(r"electrodePositions\s*\{(.*?)\n\s*\}", text, re.S)
    if block is None:
        return {}
    body = re.sub(r"//[^\n]*", "", block.group(1))
    return {m[1]: np.array([float(v) for v in m[2].split()])
            for m in re.finditer(r"(\w+)\s*\(([^)]*)\)\s*;", body)}


BUILTIN_CUBE_ELECTRODES = [(-0.5, 0.5, 0.5), (1.5, 0.5, 0.5), (1.2, 0.23, 0.61),
                           (1.35, 0.74, 0.28), (1.55, 0.41, 0.83)]
BUILTIN_SHELL_ELECTRODES = {
    "E1": (0.4, -0.8, 0.8), "E2": (0.6, 1.2, 0.4), "E3": (-0.96, 0.0, -1.28),
    "E4": (-0.6, 1.2, -1.2), "E5": (1.2, -0.96, -1.28)}


def tutorial_electrodes(case_dir):
    """Cube electrodes of a tutorial's electroProperties; the built-in set when absent."""
    path = os.path.join(GENERATORS, "..", case_dir, "constant", "electroProperties")
    found = read_electrodes(path) if os.path.exists(path) else {}
    if not found:
        print("note: %s not found; using the built-in cube electrode set" % path)
    return [tuple(v) for v in found.values()] or list(BUILTIN_CUBE_ELECTRODES)


def shell_electrodes(case_dir):
    """E1-E5 of a tutorial's sphericalShell case; the built-in set when absent."""
    path = os.path.join(GENERATORS, "..", case_dir, "sphericalShell", "constant",
                        "electroProperties")
    found = read_electrodes(path) if os.path.exists(path) else {}
    if not found:
        print("note: %s not found; using the built-in shell electrode set" % path)
    return found or {k: np.array(v) for k, v in BUILTIN_SHELL_ELECTRODES.items()}


def _tidy3d(ax, lim, azim=-58, elev=18, zoom=1.35):
    for set_lim, (lo, hi) in zip((ax.set_xlim, ax.set_ylim, ax.set_zlim), lim):
        set_lim(lo, hi)
    ax.set_box_aspect((1, 1, 1), zoom=zoom)
    ax.set_axis_off()
    ax.view_init(elev=elev, azim=azim)


def cube_electrode_axes(ax, electrodes, title=None):
    """Wireframe unit cube with the electrodes as small blue dots (matplotlib 3D)."""
    c = np.array([[x, y, z] for x in (0, 1) for y in (0, 1) for z in (0, 1)])
    for i, a in enumerate(c):
        for b in c[i + 1:]:
            if np.abs(a - b).sum() == 1:
                ax.plot(*np.array([a, b]).T, color=DOMAIN, lw=1.0, zorder=5)
    _tidy3d(ax, [(-0.75, 1.75)] * 3)
    ax.scatter(*np.array(electrodes).T, s=6, c=COHORT, alpha=0.7, edgecolors="none",
               depthshade=False, zorder=10)
    if title:
        ax.set_title(title, loc="left")


def shell_electrode_axes(ax, electrodes, r1=0.5, r2=1.0, title=None):
    """Wireframe shell spheres with named electrodes as labelled red triangles."""
    t = np.linspace(0, 2 * np.pi, 120)
    for r, colour, lw, n in ((r2, "#bbbbbb", 0.5, 13), (r1, DOMAIN, 0.6, 9)):
        for th in np.linspace(0, np.pi, n)[1:-1]:
            ax.plot(r * np.sin(th) * np.cos(t), r * np.sin(th) * np.sin(t),
                    np.full_like(t, r * np.cos(th)), color=colour, lw=lw, zorder=4)
        for ph in np.linspace(0, np.pi, n // 2, endpoint=False):
            ax.plot(r * np.sin(t) * np.cos(ph), r * np.sin(t) * np.sin(ph),
                    r * np.cos(t), color=colour, lw=lw, zorder=4)
    _tidy3d(ax, [(-1.6, 1.6)] * 3)
    xyz = np.array(list(electrodes.values()))
    ax.scatter(*xyz.T, s=26, marker="^", c=NAMED, edgecolors="k", linewidths=0.4,
               depthshade=False, zorder=20)
    ax.figure.canvas.draw()
    to_px = lambda q: ax.transData.transform(proj3d.proj_transform(*q, ax.get_proj())[:2])
    c2 = to_px((0, 0, 0))
    for name, p in electrodes.items():
        d = to_px(p) - c2
        d = d / (np.linalg.norm(d) or 1.0)
        x2, y2, _ = proj3d.proj_transform(*p, ax.get_proj())
        ax.annotate(name, xy=(x2, y2), xytext=tuple(11 * d), textcoords="offset points",
                    ha="center", va="center", fontsize=10, color=NAMED, zorder=30)
    if title:
        ax.set_title(title, loc="left")


# ---------------------------------------------------------------- pyvista ----

def plotter(settings):
    pl = pv.Plotter(off_screen=True, window_size=list(settings.window))
    pl.set_background("white")
    return pl


def cube_grid(n, lo=0.0, hi=1.0):
    """Point grid of the box [lo, hi]^3 with n points per axis."""
    return pv.ImageData(dimensions=(n, n, n), spacing=((hi - lo) / (n - 1),) * 3,
                        origin=(lo, lo, lo))


def point_array(grid):
    return grid.points[:, 0], grid.points[:, 1], grid.points[:, 2]


def camera(pl, centre, distance, direction=(1.0, -1.15, 0.8)):
    """Common oblique view of an object centred at `centre`."""
    d = np.asarray(direction, float)
    d = d / np.linalg.norm(d)
    pl.camera_position = [tuple(np.asarray(centre) + distance * d), tuple(centre),
                          (0, 0, 1)]


def add_box(pl, lo=0.0, hi=1.0, color="black"):
    pl.add_mesh(pv.Box(bounds=(lo, hi, lo, hi, lo, hi)), style="wireframe",
                color=color, line_width=2)


def add_axes_arrows(pl, origin=(0, 0, 0), length=0.35, color="#9a9a9a"):
    for i in range(3):
        v = np.zeros(3)
        v[i] = 1.0
        pl.add_mesh(pv.Arrow(start=origin, direction=v, scale=length, tip_length=0.25,
                             tip_radius=0.06, shaft_radius=0.025), color=color)


def render(pl, crop=True):
    """Image array of the plotter; the white margin is cropped."""
    img = pl.screenshot(return_img=True)
    pl.close()
    if crop:
        mask = np.any(img < 250, axis=2)
        rows, cols = np.where(mask)
        if rows.size:
            pad = 6
            img = img[max(rows.min() - pad, 0):rows.max() + pad,
                      max(cols.min() - pad, 0):cols.max() + pad]
    return img


def panel(fig, ax, img, title=None, label=None, clim=None, cmap=CMAP, tag=None):
    """Place a rendered image; with `label` and `clim` add a mathtext colourbar."""
    ax.imshow(img)
    ax.set_axis_off()
    if title:
        ax.set_title(title, pad=3)
    if tag:
        ax.text(0.0, 1.0, tag, transform=ax.transAxes, va="top", ha="left",
                fontweight="bold")
    if label is not None and clim is not None:
        sm = cm.ScalarMappable(norm=Normalize(*clim), cmap=cmap)
        cb = fig.colorbar(sm, ax=ax, orientation="horizontal", fraction=0.05,
                          pad=0.03, shrink=0.8)
        cb.set_label(label)
        cb.locator = MaxNLocator(5)
        cb.update_ticks()
        cb.ax.tick_params(labelsize=10)
    return ax


def legend(ax, entries, ncol=3, inside=None):
    """Legend of (label, colour[, marker]) entries; `inside` is an axes location
    (default: below the axes).  A marker draws a dot instead of a bar."""
    handles = []
    for entry in entries:
        label, color = entry[0], entry[1]
        if len(entry) > 2:
            handles.append(Line2D([0], [0], color=color, marker=entry[2], ls="",
                                  ms=7, label=label))
        else:
            handles.append(Line2D([0], [0], color=color, lw=3, label=label))
    if inside:
        ax.legend(handles=handles, loc=inside, ncol=ncol, frameon=False, fontsize=9,
                  handlelength=1.2)
    else:
        ax.legend(handles=handles, loc="upper center", ncol=ncol, frameon=False,
                  fontsize=10, bbox_to_anchor=(0.5, 0.0), handlelength=1.2,
                  columnspacing=1.0)


def colourbar_2d(fig, ax, mappable, label):
    cb = fig.colorbar(mappable, ax=ax, orientation="vertical", fraction=0.05,
                      pad=0.04, shrink=0.8)
    cb.set_label(label)
    cb.ax.tick_params(labelsize=10)
    return cb


def symmetric(values):
    m = float(np.nanmax(np.abs(values)))
    return (-m, m)


def span(values):
    return (float(np.nanmin(values)), float(np.nanmax(values)))


def save(fig, outdir, name, settings):
    os.makedirs(outdir, exist_ok=True)
    path = os.path.join(outdir, name)
    fig.savefig(path, dpi=settings.dpi, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print("wrote", path)
    return path


# ----------------------------------------------------------- spherical shell ----

def shell_grid(r1, r2, settings, phi=(0.5 * np.pi, 2.0 * np.pi)):
    """Structured grid of the shell r1 <= |x| <= r2 restricted to azimuth
    phi = (phi0, phi1); the removed wedge exposes the interior on two cut faces."""
    nr, nth, nph = settings.n_shell
    r = np.linspace(r1, r2, nr)
    th = np.linspace(0.0, np.pi, nth)
    ph = np.linspace(phi[0], phi[1], nph)
    R, TH, PH = np.meshgrid(r, th, ph, indexing="ij")
    return pv.StructuredGrid(R * np.sin(TH) * np.cos(PH), R * np.sin(TH) * np.sin(PH),
                             R * np.cos(TH))


def shell_surface(grid, scalars, name):
    grid.point_data[name] = np.asarray(scalars).ravel(order="F")
    return grid.extract_surface(algorithm="dataset_surface")


def principal_frame(sigma):
    """Unit eigenvectors of sigma (3, 3, N) ordered by decreasing eigenvalue:
    fibre, sheet, sheet normal (each (N, 3))."""
    s = np.moveaxis(np.asarray(sigma), (0, 1), (-2, -1))
    _, v = np.linalg.eigh(s)
    return v[..., :, 2], v[..., :, 1], v[..., :, 0]


def depth_shades(n_radii, n_per_radius):
    """One colour per frame point: endocardium light to epicardium dark."""
    lo, hi = FIBRE_DEPTH
    shades = [tuple(v) for v in np.linspace(lo, hi, n_radii)]
    return [c for c in shades for _ in range(n_per_radius)]


def add_frame_glyphs(pl, pts, frame, length):
    """The orthonormal frame at each point as three short arrows of equal length
    centred on the point; frame is a list of (vectors (N, 3), colour or a list of
    one colour per point)."""
    for vec, color in frame:
        colors = color if isinstance(color, list) else [color] * len(pts)
        for point, d, c in zip(pts, vec, colors):
            pl.add_mesh(pv.Arrow(start=point - 0.5 * length * d, direction=d,
                                 scale=length, tip_length=0.3, tip_radius=0.13,
                                 shaft_radius=0.05), color=c)


def shell_cutaway(settings, ref, scalar, clim=None, glyphs=None, direction=(1, 1, 0.7),
                  phi=(0.5 * np.pi, 2.0 * np.pi)):
    """Cut-away image of the shell coloured by scalar(x, y, z)."""
    p = ref.PARAMS["default"]
    grid = shell_grid(p["r1"], p["r2"], settings, phi)
    surf = shell_surface(grid, scalar(*point_array(grid)), "f")
    pl = plotter(settings)
    if clim is None:
        pl.add_mesh(surf, color=MESH_FACE, smooth_shading=True)
    else:
        pl.add_mesh(surf, scalars="f", cmap=CMAP, clim=clim, show_scalar_bar=False,
                    smooth_shading=True, lighting=False)
    if glyphs:
        glyphs(pl)
    camera(pl, (0, 0, 0), 4.4, direction)
    return render(pl), surf


def frame_view(vectors, centre, max_angle=25.0, samples=8000, seed=1):
    """View direction within max_angle degrees of `centre` (a direction looking into
    the removed wedge onto the cut face) maximising the shortest projected length of
    any frame vector at any site."""
    rng = np.random.default_rng(seed)
    d = rng.normal(size=(samples, 3))
    d /= np.linalg.norm(d, axis=1)[:, None]
    d = d[d @ centre >= np.cos(np.radians(max_angle))]
    all_vec = np.concatenate(vectors, axis=0)
    proj = np.sqrt(np.maximum(0.0, 1.0 - (d @ all_vec.T) ** 2)).min(axis=1)
    return d[np.argmax(proj)], float(proj.max())


def screen_basis(direction):
    """Right and up vectors of a camera that looks along -direction with up = z."""
    forward = -np.asarray(direction, float)
    forward /= np.linalg.norm(forward)
    right = np.cross(forward, (0.0, 0.0, 1.0))
    right /= np.linalg.norm(right)
    return right, np.cross(right, forward)


def glyph_separation(pts, vectors, length, direction, samples=24):
    """Smallest screen-space distance between arrows of different sites: every
    arrow (a segment from -L/2 to +L/2 about its site) is projected along the view
    direction and sampled; the distance is taken over all pairs of arrows that
    belong to different sites."""
    right, up = screen_basis(direction)
    u = np.linspace(-0.5, 0.5, samples)
    sets = []
    for i, point in enumerate(pts):
        arrows = []
        for vec in vectors:
            xyz = point[None, :] + u[:, None] * length * vec[i][None, :]
            arrows.append(np.stack([xyz @ right, xyz @ up], axis=1))
        sets.append(np.concatenate(arrows, axis=0))
    best = np.inf
    for i in range(len(sets)):
        for j in range(i + 1, len(sets)):
            diff = sets[i][:, None, :] - sets[j][None, :, :]
            best = min(best, float(np.sqrt((diff ** 2).sum(axis=2)).min()))
    return best


def shell_frame_glyphs(ref, p, phi=0.5 * np.pi, length=0.16, margin=0.05,
                       n_theta=2, n_candidates=300, seed=3):
    """Draw the frame (fibre, sheet, sheet normal) of sigma on the cut face phi as
    arrows of equal length centred on 3 x n_theta sites (three radii, n_theta polar
    angles each), and carry the view direction in `.direction`.  The sites are the
    candidate, out of a seeded search, whose arrows are best separated on screen.
    Asserted: the vectors are orthonormal; both ends of every arrow lie inside
    r1 < |x| < r2; arrows of different sites are more than `margin` apart in the
    image; every frame vector keeps a visible projected length."""
    half = 0.5 * length
    pad = half + 0.03
    radii = (p["r1"] + pad, 0.5 * (p["r1"] + p["r2"]), p["r2"] - pad)
    face = np.array([np.sin(phi), -np.cos(phi), 0.0])
    along = np.array([np.cos(phi), np.sin(phi), 0.0])
    centre = 0.94 * face + 0.34 * along + np.array([0.0, 0.0, 0.3])
    centre /= np.linalg.norm(centre)
    rng = np.random.default_rng(seed)
    best = None
    for _ in range(n_candidates):
        theta = np.sort(rng.uniform(0.5, np.pi - 0.5, size=(3, n_theta)), axis=1)
        pts = np.array([(r * np.sin(t) * np.cos(phi), r * np.sin(t) * np.sin(phi),
                         r * np.cos(t)) for r, row in zip(radii, theta) for t in row])
        sigma = ref.field("sigma", pts[:, 0], pts[:, 1], pts[:, 2])
        fibre, sheet, normal = principal_frame(sigma)
        direction, shortest = frame_view([fibre, sheet, normal], centre)
        if shortest < 0.5:
            continue
        sep = glyph_separation(pts, (fibre, sheet, normal), length, direction)
        if best is None or sep > best[0]:
            best = (sep, pts, fibre, sheet, normal, direction, shortest)
    assert best is not None, "no candidate keeps every frame vector visible"
    sep, pts, fibre, sheet, normal, direction, shortest = best
    assert np.abs((fibre * sheet).sum(1)).max() < 1e-12
    assert np.abs((fibre * normal).sum(1)).max() < 1e-12
    assert np.abs((sheet * normal).sum(1)).max() < 1e-12
    for vec in (fibre, sheet, normal):
        for sign in (-1.0, 1.0):
            r = np.linalg.norm(pts + sign * half * vec, axis=1)
            assert np.all(r > p["r1"]) and np.all(r < p["r2"]), "arrow leaves the wall"
    assert sep > margin, "glyphs overlap on screen (separation %.3f)" % sep

    def draw(pl):
        add_frame_glyphs(pl, pts, [(fibre, depth_shades(3, n_theta)),
                                   (sheet, SHEET), (normal, NORMAL)], length)
    draw.direction = tuple(direction)
    draw.shortest_projection = shortest
    draw.separation = sep
    return draw


# ------------------------------------------------------------------ meshes ----

def mesh_figure(name, path, outdir, settings, title):
    """Surface-with-edges render of constant/polyMesh of an author-supplied case."""
    if not os.path.isdir(os.path.join(path, "constant", "polyMesh")):
        raise FileNotFoundError(
            "mesh case %r: %s has no constant/polyMesh" % (name, path))
    tmp = tempfile.mkdtemp(prefix="mesh_case_")
    try:
        os.symlink(os.path.abspath(os.path.join(path, "constant")),
                   os.path.join(tmp, "constant"))
        system = os.path.join(path, "system")
        if os.path.isdir(system):
            os.symlink(os.path.abspath(system), os.path.join(tmp, "system"))
        foam = os.path.join(tmp, "case.foam")
        open(foam, "w").close()
        reader = pv.OpenFOAMReader(foam)
        reader.skip_zero_time = True
        mesh = reader.read()["internalMesh"]
        pl = plotter(settings)
        clipped = mesh.clip(normal=(0, 1, 0), origin=mesh.center, invert=False)
        pl.add_mesh(clipped, color="#e5e5e5", show_edges=True, edge_color="#444444",
                    line_width=0.6)
        b = mesh.bounds
        camera(pl, mesh.center, 1.7 * max(b[1] - b[0], b[3] - b[2], b[5] - b[4]))
        img = render(pl)
        fig, ax = plt.subplots(figsize=(4.2, 3.8))
        panel(fig, ax, img, title="%s: %d cells" % (title, mesh.n_cells))
        return save(fig, outdir, "mesh_%s.png" % name, settings)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def cli(description, mesh_cases):
    """Argument parser shared by the figure scripts; `mesh_cases` maps the
    accepted --mesh-case names to their titles."""
    ap = argparse.ArgumentParser(description=description)
    ap.add_argument("--outdir", default=None, help="output directory (default: figures/)")
    ap.add_argument("--draft", action="store_true", help="low resolution, fast")
    ap.add_argument("--mesh-case", action="append", default=[], metavar="NAME=PATH",
                    help="meshed case for a mesh figure; NAME is one of: "
                    + ", ".join(mesh_cases))
    return ap


def run_mesh_cases(args, mesh_cases, outdir, settings):
    """Mesh figures of the cases passed with --mesh-case; exit 1 on a bad case."""
    given = {}
    for item in args.mesh_case:
        if "=" not in item:
            print("error: --mesh-case expects NAME=PATH, got %r" % item, file=sys.stderr)
            sys.exit(1)
        name, path = item.split("=", 1)
        given[name] = path
    if not given:
        print("note: mesh figures skipped. Pass --mesh-case NAME=PATH (PATH holds "
              "constant/polyMesh of a meshed case) for: %s" % ", ".join(mesh_cases))
        return
    for name, path in given.items():
        if name not in mesh_cases:
            print("note: %r is not a mesh case of this script; skipped" % name)
            continue
        try:
            mesh_figure(name, path, outdir, settings, mesh_cases[name])
        except FileNotFoundError as err:
            print("error:", err, file=sys.stderr)
            sys.exit(1)
