"""Shared style and helpers of the figure scripts.

One colormap (viridis) for every scalar field, one serif font with Computer
Modern mathtext, one set of accent colours, one camera.  Fields are rendered
with pyvista off-screen into an image and placed in a matplotlib panel, so the
colourbars and the symbols are matplotlib mathtext.  Figures show the problem
only: geometry, exact fields, conductivity frames and sources.
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

pv.OFF_SCREEN = True

HERE = os.path.dirname(os.path.abspath(__file__))
GENERATORS = os.path.abspath(os.path.join(HERE, "..", ".."))
REFERENCE_DIR = os.path.join(GENERATORS, "reference")

CMAP = "viridis"
FIBRE = "#d95f02"
SHEET = "#e7298a"
GRAPH = "#d95f02"
KERNEL = "#7570b3"
ELECTRODE = "#d95f02"
GREY = "#808080"
DPI = 300
DRAFT_DPI = 80

plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["DejaVu Serif"],
    "mathtext.fontset": "cm",
    "font.size": 9,
    "axes.titlesize": 9,
    "axes.labelsize": 9,
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


def tutorial_electrodes(case_dir):
    """Electrodes of a tutorial's electroProperties; the built-in set when absent."""
    builtin = [(-0.5, 0.5, 0.5), (1.5, 0.5, 0.5), (1.2, 0.23, 0.61),
               (1.35, 0.74, 0.28), (1.55, 0.41, 0.83)]
    path = os.path.join(GENERATORS, "..", case_dir, "constant", "electroProperties")
    if not os.path.exists(path):
        return builtin
    with open(path) as fh:
        found = re.findall(r"^\s*([ER]\d+)\s+\(([^)]*)\);", fh.read(), re.M)
    return [tuple(float(v) for v in xyz.split()) for _, xyz in found] or builtin


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
        cb.ax.tick_params(labelsize=7)
    return ax


def legend(ax, entries, ncol=3):
    handles = [Line2D([0], [0], color=c, lw=3, label=l) for l, c in entries]
    ax.legend(handles=handles, loc="upper center", ncol=ncol, frameon=False,
              fontsize=8, bbox_to_anchor=(0.5, 0.0), handlelength=1.2,
              columnspacing=1.0)


def colourbar_2d(fig, ax, mappable, label):
    cb = fig.colorbar(mappable, ax=ax, orientation="vertical", fraction=0.05,
                      pad=0.04, shrink=0.8)
    cb.set_label(label)
    cb.ax.tick_params(labelsize=7)
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


def frame_points(r1, r2, n_theta=6, phi=0.5 * np.pi):
    """Points on the cut plane azimuth=phi, on three radii."""
    radii = (1.06 * r1, 0.5 * (r1 + r2), 0.94 * r2)
    th = np.linspace(0.35, np.pi - 0.35, n_theta)
    pts = [(r * np.sin(t) * np.cos(phi), r * np.sin(t) * np.sin(phi), r * np.cos(t))
           for r in radii for t in th]
    return np.array(pts)


def add_frame_glyphs(pl, pts, fibre, sheet, length):
    """Fibre and sheet directions as thin bars centred on pts."""
    for vec, color in ((fibre, FIBRE), (sheet, SHEET)):
        cloud = pv.PolyData(pts)
        cloud["dir"] = vec
        bar = pv.Cylinder(direction=(1, 0, 0), radius=0.09, height=1.0, resolution=12)
        pl.add_mesh(cloud.glyph(orient="dir", scale=False, factor=length, geom=bar),
                    color=color, lighting=False)


def shell_cutaway(settings, ref, scalar, clim=None, glyphs=None, direction=(1, 1, 0.7),
                  phi=(0.5 * np.pi, 2.0 * np.pi)):
    """Cut-away image of the shell coloured by scalar(x, y, z)."""
    p = ref.PARAMS["default"]
    grid = shell_grid(p["r1"], p["r2"], settings, phi)
    surf = shell_surface(grid, scalar(*point_array(grid)), "f")
    pl = plotter(settings)
    if clim is None:
        pl.add_mesh(surf, color="#e0e0e0", smooth_shading=True)
    else:
        pl.add_mesh(surf, scalars="f", cmap=CMAP, clim=clim, show_scalar_bar=False,
                    smooth_shading=True, lighting=False)
    if glyphs:
        glyphs(pl)
    camera(pl, (0, 0, 0), 4.4, direction)
    return render(pl), surf


def shell_frame_glyphs(ref, p, phi=0.5 * np.pi):
    pts = frame_points(p["r1"], p["r2"], phi=phi)
    sigma = ref.field("sigma", pts[:, 0], pts[:, 1], pts[:, 2])
    fibre, sheet, _ = principal_frame(sigma)

    def draw(pl):
        add_frame_glyphs(pl, pts * 1.0, fibre, sheet, 0.2)
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
