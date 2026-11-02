#!/usr/bin/env python3
"""Figures of the lead-field manufactured solutions: the problem only.

    python figures.py [--outdir DIR] [--draft] [--mesh-case NAME=PATH ...]

Analytic figures need leadfield_reference.py (this folder) and
generators/shell_mms.py.  Mesh figures read constant/polyMesh and its
cellZones of a meshed case the author supplies:

    slab_hex, slab_tet, layered_hex, layered_tet, shell_hex, shell_tet
"""

import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "_common"))
sys.path.insert(0, HERE)

import figstyle as fs  # noqa: E402
import leadfield_reference as lf  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pyvista as pv  # noqa: E402
from matplotlib.gridspec import GridSpec  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import FancyArrowPatch, Rectangle  # noqa: E402
from mpl_toolkits.mplot3d.art3d import Poly3DCollection  # noqa: E402

MESH_CASES = {
    "slab_hex": "Slab, hex",
    "slab_tet": "Slab, tet",
    "layered_hex": "Layered slab, hex",
    "layered_tet": "Layered slab, tet",
    "shell_hex": "Curved conductor, hex",
    "shell_tet": "Curved conductor, tet",
}
XL = r"$x/(1\,\mathrm{m})$"
YL = r"$y/(1\,\mathrm{m})$"

# region colours: heart, torso layers, the bath ball
HEART = "#e8a39b"
TORSO = "#bcd4e6"
TORSO_INNER = "#9cc0dd"
TORSO_OUTER = "#dbe8f2"
POOL = "#f1d38e"
ZONES = {
    "myocardium": (HEART, "myocardium"),
    "bath": (TORSO, "bath"),
    "torsoInner": (TORSO_INNER, "torsoInner"),
    "torsoOuter": (TORSO_OUTER, "torsoOuter"),
    "bloodPool": (POOL, "bloodPool"),
}
SURFACE = fs.NAMED                   # body-surface electrodes: the ECG
PROBE = fs.KERNEL                    # interior probes
PIN = fs.PIN
SURFACE_MARK, PROBE_MARK = "^", "o"


def split_electrodes(electrodes):
    """(surface, probes): the body-surface electrodes surface1..surface8 and
    the interior probes, each {name: point}."""
    surface = {n: p for n, p in electrodes.items() if n.startswith("surface")}
    probes = {n: p for n, p in electrodes.items() if not n.startswith("surface")}
    return surface, probes


def short(name):
    """Label of an electrode: the number of a body-surface electrode, else the name."""
    return name[len("surface"):] if name.startswith("surface") else name


def electrode_legend(shape_surface=SURFACE_MARK, shape_probe=PROBE_MARK):
    return [("body-surface electrode (surface1 to surface8, labelled 1 to 8)", SURFACE,
             shape_surface),
            ("interior probe", PROBE, shape_probe),
            ("reference point (pin)", PIN, "s")]


# ------------------------------------------------------------------ slab ----

def _box_faces(lo, hi):
    (x0, y0, z0), (x1, y1, z1) = lo, hi
    v = np.array([[x0, y0, z0], [x1, y0, z0], [x1, y1, z0], [x0, y1, z0],
                  [x0, y0, z1], [x1, y0, z1], [x1, y1, z1], [x0, y1, z1]])
    idx = [(0, 1, 2, 3), (4, 5, 6, 7), (0, 1, 5, 4), (2, 3, 7, 6), (1, 2, 6, 5), (0, 3, 7, 4)]
    return [v[list(f)] for f in idx]


SLAB_REGIONS = ((-1.0, 0.0, TORSO, "bath"), (0.0, 1.0, HEART, "myocardium"),
                (1.0, 2.0, TORSO, "bath"))
LAYERED_REGIONS = ((-1.0, lf.X_LEFT, TORSO_OUTER, "torsoOuter"),
                   (lf.X_LEFT, 0.0, TORSO_INNER, "torsoInner"),
                   (0.0, 1.0, HEART, "myocardium"),
                   (1.0, lf.X_RIGHT, TORSO_INNER, "torsoInner"),
                   (lf.X_RIGHT, 2.0, TORSO_OUTER, "torsoOuter"))

# label offsets (points), chosen so that no two labels meet
SLAB_3D_OFFSET = {"surface1": (-0.05, 0.0, -0.03, "right")}
SLAB_LABEL_OFFSET = {"torsoLeft": (0, 9), "torsoLeftFar": (30, -3), "torsoRight": (0, 9),
                     "nearLeft": (-4, -14), "nearRight": (4, -14),
                     "surface1": (-9, 0), "surface2": (-9, 0), "surface3": (9, 0),
                     "surface4": (9, 0), "surface5": (0, -9), "surface6": (0, 9),
                     "surface7": (0, 9), "surface8": (0, -9)}


def _scatter3d(ax, points, colour, marker):
    e = np.array(list(points.values()))
    ax.scatter(*e.T, marker=marker, s=34, c=colour, edgecolors="k", linewidths=0.4,
               depthshade=False, zorder=20)
    for p in points.values():
        ax.plot([p[0], p[0]], [p[1], p[1]], [0, p[2]], color=colour, lw=0.6, ls=":")


def _slab_3d(ax, electrodes, pin):
    for x0, x1, colour, _ in SLAB_REGIONS:
        ax.add_collection3d(Poly3DCollection(_box_faces((x0, 0, 0), (x1, 1, 1)),
                                             facecolor=colour, edgecolor=fs.DOMAIN,
                                             linewidths=0.6, alpha=0.28))
    surface, probes = split_electrodes(electrodes)
    _scatter3d(ax, surface, SURFACE, SURFACE_MARK)
    _scatter3d(ax, probes, PROBE, PROBE_MARK)
    ax.scatter(*pin, marker="s", s=30, c=PIN, edgecolors="k", linewidths=0.4,
               depthshade=False, zorder=20)
    ax.plot([pin[0]] * 2, [pin[1]] * 2, [0, pin[2]], color=PIN, lw=0.6, ls=":")
    ax.set_xlim(-1, 2)
    ax.set_ylim(0, 1)
    ax.set_zlim(0, 1)
    ax.set_box_aspect((3, 1, 1), zoom=1.02)
    ax.view_init(elev=24, azim=-62)
    ax.set_xticks([-1, 0, 1, 2])
    ax.set_yticks([0, 1])
    ax.set_zticks([0, 1])
    ax.set_xlabel(XL, labelpad=6)
    ax.set_ylabel(YL, labelpad=2)
    ax.set_zlabel(r"$z/(1\,\mathrm{m})$", labelpad=0)
    for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
        axis.pane.set_visible(False)
    ax.grid(False)
    for name, p in electrodes.items():
        off = SLAB_3D_OFFSET.get(name, (0.0, 0.0, 0.09, "center"))
        ax.text(p[0] + off[0], p[1] + off[1], p[2] + off[2], short(name),
                fontsize=9 if name in surface else 8,
                fontweight="bold" if name in surface else "normal",
                color=SURFACE if name in surface else PROBE, ha=off[3], zorder=30)
    ax.text(pin[0], pin[1], pin[2] + 0.12, "reference point", fontsize=8, color=PIN,
            ha="left", zorder=30)


def _plan(ax, regions, electrodes, pin, labels=True):
    for x0, x1, colour, name in regions:
        ax.add_patch(Rectangle((x0, 0), x1 - x0, 1, facecolor=colour, edgecolor=fs.DOMAIN,
                               lw=0.8))
    surface, probes = split_electrodes(electrodes)
    for points, colour, marker in ((surface, SURFACE, SURFACE_MARK),
                                   (probes, PROBE, PROBE_MARK)):
        e = np.array(list(points.values()))
        ax.scatter(e[:, 0], e[:, 1], marker=marker, s=34, c=colour, edgecolors="k",
                   linewidths=0.4, zorder=5)
    ax.scatter([pin[0]], [pin[1]], marker="s", s=28, c=PIN, edgecolors="k", linewidths=0.4,
               zorder=5)
    if labels:
        for name, p in electrodes.items():
            ax.annotate(short(name), xy=(p[0], p[1]),
                        xytext=SLAB_LABEL_OFFSET.get(name, (0, 9)),
                        textcoords="offset points", ha="center", va="center",
                        fontsize=8, fontweight="bold" if name in surface else "normal",
                        color=SURFACE if name in surface else PROBE)
        ax.annotate("reference point", xy=(pin[0], pin[1]), xytext=(7, 0),
                    textcoords="offset points", ha="left", va="center", fontsize=8,
                    color=PIN)
    ax.set_xlim(-1.32, 2.32)
    ax.set_ylim(-0.15, 1.15)
    ax.set_aspect("equal")
    ax.set_xticks([-1, -0.5, 0, 0.5, 1, 1.5, 2])
    ax.set_yticks([0, 1])
    ax.set_xlabel(XL)
    ax.set_ylabel(YL)


def slab_geometry(settings, outdir):
    electrodes, pin = lf.case_points("slab")
    fig = plt.figure(figsize=(12.6, 7.6))
    gs = GridSpec(1, 2, figure=fig, top=0.36, bottom=0.08, wspace=0.18)
    ax = fig.add_axes([-0.06, 0.34, 1.12, 0.70], projection="3d", computed_zorder=False)
    _slab_3d(ax, electrodes, pin)
    ax.set_title(r"(a) slab $[-1,2]\times[0,1]^2$: torso, heart, torso", loc="left")
    fs.legend(ax, [("bath (torso slabs)", TORSO), ("myocardium (heart slab)", HEART)]
              + electrode_legend(), ncol=3, inside="upper center")
    ax = fig.add_subplot(gs[0, 0])
    _plan(ax, SLAB_REGIONS, electrodes, pin)
    for x, text in ((-0.5, "bath"), (0.5, "myocardium"), (1.5, "bath")):
        ax.text(x, 0.88, text, ha="center", fontsize=8, color="#333333")
    ax.set_title(r"(b) slab, seen along $z$: every outer face insulated", loc="left",
                 fontsize=11)
    ax = fig.add_subplot(gs[0, 1])
    _plan(ax, LAYERED_REGIONS, electrodes, pin, labels=False)
    for x0, x1, _, name in LAYERED_REGIONS:
        sig = {"torsoOuter": r"$\sigma_{b,2}=0.02$", "torsoInner": r"$\sigma_{b,1}=0.1$",
               "myocardium": r"$(\boldsymbol{\sigma}_i+\boldsymbol{\sigma}_e)_{xx}$"}[name]
        ax.text(0.5 * (x0 + x1), 0.88, sig, ha="center", fontsize=8, color="#333333")
    for x in (lf.X_LEFT, lf.X_RIGHT):
        ax.plot([x, x], [0, 1], color=fs.DOMAIN, lw=1.4)
    for x0, name, n in ((-1.0, "xMin", -1), (2.0, "xMax", 1)):
        for yy in (0.3, 0.7):
            ax.add_patch(FancyArrowPatch((x0 + 0.28 if n > 0 else x0 + 0.02, yy),
                                         (x0 + 0.02 if n > 0 else x0 - 0.26, yy),
                                         arrowstyle="-|>", mutation_scale=10,
                                         color=fs.ECG_DARK, lw=1.4, zorder=6))
        ax.text(x0 + 0.17 * n, 0.5, name, ha="center", va="center", fontsize=8,
                color=fs.ECG_DARK, rotation=90)
    ax.set_title(r"(c) layered: layers and the current $-\alpha\,\mathbf{e}_x$",
                 loc="left", fontsize=11)
    ax.legend(handles=[Line2D([0], [0], color=c, lw=3, label=l) for l, c in
                       (("torsoOuter", TORSO_OUTER), ("torsoInner", TORSO_INNER),
                        ("myocardium", HEART),
                        (r"current density $-\sigma\nabla\phi$", fs.ECG_DARK))],
              loc="upper center", bbox_to_anchor=(0.5, -0.32), ncol=4, frameon=False,
              fontsize=9, handlelength=1.2, columnspacing=1.0)
    fs.save(fig, outdir, "slab_geometry.png", settings)


def _shade(ax, regions):
    for x0, x1, colour, _ in regions:
        ax.axvspan(x0, x1, color=colour, alpha=0.45, lw=0)


def _mark(ax, electrodes, pin, func, offsets):
    """Electrodes at their x on the curve func; electrodes that share an x share
    one label."""
    surface, _ = split_electrodes(electrodes)
    groups = {}
    for name, p in electrodes.items():
        groups.setdefault(round(float(p[0]), 6), []).append(name)
    for x, names in groups.items():
        v = func(x)
        kind = names[0] in surface
        ax.plot([x], [v], SURFACE_MARK if kind else PROBE_MARK,
                color=SURFACE if kind else PROBE, mec="k", mew=0.4, ms=7, zorder=6)
        ax.annotate(", ".join(short(n) for n in names), xy=(x, v),
                    xytext=offsets.get(names[0], (0, 9)), textcoords="offset points",
                    ha="center", va="center", fontsize=8,
                    fontweight="bold" if kind else "normal",
                    color=SURFACE if kind else PROBE)
    ax.plot([pin[0]], [func(pin[0])], "s", color=PIN, mec="k", mew=0.4, ms=6, zorder=6)


def slab_potential(settings, outdir):
    electrodes, pin = lf.case_points("slab")
    x = np.linspace(-1, 2, 1201)
    xh = np.linspace(0, 1, 401)
    fig, axes = plt.subplots(2, 1, figsize=(12.6, 9.6), gridspec_kw={"hspace": 0.32})
    ax = axes[0]
    _shade(ax, SLAB_REGIONS)
    for t, colour in ((0.0, fs.ECG), (0.05, fs.ECG_LIGHT)):
        ax.plot(x, lf.gauged(lf.slab_phi_raw, x, pin[0], t), color=colour)
        ax.plot(xh, lf.slab_vm(xh, t), "--", color=fs.VM if t == 0 else fs.VM_LIGHT)
    phi0 = lambda q: float(lf.gauged(lf.slab_phi_raw, q, pin[0], 0.0))  # noqa: E731
    _mark(ax, electrodes, pin, phi0, SLAB_POTENTIAL_OFFSET)
    ax.axhline(2 * lf.K, color=fs.GREY, lw=0.6, ls=":")
    ax.text(1.18, 2 * lf.K + 0.06, r"$2\kappa$", ha="center", color=fs.GREY)
    ax.set_xlim(-1.05, 2.05)
    ax.set_ylim(-1.25, 1.85)
    ax.set_xlabel(XL)
    ax.set_ylabel(r"$\phi$, $V_m$")
    ax.set_title(r"(a) slab: $\phi$ and the heart $V_m$", loc="left", fontsize=11)
    ax.legend(handles=[Line2D([0], [0], color=fs.ECG, label=r"$\phi$, $t=0$"),
                       Line2D([0], [0], color=fs.ECG_LIGHT, label=r"$\phi$, $t=0.05$"),
                       Line2D([0], [0], color=fs.VM, ls="--", label=r"$V_m$, $t=0$"),
                       Line2D([0], [0], color=fs.VM_LIGHT, ls="--", label=r"$V_m$, $t=0.05$")],
              frameon=False, fontsize=9, loc="lower right", ncol=2)
    ax = axes[1]
    _shade(ax, LAYERED_REGIONS)
    for xi in (lf.X_LEFT, lf.X_RIGHT):
        ax.axvline(xi, color=fs.DOMAIN, lw=0.8, ls=":")
    on = lambda q: lf.gauged(lf.layered_phi_raw, q, pin[0], 0.0)  # noqa: E731
    off = lambda q: lf.gauged(lf.layered_phi_raw, q, pin[0], 0.0, alpha=0.0)  # noqa: E731
    ax.plot(x, on(x), color=fs.ECG)
    ax.plot(x, off(x), "--", color=fs.ECG_LIGHT)
    _mark(ax, electrodes, pin, lambda q: float(on(q)), LAYERED_POTENTIAL_OFFSET)
    for x0, x1, _, name in LAYERED_REGIONS:
        sig = {"torsoOuter": r"$\sigma_{b,2}$", "torsoInner": r"$\sigma_{b,1}$",
               "myocardium": r"$(\boldsymbol{\sigma}_i+\boldsymbol{\sigma}_e)_{xx}$"}[name]
        ax.text(0.5 * (x0 + x1), -1.15, sig, ha="center", fontsize=9, color="#333333")
    ax.set_xlim(-1.05, 2.05)
    ax.set_ylim(-1.3, 2.3)
    ax.set_xlabel(XL)
    ax.set_ylabel(r"$\phi$")
    ax.set_title(r"(b) layered: $\phi=\phi_{\mathrm{slab}}+\alpha R(x)$ at $t=0$", loc="left",
                 fontsize=11)
    ax.legend(handles=[Line2D([0], [0], color=fs.ECG, label=r"$\alpha=0.01$"),
                       Line2D([0], [0], color=fs.ECG_LIGHT, ls="--",
                              label=r"$\alpha=0$ (insulated)")]
              + [Line2D([0], [0], color=c, marker=m, mec="k", mew=0.4, ls="", label=l)
                 for l, c, m in electrode_legend()],
              frameon=False, fontsize=9, loc="upper left")
    fs.save(fig, outdir, "slab_potential.png", settings)


SLAB_POTENTIAL_OFFSET = {"nearLeft": (-6, -15), "nearRight": (-8, 9),
                         "torsoLeftFar": (0, -11), "torsoLeft": (0, -11),
                         "surface5": (0, 9), "surface6": (0, 9), "surface1": (0, 9),
                         "surface3": (0, 9), "surface7": (0, 9), "surface8": (0, 9),
                         "torsoRight": (0, -11)}
LAYERED_POTENTIAL_OFFSET = {"nearLeft": (-2, 9), "nearRight": (-6, 9), "torsoLeft": (0, -11),
                            "surface1": (12, -4), "surface3": (-14, 4)}


# ----------------------------------------------------------------- shell ----

R1, R2, R3 = 0.5, 1.0, 1.5
WEDGE = (-0.5 * np.pi, 0.0)          # removed azimuth range of the cut-away


def _in_wedge(p):
    az = np.arctan2(p[..., 1], p[..., 0])
    return (az > WEDGE[0]) & (az < WEDGE[1])


MARKER_CLEARANCE = 0.07             # larger than the marker radii


def _visible(p, d, clearance=MARKER_CLEARANCE):
    """No opaque solid (|x| < r2 outside the removed wedge) between p and a
    distant camera in direction d, with `clearance` around the ray so that no
    part of the marker is hidden either."""
    s = np.linspace(clearance + 0.02, 6.0, 3000)[:, None]
    q = p[None, :] + s * d[None, :]
    r = np.linalg.norm(q, axis=1)
    return not np.any((r < R2 + clearance) & ~_in_wedge(q))


def shell_view(points, nominal=(0.35, -0.45, 0.85), max_angle=35.0, samples=6000, seed=2):
    """Camera direction within max_angle of `nominal` from which every point is
    visible, keeping the largest angular margin from occlusion."""
    rng = np.random.default_rng(seed)
    nominal = np.asarray(nominal, float) / np.linalg.norm(nominal)
    d = rng.normal(size=(samples, 3))
    d /= np.linalg.norm(d, axis=1)[:, None]
    d = d[d @ nominal >= np.cos(np.radians(max_angle))]
    best = None
    for di in d:
        if all(_visible(p, di) for p in points):
            # margin: the smallest perturbation that hides a point, sampled
            jit = di[None, :] + 0.12 * rng.normal(size=(12, 3))
            jit /= np.linalg.norm(jit, axis=1)[:, None]
            score = sum(all(_visible(p, j) for p in points) for j in jit)
            score += 0.5 * float(di @ nominal)
            if best is None or score > best[0]:
                best = (score, di)
    assert best is not None, "no view shows every electrode and the reference point"
    return best[1]


def _region_surface(r_in, r_out, settings):
    grid = fs.shell_grid(max(r_in, 1e-4), r_out, settings,
                         phi=(WEDGE[1], WEDGE[0] + 2 * np.pi))
    return grid.extract_surface(algorithm="dataset_surface")


def place_labels(points, names, right, up, step=0.15, char=0.028):
    """Screen positions of labels: for each point in turn, of eight directions
    the one whose label (a row of samples as wide as the text) lies farthest
    from every marker and every label already placed, outward preferred."""
    xy = np.array([[np.dot(p, right), np.dot(p, up)] for p in points])
    angles = np.radians(np.arange(0, 360, 45))
    dirs = np.stack([np.cos(angles), np.sin(angles)], -1)
    placed = []
    out = []
    for i, (q, name) in enumerate(zip(xy, names)):
        half = 0.5 * char * len(name)
        outward = q / (np.linalg.norm(q) or 1.0)
        best = None
        for d in dirs:
            c = q + step * d + np.array([np.sign(d[0]) * half * (abs(d[0]) > 0.5), 0.0])
            row = c[None, :] + np.linspace(-half, half, 7)[:, None] * np.array([1.0, 0.0])
            others = np.delete(xy, i, axis=0)
            gap = min(np.min(np.linalg.norm(row[:, None] - others[None], axis=-1)),
                      min((np.min(np.linalg.norm(row[:, None] - r[None], axis=-1))
                           for r in placed), default=np.inf))
            score = min(gap, 0.3) + 0.02 * float(d @ outward)
            if best is None or score > best[0]:
                best = (score, c, row)
        placed.append(best[2])
        out.append(np.asarray(points[i]) + (best[1] - q)[0] * right + (best[1] - q)[1] * up)
    return np.array(out)


def shell_regions_image(settings, electrodes, pin):
    pts = list(electrodes.values()) + [pin]
    d = shell_view(pts)
    surface, probes = split_electrodes(electrodes)
    pl = fs.plotter(settings)
    pl.enable_depth_peeling()
    pl.add_mesh(_region_surface(0.0, R1, settings), color=POOL, smooth_shading=True,
                ambient=0.3)
    pl.add_mesh(_region_surface(R1, R2, settings), color=HEART, smooth_shading=True,
                ambient=0.3)
    # the outer bath whole and translucent: the surface electrodes in the
    # removed quarter stay on a drawn sphere
    pl.add_mesh(pv.Sphere(radius=R3, theta_resolution=120, phi_resolution=60),
                color=TORSO, opacity=0.16, smooth_shading=True)
    for p in surface.values():
        pl.add_mesh(pv.Sphere(radius=0.05, center=p), color=SURFACE)
    for p in probes.values():
        pl.add_mesh(pv.PlatonicSolid("octahedron", radius=0.065, center=p), color=PROBE)
    pl.add_mesh(pv.Cube(center=pin, x_length=0.09, y_length=0.09, z_length=0.09), color=PIN)
    right, up = fs.screen_basis(d)
    font = max(8, settings.window[0] // 42)
    names = list(surface) + list(probes) + ["reference point"]
    where = place_labels(list(surface.values()) + list(probes.values()) + [pin],
                         [short(n) for n in names[:-1]] + [names[-1]], right, up)
    n_s, n_p = len(surface), len(probes)
    pl.add_point_labels(where[:n_s], [short(n) for n in surface], font_size=font, bold=True,
                        shape=None, show_points=False, always_visible=True,
                        text_color=SURFACE, font_family="arial",
                        justification_horizontal="center", justification_vertical="center")
    pl.add_point_labels(where[n_s:n_s + n_p], list(probes), font_size=font, shape=None,
                        show_points=False, always_visible=True, text_color=PROBE,
                        font_family="arial", justification_horizontal="center",
                        justification_vertical="center")
    pl.add_point_labels(where[-1:], ["reference point"], font_size=font, shape=None,
                        show_points=False, always_visible=True, text_color=PIN,
                        font_family="arial", justification_horizontal="center",
                        justification_vertical="center")
    fs.camera(pl, (0, 0, 0), 6.6, tuple(d))
    return fs.render(pl)


class ShellTensor:
    """The interface fs.shell_frame_glyphs and fs.shell_cutaway expect:
    PARAMS["default"] with r1, r2 and field("sigma", x, y, z) -> (3, 3, N)."""

    def __init__(self):
        self.PARAMS = {"default": {"r1": R1, "r2": R2}}

    def field(self, name, x, y, z):
        assert name == "sigma"
        pts = np.stack([np.ravel(x), np.ravel(y), np.ravel(z)], -1)
        g = lf.shell_model().eval_xyz(1, "GiZ", pts, 0.0)
        full = np.stack([np.stack([g[:, 0], g[:, 1], g[:, 2]], 0),
                         np.stack([g[:, 1], g[:, 3], g[:, 4]], 0),
                         np.stack([g[:, 2], g[:, 4], g[:, 5]], 0)], 0)
        return full


def check_frames_coincide():
    """Gi and Ge share their eigenvectors (one fibre angle, one frame) and order
    them alike: f, s, n_s by decreasing eigenvalue."""
    m = lf.shell_model()
    rng = np.random.default_rng(4)
    u = rng.normal(size=(50, 3))
    pts = (R1 + (R2 - R1) * rng.uniform(size=(50, 1))) * u / np.linalg.norm(u, axis=1)[:, None]
    frames = []
    for name in ("GiZ", "GeZ"):
        g = m.eval_xyz(1, name, pts, 0.0)
        full = np.stack([np.stack([g[:, 0], g[:, 1], g[:, 2]], -1),
                         np.stack([g[:, 1], g[:, 3], g[:, 4]], -1),
                         np.stack([g[:, 2], g[:, 4], g[:, 5]], -1)], -2)
        _, v = np.linalg.eigh(full)
        frames.append(v)
    dots = np.abs(np.einsum("nik,nik->nk", frames[0], frames[1]))
    assert np.all(dots > 1 - 1e-10), dots.min()


def shell_regions(settings, outdir):
    electrodes, pin = lf.case_points("shell")
    check_frames_coincide()
    img_r = shell_regions_image(settings, electrodes, pin)
    ref = ShellTensor()
    glyph = fs.shell_frame_glyphs(ref, ref.PARAMS["default"])
    img_g, _ = fs.shell_cutaway(settings, ref, lambda x, y, z: 0 * x, None, glyph,
                                glyph.direction)
    fig, axes = plt.subplots(1, 2, figsize=(12.6, 6.0), gridspec_kw={"wspace": 0.08})
    fs.panel(fig, axes[0], img_r, title="regions, electrodes, probes and the reference point",
             tag="(a)")
    fs.legend(axes[0], [(r"bloodPool $|\mathbf{x}|<0.5$", POOL),
                        (r"myocardium $0.5<|\mathbf{x}|<1$", HEART),
                        (r"bath $1<|\mathbf{x}|<1.5$", TORSO),
                        ("body-surface electrode (labelled 1 to 8)", SURFACE, "o"),
                        ("interior probe", PROBE, "D"), ("reference point", PIN, "s")],
              ncol=2)
    fs.panel(fig, axes[1], img_g, title=r"frame of $\boldsymbol{\sigma}_i$ and $\boldsymbol{\sigma}_e$ in the wall", tag="(b)")
    fs.legend(axes[1], [("fibre (light: endo, dark: epi)", fs.FIBRE),
                        ("sheet", fs.SHEET), ("sheet normal", fs.NORMAL)], ncol=1)
    fs.save(fig, outdir, "shell_regions.png", settings)


def _plane(n):
    s = np.linspace(-R3, R3, n)
    X, Y = np.meshgrid(s, s, indexing="ij")
    pts = np.stack([X, Y, np.zeros_like(X)], -1)
    return X, Y, pts


def _circles(ax):
    t = np.linspace(0, 2 * np.pi, 300)
    for r in (R1, R2, R3):
        ax.plot(r * np.cos(t), r * np.sin(t), color="k" if r < R3 else fs.DOMAIN, lw=0.8)


def _image(fig, ax, X, Y, V, label, title, clim=None, cmap=fs.CMAP):
    clim = clim or fs.span(V[np.isfinite(V)])
    im = ax.pcolormesh(X, Y, np.ma.masked_invalid(V), cmap=cmap, vmin=clim[0], vmax=clim[1],
                       shading="auto", rasterized=True)
    _circles(ax)
    ax.set_aspect("equal")
    ax.set_xlim(-1.55, 1.55)
    ax.set_ylim(-1.55, 1.55)
    ax.set_xlabel(XL)
    ax.set_ylabel(YL)
    ax.set_title(title, loc="left", fontsize=11)
    fs.colourbar_2d(fig, ax, im, label)


def shell_fields(settings, outdir):
    electrodes, pin = lf.case_points("shell")
    n = 121 if settings.draft else 401
    X, Y, pts = _plane(n)
    flat = pts.reshape(-1, 3)
    r = np.linalg.norm(flat, axis=1)
    outside = r > R3
    phi = lf.shell_phi(flat, 0.0, pin)
    phi[outside] = np.nan
    zones = lf.shell_zone(flat)
    sig = lf.shell_eval("sigmaB", flat, 0.0, only=(0, 2))
    sig[outside] = np.nan
    vm = lf.shell_eval("Vm", flat, 0.0, only=(1,))
    fig, axes = plt.subplots(1, 3, figsize=(14.6, 4.4), gridspec_kw={"wspace": 0.62})
    _image(fig, axes[0], X, Y, phi.reshape(X.shape), r"$\phi$",
           r"(a) exact $\phi$ on $z=0$, $t=0$")
    axes[0].plot([pin[0]], [pin[1]], "s", color=PIN, mec="k", mew=0.4, ms=6)
    _image(fig, axes[1], X, Y, sig.reshape(X.shape), r"$\sigma_b$",
           r"(b) $\sigma_b$ on $z=0$ (torso zones)")
    _image(fig, axes[2], X, Y, vm.reshape(X.shape), r"$V_m$",
           r"(c) exact $V_m$ on $z=0$, $t=0$ (heart)")
    assert np.all(zones[np.isfinite(vm)] == 1)
    fs.save(fig, outdir, "shell_fields.png", settings)


# ------------------------------------------------------------------ meshes ----

def mesh_zone_figure(name, path, outdir, settings, title):
    """Cut-away of constant/polyMesh coloured by cellZone, cell edges drawn."""
    data = lf.read_polymesh(path, zones=True)
    mesh = data["internalMesh"]
    zones = data["zones"]["cellZones"]
    mid = np.asarray(mesh.center)
    side = np.sign(np.asarray(fs.MESH_VIEW))
    pl = fs.plotter(settings)
    entries = []
    total = 0
    order = list(ZONES)
    for zone in sorted(zones.keys(), key=lambda z: order.index(z) if z in order else 99):
        grid = zones[zone]
        total += grid.n_cells
        centres = grid.cell_centers().points
        removed = np.all((centres - mid) * side > 0.0, axis=1)
        colour, label = ZONES.get(zone, ("#dfe3e6", zone))
        pl.add_mesh(grid.extract_cells(np.flatnonzero(~removed)), color=colour,
                    show_edges=True, edge_color="#6b757d", line_width=0.5,
                    smooth_shading=False, ambient=0.35, diffuse=0.7, specular=0.0)
        entries.append(("%s: %d" % (label, grid.n_cells), colour))
    assert total == mesh.n_cells, (name, total, mesh.n_cells)
    b = mesh.bounds
    size = max(b[1] - b[0], b[3] - b[2], b[5] - b[4])
    camera_dir = fs.MESH_VIEW if "shell" in name else (0.55, -1.15, 0.8)
    fs.camera(pl, mesh.center, 2.3 * size, camera_dir)
    img = fs.render(pl)
    fig, ax = plt.subplots(figsize=(5.4, 4.4) if "shell" in name else (6.4, 3.6))
    fs.panel(fig, ax, img, title="%s: %d cells" % (title, mesh.n_cells))
    fs.legend(ax, entries, ncol=1 if "shell" in name else 3)
    return fs.save(fig, outdir, "mesh_%s.png" % name, settings)


def run_mesh_cases(args, outdir, settings):
    given = {}
    for item in args.mesh_case:
        if "=" not in item:
            print("error: --mesh-case expects NAME=PATH, got %r" % item, file=sys.stderr)
            sys.exit(1)
        key, path = item.split("=", 1)
        given[key] = path
    if not given:
        print("note: mesh figures skipped. Pass --mesh-case NAME=PATH (PATH holds "
              "constant/polyMesh of a meshed case) for: %s" % ", ".join(MESH_CASES))
        return
    for key, path in given.items():
        if key not in MESH_CASES:
            print("note: %r is not a mesh case of this script; skipped" % key)
            continue
        try:
            mesh_zone_figure(key, path, outdir, settings, MESH_CASES[key])
        except FileNotFoundError as err:
            print("error:", err, file=sys.stderr)
            sys.exit(1)


def main():
    ap = fs.cli(__doc__, MESH_CASES)
    ap.add_argument("--only", default=None, help="comma-separated figure builders to run")
    args = ap.parse_args()
    settings = fs.Settings(args.draft)
    outdir = args.outdir or os.path.join(HERE, "figures")
    builders = (slab_geometry, slab_potential, shell_regions, shell_fields)
    only = set(args.only.split(",")) if args.only else None
    for build in builders:
        if only is None or build.__name__ in only:
            build(settings, outdir)
    run_mesh_cases(args, outdir, settings)


if __name__ == "__main__":
    main()
