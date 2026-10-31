#!/usr/bin/env python3
"""Figures of the bidomain manufactured solutions: the problem only.

    python figures.py [--outdir DIR] [--draft] [--mesh-case NAME=PATH ...]

The analytic figures use bidomain_reference.py (FDA cube through
reference/fda_bidomain.py, sealed shell through shell_bidomain_mms.py).  Mesh
figures read constant/polyMesh of a meshed case the author supplies:

    fda_cube_hex, fda_cube_tet, shell_hex, shell_tet
"""

import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "_common"))
sys.path.insert(0, HERE)

import figstyle as fs  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pyvista as pv  # noqa: E402
import vtk  # noqa: E402
from matplotlib.gridspec import GridSpec  # noqa: E402

import bidomain_reference as br  # noqa: E402

MESH_CASES = {
    "fda_cube_hex": "FDA cube, hex",
    "fda_cube_tet": "FDA cube, tet",
    "shell_hex": "Shell, hex",
    "shell_tet": "Shell, tet",
}
CUBE_CENTRE = (0.5, 0.5, 0.5)
XL = r"$x/(1\,\mathrm{m})$"
SLICE_TITLE = r"$x=0.3,\ y=0.2,\ z=0.4$"
PHIE = "#2ca02c"
PHII = "#d62728"


def cube_slices(settings, func, clim, offsets=(0.3, 0.2, 0.4), marker=None):
    grid = fs.cube_grid(settings.n_cube)
    grid.point_data["f"] = func(*fs.point_array(grid))
    pl = fs.plotter(settings)
    pl.add_mesh(grid.slice_orthogonal(x=offsets[0], y=offsets[1], z=offsets[2]),
                scalars="f", cmap=fs.CMAP, clim=clim, show_scalar_bar=False,
                lighting=False)
    fs.add_box(pl, color=fs.DOMAIN)
    fs.camera(pl, CUBE_CENTRE, 3.4)
    if marker is None:
        return fs.render(pl)
    pl.add_mesh(pv.Sphere(radius=0.035, center=marker), color="black")
    img = pl.screenshot(return_img=True)
    return crop_with_point(img, pl, marker)


def crop_with_point(img, pl, point, pad=6):
    """Crop the white margin as fs.render does and return the image with the
    pixel position (column, row) of the world point `point` in it."""
    coord = vtk.vtkCoordinate()
    coord.SetCoordinateSystemToWorld()
    coord.SetValue(*point)
    px, py = coord.GetComputedDoubleDisplayValue(pl.renderer)
    scale = img.shape[0] / pl.window_size[1]
    col, row = px * scale, img.shape[0] - py * scale
    pl.close()
    rows, cols = np.where(np.any(img < 250, axis=2))
    r0, c0 = max(rows.min() - pad, 0), max(cols.min() - pad, 0)
    img = img[r0:rows.max() + pad, c0:cols.max() + pad]
    return img, (col - c0, row - r0)


def fda_cube(settings, outdir):
    d = 3
    vm = lambda x, y, z: br.Vm(d, x, y, z, 0.0)
    pe = lambda x, y, z: br.phiE(d, x, y, z, 0.0)
    grid = fs.cube_grid(settings.n_cube)
    e_clim = fs.span(pe(*fs.point_array(grid)))
    img_v = cube_slices(settings, vm, (-1.0, 1.0))
    img_e, ref_px = cube_slices(settings, pe, e_clim, marker=br.PHIE_REF_POINT)
    fig = plt.figure(figsize=(12.8, 4.6))
    gs = GridSpec(2, 3, width_ratios=[1.25, 1.25, 1], figure=fig, wspace=0.25, hspace=0.8)
    fs.panel(fig, fig.add_subplot(gs[:, 0]), img_v, title=r"$V_m$ at $t=0$, $d=3$",
             label=r"$V_m$", clim=(-1, 1), tag="(a)")
    ax = fig.add_subplot(gs[:, 1])
    fs.panel(fig, ax, img_e, title=r"$\phi_e$ at $t=0$, $d=3$", label=r"$\phi_e$",
             clim=e_clim, tag="(b)")
    ax.annotate(r"$\mathbf{x}_{\mathrm{ref}}$, $\phi_e=0$", xy=ref_px,
                xytext=(ref_px[0] + 0.12 * img_e.shape[1], ref_px[1] + 0.1 * img_e.shape[0]),
                fontsize=10, va="center",
                arrowprops=dict(arrowstyle="->", color="black", lw=0.9, shrinkB=4))
    x = np.linspace(0, 1, 200)
    z0 = 0 * x
    ax = fig.add_subplot(gs[0, 2])
    ax.plot(x, br.Vm(1, x, z0, z0, 0.0), color=fs.VM, label=r"$V_m$")
    ax.plot(x, br.phiE(1, x, z0, z0, 0.0), color=PHIE, label=r"$\phi_e$")
    ax.plot(x, br.phiI(1, x, z0, z0, 0.0), color=PHII, label=r"$\phi_i$")
    ax.axvline(br.PHIE_REF_POINT[0], color=fs.GREY, lw=0.8, ls=":")
    ax.set_xlabel(XL)
    ax.set_title(r"(c) $d=1$, $t=0$", loc="left")
    ax.legend(frameon=False, fontsize=9, ncol=3, loc="upper center",
              bbox_to_anchor=(0.5, 1.02))
    ax.set_ylim(-1.2, 1.9)
    ax = fig.add_subplot(gs[1, 2])
    idx = np.arange(3)
    ax.bar(idx - 0.18, br.GI, width=0.36, color=fs.VM, label=r"$G_{i,mm}$")
    ax.bar(idx + 0.18, br.GE, width=0.36, color=PHIE, label=r"$G_{e,mm}$")
    ax.set_xticks(idx, [r"$xx$", r"$yy$", r"$zz$"])
    ax.set_ylabel(r"S m$^{-1}$")
    ax.set_title(r"(d) diagonal $\mathbf{G}_i$, $\mathbf{G}_e$", loc="left")
    ax.legend(frameon=False, fontsize=9, loc="upper right")
    ax.text(0.02, 0.95, r"$\frac{G_{i,mm}}{G_{i,mm}+G_{e,mm}}=\kappa$", transform=ax.transAxes,
            va="top", fontsize=10)
    ax.set_ylim(0, 0.17)
    fs.save(fig, outdir, "fda_bidomain_cube.png", settings)


def shell_fields(settings, outdir):
    ref = br.SHELL
    p = ref.PARAMS["default"]
    vm = lambda x, y, z: ref.field("Vm", x, y, z)
    mean_e = br.shell_volume_mean("Fe", *((10, 10, 20) if settings.draft else (40, 40, 80)))
    pe = lambda x, y, z: ref.field("Fe", x, y, z) - mean_e
    probe = fs.shell_grid(p["r1"], p["r2"], settings)
    v_clim = fs.span(vm(*fs.point_array(probe)))
    e_clim = fs.span(pe(*fs.point_array(probe)))
    img_v, _ = fs.shell_cutaway(settings, ref, vm, v_clim)
    img_e, _ = fs.shell_cutaway(settings, ref, pe, e_clim)
    glyph = fs.shell_frame_glyphs(ref, p)
    img_g, _ = fs.shell_cutaway(settings, ref, vm, None, glyph, glyph.direction)
    fig, axes = plt.subplots(1, 3, figsize=(12.6, 4.6))
    fs.panel(fig, axes[0], img_v, title=r"$V_m$ at $t=0$", label=r"$V_m$", clim=v_clim,
             tag="(a)")
    fs.panel(fig, axes[1], img_e, title=r"$\phi_e$ at $t=0$, zero volume mean",
             label=r"$\phi_e$", clim=e_clim, tag="(b)")
    fs.panel(fig, axes[2], img_g, title=r"frame of $\mathbf{G}_i$ and $\mathbf{G}_e$",
             tag="(c)")
    fs.legend(axes[2], [(r"$\mathbf{f}$: $\sigma_{i,f}=%g$, $\sigma_{e,f}=%g$"
                         % (p["gil"], p["gel"]), fs.FIBRE),
                        (r"$\mathbf{s}$: $\sigma_{i,s}=%g$, $\sigma_{e,s}=%g$"
                         % (p["git"], p["get"]), fs.SHEET),
                        (r"$\mathbf{n}_s$: $\sigma_{i,n}=%g$, $\sigma_{e,n}=%g$"
                         % (p["gin"], p["gen"]), fs.NORMAL)], ncol=1)
    fs.save(fig, outdir, "shell_bidomain.png", settings)


def wall_profiles(ax, radius, sign, title):
    """Normal derivative and conormal flux of Fi and Fe on the great circle of the
    sphere in the plane y = 0."""
    a = np.linspace(0, 2 * np.pi, 361)
    pts = np.stack([radius * np.sin(a), 0 * a, radius * np.cos(a)], -1)
    n = sign * pts / radius
    out = {}
    for F, G in (("Fi", "Gi"), ("Fe", "Ge")):
        g = br.sbm.evaluate("grad" + F, pts)
        T = br.sbm._full(br.sbm.evaluate(G, pts))
        out[F] = (np.sum(n * g, -1), np.einsum("pi,pij,pj->p", n, T, g))
    deg = np.degrees(a)
    ax.plot(deg, out["Fi"][0], color=fs.VM, label=r"$\mathbf{n}\cdot\nabla F_i$")
    ax.plot(deg, out["Fe"][0], color=PHIE, label=r"$\mathbf{n}\cdot\nabla F_e$")
    ax.plot(deg, out["Fi"][1], color=fs.VM, ls="--", lw=2.2,
            label=r"$\mathbf{n}\cdot\mathbf{G}_i\nabla F_i$")
    ax.plot(deg, out["Fe"][1], color=PHIE, ls=":", lw=2.2,
            label=r"$\mathbf{n}\cdot\mathbf{G}_e\nabla F_e$")
    ax.axhline(0, color=fs.GREY, lw=0.6)
    ax.set_xlim(0, 360)
    ax.set_xticks([0, 90, 180, 270, 360])
    ax.set_xlabel(r"angle from $+z$ in the plane $y=0$ (deg)")
    ax.set_title(title, loc="left", fontsize=11)
    return max(np.abs(out["Fi"][1]).max(), np.abs(out["Fe"][1]).max())


def shell_sources(settings, outdir):
    ref = br.SHELL
    p = ref.PARAMS["default"]
    se = lambda x, y, z: ref.field("divGiGradFi", x, y, z) + ref.field("divGeGradFe", x, y, z)
    di = lambda x, y, z: -ref.field("divGiGradFi", x, y, z)
    probe = fs.shell_grid(p["r1"], p["r2"], settings)
    s_clim = fs.symmetric(se(*fs.point_array(probe)))
    d_clim = fs.symmetric(di(*fs.point_array(probe)))
    img_s, _ = fs.shell_cutaway(settings, ref, se, s_clim)
    img_d, _ = fs.shell_cutaway(settings, ref, di, d_clim)
    fig = plt.figure(figsize=(12.6, 9.4))
    gs = GridSpec(2, 2, figure=fig, hspace=0.32, wspace=0.25, height_ratios=[1.25, 1])
    fs.panel(fig, fig.add_subplot(gs[0, 0]), img_s, title=r"$s_e$ at $t=0$",
             label=r"$\nabla\!\cdot(\mathbf{G}_i\nabla F_i)+\nabla\!\cdot(\mathbf{G}_e\nabla F_e)$",
             clim=s_clim, tag="(a)")
    fs.panel(fig, fig.add_subplot(gs[0, 1]), img_d, title=r"diffusion term of $I_{\mathrm{app}}$ at $t=0$",
             label=r"$-\nabla\!\cdot(\mathbf{G}_i\nabla F_i)$", clim=d_clim, tag="(b)")
    ax1 = fig.add_subplot(gs[1, 0])
    ax2 = fig.add_subplot(gs[1, 1], sharey=ax1)
    m1 = wall_profiles(ax1, p["r1"], -1, r"(c) inner sphere $|\mathbf{x}|=r_1$ (endocardium)")
    m2 = wall_profiles(ax2, p["r2"], 1, r"(d) outer sphere $|\mathbf{x}|=r_2$ (epicardium)")
    assert max(m1, m2) < 1e-12, "conormal flux on a wall"
    ax1.legend(frameon=False, fontsize=9, ncol=2, loc="upper center",
               bbox_to_anchor=(1.1, -0.22))
    fs.save(fig, outdir, "shell_bidomain_sources.png", settings)


def main():
    ap = fs.cli(__doc__, MESH_CASES)
    args = ap.parse_args()
    settings = fs.Settings(args.draft)
    outdir = args.outdir or os.path.join(HERE, "figures")
    for build in (fda_cube, shell_fields, shell_sources):
        build(settings, outdir)
    fs.run_mesh_cases(args, MESH_CASES, outdir, settings)


if __name__ == "__main__":
    main()
