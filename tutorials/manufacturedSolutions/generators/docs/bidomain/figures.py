#!/usr/bin/env python3
"""Figures of the bidomain manufactured solutions: the problem only.

    python figures.py [--outdir DIR] [--draft] [--mesh-case NAME=PATH ...]

The analytic figures use bidomain_reference.py (FDA cube through
reference/fda_bidomain.py, sealed shell through shell_bidomain_mms.py), the
FDA bath slab through reference/fda_bath_bidomain.py and the bath shell through
bath_shell_bidomain_mms.py.  Mesh
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
import bath_shell_bidomain_mms as bsb  # noqa: E402

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
    ax.bar(idx - 0.18, br.GI, width=0.36, color=fs.VM, label=r"$\sigma_{i,mm}$")
    ax.bar(idx + 0.18, br.GE, width=0.36, color=PHIE, label=r"$\sigma_{e,mm}$")
    ax.set_xticks(idx, [r"$xx$", r"$yy$", r"$zz$"])
    ax.set_ylabel(r"S m$^{-1}$")
    ax.set_title(r"(d) diagonal $\boldsymbol{\sigma}_i$, $\boldsymbol{\sigma}_e$", loc="left")
    ax.legend(frameon=False, fontsize=9, loc="upper right")
    ax.text(0.02, 0.95, r"$\frac{\sigma_{i,mm}}{\sigma_{i,mm}+\sigma_{e,mm}}=\kappa$", transform=ax.transAxes,
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
    fs.panel(fig, axes[2], img_g, title=r"frame of $\boldsymbol{\sigma}_i$ and $\boldsymbol{\sigma}_e$",
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
            label=r"$\mathbf{n}\cdot\boldsymbol{\sigma}_i\nabla F_i$")
    ax.plot(deg, out["Fe"][1], color=PHIE, ls=":", lw=2.2,
            label=r"$\mathbf{n}\cdot\boldsymbol{\sigma}_e\nabla F_e$")
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
             label=r"$\nabla\!\cdot(\boldsymbol{\sigma}_i\nabla F_i)+\nabla\!\cdot(\boldsymbol{\sigma}_e\nabla F_e)$",
             clim=s_clim, tag="(a)")
    fs.panel(fig, fig.add_subplot(gs[0, 1]), img_d, title=r"diffusion term of $I_{\mathrm{app}}$ at $t=0$",
             label=r"$-\nabla\!\cdot(\boldsymbol{\sigma}_i\nabla F_i)$", clim=d_clim, tag="(b)")
    ax1 = fig.add_subplot(gs[1, 0])
    ax2 = fig.add_subplot(gs[1, 1], sharey=ax1)
    m1 = wall_profiles(ax1, p["r1"], -1, r"(c) inner sphere $|\mathbf{x}|=r_1$ (endocardium)")
    m2 = wall_profiles(ax2, p["r2"], 1, r"(d) outer sphere $|\mathbf{x}|=r_2$ (epicardium)")
    assert max(m1, m2) < 1e-12, "conormal flux on a wall"
    ax1.legend(frameon=False, fontsize=9, ncol=2, loc="upper center",
               bbox_to_anchor=(1.1, -0.22))
    fs.save(fig, outdir, "shell_bidomain_sources.png", settings)


def _heart_bath_axes(ax):
    """Dotted lines at the two heart walls and the slab names on the x axis."""
    for wall in (0.0, 1.0):
        ax.axvline(wall, color=fs.GREY, lw=0.8, ls=":")
    ax.set_xlim(-1, 2)
    ax.set_xlabel(XL)


def bath_slab(settings, outdir):
    """FDA bath slab: geometry, potentials and current densities at t = 0."""
    ref = fs.load_reference("fda_bath_bidomain")
    par = (ref.PARAMS["K"], ref.PARAMS["ALPHA"], ref.PARAMS["SE"])
    gi_xx, se, sb = br.GI[0], ref.PARAMS["SE"], ref.PARAMS["SIGMA_B"]
    alpha = ref.PARAMS["ALPHA"]
    x = np.linspace(-1, 2, 1201)
    xh = x[(x >= 0) & (x <= 1)]
    z0, zh = 0 * x, 0 * xh
    pe = ref.phiE_d1(x, z0, z0, 0.0, *par)
    vm = ref.Vm_d1(xh, zh, zh, 0.0, *par)
    pi = ref.phiI_d1(xh, zh, zh, 0.0, *par)
    h = 1e-7

    def slope(f, xs, lo, hi):
        """Central difference of f inside one slab (lo, hi): the points are kept
        2h inside it, since phi_e has a kink on each heart wall."""
        xc = np.clip(xs, lo + 2 * h, hi - 2 * h)
        zs = 0 * xc
        return (f(xc + h, zs, zs, 0.0, *par) - f(xc - h, zs, zs, 0.0, *par)) / (2 * h)

    left, right = x < 0, x > 1
    heart = ~left & ~right
    j_e = np.concatenate([sb * slope(ref.phiE_d1, x[left], -1, 0),
                          se * slope(ref.phiE_d1, x[heart], 0, 1),
                          sb * slope(ref.phiE_d1, x[right], 1, 2)])
    j_i = gi_xx * slope(ref.phiI_d1, xh, 0, 1)
    assert abs(j_i[0]) < 1e-4 * alpha and abs(j_i[-1]) < 1e-4 * alpha, "intracellular current on a wall"
    assert np.allclose(j_e[left], alpha, atol=1e-6) and np.allclose(j_e[right], alpha, atol=1e-6)
    assert abs(j_e[heart][0] - alpha) < 1e-4 * alpha and abs(j_e[heart][-1] - alpha) < 1e-4 * alpha

    fig, axes = plt.subplots(1, 3, figsize=(12.8, 3.9), gridspec_kw={"width_ratios": [1.15, 1, 1]})
    ax = axes[0]
    for lo, hi, colour, name in ((-1, 0, fs.MESH_FACE, "bath"), (0, 1, "#b7cde3", "myocardium"),
                                 (1, 2, fs.MESH_FACE, "bath")):
        ax.add_patch(plt.Rectangle((lo, 0), hi - lo, 1, facecolor=colour, edgecolor=fs.MESH_EDGE, lw=1.0))
        ax.text((lo + hi) / 2, 0.62, name, ha="center", va="center", fontsize=10)
    ax.text(-0.5, 0.36, r"$\sigma_b$", ha="center", va="center", fontsize=10)
    ax.text(0.5, -0.05, r"$\sigma_b=\sigma_{e,xx}/2$", ha="center", va="top", fontsize=10)
    ax.text(1.5, 0.36, r"$\sigma_b$", ha="center", va="center", fontsize=10)
    ax.text(0.5, 0.36, r"$\boldsymbol{\sigma}_i,\ \boldsymbol{\sigma}_e$", ha="center", va="center", fontsize=10)
    ax.plot([-1.0, -1.0], [0, 1], color="black", lw=3.0, solid_capstyle="butt")
    ax.text(-1.0, 1.1, r"grounded, $\phi_e=0$", ha="left", va="bottom", fontsize=10)
    ax.annotate("", xy=(2.0, 0.5), xytext=(2.4, 0.5),
                arrowprops=dict(arrowstyle="-|>", lw=1.8, color="black", shrinkA=0, shrinkB=0))
    ax.text(2.0, 1.1, r"current $\alpha$ in", ha="right", va="bottom", fontsize=10)
    ax.set_xlim(-1.1, 2.5)
    ax.set_ylim(-0.25, 1.35)
    ax.set_aspect("equal")
    ax.set_axis_off()
    ax.set_title("(a) heart between two bath slabs", loc="left", fontsize=11)
    ax = axes[1]
    ax.plot(xh, vm, color=fs.VM, label=r"$V_m$")
    ax.plot(x, pe, color=PHIE, label=r"$\phi_e$")
    ax.plot(xh, pi, color=PHII, label=r"$\phi_i$")
    _heart_bath_axes(ax)
    ax.set_title(r"(b) potentials at $t=0$", loc="left", fontsize=11)
    ax.legend(frameon=False, fontsize=9, ncol=3, loc="upper center", bbox_to_anchor=(0.5, -0.2))
    ax = axes[2]
    ax.plot(x, j_e / alpha, color=PHIE, label=r"$\sigma_e\,\partial_x\phi_e$ (bath: $\sigma_b\,\partial_x\phi_e$)")
    ax.plot(xh, j_i / alpha, color=PHII, label=r"$\sigma_{i,xx}\,\partial_x\phi_i$")
    ax.axhline(0, color=fs.GREY, lw=0.6)
    _heart_bath_axes(ax)
    ax.set_ylabel(r"current density$/\alpha$")
    ax.set_title(r"(c) current densities at $t=0$", loc="left", fontsize=11)
    ax.legend(frameon=False, fontsize=9, loc="upper center", bbox_to_anchor=(0.5, -0.2))
    fig.subplots_adjust(wspace=0.4)
    fs.save(fig, outdir, "fda_bath_slab.png", settings)


def _t6(components):
    """Full tensor of the six stored components (xx, xy, xz, yy, yz, zz)."""
    c = components
    return np.stack([np.stack([c[..., 0], c[..., 1], c[..., 2]], -1),
                     np.stack([c[..., 1], c[..., 3], c[..., 4]], -1),
                     np.stack([c[..., 2], c[..., 4], c[..., 5]], -1)], -2)


def bath_shell(settings, outdir):
    """Bath shell: phi_e and V_m in the plane z = 0 and the wall currents."""
    params = bsb.default_params()
    model = bsb.build("insulated", params)
    r1, r2, r3 = params["r1"], params["r2"], params["r3"]
    n = 161 if settings.draft else 481
    s = np.linspace(-r3, r3, n)
    X, Y = np.meshgrid(s, s)
    pts = np.stack([X, Y, 0 * X], -1)
    r = np.hypot(X, Y)
    phie = np.full(r.shape, np.nan)
    vm = np.full(r.shape, np.nan)
    for zone, mask in ((0, r < r1), (1, (r >= r1) & (r <= r2)), (2, (r > r2) & (r <= r3))):
        phie[mask] = model.eval_xyz(zone, "phiE", pts[mask], 0.0)
        if zone == 1:
            vm[mask] = model.eval_xyz(1, "Vm", pts[mask], 0.0)
    a = np.linspace(0, 2 * np.pi, 361)
    flux = {}
    for j, sign in ((1, -1.0), (2, 1.0)):
        radius = params["r%d" % j]
        ring = np.stack([radius * np.cos(a), radius * np.sin(a), 0 * a], -1)
        nrm = sign * ring / radius
        gi = _t6(model.eval_xyz(1, "GiZ", ring, 0.0))
        g_vm = model.eval_xyz(1, "gradVm", ring, 0.0)
        g_pi = model.eval_xyz(1, "gradPhiI", ring, 0.0)
        flux[j] = (np.einsum("pi,pij,pj->p", nrm, gi, g_vm), np.einsum("pi,pij,pj->p", nrm, gi, g_pi))
        assert np.abs(flux[j][1]).max() < 1e-12, "intracellular conormal flux on a wall"
    fig, axes = plt.subplots(1, 3, figsize=(14.0, 4.6), constrained_layout=True,
                             gridspec_kw={"width_ratios": [1, 1, 1.25]})
    for ax, field, title, label, tag in ((axes[0], phie, r"$\phi_e$ at $t=0$", r"$\phi_e$", "(a)"),
                                         (axes[1], vm, r"$V_m$ at $t=0$", r"$V_m$", "(b)")):
        im = ax.imshow(field, origin="lower", extent=(-r3, r3, -r3, r3), cmap=fs.CMAP)
        for radius in (r1, r2):
            ax.add_patch(plt.Circle((0, 0), radius, fill=False, ls="--", lw=0.9, color="white"))
        ax.set_title("%s %s, plane $z=0$" % (tag, title), loc="left", fontsize=11)
        ax.set_xlabel(r"$x/(1\,\mathrm{m})$")
        if tag == "(a)":
            ax.set_ylabel(r"$y/(1\,\mathrm{m})$")
        fs.colourbar_2d(fig, ax, im, label)
    axes[2].set_ylabel("conormal flux")
    ax = axes[2]
    deg = np.degrees(a)
    ax.plot(deg, flux[1][0], color=fs.VM, label=r"$r_1$: $\mathbf{n}\cdot\boldsymbol{\sigma}_i\nabla V_m$")
    ax.plot(deg, flux[2][0], color=PHII, label=r"$r_2$: $\mathbf{n}\cdot\boldsymbol{\sigma}_i\nabla V_m$")
    ax.plot(deg, flux[1][1], color="black", ls="--", lw=1.4, label=r"$\mathbf{n}\cdot\boldsymbol{\sigma}_i\nabla\phi_i$ (both)")
    ax.axhline(0, color=fs.GREY, lw=0.6)
    ax.set_xlim(0, 360)
    ax.set_xticks([0, 90, 180, 270, 360])
    ax.set_xlabel(r"angle from $+x$ in the plane $z=0$ (deg)")
    ax.set_title("(c) intracellular wall current", loc="left", fontsize=11)
    ax.legend(frameon=False, fontsize=9, loc="upper center", bbox_to_anchor=(0.5, -0.2))
    fs.save(fig, outdir, "bath_shell.png", settings)


def main():
    ap = fs.cli(__doc__, MESH_CASES)
    args = ap.parse_args()
    settings = fs.Settings(args.draft)
    outdir = args.outdir or os.path.join(HERE, "figures")
    for build in (fda_cube, shell_fields, shell_sources, bath_slab, bath_shell):
        build(settings, outdir)
    fs.run_mesh_cases(args, MESH_CASES, outdir, settings)


if __name__ == "__main__":
    main()
