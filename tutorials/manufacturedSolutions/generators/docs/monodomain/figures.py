#!/usr/bin/env python3
"""Figures of the monodomain manufactured solutions: the problem only.

    python figures.py [--outdir DIR] [--draft] [--mesh-case NAME=PATH ...]

Analytic figures need only the NumPy modules in generators/reference/.  Mesh
figures read constant/polyMesh of a meshed case the author supplies:

    fda_cube_hex, fda_cube_tet, rotated_cube_hex, rotated_cube_tet,
    shell_hex, shell_tet, graph_cube_hex, graph_cube_tet
"""

import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "_common"))

import figstyle as fs  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pyvista as pv  # noqa: E402
from matplotlib.gridspec import GridSpec  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

MESH_CASES = {
    "fda_cube_hex": "FDA cube, hex",
    "fda_cube_tet": "FDA cube, tet",
    "rotated_cube_hex": "Rotated cube, hex",
    "rotated_cube_tet": "Rotated cube, tet",
    "shell_hex": "Shell, hex",
    "shell_tet": "Shell, tet",
    "graph_cube_hex": "Graph cube, hex",
    "graph_cube_tet": "Graph cube, tet",
}
CUBE_CENTRE = (0.5, 0.5, 0.5)
XL = r"$x/(1\,\mathrm{m})$"
YL = r"$y/(1\,\mathrm{m})$"
ZL = r"$z/(1\,\mathrm{m})$"
TL = r"$t/(1\,\mathrm{s})$"
SLICE_TITLE = r"$x=0.3,\ y=0.2,\ z=0.4$"


def cube_slices(settings, func, clim, offsets=(0.3, 0.2, 0.4)):
    """Image of three orthogonal slices of func(x, y, z) on the unit cube."""
    grid = fs.cube_grid(settings.n_cube)
    grid.point_data["f"] = func(*fs.point_array(grid))
    pl = fs.plotter(settings)
    pl.add_mesh(grid.slice_orthogonal(x=offsets[0], y=offsets[1], z=offsets[2]),
                scalars="f", cmap=fs.CMAP, clim=clim, show_scalar_bar=False,
                lighting=False)
    fs.add_box(pl, color=fs.DOMAIN)
    fs.camera(pl, CUBE_CENTRE, 3.4)
    return fs.render(pl)


def fda_cube(settings, outdir):
    ref = fs.load_reference("fda_monodomain")
    img = cube_slices(settings, lambda x, y, z: ref.F_d3(x, y, z, 0.0), (-1.0, 1.0))
    fig = plt.figure(figsize=(8.6, 4.2))
    gs = GridSpec(2, 2, width_ratios=[1.5, 1], figure=fig, wspace=0.3, hspace=0.75)
    fs.panel(fig, fig.add_subplot(gs[:, 0]), img, title=r"$F_3$ on " + SLICE_TITLE,
             label=r"$F_3$", clim=(-1, 1), tag="(a)")
    x = np.linspace(0, 1, 200)
    ax = fig.add_subplot(gs[0, 1])
    ax.plot(x, ref.F_d1(x, 0 * x, 0 * x, 0.0), color=fs.VM)
    ax.set_xlabel(XL)
    ax.set_ylabel(r"$F_1$", fontsize=14)
    ax.set_title(r"(b) $d=1$", loc="left")
    xx, yy = np.meshgrid(x, x, indexing="ij")
    ax = fig.add_subplot(gs[1, 1])
    im = ax.pcolormesh(xx, yy, ref.F_d2(xx, yy, 0 * xx, 0.0), cmap=fs.CMAP, vmin=-1,
                       vmax=1, shading="auto", rasterized=True)
    ax.set_aspect("equal")
    ax.set_xlabel(XL)
    ax.set_ylabel(YL)
    ax.set_title(r"(c) $F_2$, $d=2$", loc="left")
    cb = fs.colourbar_2d(fig, ax, im, r"$F_2$")
    cb.set_ticks([-1, 0, 1])
    fs.save(fig, outdir, "fda_cube.png", settings)


def best_view(vectors, samples=4000, seed=0):
    """Unit view direction (positive z) maximising the smallest projected length of
    the given unit vectors, so that none is seen end-on."""
    rng = np.random.default_rng(seed)
    d = rng.normal(size=(samples, 3))
    d[:, 2] = np.abs(d[:, 2]) + 0.2
    d /= np.linalg.norm(d, axis=1)[:, None]
    v = np.array(vectors)
    proj = np.sqrt(np.maximum(0.0, 1.0 - (d @ v.T) ** 2))
    return d[np.argmax(proj.min(axis=1))]


def ellipsoid_glyph(pl, K):
    """Surface {K u : |u| = 1}, its principal axes (common length) and the cube
    axes.  Returns the eigenvalues, eigenvectors and a good view direction."""
    w, v = np.linalg.eigh(K)
    order = np.argsort(w)[::-1]
    w, v = w[order], v[:, order]
    sphere = pv.Sphere(radius=1.0, theta_resolution=60, phi_resolution=60)
    sphere.points = (K @ sphere.points.T).T / w[0]
    pl.add_mesh(sphere, color=fs.MESH_FACE, opacity=0.45, smooth_shading=True)
    for i, color in zip(range(3), (fs.FIBRE, fs.SHEET, fs.NORMAL)):
        for sign in (1, -1):
            pl.add_mesh(pv.Arrow(start=(0, 0, 0), direction=sign * v[:, i], scale=1.3,
                                 tip_length=0.2, tip_radius=0.06, shaft_radius=0.025),
                        color=color)
    fs.add_axes_arrows(pl, length=0.9, color="#222222")
    pl.add_point_labels(np.array([[1.12, 0, 0], [0, 1.12, 0], [0, 0, 1.12]]),
                        ["x", "y", "z"], font_size=26, shape=None, show_points=False,
                        text_color="#444444", font_family="arial")
    axes = [v[:, i] for i in range(3)] + [np.eye(3)[i] for i in range(3)]
    return w, v, best_view(axes)


def rotated_cube(settings, outdir):
    ref = fs.load_reference("rotated_monodomain")
    K = np.array(ref.PARAMS["K"], float)
    vm = lambda x, y, z: ref.Vm(x, y, z, 0.0)
    src = lambda x, y, z: ref.S(x, y, z, 0.0)
    grid = fs.cube_grid(settings.n_cube)
    s_clim = fs.symmetric(src(*fs.point_array(grid)))
    img_v = cube_slices(settings, vm, (0.0, 1.0))
    img_s = cube_slices(settings, src, s_clim)
    pl = fs.plotter(settings)
    w, v, view = ellipsoid_glyph(pl, K)
    fs.camera(pl, (0, 0, 0), 5.5, tuple(view))
    img_k = fs.render(pl)
    fig, axes = plt.subplots(1, 3, figsize=(12.6, 4.8))
    fs.panel(fig, axes[0], img_v, title=r"$V_m$ at $t=0$", label=r"$V_m$",
             clim=(0, 1), tag="(a)")
    fs.panel(fig, axes[1], img_k, title=r"axes of $\boldsymbol{\sigma}$ against $x,y,z$",
             tag="(b)")
    fs.legend(axes[1], [(r"$\mathbf{e}_%d$, $\lambda_%d=%.4f$" % (i + 1, i + 1, w[i]), c)
                        for i, c in enumerate((fs.FIBRE, fs.SHEET, fs.NORMAL))]
                       + [(r"cube axes $x,y,z$", "#222222")], ncol=2)
    fs.panel(fig, axes[2], img_s, title=r"$S$ at $t=0$", label=r"$S$", clim=s_clim,
             tag="(c)")
    fig.suptitle("slices at " + SLICE_TITLE, y=0.03, fontsize=10)
    fs.save(fig, outdir, "rotated_cube.png", settings)


def shell_monodomain(settings, outdir):
    ref = fs.load_reference("shell_monodomain")
    p = ref.PARAMS["default"]
    F = lambda x, y, z: ref.field("F", x, y, z)
    div = lambda x, y, z: -ref.field("divSigmaGradF", x, y, z)
    probe = fs.shell_grid(p["r1"], p["r2"], settings)
    f_clim = fs.span(F(*fs.point_array(probe)))
    d_clim = fs.symmetric(div(*fs.point_array(probe)))
    img_v, _ = fs.shell_cutaway(settings, ref, F, f_clim)
    glyph = fs.shell_frame_glyphs(ref, p)
    img_g, _ = fs.shell_cutaway(settings, ref, F, None, glyph, glyph.direction)
    img_s, _ = fs.shell_cutaway(settings, ref, div, d_clim)
    fig, axes = plt.subplots(1, 3, figsize=(12.6, 4.6))
    fs.panel(fig, axes[0], img_v, title=r"$V_m$ at $t=0$", label=r"$V_m$", clim=f_clim,
             tag="(a)")
    fs.panel(fig, axes[1], img_g, title=r"frame of $\boldsymbol{\sigma}$ in the wall",
             tag="(b)")
    fs.legend(axes[1], [("fibre (light: endo, dark: epi)", fs.FIBRE),
                        ("sheet", fs.SHEET), ("sheet normal", fs.NORMAL)], ncol=1)
    fs.panel(fig, axes[2], img_s, title=r"diffusion term of $S$ at $t=0$",
             label=r"$-\nabla\!\cdot(\boldsymbol{\sigma}\nabla F)$", clim=d_clim, tag="(c)")
    fs.save(fig, outdir, "shell_monodomain.png", settings)


def pseudo_ecg(settings, outdir):
    ref = fs.load_reference("pseudo_ecg")
    electrodes = fs.tutorial_electrodes("monodomainPseudoECG")
    shell = fs.shell_electrodes("monodomainPseudoECG")
    sig_fda = np.diag(ref.PARAMS["SIGMA_FDA"])
    sig_rot = np.array(ref.PARAMS["SIGMA_ROTATED"], float)
    order = settings.ecg_order if settings.draft else settings.ecg_order // 2 + 4
    phi = {}
    for case, sig in (("d3", sig_fda), ("rotated", sig_rot)):
        phi[case] = np.array([ref.spatial_reference(case, sig, e, order)
                              for e in electrodes])
    e = np.array(electrodes)
    fig = plt.figure(figsize=(9.6, 8.4))
    gs = GridSpec(2, 2, figure=fig, hspace=0.12, wspace=0.28)
    fs.cube_electrode_axes(fig.add_subplot(gs[0, 0], projection="3d", computed_zorder=False),
                           electrodes, "(a) unit cube: %d electrodes" % len(electrodes))
    fs.shell_electrode_axes(fig.add_subplot(gs[0, 1], projection="3d", computed_zorder=False),
                            shell, title=r"(b) sealed shell: $\mathrm{E1}$ to $\mathrm{E5}$")
    ax = fig.add_subplot(gs[1, 0])
    dist = np.linalg.norm(e - 0.5, axis=1)
    ax.semilogy(dist, np.abs(phi["d3"]), "D", color=fs.ECG, ms=4, label=r"FDA $d=3$")
    ax.semilogy(dist, np.abs(phi["rotated"]), "D", mfc="none", color=fs.ECG, ms=4,
                label="rotated")
    ax.set_ylim(1e-6, 0.5)
    ax.set_xlabel(r"$|\mathbf{e}-\mathbf{x}_c|/(1\,\mathrm{m})$")
    ax.set_ylabel(r"$|\Phi_e(\mathbf{e},0)|$")
    ax.set_title(r"(c) $\Phi_e$ against distance from the cube centre $\mathbf{x}_c$",
                 loc="left", fontsize=11)
    ax.legend(frameon=False)
    ax = fig.add_subplot(gs[1, 1])
    t = np.linspace(0.0, 2.0, 100)
    top = np.argsort(-np.abs(phi["d3"]))[:3]
    for k, color, ls in zip(top, (fs.ECG, fs.ECG_LIGHT, fs.ECG_DARK), ("-", "--", ":")):
        ax.plot(t, ref.time_factor(t) * phi["d3"][k], color=color, ls=ls,
                label=r"$\mathbf{e}=(%.2f,\,%.2f,\,%.2f)$" % tuple(e[k]))
    ax.set_xlabel(TL)
    ax.set_ylabel(r"$\Phi_e(\mathbf{e},t)$")
    ax.set_title(r"(d) $\Phi_e(\mathbf{e},t)=\sqrt{1+t}\,\Phi_e(\mathbf{e},0)$",
                 loc="left", fontsize=11)
    ax.margins(y=0.3)
    ax.legend(frameon=False, fontsize=9, loc="upper center",
              title=r"three largest $|\Phi_e(\mathbf{e},0)|$, FDA $d=3$", title_fontsize=9)
    fs.save(fig, outdir, "pseudo_ecg.png", settings)


def half_ball(centre, radius, side):
    """Half of the ball of radius `radius` at `centre` inside the cube: the part with
    (x - centre_x) * side >= 0."""
    ball = pv.Sphere(radius=radius, center=centre, theta_resolution=48, phi_resolution=48)
    half = ball.clip(normal=(side, 0, 0), origin=centre, invert=False)
    if half.n_points and np.all(side * (half.points[:, 0] - centre[0]) > -1e-9):
        return half
    return ball.clip(normal=(side, 0, 0), origin=centre, invert=True)


def graph_geometry(settings, outdir):
    ref = fs.load_reference("graph_monodomain")
    cpl = fs.load_reference("coupled_1d3d")
    root = np.array(ref.PARAMS["ROOT"])
    terminals = [np.array(tt) for tt in cpl.TERMINALS]
    zc = terminals[0][2]
    grid = fs.cube_grid(settings.n_cube)
    x, y, z = fs.point_array(grid)
    grid.point_data["f"] = cpl.Vm3D(x, y, z, 0.0)
    clim = fs.span(grid["f"])
    pl = fs.plotter(settings)
    pl.add_mesh(grid.slice(normal=(0, 0, 1), origin=(0.5, 0.5, zc)), scalars="f",
                cmap=fs.CMAP, clim=clim, show_scalar_bar=False, lighting=False)
    fs.add_box(pl, color=fs.DOMAIN)
    pl.add_mesh(pv.Line(terminals[0], terminals[1]).tube(radius=0.012), color=fs.GRAPH)
    pl.add_mesh(pv.Sphere(radius=0.03, center=root), color="black")
    for tt, side in zip(terminals, (1, -1)):
        pl.add_mesh(pv.Sphere(radius=0.028, center=tt), color=fs.GRAPH)
        pl.add_mesh(half_ball(tt, cpl.RADIUS, side), color=fs.KERNEL, opacity=0.55)
    fs.camera(pl, CUBE_CENTRE, 3.4, (0.8, -1.3, 1.0))
    img = fs.render(pl)
    fig, ax = plt.subplots(figsize=(6.0, 5.2))
    fs.panel(fig, ax, img, title=r"graph, kernels and $V_m$ on $z=1/3$",
             label=r"$V_m$ at $t=0$", clim=clim)
    fs.legend(ax, [("graph", fs.GRAPH),
                   (r"root $(\frac{1}{2},\frac{1}{6},\frac{1}{3})$", "black", "o"),
                   ("terminals", fs.GRAPH, "o"),
                   (r"kernel half-balls, radius $0.11$", fs.KERNEL)],
              ncol=1, inside="upper left")
    fs.save(fig, outdir, "graph_1d3d_geometry.png", settings)


def graph_fields(settings, outdir):
    ref = fs.load_reference("graph_monodomain")
    cpl = fs.load_reference("coupled_1d3d")
    t0, t1 = cpl.TERMINALS
    fig, axes = plt.subplots(1, 3, figsize=(13.6, 4.6),
                             gridspec_kw={"wspace": 0.5, "width_ratios": [1, 1, 1.15]})
    ax = axes[0]
    x = np.linspace(0, 1, 200)
    for t, color in ((0.0, fs.VM), (1.0, fs.VM_LIGHT)):
        ax.plot(x, ref.V1D(x, t), color=color)
        ax.plot(x, cpl.Vm3D(x, t0[1], t0[2], t), "--", color=color)
    ax.axvline(0.5, color=fs.GREY, lw=0.8, ls=":")
    ax.set_ylim(-1.7, 3.3)
    ax.text(0.51, -1.6, "root", color=fs.GREY, va="bottom")
    ax.set_xlabel(XL + " along the graph")
    ax.set_ylabel(r"$V_{1D}$, $V_m$")
    ax.set_title(r"(a) $V_{1D}$ and $V_m$ on the graph line", loc="left", fontsize=11)
    ax.legend(handles=[
        Line2D([0], [0], color=fs.NEUTRAL, label=r"$V_{1D}$"),
        Line2D([0], [0], color=fs.NEUTRAL, ls="--", label=r"$V_m(x,\frac{1}{6},\frac{1}{3})$"),
        Line2D([0], [0], color=fs.VM, label=r"$t/(1\,\mathrm{s})=0$"),
        Line2D([0], [0], color=fs.VM_LIGHT, label=r"$t/(1\,\mathrm{s})=1$")],
        frameon=False, fontsize=9, loc="upper center", ncol=2)
    ax = axes[1]
    r = np.linspace(0, 1.3 * cpl.RADIUS, 200)
    ax.plot(r, np.maximum(0.0, 1.0 - r / cpl.RADIUS), color=fs.NEUTRAL)
    ax.axvline(cpl.RADIUS, color=fs.GREY, lw=0.8, ls=":")
    ax.set_xlabel(r"distance from the terminal$/(1\,\mathrm{m})$")
    ax.set_ylabel(r"kernel weight $w_k$")
    ax.set_title("(b) linear kernel", loc="left", fontsize=11)
    ax = axes[2]
    n_quad = 6 if settings.draft else 16
    npts = 20 if settings.draft else 80
    xs = np.linspace(0, 1, npts)
    xx, zz = np.meshgrid(xs, xs, indexing="ij")
    pts = np.stack([xx.ravel(), np.full(xx.size, t0[1]), zz.ravel()], axis=1)
    vals = cpl.continuum_junction_source(pts, 0.0, n=n_quad).reshape(xx.shape)
    clim = fs.symmetric(vals)
    im = ax.pcolormesh(xx, zz, vals, cmap=fs.CMAP, vmin=clim[0], vmax=clim[1],
                       shading="auto", rasterized=True)
    ax.plot([0, 1], [t0[2], t0[2]], color=fs.GRAPH, lw=2)
    for tt, side in ((t0, 1), (t1, -1)):
        ang = np.linspace(-np.pi / 2, np.pi / 2, 60)
        ax.plot(tt[0] + side * cpl.RADIUS * np.cos(ang), tt[2] + cpl.RADIUS * np.sin(ang),
                color="white", lw=1)
    ax.set_aspect("equal")
    ax.set_xlabel(XL)
    ax.set_ylabel(ZL)
    ax.set_title(r"(c) $S_{\mathrm{pvj}}$ on $y=1/6$ at $t=0$", loc="left", fontsize=11)
    fs.colourbar_2d(fig, ax, im, r"$S_{\mathrm{pvj}}$")
    fs.save(fig, outdir, "graph_1d3d_fields.png", settings)


def main():
    ap = fs.cli(__doc__, MESH_CASES)
    args = ap.parse_args()
    settings = fs.Settings(args.draft)
    outdir = args.outdir or os.path.join(HERE, "figures")
    for build in (fda_cube, rotated_cube, shell_monodomain, pseudo_ecg,
                  graph_geometry, graph_fields):
        build(settings, outdir)
    fs.run_mesh_cases(args, MESH_CASES, outdir, settings)


if __name__ == "__main__":
    main()
