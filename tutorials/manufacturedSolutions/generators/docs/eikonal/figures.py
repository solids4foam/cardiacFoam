#!/usr/bin/env python3
"""Figures of the eikonal manufactured solutions: the problem only.

    python figures.py [--outdir DIR] [--draft] [--mesh-case NAME=PATH ...]

Analytic figures need only the NumPy modules in generators/reference/.  Mesh
figures read constant/polyMesh of a meshed case the author supplies:

    cartesian_cube_hex, cartesian_cube_tet, shell_eikonal_hex, shell_eikonal_tet
"""

import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "_common"))

import figstyle as fs  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pyvista as pv  # noqa: E402

MESH_CASES = {
    "cartesian_cube_hex": "Cartesian eikonal cube, hex",
    "cartesian_cube_tet": "Cartesian eikonal cube, tet",
    "shell_eikonal_hex": "Eikonal shell, hex",
    "shell_eikonal_tet": "Eikonal shell, tet",
}
CUBE_CENTRE = (0.5, 0.5, 0.5)
OFFSETS = (0.3, 0.2, 0.4)


def slices_image(settings, func, clim):
    grid = fs.cube_grid(settings.n_cube)
    grid.point_data["f"] = func(*fs.point_array(grid))
    pl = fs.plotter(settings)
    pl.add_mesh(grid.slice_orthogonal(x=OFFSETS[0], y=OFFSETS[1], z=OFFSETS[2]),
                scalars="f", cmap=fs.CMAP, clim=clim, show_scalar_bar=False,
                lighting=False)
    fs.add_box(pl)
    fs.camera(pl, CUBE_CENTRE, 3.4)
    return fs.render(pl)


def cartesian_eikonal(settings, outdir):
    ref = fs.load_reference("cartesian_eikonal")
    sigma = ref.PARAMS["SIGMA_AXIS"]
    chi, cm = ref.PARAMS["CHI"], ref.PARAMS["CM"]
    psi = lambda x, y, z: ref.tau_d3(x, y, z)
    src = lambda x, y, z: ref.source_d3(x, y, z, 1.0, chi, cm, *sigma)
    grid = fs.cube_grid(settings.n_cube)
    pts = fs.point_array(grid)
    p_clim, s_clim = fs.span(psi(*pts)), fs.symmetric(src(*pts))
    fig, axes = plt.subplots(1, 3, figsize=(10.5, 3.9),
                             gridspec_kw={"width_ratios": [1, 1, 0.8]})
    fs.panel(fig, axes[0], slices_image(settings, psi, p_clim),
             title=r"$\psi=\exp(\mathbf{k}_3\cdot\mathbf{x})$", label=r"$\psi$",
             clim=p_clim, tag="(a)")
    fs.panel(fig, axes[1], slices_image(settings, src, s_clim),
             title=r"source $S_\psi$ at $c_0=1$", label=r"$S_\psi$", clim=s_clim,
             tag="(b)")
    ax = axes[2]
    x = np.linspace(0, 1, 200)
    for d, color in ((1, "black"), (2, fs.ELECTRODE), (3, fs.KERNEL)):
        zero = 0 * x
        args = (x, x if d >= 2 else zero, x if d == 3 else zero)
        ax.plot(x, getattr(ref, "tau_d%d" % d)(*args), color=color,
                label=r"$\psi$, $d=%d$" % d)
    ax.set_xlabel(r"$s$ along $x=y=z=s$ (up to $d$ coordinates)")
    ax.set_ylabel(r"$\psi$")
    ax.set_title("(c) $d=1,2,3$", loc="left")
    ax.legend(frameon=False, fontsize=8)
    fs.save(fig, outdir, "cartesian_eikonal.png", settings)


def template_u(ref, kind, s):
    """U(s) = integral of U' from s[0] (trapezoidal), U'(s) from the module."""
    if kind == "single_event":
        du = ref.template_derivative_single_event(s, 64.0, 4.0)
    else:
        du = ref.template_derivative_sine(s)
    return np.concatenate([[0.0], np.cumsum(0.5 * (du[1:] + du[:-1]) * np.diff(s))])


def eikonal_ecg(settings, outdir):
    ref = fs.load_reference("eikonal_ecg")
    sigma = np.array(ref.PARAMS["SIGMA_ROTATED"], float)
    electrodes = fs.tutorial_electrodes("eikonalECG")[:3]
    order = 12 if settings.draft else 32
    colors = ("black", fs.ELECTRODE, fs.KERNEL)
    fig, axes = plt.subplots(2, 2, figsize=(8.4, 5.6), gridspec_kw={"hspace": 0.5,
                                                                      "wspace": 0.3})
    specs = (("single_event", np.linspace(-10, 90, 600), np.linspace(0, 90, settings.n_time)),
             ("sine", np.linspace(-0.5, 2.0, 600), np.linspace(0, 3.0, settings.n_time)))
    for col, (kind, s_axis, t_axis) in enumerate(specs):
        name = "single event" if kind == "single_event" else "sine"
        ax = axes[0, col]
        ax.plot(s_axis, template_u(ref, kind, s_axis), color="black")
        ax.set_xlabel(r"$s$")
        ax.set_ylabel(r"$U(s)$")
        ax.set_title("(%s) template U, %s" % ("ab"[col], name), loc="left")
        ax = axes[1, col]
        for k, (e, color) in enumerate(zip(electrodes, colors)):
            phi = [ref.spatial_reference(3, kind, sigma, e, t, order) for t in t_axis]
            ax.plot(t_axis, phi, color=color, label=r"$\mathbf{e}_%d$" % (k + 1))
        ax.set_xlabel(r"$t$")
        ax.set_ylabel(r"$\phi(\mathbf{e},t)$")
        ax.set_title("(%s) exact eikonal ECG, %s" % ("cd"[col], name), loc="left")
        ax.margins(y=0.3)
        ax.legend(frameon=False, fontsize=8, ncol=3, loc="upper center")
    fs.save(fig, outdir, "eikonal_ecg.png", settings)


def shell_params():
    ref = fs.load_reference("shell_eikonal")
    return ref, ref.PARAMS["default"]


def shell_eikonal(settings, outdir):
    ref, p = shell_params()
    xs = np.array([p["xs0"], p["xs1"], p["xs2"]])
    phi_s = float(np.arctan2(xs[1], xs[0]))
    phi = (phi_s, phi_s + 1.5 * np.pi)
    direction = (np.cos(phi_s - np.pi / 4), np.sin(phi_s - np.pi / 4), 0.7)
    psi = lambda x, y, z: ref.field("psi", x, y, z)
    src = lambda x, y, z: ref.field("source", x, y, z)
    probe = fs.shell_grid(p["r1"], p["r2"], settings, phi)
    pts = fs.point_array(probe)
    p_clim, s_clim = fs.span(psi(*pts)), fs.symmetric(src(*pts))

    def pin(pl):
        pl.add_mesh(pv.Sphere(radius=p["rPin"], center=xs, theta_resolution=40,
                              phi_resolution=40), color=fs.ELECTRODE, opacity=0.55)
        pl.add_mesh(pv.Sphere(radius=0.02, center=xs), color="black")

    img_p, _ = fs.shell_cutaway(settings, ref, psi, p_clim, pin, direction, phi)
    glyph = fs.shell_frame_glyphs(ref, p, phi_s)
    img_g, _ = fs.shell_cutaway(settings, ref, psi, None, glyph, direction, phi)
    img_s, _ = fs.shell_cutaway(settings, ref, src, s_clim, None, direction, phi)
    fig, axes = plt.subplots(1, 3, figsize=(10.5, 3.9))
    fs.panel(fig, axes[0], img_p, title=r"$\psi$ at $c_0=%g$; pinned ball $r_{\mathrm{Pin}}=%g$ (orange)" % (p["c0"], p["rPin"]),
             label=r"$\psi$", clim=p_clim, tag="(a)")
    fs.panel(fig, axes[1], img_g, title=r"frame of $\boldsymbol{\sigma}$ through the wall",
             tag="(b)")
    fs.legend(axes[1], [("fibre", fs.FIBRE), ("sheet", fs.SHEET)], ncol=2)
    fs.panel(fig, axes[2], img_s, title=r"source $S_\psi$", label=r"$S_\psi$", clim=s_clim,
             tag="(c)")
    fs.save(fig, outdir, "shell_eikonal.png", settings)


def shell_eikonal_c0(settings, outdir):
    ref, p = shell_params()
    xs = np.array([p["xs0"], p["xs1"], p["xs2"]])
    phi_s = float(np.arctan2(xs[1], xs[0]))
    phi = (phi_s, phi_s + 1.5 * np.pi)
    direction = (np.cos(phi_s - np.pi / 4), np.sin(phi_s - np.pi / 4), 0.7)
    probe = fs.shell_grid(p["r1"], p["r2"], settings, phi)
    pts = fs.point_array(probe)
    fig, axes = plt.subplots(1, 2, figsize=(7.4, 3.9))
    for ax, c0, tag in zip(axes, (1.0, 5.0), ("(a)", "(b)")):
        func = lambda x, y, z, c0=c0: ref.field("psi", x, y, z, c0=c0)
        clim = fs.span(func(*pts))
        img, _ = fs.shell_cutaway(settings, ref, func, clim, None, direction, phi)
        fs.panel(fig, ax, img, title=r"$\psi$ at $c_0=%g$" % c0, label=r"$\psi$",
                 clim=clim, tag=tag)
    fs.save(fig, outdir, "shell_eikonal_c0.png", settings)


def main():
    ap = fs.cli(__doc__, MESH_CASES)
    args = ap.parse_args()
    settings = fs.Settings(args.draft)
    outdir = args.outdir or os.path.join(HERE, "figures")
    for build in (cartesian_eikonal, eikonal_ecg, shell_eikonal, shell_eikonal_c0):
        build(settings, outdir)
    fs.run_mesh_cases(args, MESH_CASES, outdir, settings)


if __name__ == "__main__":
    main()
