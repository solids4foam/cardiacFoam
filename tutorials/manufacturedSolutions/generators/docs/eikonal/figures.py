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
from matplotlib.gridspec import GridSpec  # noqa: E402

MESH_CASES = {
    "cartesian_cube_hex": "Cartesian eikonal cube, hex",
    "cartesian_cube_tet": "Cartesian eikonal cube, tet",
    "shell_eikonal_hex": "Eikonal shell, hex",
    "shell_eikonal_tet": "Eikonal shell, tet",
}
CUBE_CENTRE = (0.5, 0.5, 0.5)
OFFSETS = (0.3, 0.2, 0.4)
XL = r"$x/(1\,\mathrm{m})$"
TL = r"$t/(1\,\mathrm{s})$"
SL = r"$s/(1\,\mathrm{s})$"
SLICE_TITLE = r"$x=0.3,\ y=0.2,\ z=0.4$"


def slices_image(settings, func, clim):
    grid = fs.cube_grid(settings.n_cube)
    grid.point_data["f"] = func(*fs.point_array(grid))
    pl = fs.plotter(settings)
    pl.add_mesh(grid.slice_orthogonal(x=OFFSETS[0], y=OFFSETS[1], z=OFFSETS[2]),
                scalars="f", cmap=fs.CMAP, clim=clim, show_scalar_bar=False,
                lighting=False)
    fs.add_box(pl, color=fs.DOMAIN)
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
    fig, axes = plt.subplots(1, 3, figsize=(13.4, 4.8),
                             gridspec_kw={"width_ratios": [1, 1, 0.85], "wspace": 0.25})
    fs.panel(fig, axes[0], slices_image(settings, psi, p_clim),
             title=r"$\psi=\exp(\mathbf{k}_3\cdot\mathbf{x})$", label=r"$\psi$",
             clim=p_clim, tag="(a)")
    fs.panel(fig, axes[1], slices_image(settings, src, s_clim),
             title=r"$S_\psi$ at $c_0=1$", label=r"$S_\psi$", clim=s_clim, tag="(b)")
    fig.text(0.37, 0.02, "slices at " + SLICE_TITLE, ha="center", fontsize=10)
    ax = axes[2]
    x = np.linspace(0, 1, 200)
    zero = 0 * x
    for d, color in ((1, fs.PSI_LIGHT), (2, fs.PSI), (3, fs.PSI_DARK)):
        args = (x, x if d >= 2 else zero, x if d == 3 else zero)
        ax.plot(x, getattr(ref, "tau_d%d" % d)(*args), color=color,
                label=r"$d=%d$" % d)
    ax.set_xlabel(XL + r", with $y=x$ ($d\geq 2$), $z=x$ ($d=3$)", fontsize=10)
    ax.set_ylabel(r"$\psi$")
    ax.set_title(r"(c) $\psi$ for $d=1,2,3$", loc="left", fontsize=11)
    ax.legend(frameon=False)
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
    electrodes = fs.tutorial_electrodes("eikonalECG")
    shell = fs.shell_electrodes("eikonalECG")
    order = 12 if settings.draft else 32
    fig = plt.figure(figsize=(9.6, 10.4))
    gs = GridSpec(3, 2, figure=fig, hspace=0.4, wspace=0.3,
                  height_ratios=[0.85, 0.8, 0.8])
    fs.cube_electrode_axes(fig.add_subplot(gs[0, 0], projection="3d", computed_zorder=False),
                           electrodes, "(a) unit cube: %d electrodes" % len(electrodes))
    fs.shell_electrode_axes(fig.add_subplot(gs[0, 1], projection="3d", computed_zorder=False),
                            shell, title=r"(b) sealed shell: $\mathrm{E1}$ to $\mathrm{E5}$")
    specs = (("single_event", np.linspace(-10, 90, 600), np.linspace(0, 90, settings.n_time),
              "single event (duration $64$, rise $4$)", "single event"),
             ("sine", np.linspace(-0.5, 2.0, 600), np.linspace(0, 3.0, settings.n_time),
              "sine", "sine"))
    styles = ((fs.ECG, "-"), (fs.ECG_LIGHT, "--"), (fs.ECG_DARK, ":"))
    for col, (kind, s_axis, t_axis, long_name, short) in enumerate(specs):
        ax = fig.add_subplot(gs[1, col])
        ax.plot(s_axis, template_u(ref, kind, s_axis), color=fs.VM)
        ax.set_xlabel(SL)
        ax.set_ylabel(r"$U(s)$")
        ax.set_title("(%s) template $U$: %s" % ("cd"[col], long_name), loc="left",
                     fontsize=11)
        ax = fig.add_subplot(gs[2, col])
        for e, (color, ls) in zip(electrodes[:3], styles):
            phi = [ref.spatial_reference(3, kind, sigma, e, t, order) for t in t_axis]
            ax.plot(t_axis, phi, color=color, ls=ls,
                    label=r"$\mathbf{e}=(%.2f,\,%.2f,\,%.2f)$" % tuple(e))
        ax.set_xlabel(TL)
        ax.set_ylabel(r"$\Phi_e(\mathbf{e},t)$")
        ax.set_title("(%s) exact eikonal ECG, %s" % ("ef"[col], short), loc="left",
                     fontsize=11)
        ax.margins(y=0.6)
        ax.legend(frameon=False, fontsize=8, loc="upper center")
    fs.save(fig, outdir, "eikonal_ecg.png", settings)


def shell_params():
    ref = fs.load_reference("shell_eikonal")
    return ref, ref.PARAMS["default"]


def shell_geometry(p):
    xs = np.array([p["xs0"], p["xs1"], p["xs2"]])
    phi_s = float(np.arctan2(xs[1], xs[0]))
    phi = (phi_s, phi_s + 1.5 * np.pi)
    direction = (np.cos(phi_s - np.pi / 4), np.sin(phi_s - np.pi / 4), 0.7)
    return xs, phi_s, phi, direction


def shell_eikonal(settings, outdir):
    ref, p = shell_params()
    xs, phi_s, phi, direction = shell_geometry(p)
    psi = lambda x, y, z: ref.field("psi", x, y, z)
    src = lambda x, y, z: ref.field("source", x, y, z)
    probe = fs.shell_grid(p["r1"], p["r2"], settings, phi)
    pts = fs.point_array(probe)
    p_clim, s_clim = fs.span(psi(*pts)), fs.symmetric(src(*pts))

    def pin(pl):
        pl.add_mesh(pv.Sphere(radius=p["rPin"], center=xs, theta_resolution=40,
                              phi_resolution=40), color=fs.PIN, opacity=0.25)
        pl.add_mesh(pv.Sphere(radius=0.035, center=xs), color="black")

    img_p, _ = fs.shell_cutaway(settings, ref, psi, p_clim, pin, direction, phi)
    glyph = fs.shell_frame_glyphs(ref, p, phi_s)
    img_g, _ = fs.shell_cutaway(settings, ref, psi, None, glyph, glyph.direction, phi)
    img_s, _ = fs.shell_cutaway(settings, ref, src, s_clim, None, direction, phi)
    fig, axes = plt.subplots(1, 3, figsize=(13.2, 4.8))
    fs.panel(fig, axes[0], img_p, title=r"$\psi$ at $c_0=%g$" % p["c0"], label=r"$\psi$",
             clim=p_clim, tag="(a)")
    fs.legend(axes[0], [("stimulus point $\\mathbf{x}_s$", "black", "o"),
                        (r"pinned ball, $r_{\mathrm{Pin}}=%g$" % p["rPin"], fs.PIN)],
              ncol=1, inside="lower left")
    fs.panel(fig, axes[1], img_g, title=r"frame of $\boldsymbol{\sigma}$ in the wall",
             tag="(b)")
    fs.legend(axes[1], [("fibre (light: endo, dark: epi)", fs.FIBRE),
                        ("sheet", fs.SHEET), ("sheet normal", fs.NORMAL)], ncol=1)
    fs.panel(fig, axes[2], img_s, title=r"$S_\psi$ at $c_0=%g$" % p["c0"],
             label=r"$S_\psi$", clim=s_clim, tag="(c)")
    fs.save(fig, outdir, "shell_eikonal.png", settings)


def shell_eikonal_c0(settings, outdir):
    ref, p = shell_params()
    _, _, phi, direction = shell_geometry(p)
    probe = fs.shell_grid(p["r1"], p["r2"], settings, phi)
    pts = fs.point_array(probe)
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 4.8))
    for ax, c0, tag in zip(axes, (1.0, 5.0), ("(a)", "(b)")):
        func = lambda x, y, z, c0=c0: ref.field("psi", x, y, z, c0=c0)
        clim = fs.span(func(*pts))
        img, _ = fs.shell_cutaway(settings, ref, func, clim, None, direction, phi)
        fs.panel(fig, ax, img, title=r"$\psi$ at $c_0=%g$ (own colour scale)" % c0,
                 label=r"$\psi$", clim=clim, tag=tag)
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
