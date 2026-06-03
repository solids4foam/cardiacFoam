"""Convergence plot for the manufactured eikonal + ECG verification.

Usage:
    python plot_convergence.py [postProcessing_dir]

Reads the CSV files produced by post_processing_manufactured_eikonal_ecg.py
and saves convergence_plot.pdf (and .png) in the same output directory.
"""

from __future__ import annotations

import csv
import math
import re
import sys
from pathlib import Path

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except ImportError:
    print("matplotlib is required: pip install matplotlib", file=sys.stderr)
    raise SystemExit(1)


# ── colour / marker scheme ────────────────────────────────────────────────────

DIM_COLOUR = {"1D": "#1f77b4", "2D": "#ff7f0e", "3D": "#2ca02c"}
NORM_STYLE = {
    "L1":   {"ls": "-",  "marker": "o"},
    "L2":   {"ls": "--", "marker": "s"},
    "Linf": {"ls": ":",  "marker": "^"},
}


# ── CSV helpers ───────────────────────────────────────────────────────────────

def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def _f(val: str) -> float:
    try:
        return float(val)
    except (TypeError, ValueError):
        return math.nan


# ── reference slope (order in terms of N, so O(h^p) = O(N^-p)) ───────────────

def _reference_slope(ax, n_vals, anchor_n, anchor_err, order, label, colour="grey"):
    xs = sorted(n_vals)
    ys = [anchor_err * (x / anchor_n) ** order for x in xs]
    ax.plot(xs, ys, ls=(0, (4, 4)), lw=1.0, color=colour, zorder=0)
    ax.text(xs[-1] * 1.06, ys[-1], label, fontsize=7, color=colour, va="center")


# ── load data ─────────────────────────────────────────────────────────────────

def _load_pimple_iters(output_dir: Path) -> dict[tuple[str, int], int]:
    iters: dict[tuple[str, int], int] = {}
    log_root = output_dir / "logs"
    if not log_root.is_dir():
        return iters
    pat           = re.compile(r"(?P<dim>\dD)_(?P<n>\d+)_cells")
    converged_pat = re.compile(r"PIMPLE: converged in (\d+) iterations")
    not_conv_pat  = re.compile(r"PIMPLE: not converged within (\d+) iterations")
    for log_dir in sorted(log_root.iterdir()):
        m = pat.search(log_dir.name)
        if not m:
            continue
        log = log_dir / "log.cardiacFoam"
        if not log.is_file():
            continue
        text = log.read_text(errors="ignore")
        mc = converged_pat.search(text)
        mn = not_conv_pat.search(text)
        if mc:
            iters[(m.group("dim"), int(m.group("n")))] = int(mc.group(1))
        elif mn:
            iters[(m.group("dim"), int(m.group("n")))] = int(mn.group(1))
    return iters


def _load(output_dir: Path):
    activation = _read_csv(output_dir / "manufacturedEikonalActivationSummary.csv")
    ecg_agg    = _read_csv(output_dir / "manufacturedEikonalECGAggregateSummary.csv")
    pimple     = _load_pimple_iters(output_dir)
    return activation, ecg_agg, pimple


# ── main plot ─────────────────────────────────────────────────────────────────

def plot(output_dir: Path) -> None:
    activation, ecg_agg, pimple_iters = _load(output_dir)

    dims = sorted({r["Dimension"] for r in activation if r.get("N", "?") != "?"})

    fig, axes = plt.subplots(1, 3, figsize=(13, 4.5))
    ax_act, ax_ecg, ax_pimple = axes

    for dim in dims:
        colour = DIM_COLOUR.get(dim, "black")

        # ── Activation-time error ─────────────────────────────────────────
        act_rows = sorted(
            [r for r in activation if r["Dimension"] == dim and r.get("N", "?") != "?"],
            key=lambda r: int(r["N"]),
        )
        if act_rows:
            for norm_name, style in NORM_STYLE.items():
                ns, errs = [], []
                for r in act_rows:
                    n   = int(r["N"])
                    err = _f(r.get(f"activation_{norm_name}", ""))
                    if not math.isnan(err):
                        ns.append(n)
                        errs.append(err)
                if ns:
                    ax_act.plot(
                        ns, errs,
                        color=colour, ls=style["ls"], marker=style["marker"],
                        ms=5, label=f"{dim} {norm_name}",
                    )

            # reference slopes anchored at coarsest point
            pairs = [
                (int(r["N"]), _f(r.get("activation_L1", "")))
                for r in act_rows
                if not math.isnan(_f(r.get("activation_L1", "")))
            ]
            if pairs:
                n0, e0 = pairs[0]
                all_n  = [p[0] for p in pairs]
                _reference_slope(ax_act, all_n, n0, e0, -2.0, "O(N⁻²)")
                _reference_slope(ax_act, all_n, n0, e0, -1.0, "O(N⁻¹)", colour="#aaaaaa")

        # ── ECG max / mean Linf error ─────────────────────────────────────
        ecg_rows = sorted(
            [r for r in ecg_agg if r["Dimension"] == dim and r.get("N", "?") != "?"],
            key=lambda r: int(r["N"]),
        )
        if ecg_rows:
            ns_max, errs_max, ns_mean, errs_mean = [], [], [], []
            for r in ecg_rows:
                n        = int(r["N"])
                err_max  = _f(r.get("max_Linf_err_ref",  ""))
                err_mean = _f(r.get("mean_Linf_err_ref", ""))
                if not math.isnan(err_max):
                    ns_max.append(n)
                    errs_max.append(err_max)
                if not math.isnan(err_mean):
                    ns_mean.append(n)
                    errs_mean.append(err_mean)
            if ns_max:
                ax_ecg.plot(ns_max,  errs_max,  color=colour, ls="-",  marker="o", ms=5, label=f"{dim} max Linf")
                ax_ecg.plot(ns_mean, errs_mean, color=colour, ls="--", marker="s", ms=5, label=f"{dim} mean Linf")

        # ── PIMPLE iteration count ────────────────────────────────────────
        pimple_ns   = []
        pimple_vals = []
        for r in act_rows:
            n = int(r["N"])
            it = pimple_iters.get((dim, n))
            if it is not None:
                pimple_ns.append(n)
                pimple_vals.append(it)
        if pimple_ns:
            ax_pimple.plot(pimple_ns, pimple_vals, color=colour, ls="-", marker="o", ms=5, label=dim)

    # ── axis formatting ───────────────────────────────────────────────────────

    for ax in (ax_act, ax_ecg, ax_pimple):
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("N  (cells)")
        ax.grid(True, which="both", ls=":", lw=0.4)

    ax_act.set_ylabel("error")
    ax_act.set_title("Activation time error vs N")
    ax_act.legend(fontsize=7, ncol=1)

    ax_ecg.set_ylabel("error")
    ax_ecg.set_title("ECG Linf error vs N\n(numeric vs exact manufactured reference)")
    ax_ecg.legend(fontsize=7)

    ax_pimple.set_ylabel("PIMPLE iterations")
    ax_pimple.set_title("Nonlinear solver iterations vs N")
    ax_pimple.legend(fontsize=7)

    fig.suptitle("Manufactured eikonal ECG — convergence study", fontsize=11)
    fig.tight_layout()

    out_pdf = output_dir / "convergence_plot.pdf"
    out_png = output_dir / "convergence_plot.png"
    fig.savefig(out_pdf, bbox_inches="tight")
    fig.savefig(out_png, bbox_inches="tight", dpi=150)
    print(f"Saved: {out_pdf}")
    print(f"Saved: {out_png}")
    plt.close(fig)


def plot_quadrature(output_dir: Path) -> None:
    """Separate figure: quadrature delta vs N for each check order.

    Shows |ECG(q_check) - ECG(q_ref)| aggregated over electrodes.
    Flat lines = quadrature accuracy independent of mesh (correct).
    Lower lines = higher q is closer to the reference (quadrature converging).
    """
    quad_csv = output_dir / "manufacturedEikonalECGQuadratureSummary.csv"
    if not quad_csv.is_file():
        print(f"Quadrature summary not found, skipping: {quad_csv}")
        return

    rows = _read_csv(quad_csv)
    if not rows:
        return

    # discover q orders from column names
    q_pat = re.compile(r"max_Linf_delta_q(\d+)_ref")
    q_checks = sorted(
        {int(m.group(1)) for r in rows for k in r if (m := q_pat.fullmatch(k))}
    )
    if not q_checks:
        return

    dims = sorted({r["Dimension"] for r in rows if r.get("N", "?") != "?"})

    # colour for q orders (darker = higher order)
    import colorsys
    def _q_colour(q, q_list):
        idx = q_list.index(q) / max(len(q_list) - 1, 1)
        h, s, v = 0.60, 0.8, 0.3 + 0.6 * idx   # blue family, light→dark
        r, g, b = colorsys.hsv_to_rgb(h, s, v)
        return (r, g, b)

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    ax_max, ax_mean = axes

    for dim in dims:
        dim_rows = sorted(
            [r for r in rows if r["Dimension"] == dim and r.get("N", "?") != "?"],
            key=lambda r: int(r["N"]),
        )
        if not dim_rows:
            continue
        dim_marker = {"1D": "o", "2D": "s", "3D": "^"}.get(dim, "o")

        for q in q_checks:
            colour = _q_colour(q, q_checks)
            ns, maxs, means = [], [], []
            for r in dim_rows:
                n = int(r["N"])
                mx = _f(r.get(f"max_Linf_delta_q{q}_ref", ""))
                mn = _f(r.get(f"mean_Linf_delta_q{q}_ref", ""))
                if not math.isnan(mx):
                    ns.append(n); maxs.append(mx); means.append(mn)

            lbl = f"{dim} q={q}"
            ax_max.plot(ns, maxs, color=colour, marker=dim_marker, ms=5,
                        ls="-", label=lbl)
            ax_mean.plot(ns, means, color=colour, marker=dim_marker, ms=5,
                         ls="-", label=lbl)

    for ax, title in (
        (ax_max,  "Max electrode |ECG(q) − ECG(q=96)| vs N"),
        (ax_mean, "Mean electrode |ECG(q) − ECG(q=96)| vs N"),
    ):
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("N  (cells)")
        ax.set_ylabel("|Δ ECG|  (quadrature delta)")
        ax.set_title(title)
        ax.legend(fontsize=7, ncol=2)
        ax.grid(True, which="both", ls=":", lw=0.4)

    fig.suptitle(
        "Manufactured eikonal ECG — quadrature convergence\n"
        r"Template voltage: $U(s)=\sin(2\pi s)$,  $s = t - \tau(\mathbf{x})$",
        fontsize=10,
    )
    fig.tight_layout()

    out_pdf = output_dir / "quadrature_plot.pdf"
    out_png = output_dir / "quadrature_plot.png"
    fig.savefig(out_pdf, bbox_inches="tight")
    fig.savefig(out_png, bbox_inches="tight", dpi=150)
    print(f"Saved: {out_pdf}")
    print(f"Saved: {out_png}")
    plt.close(fig)


def main(argv: list[str]) -> int:
    output_dir = Path(argv[1]) if len(argv) > 1 else Path("postProcessing")
    plot(output_dir)
    plot_quadrature(output_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
