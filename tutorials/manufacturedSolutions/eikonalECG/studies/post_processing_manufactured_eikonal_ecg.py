"""Post-processing helpers for manufactured eikonal ECG runs."""

from __future__ import annotations

import csv
import math
from pathlib import Path
import re
import sys

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except ImportError:
    print("matplotlib is required: pip install matplotlib", file=sys.stderr)
    raise SystemExit(1)


CASE_PATTERN = re.compile(r"(?P<dimension>\dD)_(?P<cells>\d+)_cells")


def _case_metadata(path: Path) -> dict[str, str]:
    match = CASE_PATTERN.search(path.name)
    if match:
        return {
            "Dimension": match.group("dimension"),
            "N": match.group("cells"),
        }

    for parent in path.parents:
        match = CASE_PATTERN.search(parent.name)
        if match:
            return {
                "Dimension": match.group("dimension"),
                "N": match.group("cells"),
            }
    return {"Dimension": "unknown", "N": "unknown"}


def _as_float(value: str) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return math.nan


def _safe_rate(coarse_error: float, fine_error: float, coarse_n: int, fine_n: int) -> float:
    if (
        coarse_error <= 0.0
        or fine_error <= 0.0
        or coarse_n <= 0
        or fine_n <= coarse_n
        or not math.isfinite(coarse_error)
        or not math.isfinite(fine_error)
    ):
        return math.nan

    return math.log(coarse_error/fine_error)/math.log(float(fine_n)/float(coarse_n))


def _format_float(value: float) -> str:
    if not math.isfinite(value):
        return ""
    return f"{value:.8e}"


def _parse_activation_summary(path: Path) -> dict[str, str]:
    row: dict[str, str] = {}
    for line in path.read_text(encoding="utf-8", errors="ignore").splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        tokens = stripped.split()
        if tokens[0] == "activationTime" and len(tokens) >= 4:
            row.update(
                {
                    "activation_L1": tokens[1],
                    "activation_L2": tokens[2],
                    "activation_Linf": tokens[3],
                }
            )
    return row


def _parse_ecg_summary(path: Path) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    header: list[str] | None = None
    samples = ""
    tau0 = ""
    grad_tau = ""
    planar_fit_linf = ""

    for line in path.read_text(encoding="utf-8", errors="ignore").splitlines():
        stripped = line.strip()
        if not stripped:
            continue
        tokens = stripped.split()
        if tokens[0] == "samples" and len(tokens) >= 2:
            samples = tokens[1]
        elif tokens[0] == "tau0" and len(tokens) >= 2:
            tau0 = tokens[1]
        elif tokens[0] == "gradTau" and len(tokens) >= 2:
            grad_tau = " ".join(tokens[1:])
        elif tokens[0] == "planarFitLinf" and len(tokens) >= 2:
            planar_fit_linf = tokens[1]
        elif tokens[0] == "Electrode":
            header = tokens
        elif header and len(tokens) == len(header):
            row = dict(zip(header, tokens))
            row["samples"] = samples
            row["tau0"] = tau0
            row["gradTau"] = grad_tau
            row["planarFitLinf"] = planar_fit_linf
            rows.append(row)

    return rows


def _activation_convergence_rows(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    rates: list[dict[str, str]] = []
    grouped: dict[str, list[dict[str, str]]] = {}

    for row in rows:
        if row.get("N", "unknown") == "unknown":
            continue
        grouped.setdefault(row["Dimension"], []).append(row)

    for dimension, group in sorted(grouped.items()):
        ordered = sorted(group, key=lambda row: int(row["N"]))
        for coarse, fine in zip(ordered, ordered[1:]):
            coarse_n = int(coarse["N"])
            fine_n = int(fine["N"])
            rate_row = {
                "Dimension": dimension,
                "N_coarse": str(coarse_n),
                "N_fine": str(fine_n),
            }

            for field in ("activation_L1", "activation_L2", "activation_Linf"):
                rate_row[f"rate_{field}"] = _format_float(
                    _safe_rate
                    (
                        _as_float(coarse.get(field, "")),
                        _as_float(fine.get(field, "")),
                        coarse_n,
                        fine_n,
                    )
                )

            rates.append(rate_row)

    return rates


def _ecg_convergence_rows(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    rates: list[dict[str, str]] = []
    grouped: dict[tuple[str, str], list[dict[str, str]]] = {}

    for row in rows:
        if row.get("N", "unknown") == "unknown":
            continue
        key = (row["Dimension"], row["Electrode"])
        grouped.setdefault(key, []).append(row)

    for (dimension, electrode), group in sorted(grouped.items()):
        ordered = sorted(group, key=lambda row: int(row["N"]))
        for coarse, fine in zip(ordered, ordered[1:]):
            coarse_n = int(coarse["N"])
            fine_n = int(fine["N"])
            rates.append(
                {
                    "Dimension": dimension,
                    "Electrode": electrode,
                    "N_coarse": str(coarse_n),
                    "N_fine": str(fine_n),
                    "rate_L1_err_ref": _format_float(
                        _safe_rate
                        (
                            _as_float(coarse.get("L1_err_ref", "")),
                            _as_float(fine.get("L1_err_ref", "")),
                            coarse_n,
                            fine_n,
                        )
                    ),
                    "rate_L2_err_ref": _format_float(
                        _safe_rate
                        (
                            _as_float(coarse.get("L2_err_ref", "")),
                            _as_float(fine.get("L2_err_ref", "")),
                            coarse_n,
                            fine_n,
                        )
                    ),
                    "rate_Linf_err_ref": _format_float(
                        _safe_rate
                        (
                            _as_float(coarse.get("Linf_err_ref", "")),
                            _as_float(fine.get("Linf_err_ref", "")),
                            coarse_n,
                            fine_n,
                        )
                    ),
                }
            )

    return rates


_DELTA_Q_PAT = re.compile(r"Linf_delta_q(\d+)_ref")


def _q_checks_from_rows(rows: list[dict[str, str]]) -> list[int]:
    """Extract unique check quadrature orders from row column names."""
    orders: set[int] = set()
    for row in rows:
        for key in row:
            m = _DELTA_Q_PAT.fullmatch(key)
            if m:
                orders.add(int(m.group(1)))
    return sorted(orders)


def _ecg_quadrature_aggregate_rows(
    rows: list[dict[str, str]],
) -> tuple[list[dict[str, str]], list[int]]:
    """Aggregate |ECG(q_check) - ECG(q_ref)| delta columns over electrodes."""
    q_checks = _q_checks_from_rows(rows)
    if not q_checks:
        return [], []

    aggregated: list[dict[str, str]] = []
    grouped: dict[tuple[str, str], list[dict[str, str]]] = {}

    for row in rows:
        if row.get("N", "unknown") == "unknown":
            continue
        grouped.setdefault((row["Dimension"], row["N"]), []).append(row)

    for (dimension, n_value), group in sorted(
        grouped.items(),
        key=lambda item: (item[0][0], int(item[0][1])),
    ):
        summary = {"Dimension": dimension, "N": n_value}
        for q in q_checks:
            field = f"Linf_delta_q{q}_ref"
            values = [_as_float(row.get(field, "")) for row in group]
            values = [v for v in values if math.isfinite(v)]
            if values:
                summary[f"max_{field}"] = _format_float(max(values))
                summary[f"mean_{field}"] = _format_float(sum(values) / len(values))
            else:
                summary[f"max_{field}"] = ""
                summary[f"mean_{field}"] = ""
        aggregated.append(summary)

    return aggregated, q_checks


def _ecg_aggregate_rows(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    aggregated: list[dict[str, str]] = []
    grouped: dict[tuple[str, str], list[dict[str, str]]] = {}

    for row in rows:
        if row.get("N", "unknown") == "unknown":
            continue
        grouped.setdefault((row["Dimension"], row["N"]), []).append(row)

    for (dimension, n_value), group in sorted(
        grouped.items(),
        key=lambda item: (item[0][0], int(item[0][1])),
    ):
        summary = {"Dimension": dimension, "N": n_value}

        for field in ("L1_err_ref", "L2_err_ref", "Linf_err_ref", "planarFitLinf"):
            values = [_as_float(row.get(field, "")) for row in group]
            values = [value for value in values if math.isfinite(value)]
            if values:
                summary[f"max_{field}"] = _format_float(max(values))
                summary[f"mean_{field}"] = _format_float(sum(values)/len(values))
            else:
                summary[f"max_{field}"] = ""
                summary[f"mean_{field}"] = ""

        aggregated.append(summary)

    return aggregated


def _write_report(
    path: Path,
    activation_rows: list[dict[str, str]],
    activation_rate_rows: list[dict[str, str]],
    ecg_aggregate_rows: list[dict[str, str]],
    ecg_rate_rows: list[dict[str, str]],
) -> None:
    lines: list[str] = [
        "# Manufactured eikonal ECG post-processing",
        "",
        "## Activation-time error",
        "",
        "| Dimension | N | L1 | L2 | Linf |",
        "|---|---:|---:|---:|---:|",
    ]

    for row in sorted(
        activation_rows,
        key=lambda item: (item.get("Dimension", ""), int(item.get("N", "0"))),
    ):
        if row.get("N", "unknown") == "unknown":
            continue
        lines.append(
            "| {Dimension} | {N} | {activation_L1} | {activation_L2} | {activation_Linf} |".format(**row)
        )

    lines.extend(
        [
            "",
            "## Activation-time convergence rates",
            "",
            "| Dimension | N coarse | N fine | rate L1 | rate L2 | rate Linf |",
            "|---|---:|---:|---:|---:|---:|",
        ]
    )

    for row in activation_rate_rows:
        lines.append(
            "| {Dimension} | {N_coarse} | {N_fine} | {rate_activation_L1} | {rate_activation_L2} | {rate_activation_Linf} |".format(**row)
        )

    lines.extend(
        [
            "",
            "## ECG total error summary",
            "",
            "ECG errors are aggregated over electrodes. They compare numerical ECG computed from numerical psi against the manufactured reference ECG from exact psi.",
            "",
            "| Dimension | N | max Linf ECG | mean Linf ECG | max planar-fit Linf |",
            "|---|---:|---:|---:|---:|",
        ]
    )

    for row in ecg_aggregate_rows:
        lines.append(
            "| {Dimension} | {N} | {max_Linf_err_ref} | {mean_Linf_err_ref} | {max_planarFitLinf} |".format(**row)
        )

    lines.extend(
        [
            "",
            "## ECG convergence rates",
            "",
            "| Dimension | Electrode | N coarse | N fine | rate L1 | rate L2 | rate Linf |",
            "|---|---|---:|---:|---:|---:|---:|",
        ]
    )

    for row in ecg_rate_rows:
        lines.append(
            "| {Dimension} | {Electrode} | {N_coarse} | {N_fine} | {rate_L1_err_ref} | {rate_L2_err_ref} | {rate_Linf_err_ref} |".format(**row)
        )

    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _write_csv(path: Path, fieldnames: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in fieldnames})


def postprocess(output_dir: str | Path) -> None:
    root = Path(output_dir)
    activation_rows: list[dict[str, str]] = []
    ecg_rows: list[dict[str, str]] = []

    for summary in root.rglob("*manufacturedEikonalActivationTime.dat"):
        row = _case_metadata(summary)
        row.update(_parse_activation_summary(summary))
        row["source"] = str(summary)
        activation_rows.append(row)

    for summary in root.rglob("*manufacturedEikonalECGSummary.dat"):
        metadata = _case_metadata(summary)
        for row in _parse_ecg_summary(summary):
            row.update(metadata)
            row["source"] = str(summary)
            ecg_rows.append(row)

    activation_rows = [row for row in activation_rows if row.get("N", "unknown") != "unknown"]
    ecg_rows = [row for row in ecg_rows if row.get("N", "unknown") != "unknown"]

    activation_rows.sort(
        key=lambda row: (row.get("Dimension", ""), int(row.get("N", "0")))
    )
    ecg_rows.sort(
        key=lambda row:
        (
            row.get("Dimension", ""),
            int(row.get("N", "0")),
            row.get("Electrode", ""),
        )
    )

    activation_rate_rows = _activation_convergence_rows(activation_rows)
    ecg_rate_rows = _ecg_convergence_rows(ecg_rows)
    ecg_aggregate_rows = _ecg_aggregate_rows(ecg_rows)
    ecg_quadrature_rows, q_checks = _ecg_quadrature_aggregate_rows(ecg_rows)

    _write_csv(
        root / "manufacturedEikonalActivationSummary.csv",
        [
            "Dimension",
            "N",
            "activation_L1",
            "activation_L2",
            "activation_Linf",
            "source",
        ],
        activation_rows,
    )

    _write_csv(
        root / "manufacturedEikonalActivationConvergenceRates.csv",
        [
            "Dimension",
            "N_coarse",
            "N_fine",
            "rate_activation_L1",
            "rate_activation_L2",
            "rate_activation_Linf",
        ],
        activation_rate_rows,
    )

    _write_csv(
        root / "manufacturedEikonalECGSummary.csv",
        [
            "Dimension",
            "N",
            "Electrode",
            "samples",
            "tau0",
            "gradTau",
            "planarFitLinf",
            "L1_err_ref",
            "L2_err_ref",
            "Linf_err_ref",
            "source",
        ],
        ecg_rows,
    )

    _write_csv(
        root / "manufacturedEikonalECGAggregateSummary.csv",
        [
            "Dimension",
            "N",
            "max_L1_err_ref",
            "mean_L1_err_ref",
            "max_L2_err_ref",
            "mean_L2_err_ref",
            "max_Linf_err_ref",
            "mean_Linf_err_ref",
            "max_planarFitLinf",
            "mean_planarFitLinf",
        ],
        ecg_aggregate_rows,
    )

    _write_csv(
        root / "manufacturedEikonalECGConvergenceRates.csv",
        [
            "Dimension",
            "Electrode",
            "N_coarse",
            "N_fine",
            "rate_L1_err_ref",
            "rate_L2_err_ref",
            "rate_Linf_err_ref",
        ],
        ecg_rate_rows,
    )

    if ecg_quadrature_rows:
        q_delta_fields = ["Dimension", "N"]
        for q in q_checks:
            q_delta_fields += [f"max_Linf_delta_q{q}_ref", f"mean_Linf_delta_q{q}_ref"]
        _write_csv(
            root / "manufacturedEikonalECGQuadratureSummary.csv",
            q_delta_fields,
            ecg_quadrature_rows,
        )

    _write_report(
        root / "manufacturedEikonalECGReport.md",
        activation_rows,
        activation_rate_rows,
        ecg_aggregate_rows,
        ecg_rate_rows,
    )


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



def run_postprocessing(
    *,
    output_dir: str,
    setup_root: str | None = None,
    **_: object,
) -> list[dict]:
    """driverFOAM post-processing entry point (called for 'all' runs)."""
    root = Path(output_dir)
    postprocess(root)

    artifacts: list[dict] = [
        {
            "path": str(root / "manufacturedEikonalActivationSummary.csv"),
            "label": "Eikonal activation-time error table",
            "kind": "table",
            "format": "csv",
        },
        {
            "path": str(root / "manufacturedEikonalActivationConvergenceRates.csv"),
            "label": "Eikonal activation-time convergence rates",
            "kind": "table",
            "format": "csv",
        },
        {
            "path": str(root / "manufacturedEikonalECGAggregateSummary.csv"),
            "label": "Eikonal ECG aggregate error table",
            "kind": "table",
            "format": "csv",
        },
        {
            "path": str(root / "manufacturedEikonalECGConvergenceRates.csv"),
            "label": "Eikonal ECG convergence rates",
            "kind": "table",
            "format": "csv",
        },
        {
            "path": str(root / "manufacturedEikonalECGReport.md"),
            "label": "Eikonal ECG verification report",
            "kind": "report",
            "format": "markdown",
        },
        {
            "path": str(root / "manufacturedEikonalECGQuadratureSummary.csv"),
            "label": "Eikonal ECG quadrature delta table",
            "kind": "table",
            "format": "csv",
        },
    ]

    try:
        plot(root)
        plot_quadrature(root)
        for stem, label in (
            ("convergence_plot", "Convergence plot"),
            ("quadrature_plot",  "Quadrature convergence plot"),
        ):
            for suffix in ("pdf", "png"):
                p = root / f"{stem}.{suffix}"
                if p.exists():
                    artifacts.append(
                        {
                            "path": str(p),
                            "label": f"{label} ({suffix})",
                            "kind": "plot",
                            "format": suffix,
                        }
                    )
    except Exception as exc:  # noqa: BLE001
        print(f"Warning: plots skipped ({exc})")

    return artifacts


def main(argv: list[str]) -> int:
    output_dir = Path(argv[1]) if len(argv) > 1 else Path("postProcessing")
    postprocess(output_dir)
    try:
        plot(output_dir)
        plot_quadrature(output_dir)
    except Exception as exc:
        print(f"Warning: plots skipped ({exc})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
