"""Compare the pig/Gaur and human/TWORLD single-cell traces.

The script consumes the per-case directories produced by the focused driverFOAM
sweep. It creates an audience-facing six-panel waveform figure (Vm, Ca, ICaL,
SR fluxes, repolarisation currents, and Land--Niederer tension), a rate-
dependence figure, and a machine-readable metrics table. Model-specific names
are mapped explicitly because Gaur and TWORLD do not use the same names for
SR release and SERCA uptake.
"""
from __future__ import annotations

import argparse
import csv
import re
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import LinearLocator, ScalarFormatter


MODEL_INFO = {
    "Gaur": {
        "species": "Pig",
        "tissue": "myocyte",
        "jrel": "AV_Jrel",
        "jup": "AV_Jup",
        "ito": "AV_ITo_ITo",
    },
    "TWorld": {
        "species": "Human",
        "tissue": "endocardialCells",
        "jrel": "AV_J_SRCarel",
        "jup": "AV_Jserca",
        "ito": "AV_Ito_Ito",
    },
}
_CL_RE = re.compile(r"_S1_(\d+(?:p\d+)?)")


@dataclass
class Trace:
    model: str
    cl_ms: float
    data: dict[str, np.ndarray]
    tension: dict[str, np.ndarray]

    @property
    def species(self) -> str:
        return MODEL_INFO[self.model]["species"]


def _read_table(path: Path) -> dict[str, np.ndarray]:
    lines = path.read_text().splitlines()
    if not lines:
        raise ValueError(f"empty trace: {path}")
    names = lines[0].split()
    values = np.loadtxt(path, skiprows=1, ndmin=2)
    if values.shape[1] != len(names):
        raise ValueError(f"header/data mismatch in {path}")
    return {name: values[:, i] for i, name in enumerate(names)}


def _cl_from_name(path: Path) -> float:
    match = _CL_RE.search(path.name)
    if not match:
        raise ValueError(f"cannot determine pacing CL from {path.name}")
    return float(match.group(1).replace("p", "."))


def _find_column(data: dict[str, np.ndarray], *names: str) -> np.ndarray:
    for name in names:
        if name in data:
            return data[name]
    raise KeyError(f"none of {names!r} found; available columns: {sorted(data)}")


def _load_traces(root: Path) -> list[Trace]:
    traces: list[Trace] = []
    for path in sorted(root.glob("**/postProcessing/*.txt")):
        if path.name.endswith("_Ta.txt"):
            continue
        model = next((name for name in MODEL_INFO if path.name.startswith(name + "_")), None)
        if model is None:
            continue
        data = _read_table(path)
        tension_path = path.with_name(path.stem + "_Ta.txt")
        tension = _read_table(tension_path) if tension_path.exists() else {}
        # Keep the common semantic names in one table while preserving the
        # native columns for the plot labels and auditability.
        info = MODEL_INFO[model]
        data["Vm"] = _find_column(data, "Vm")
        data["cai"] = _find_column(data, "cai")
        data["ICaL"] = _find_column(data, "AV_ICaL_ICaL")
        data["IKr"] = _find_column(data, "AV_IKr_IKr")
        data["IK1"] = _find_column(data, "AV_IK1_IK1")
        data["Ito"] = _find_column(data, info["ito"], "AV_Ito_Ito", "AV_ITo_ITo")
        data["Jrel"] = _find_column(data, info["jrel"])
        data["Jup"] = _find_column(data, info["jup"])
        traces.append(Trace(model, _cl_from_name(path), data, tension))
    if not traces:
        raise FileNotFoundError(f"no Gaur/TWorld traces found below {root}")
    return traces


def _last_beat(trace: Trace) -> tuple[np.ndarray, np.ndarray]:
    time_ms = trace.data["time"] * 1000.0
    vm = trace.data["Vm"]
    candidates = np.flatnonzero((vm[1:-1] > vm[:-2]) & (vm[1:-1] >= vm[2:]) & (vm[1:-1] > -20.0)) + 1
    end_window_ms = max(250.0, 0.95 * trace.cl_ms)
    peak_i = int(np.argmax(vm))
    if candidates.size:
        # The generated protocol ends shortly after the last stimulus.  Use
        # the latest activation that has a full post-peak window, rather than
        # reporting APD90 as unavailable merely because the final beat is
        # truncated by endTime.
        peak_i = int(candidates[-1])
        for candidate in reversed(candidates):
            if time_ms[-1] - time_ms[candidate] >= end_window_ms:
                peak_i = int(candidate)
                break
    start = time_ms[peak_i] - 80.0
    end = time_ms[peak_i] + end_window_ms
    mask = (time_ms >= start) & (time_ms <= end)
    return time_ms[mask] - time_ms[peak_i], mask


def _metrics(trace: Trace) -> dict[str, object]:
    rel_time, mask = _last_beat(trace)
    vm = trace.data["Vm"][mask]
    cai = trace.data["cai"][mask]
    peak_i = int(np.argmax(vm))
    baseline_vm = float(np.median(vm[rel_time < -20.0])) if np.any(rel_time < -20.0) else float(vm[0])
    baseline_ca = float(np.median(cai[rel_time < -20.0])) if np.any(rel_time < -20.0) else float(cai[0])
    peak_vm = float(np.max(vm))
    threshold = baseline_vm + 0.10 * (peak_vm - baseline_vm)
    after_peak = np.flatnonzero(vm[peak_i:] <= threshold)
    apd90 = "N/A"
    if after_peak.size and peak_i > 0:
        repol_i = peak_i + int(after_peak[0])
        apd90 = float(rel_time[repol_i])
    peak_ca = float(np.max(cai))
    row = {
        "species": trace.species,
        "model": trace.model,
        "pacing_CL_ms": trace.cl_ms,
        "APD90_ms": apd90,
        "resting_Vm_mV": baseline_vm,
        "peak_Vm_mV": peak_vm,
        "peak_Ca_mM": peak_ca,
        "Ca_transient_amplitude_mM": peak_ca - baseline_ca,
        "peak_ICaL": float(np.min(trace.data["ICaL"][mask])),
        "peak_Jrel": float(np.max(trace.data["Jrel"][mask])),
        "peak_Jup": float(np.max(trace.data["Jup"][mask])),
        "peak_IKr": float(np.max(trace.data["IKr"][mask])),
        "peak_IK1": float(np.max(trace.data["IK1"][mask])),
        "peak_Ito": float(np.max(trace.data["Ito"][mask])),
    }
    if trace.tension:
        row["peak_Ta_kPa"] = float(np.max(_find_column(trace.tension, "AV_Ta")))
    else:
        row["peak_Ta_kPa"] = "N/A"
    return row


def _plot_waveforms(traces: list[Trace], out: Path, cl_ms: float = 1000.0) -> None:
    selected = [t for t in traces if np.isclose(t.cl_ms, cl_ms)]
    if len(selected) < 2:
        selected = sorted(traces, key=lambda t: (abs(t.cl_ms - cl_ms), t.model))[:2]
    fig, axes = plt.subplots(3, 2, figsize=(12, 10), sharex=False, constrained_layout=True)
    panels = [
        ("Vm", "Vm (mV)", "Membrane potential"),
        ("cai", "[Ca²⁺]i (mM)", "Cytosolic calcium"),
        ("ICaL", "ICaL", "L-type calcium current"),
        ("Jrel", "SR flux", "SR release and SERCA uptake"),
        ("IKr", "Current", "Repolarisation currents"),
        ("Ta", "Ta (kPa)", "Land–Niederer active tension"),
    ]
    colors = {"Gaur": "#b33c3c", "TWorld": "#2367a8"}
    for ax, (key, ylabel, title) in zip(axes.flat, panels):
        for trace in selected:
            rel_time, mask = _last_beat(trace)
            if key == "Ta":
                if not trace.tension:
                    continue
                y = _find_column(trace.tension, "AV_Ta")
                t_tension = trace.tension["time"] * 1000.0
                ionic_time_ms = trace.data["time"] * 1000.0
                peak_time_ms = ionic_time_ms[np.flatnonzero(mask)[int(np.argmax(trace.data["Vm"][mask]))]]
                tension_mask = (t_tension >= peak_time_ms - 80.0) & (t_tension <= peak_time_ms + max(250.0, 0.95 * trace.cl_ms))
                rel_time = t_tension[tension_mask] - peak_time_ms
                y = y[tension_mask]
            elif key == "Jrel":
                ax.plot(rel_time, trace.data["Jrel"][mask], color=colors[trace.model], label=f"{trace.species} / {trace.model} Jrel")
                ax.plot(rel_time, trace.data["Jup"][mask], color=colors[trace.model], linestyle="--", label=f"{trace.species} / {trace.model} Jup")
                continue
            elif key == "IKr":
                styles = {"IKr": "-", "IK1": "--", "Ito": ":"}
                for current in ("IKr", "IK1", "Ito"):
                    ax.plot(
                        rel_time,
                        trace.data[current][mask],
                        color=colors[trace.model],
                        linestyle=styles[current],
                        label=f"{trace.species} / {trace.model} {current}",
                    )
                continue
            else:
                y = trace.data[key][mask]
            ax.plot(rel_time, y, color=colors[trace.model], label=f"{trace.species} / {trace.model}")
        ax.set_title(title)
        ax.set_xlabel("Time from latest activation (ms)")
        ax.set_ylabel(ylabel)
        ax.grid(alpha=0.25)
        ax.legend(fontsize=8)
    fig.suptitle(f"Pig/Gaur versus human/TWORLD single-cell comparison (CL = {cl_ms:g} ms)")
    fig.savefig(out, dpi=220)
    plt.close(fig)


def _plot_core_waveforms(traces: list[Trace], out: Path, cl_ms: float = 1000.0) -> None:
    """Plot voltage/calcium and tension as a compact 2x1 figure."""
    selected = [t for t in traces if np.isclose(t.cl_ms, cl_ms)]
    if len(selected) < 2:
        selected = sorted(traces, key=lambda t: (abs(t.cl_ms - cl_ms), t.model))[:2]

    background = "#42516F"
    foreground = "#FFFFFF"
    fig, axes = plt.subplots(
        2,
        1,
        figsize=(10, 7.5),
        sharex=True,
        constrained_layout=True,
        facecolor=background,
    )
    calcium_ax = axes[0].twinx()
    colors = {"Gaur": "#FF8A65", "TWorld": "#90EE90"}
    calcium_styles = {"Gaur": "--", "TWorld": ":"}
    for trace in selected:
        rel_time, mask = _last_beat(trace)
        label = f"{trace.species} / {trace.model}"
        axes[0].plot(rel_time, trace.data["Vm"][mask], color=colors[trace.model], label=f"{label} Vm")
        calcium_ax.plot(
            rel_time,
            trace.data["cai"][mask],
            color=colors[trace.model],
            linestyle=calcium_styles[trace.model],
            alpha=0.30,
            label=f"{label} Ca²⁺",
        )

        if trace.tension:
            tension = _find_column(trace.tension, "AV_Ta")
            tension_time_ms = trace.tension["time"] * 1000.0
            ionic_time_ms = trace.data["time"] * 1000.0
            peak_time_ms = ionic_time_ms[np.flatnonzero(mask)[int(np.argmax(trace.data["Vm"][mask]))]]
            tension_mask = (
                (tension_time_ms >= peak_time_ms - 80.0)
                & (tension_time_ms <= peak_time_ms + max(250.0, 0.95 * trace.cl_ms))
            )
            axes[1].plot(
                tension_time_ms[tension_mask] - peak_time_ms,
                tension[tension_mask],
                color=colors[trace.model],
                linestyle="-",
                label=f"{label} Ta",
            )

    axes[0].set_ylabel("Vm (mV)")
    calcium_ax.set_ylabel("[Ca²⁺]i (mM)")
    axes[1].set_ylabel("Ta (kPa)")
    axes[1].set_xlabel("Time from latest activation (ms)")
    axes[0].set_xlim(right=700.0)
    for axis in (axes[0], calcium_ax, axes[1]):
        axis.set_facecolor(background)
        axis.tick_params(axis="both", colors=foreground, labelcolor=foreground)
        axis.xaxis.label.set_color(foreground)
        axis.yaxis.label.set_color(foreground)
        for spine in axis.spines.values():
            spine.set_color(foreground)
        axis.yaxis.set_major_locator(LinearLocator(2))
    calcium_formatter = ScalarFormatter(useMathText=True)
    calcium_formatter.set_scientific(True)
    calcium_formatter.set_powerlimits((0, 0))
    calcium_ax.yaxis.set_major_formatter(calcium_formatter)
    for ax in axes:
        ax.grid(color=foreground, alpha=0.16)
    voltage_handles, voltage_labels = axes[0].get_legend_handles_labels()
    calcium_handles, calcium_labels = calcium_ax.get_legend_handles_labels()
    axes[0].legend(
        voltage_handles + calcium_handles,
        voltage_labels + calcium_labels,
        fontsize=9,
        facecolor=background,
        edgecolor=foreground,
        labelcolor=foreground,
    )
    axes[1].legend(fontsize=9, facecolor=background, edgecolor=foreground, labelcolor=foreground)
    fig.savefig(out, dpi=220, facecolor=background, edgecolor="none")
    plt.close(fig)


def _plot_rate(rows: list[dict[str, object]], out: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), constrained_layout=True)
    colors = {"Pig": "#b33c3c", "Human": "#2367a8"}
    for species in ("Pig", "Human"):
        values = [r for r in rows if r["species"] == species and isinstance(r["APD90_ms"], float)]
        values.sort(key=lambda r: float(r["pacing_CL_ms"]))
        x = [float(r["pacing_CL_ms"]) for r in values]
        axes[0].plot(x, [float(r["APD90_ms"]) for r in values], "o-", color=colors[species], label=species)
        axes[1].plot(x, [float(r["peak_Ca_mM"]) for r in values], "o-", color=colors[species], label=species)
    for ax, ylabel, title in zip(axes, ("APD90 (ms)", "Peak [Ca²⁺]i (mM)"), ("Rate dependence of repolarisation", "Rate dependence of calcium")):
        ax.set_title(title)
        ax.set_xlabel("Pacing cycle length (ms)")
        ax.set_ylabel(ylabel)
        ax.grid(alpha=0.25)
        ax.legend()
    fig.suptitle("Pig/Gaur versus human/TWORLD rate dependence")
    fig.savefig(out, dpi=220)
    plt.close(fig)


def run_postprocessing(
    *,
    output_dir: str,
    setup_root: str | None = None,
    input_dir: str | None = None,
    **_: object,
) -> list[dict]:
    """Generate comparison figures and a metrics CSV from sweep case outputs."""
    del setup_root
    root = Path(input_dir) if input_dir is not None else Path(output_dir)
    destination = Path(output_dir)
    destination.mkdir(parents=True, exist_ok=True)
    traces = _load_traces(root)
    rows = [_metrics(trace) for trace in traces]
    columns = list(rows[0])
    with (destination / "species_comparison_metrics.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    _plot_core_waveforms(traces, destination / "species_comparison_waveforms.png")
    _plot_waveforms(traces, destination / "species_comparison_all_variables.png")
    _plot_rate(rows, destination / "species_comparison_rate_dependence.png")
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, required=True, help="run directory's cases directory")
    parser.add_argument("--output-dir", type=Path, help="run directory for figures and metrics (default: input directory)")
    args = parser.parse_args()
    destination = args.output_dir or args.input_dir
    rows = run_postprocessing(output_dir=str(destination), input_dir=str(args.input_dir))
    print(f"Wrote comparison outputs for {len(rows)} traces to {destination}")


if __name__ == "__main__":
    main()
