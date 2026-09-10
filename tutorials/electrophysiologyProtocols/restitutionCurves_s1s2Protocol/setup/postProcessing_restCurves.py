from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
import re
import csv


# ==============================
# PARAMETERS (tune here)
# ==============================

DVDT_THRESHOLD = 5000.0     # mV/s, upstroke detection
BASELINE       = 0.002       # baseline window before upstroke
PEAK_SEARCH    = 2          # window to find AP peak


# ==============================
# BASIC UTILITIES
# ==============================


def extract_S2_value(path):
    """
    Extract numeric S2 value from filename.
    Example matches: S2_300, S2-300, S2_300ms
    """
    match = re.search(r"S2[_-]?(\d+)", path.name)
    if match:
        return int(match.group(1))
    else:
        return np.inf   # puts malformed files at the end


def load_trace(filename):
    data = np.genfromtxt(
        filename,
        names=True,       # reads header
        delimiter=None    # auto-detect whitespace
    )

    time = data["time"]
    vm   = data["Vm"]

    return time, vm


def safe_dvdt(vm, time):
    dvdt = np.zeros_like(vm)
    dvdt[1:-1] = (vm[2:] - vm[:-2]) / (time[2:] - time[:-2])
    dvdt[0]  = dvdt[1]
    dvdt[-1] = dvdt[-2]
    return dvdt

def detect_beats(time, vm):

    dvdt = safe_dvdt(vm, time)

    baseline_s = BASELINE
    peak_s     = PEAK_SEARCH

    candidates = []
    for i in range(1, len(dvdt)):
        if dvdt[i-1] < DVDT_THRESHOLD <= dvdt[i]:
            frac = (DVDT_THRESHOLD - dvdt[i-1]) / (dvdt[i] - dvdt[i-1])
            if 0.0 <= frac <= 1.0:
                candidates.append(time[i-1] + frac * (time[i] - time[i-1]))

    beats = []

    for k, t_up in enumerate(candidates):

        # ---- baseline
        mask_base = (time >= t_up - baseline_s) & (time < t_up)
        if not np.any(mask_base):
            continue
        v_base = np.min(vm[mask_base])

        # ---- peak
        mask_peak = (time >= t_up) & (time <= t_up + peak_s)
        if not np.any(mask_peak):
            continue

        seg_t = time[mask_peak]
        seg_v = vm[mask_peak]
        v_peak = np.max(seg_v)

        # ---- APD90
        v90 = v_base + 0.1 * (v_peak - v_base)

        j_peak = np.argmax(seg_v)
        t_rep = np.nan

        for j in range(j_peak + 1, len(seg_v)):
            if seg_v[j-1] > v90 >= seg_v[j]:
                frac = (v90 - seg_v[j-1]) / (seg_v[j] - seg_v[j-1])
                t_rep = seg_t[j-1] + frac * (seg_t[j] - seg_t[j-1])
                break

        repol_found = np.isfinite(t_rep)

        beats.append({
            "t_up": float(t_up),
            "v_base": float(v_base),
            "v_peak": float(v_peak),
            "v90": float(v90),
            "t_repol90": float(t_rep),
            "repol_found": repol_found
        })
    return beats


def detect_beats(time, vm):

    dvdt = safe_dvdt(vm, time)

    baseline_s = BASELINE
    peak_s     = PEAK_SEARCH

    candidates = []
    last_t = None
    for i in range(1, len(dvdt)):
        if dvdt[i-1] < DVDT_THRESHOLD <= dvdt[i]:
            frac = (DVDT_THRESHOLD - dvdt[i-1]) / (dvdt[i] - dvdt[i-1])
            if 0.0 <= frac <= 1.0:
                t_cand = time[i-1] + frac * (time[i] - time[i-1])
                if last_t is None or (t_cand - last_t) > 0.05:
                    candidates.append(t_cand)
                    last_t = t_cand

    # ---- DEBUG: list all dv/dt candidates
    print("\n=== DV/DT CANDIDATES ===")
    for k, t in enumerate(candidates):
        print(f"[{k}] t_up = {t:.6f}")
    print("========================\n")

    beats = []

    for k, t_up in enumerate(candidates):

        print("\n" + "-" * 70)
        print(f"BEAT {k}")
        print(f"t_up = {t_up:.6f}")
        if k < len(candidates) - 1:
            print(f"next t_up = {candidates[k+1]:.6f}")
        else:
            print("next t_up = NONE")

        # ---- baseline
        mask_base = (time >= t_up - baseline_s) & (time < t_up)
        if not np.any(mask_base):
            print("❌ No baseline samples")
            continue
        v_base = np.min(vm[mask_base])

        # ---- peak / repolarization window
        t_end = t_up + peak_s
        if k < len(candidates) - 1:
            t_end = min(t_end, candidates[k+1])

        mask_peak = (time >= t_up) & (time <= t_end)
        if not np.any(mask_peak):
            print("❌ No peak samples")
            continue

        seg_t = time[mask_peak]
        seg_v = vm[mask_peak]
        v_peak = np.max(seg_v)

        # ---- Reject fake beats (stimulus artifact only)
        # We ensure the cell actually captured by checking the voltage 3 ms after t_up (i.e. after a 2ms stimulus + 1ms grace)
        t_check = t_up + 0.003
        mask_check = (seg_t >= t_check)
        if np.any(mask_check):
            v_check = seg_v[mask_check][0]
            if v_check < -10.0:
                print(f"❌ Voltage collapsed after stimulus ({v_check:.3f} mV at +3ms). Likely just an artifact.")
                continue
        else:
            if v_peak < -10.0:
                print(f"❌ Peak voltage too low ({v_peak:.3f} mV). Likely just a stimulus artifact.")
                continue

        # ---- APD90, APD70, APD50
        v90 = v_base + 0.1 * (v_peak - v_base)
        v70 = v_base + 0.3 * (v_peak - v_base)
        v50 = v_base + 0.5 * (v_peak - v_base)
        j_peak = np.argmax(seg_v)

        print(f"v_base = {v_base:.3f}, v_peak = {v_peak:.3f}, v90 = {v90:.3f}, v70 = {v70:.3f}, v50 = {v50:.3f}")
        print(f"j_peak index = {j_peak}, t_peak = {seg_t[j_peak]:.6f}")

        t_rep90 = np.nan
        t_rep70 = np.nan
        t_rep50 = np.nan

        for j in range(j_peak + 1, len(seg_v)):
            if np.isnan(t_rep50) and seg_v[j-1] > v50 >= seg_v[j]:
                frac = (v50 - seg_v[j-1]) / (seg_v[j] - seg_v[j-1])
                t_rep50 = seg_t[j-1] + frac * (seg_t[j] - seg_t[j-1])

            if np.isnan(t_rep70) and seg_v[j-1] > v70 >= seg_v[j]:
                frac = (v70 - seg_v[j-1]) / (seg_v[j] - seg_v[j-1])
                t_rep70 = seg_t[j-1] + frac * (seg_t[j] - seg_t[j-1])

            if np.isnan(t_rep90) and seg_v[j-1] > v90 >= seg_v[j]:
                frac = (v90 - seg_v[j-1]) / (seg_v[j] - seg_v[j-1])
                t_rep90 = seg_t[j-1] + frac * (seg_t[j] - seg_t[j-1])
                print(f"✅ REPOL FOUND at t = {t_rep90:.6f}")
                break

            if k < len(candidates) - 1 and seg_t[j] >= candidates[k+1]:
                print("⚠️ PASSED NEXT t_up WITHOUT REPOL")
                break

        repol_found = np.isfinite(t_rep90)

        print(f"RESULT: repol_found = {repol_found}")

        beats.append({
            "t_up": float(t_up),
            "v_base": float(v_base),
            "v_peak": float(v_peak),
            "v90": float(v90),
            "v70": float(v70),
            "v50": float(v50),
            "t_repol90": float(t_rep90),
            "t_repol70": float(t_rep70),
            "t_repol50": float(t_rep50),
            "repol_found": repol_found
        })

    return beats







# APD / DI
# ==============================
def get_s1_s2_beats(beats, filepath=None, config=None):
    if len(beats) < 2:
        return None, None

    if filepath is None:
        return beats[-2], beats[-1]

    name = filepath.name
    m_s1 = re.search(r"S1_(\d+)", name)

    if config:
        s1_val = config.get("s1_interval_ms", 1000)
        n_s1 = config.get("n_s1", 10)
    else:
        s1_val = int(m_s1.group(1)) if m_s1 else 1000
        n_s1 = 10

    m_s2 = re.search(r"S2_(\d+)", name)
    s2_val = int(m_s2.group(1)) if m_s2 else 250

    # The S1 train ends and the first S2 happens based on the config intervals
    s1_target_time = (n_s1 * s1_val) / 1000.0

    s1_beat = None
    s2_beat = None

    for b in beats:
        if abs(b["t_up"] - s1_target_time) < 0.1:
            s1_beat = b
            break

    if s1_beat is None:
        return None, None

    idx = beats.index(s1_beat)
    if idx + 1 < len(beats):
        potential_s2 = beats[idx + 1]

        # We must ensure this is actually the S2 beat.
        # If the true S2 beat failed, this might accidentally be the S3 beat!
        actual_interval_ms = (potential_s2["t_up"] - s1_beat["t_up"]) * 1000.0

        if abs(actual_interval_ms - s2_val) < 2.0:
            s2_beat = potential_s2
        else:
            # The next successful beat was NOT the S2 beat, meaning S2 failed!
            return None, None
    else:
        return None, None

    return s1_beat, s2_beat

def compute_apd_di(beats, filepath=None, config=None):
    s1, s2 = get_s1_s2_beats(beats, filepath, config)
    if s1 is None or s2 is None:
        return None

    res = {}

    # We don't strictly require valid APD90 for APD70/50 to be valid!
    # But we calculate what we can.

    for level in [90, 70, 50]:
        t_rep_s1 = s1.get(f"t_repol{level}", np.nan)
        t_rep_s2 = s2.get(f"t_repol{level}", np.nan)

        if np.isfinite(t_rep_s1) and np.isfinite(t_rep_s2) and t_rep_s1 < s2["t_up"]:
            di = s2["t_up"] - t_rep_s1
            apd = t_rep_s2 - s2["t_up"]
            if di > 0 and apd > 0:
                res[f"DI{level}"] = di
                res[f"APD{level}"] = apd

    if not res:
        return None

    return res



def plot_trace(time, vm, beats, filepath=None, savepath=None, config=None):
    plt.figure(figsize=(11, 4))
    plt.plot(time, vm, color="black", lw=1.2, label="Vm")

    # ---- mark all upstrokes and repolarizations
    for b in beats:
        plt.axvline(
            b["t_up"],
            color="navy",
            lw=1,
            alpha=0.6

            )
        if np.isfinite(b.get("t_repol90", np.nan)):
            plt.axvline(
                b["t_repol90"],
                color="navy",
                ls=":",
                lw=1,
                alpha=0.6
            )
    # ---- annotate last S1–S2 pair
    s1, s2 = get_s1_s2_beats(beats, filepath, config)
    if s1 is not None and s2 is not None:

        # Find best available repolarization level
        best_level = None
        for level in [90, 70, 50]:
            t_rep = s1.get(f"t_repol{level}", np.nan)
            if np.isfinite(t_rep) and t_rep < s2["t_up"]:
                best_level = level
                break

        if best_level is not None:
            t_rep_s1 = s1[f"t_repol{best_level}"]
            t_rep_s2 = s2[f"t_repol{best_level}"]

            # ---- DI
            plt.axvspan(
                t_rep_s1,
                s2["t_up"],
                color="tab:blue",
                alpha=0.5,
                label=f"DI{best_level}"
            )

            # ---- APD
            if np.isfinite(t_rep_s2):
                plt.axvspan(
                    s2["t_up"],
                    t_rep_s2,
                    color="tab:red",
                    alpha=0.5,
                    label=f"APD{best_level}"
                )

            di_ms  = (s2["t_up"] - t_rep_s1) * 1e3
            apd_ms = (t_rep_s2 - s2["t_up"]) * 1e3

            text = f"DI{best_level} = {di_ms:.1f} ms\nAPD{best_level} = {apd_ms:.1f} ms"

        else:
            text = "EXTREME EARLY DEPOLARIZATION"

            plt.axvline(
                s2["t_up"],
                color="orange",
                lw=2.5,
                label="No repol > 50%"
            )


        # ---- unified text box (always shown)
        plt.text(
            0.1, 0.8,
            text,
            transform=plt.gca().transAxes,
            fontsize=11,
            va="top",
            bbox=dict(boxstyle="round", fc="white", alpha=0.9)
        )


    plt.xlabel("Time (s)")
    plt.ylabel("Vm (mV)")
    plt.title("Action potential with S1–S2 DI and APD90")
    plt.legend(frameon=False)
    plt.tight_layout()

    if savepath is None:
        plt.show()
    else:
        plt.savefig(savepath, dpi=300)
        plt.close()


def postprocess_one_ionic_model(
    base_dir: Path,
    output_folder: str,
    ionic_model: str,
    tissues: list[str],
    show_plot: bool = False,
    config: dict = None
):


    output_dir = base_dir / output_folder / ionic_model
    print(f"DEBUG: base_dir={base_dir}, output_folder={output_folder}, output_dir={output_dir}")
    if not output_dir.exists():
        print(f"❌ Output folder not found: {output_dir}")
        return
    print(f"\n📊 Processing ionic model: {ionic_model}")

    plt.figure()
    data_rows = []
    input_dir = output_dir

    all_restitution_data = []

    for tissue in tissues:
        # Find all files for this ionic model + tissue
        file_pattern = f"*{ionic_model}*{tissue}*.txt"
        files = list(input_dir.rglob(file_pattern))

        if not files:
            print(f"No files found for {ionic_model} / {tissue}")
            continue

        for f in files:
            time, vm = load_trace(f)
            if len(time) < 2:
                continue

            beats = detect_beats(time, vm)
            # Attach filename to beats once
            for b in beats:
                b["file"] = str(f)

            tissue_dir = output_dir / tissue
            tissue_dir.mkdir(parents=True, exist_ok=True)
            plot_path = tissue_dir / f"{f.stem}_annotated.png"
            plot_trace(time, vm, beats, filepath=f, savepath=plot_path, config=config)

            res = compute_apd_di(beats, f, config=config)
            if res is None:
                continue

            row = {"tissue": tissue}

            if "DI90" in res:
                row["DI90_ms"] = res["DI90"] * 1e3
                row["APD90_ms"] = res["APD90"] * 1e3
            if "DI70" in res:
                row["DI70_ms"] = res["DI70"] * 1e3
                row["APD70_ms"] = res["APD70"] * 1e3
            if "DI50" in res:
                row["DI50_ms"] = res["DI50"] * 1e3
                row["APD50_ms"] = res["APD50"] * 1e3

            all_restitution_data.append(row)

    if not all_restitution_data:
        print(f"⚠️ No valid restitution points for {ionic_model}")
        return

    df = pd.DataFrame(all_restitution_data)

    # Save single CSV for the ionic model
    csv_path = output_dir / f"{ionic_model}_restitution.csv"
    if csv_path.exists():
        csv_path.unlink()
    df.to_csv(csv_path, index=False)

    # Plot combined restitution curves
    plt.figure(figsize=(10, 6))

    for tissue in tissues:
        tissue_df = df[df["tissue"] == tissue]
        if tissue_df.empty:
            continue

        if "DI90_ms" in tissue_df.columns:
            valid90 = tissue_df.dropna(subset=["DI90_ms", "APD90_ms"]).sort_values("DI90_ms")
            if not valid90.empty:
                plt.plot(valid90["DI90_ms"], valid90["APD90_ms"], marker='o', label=f'{tissue} APD90')

        if "DI70_ms" in tissue_df.columns:
            valid70 = tissue_df.dropna(subset=["DI70_ms", "APD70_ms"]).sort_values("DI70_ms")
            if not valid70.empty:
                plt.plot(valid70["DI70_ms"], valid70["APD70_ms"], marker='s', label=f'{tissue} APD70')

        if "DI50_ms" in tissue_df.columns:
            valid50 = tissue_df.dropna(subset=["DI50_ms", "APD50_ms"]).sort_values("DI50_ms")
            if not valid50.empty:
                plt.plot(valid50["DI50_ms"], valid50["APD50_ms"], marker='^', label=f'{tissue} APD50')

    plt.xlabel("Diastolic Interval (ms)")
    plt.ylabel("Action Potential Duration (ms)")
    plt.title(f"Restitution Curve - {ionic_model}")
    plt.grid(True, ls="--", alpha=0.6)
    plt.legend()
    plt.tight_layout()

    plot_path = output_dir / f"{ionic_model}_restitution.png"
    plt.savefig(plot_path, dpi=300)
    if show_plot:
        plt.show()
    plt.close()

    # Plot APD90 only restitution curves
    plt.figure(figsize=(10, 6))

    for tissue in tissues:
        tissue_df = df[df["tissue"] == tissue]
        if tissue_df.empty:
            continue

        if "DI90_ms" in tissue_df.columns:
            valid90 = tissue_df.dropna(subset=["DI90_ms", "APD90_ms"]).sort_values("DI90_ms")
            if not valid90.empty:
                plt.plot(valid90["DI90_ms"], valid90["APD90_ms"], marker='o', label=f'{tissue} APD90')

    plt.xlabel("Diastolic Interval (ms)")
    plt.ylabel("Action Potential Duration (ms)")
    plt.title(f"Restitution Curve (APD90) - {ionic_model}")
    plt.grid(True, ls="--", alpha=0.6)
    plt.legend()
    plt.tight_layout()

    plot_path_apd90 = output_dir / f"{ionic_model}_restitution_APD90.png"
    plt.savefig(plot_path_apd90, dpi=300)
    if show_plot:
        plt.show()
    plt.close()

    print(f"✅ Saved plots and data for {ionic_model}")



def restitution_curves(
    base_dir: Path,
    output_folder: str,
    ionic_model: str,
    tissue_types: list[str],
    show_plots: bool = False,
):

        postprocess_one_ionic_model(
            base_dir=base_dir,
            output_folder=output_folder,
            ionic_model= ionic_model,
            tissues=tissue_types,
            show_plot = show_plots,
        )



def run_postprocessing(
    *,
    output_dir: str,
    setup_root: str | None = None,
    **kwargs,
) -> list:
    """run_postprocessing entry point matching openfoam_driver.postprocessing's
    PostprocessingProtocol shape (output_dir, setup_root, **kwargs) -> list[dict].

    Not currently invoked automatically -- driverFOAM has no post-DAG hook
    calling this (openfoam_driver.postprocessing.driver, which used to wire
    tutorial postprocessing functions into the run engine, was removed
    2026-08-18 after being found unreachable). Run manually against a
    completed sweep's output_dir until a replacement hand-off exists.
    Expected kwargs:
        ionic_models  (list[str])        - Models to post-process.
        tissue_map    (dict[str, list])  - Tissue types per ionic model.
        show_plots    (bool)             - Whether to display plots interactively.
    """
    ionic_models: list[str] = kwargs.get("ionic_models", [])
    tissue_map: dict[str, list[str]] = kwargs.get("tissue_map", {})
    show_plots: bool = kwargs.get("show_plots", False)

    output_dir_path = Path(output_dir)
    artifacts: list[dict] = []

    for model in ionic_models:
        tissues = list(tissue_map.get(model, []))
        postprocess_one_ionic_model(
            base_dir=output_dir_path.parent,
            output_folder=output_dir_path.name,
            ionic_model=model,
            tissues=tissues,
            show_plot=show_plots,
            config=kwargs
        )
        model_dir = output_dir_path / model
        fig_path = model_dir / f"{model}_restitution.png"
        fig_path_apd90 = model_dir / f"{model}_restitution_APD90.png"
        csv_path = model_dir / f"{model}_restitution.csv"

        if fig_path.exists():
            artifacts.append({
                "path": f"{model}/{fig_path.name}",
                "label": f"{model} restitution curve (All)",
                "kind": "plot",
                "format": "png",
            })
        if fig_path_apd90.exists():
            artifacts.append({
                "path": f"{model}/{fig_path_apd90.name}",
                "label": f"{model} restitution curve (APD90 only)",
                "kind": "plot",
                "format": "png",
            })
        if csv_path.exists():
            artifacts.append({
                "path": f"{model}/{csv_path.name}",
                "label": f"{model} restitution data",
                "kind": "data",
                "format": "csv",
            })

    return artifacts


if __name__ == "__main__":
    raise SystemExit(
        "This module must be imported and called from mainS1-S2protocol.py"
    )