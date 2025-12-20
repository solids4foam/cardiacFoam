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
    for i in range(1, len(dvdt)):
        if dvdt[i-1] < DVDT_THRESHOLD <= dvdt[i]:
            frac = (DVDT_THRESHOLD - dvdt[i-1]) / (dvdt[i] - dvdt[i-1])
            if 0.0 <= frac <= 1.0:
                candidates.append(time[i-1] + frac * (time[i] - time[i-1]))

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

        # ---- APD90
        v90 = v_base + 0.1 * (v_peak - v_base)
        j_peak = np.argmax(seg_v)

        print(f"v_base = {v_base:.3f}, v_peak = {v_peak:.3f}, v90 = {v90:.3f}")
        print(f"j_peak index = {j_peak}, t_peak = {seg_t[j_peak]:.6f}")

        t_rep = np.nan

        for j in range(j_peak + 1, len(seg_v)):
            if seg_v[j-1] > v90 >= seg_v[j]:
                frac = (v90 - seg_v[j-1]) / (seg_v[j] - seg_v[j-1])
                t_rep = seg_t[j-1] + frac * (seg_t[j] - seg_t[j-1])
                print(f"✅ REPOL FOUND at t = {t_rep:.6f}")
                break

            if k < len(candidates) - 1 and seg_t[j] >= candidates[k+1]:
                print("⚠️ PASSED NEXT t_up WITHOUT REPOL")
                break

        repol_found = np.isfinite(t_rep)

        print(f"RESULT: repol_found = {repol_found}")

        beats.append({
            "t_up": float(t_up),
            "v_base": float(v_base),
            "v_peak": float(v_peak),
            "v90": float(v90),
            "t_repol90": float(t_rep),
            "repol_found": repol_found
        })

    return beats

def classify_early_strokes(beats):
    """
    A beat is EARLY if repolarization is interrupted
    by a subsequent upstroke.
    """

    early_beats = []
    n = len(beats)
    if n == 0:
        return early_beats

    for b in beats:
        b["early"] = False
        b["valid_apd"] = b.get("repol_found", False)
        b["valid_di"] = False

    for i in range(n - 1):
        curr = beats[i]
        nextb = beats[i + 1]

        early = False
        reason = ""

        if curr["repol_found"]:
            if nextb["t_up"] < curr["t_repol90"]:
                early = True
                reason = "next upstroke before APD90"
        else:
            early = True
            reason = "repolarization not completed before next upstroke"

        if early:
            curr["early"] = True
            curr["valid_apd"] = False
            curr["valid_di"] = False
            early_beats.append(curr)

            print("\n" + "=" * 80)
            print(f"⚠️ EARLY STROKE DETECTED (beat {i})")
            print(f"Reason: {reason}")
            print("-" * 80)

            print(f"t_up (curr)      = {curr['t_up']:.6f}")
            print(f"t_repol90        = {curr['t_repol90']}")
            print(f"t_up (next)      = {nextb['t_up']:.6f}")

            if curr["repol_found"]:
                print(f"APD90            = {curr['t_repol90'] - curr['t_up']:.6f}")
                print(f"DI (raw)         = {nextb['t_up'] - curr['t_repol90']:.6f}")
            else:
                print("APD90            = NOT FOUND (INTERRUPTED)")
                print(f"DI (raw)         = {nextb['t_up'] - curr['t_up']:.6f}")

            print("=" * 80)

    # DI belongs to the second beat ONLY if first beat was normal
    for i in range(1, n):
        if not beats[i - 1]["early"] and beats[i].get("repol_found", False):
            beats[i]["valid_di"] = True

    return early_beats




# APD / DI
# ==============================
def compute_apd_di(beats):

    if len(beats) < 2:
        return None, None

    s1 = beats[-2]
    s2 = beats[-1]

    # ---- EARLY STROKE: reject if S1 was interrupted
    if s1.get("early", False):
        return None, None

    if not (
        s1.get("valid_apd", False) and
        s2.get("valid_apd", False) and
        s2.get("valid_di",  False)
    ):
        return None, None

    di  = s2["t_up"]      - s1["t_repol90"]
    apd = s2["t_repol90"] - s2["t_up"]

    if di <= 0 or apd <= 0:
        return None, None

    return apd, di



def plot_trace(time, vm, beats, savepath=None):
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
    if len(beats) >= 2:
        s1 = beats[-2]
        s2 = beats[-1]

        valid_pair = (
            not s1.get("early", False) and
            s1.get("repol_found", False) and
            s1["t_repol90"] < s2["t_up"]
        )

        if valid_pair:
            # ---- DI
            plt.axvspan(
                s1["t_repol90"],
                s2["t_up"],
                color="tab:blue",
                alpha=0.5,
                label="DI"
            )

            # ---- APD90 (S2)
            plt.axvspan(
                s2["t_up"],
                s2["t_repol90"],
                color="tab:red",
                alpha=0.5,
                label="APD90"
            )

            di_ms  = (s2["t_up"] - s1["t_repol90"]) * 1e3
            apd_ms = (s2["t_repol90"] - s2["t_up"]) * 1e3

            text = f"DI = {di_ms:.1f} ms\nAPD90 = {apd_ms:.1f} ms"

        else:
            text = "EARLY DEPOLARIZATION"

            plt.axvline(
            s2["t_up"],
            color="orange",
            lw=2.5,
            label="Early depolarization"
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
    show_plot: bool = False
):
    
    
    output_dir = base_dir / output_folder
    if not output_dir.exists():
        print(f"❌ Output folder not found: {output_dir}")
        return
    print(f"\n📊 Processing ionic model: {ionic_model}")

    plt.figure()
    data_rows = []
    early_files = set()  

    for tissue in tissues:
        pattern = f"{ionic_model}_{tissue}_*.txt"
        files = sorted(
                output_dir.glob(pattern),
                key=extract_S2_value
                      )
        if not files:
            print(f"⚠️ No files for {ionic_model} / {tissue}")
            continue

        di_ms = []
        apd_ms = []

        for f in files:
            time, vm = load_trace(f)
            
            beats = detect_beats(time, vm)
            # Attach filename to beats once
            for b in beats:
                b["file"] = str(f)

            early_beats = classify_early_strokes(beats)
            for b in early_beats:
                early_files.add(b["file"])

            if show_plot:
                plot_trace(time, vm, beats) 
            
            apd, di = compute_apd_di(beats)
            if apd is None or di is None:
                continue

            di_ms.append(di * 1e3)
            apd_ms.append(apd * 1e3)

        if not di_ms:
            print(f"⚠️ No valid restitution points for {ionic_model} / {tissue}")
            continue

        for di, apd in zip(di_ms, apd_ms):
            data_rows.append({
                "tissue": tissue,
                "DI_ms": di,
                "APD90_ms": apd,
                            })
        
        plt.plot(di_ms, apd_ms, marker="o", label=tissue)

    plt.xlabel("DI (ms)")
    plt.ylabel("APD90 (ms)")
    plt.title(f"S1–S2 Restitution – {ionic_model}")
    plt.legend()
    plt.grid(True)

    fig_path = output_dir / f"{ionic_model}_restitution.png"
    plt.savefig(fig_path, dpi=300)
    if show_plot:
        plt.show() 
    plt.close()
    csv_path = output_dir / f"{ionic_model}_restitution.csv"

    with csv_path.open("w", newline="") as csvfile:
        writer = csv.DictWriter(
            csvfile,
            fieldnames=["tissue", "DI_ms", "APD90_ms"]
        )
        writer.writeheader()
        writer.writerows(data_rows)


    if early_files:
        print(
        f"\n⚠️ Early strokes detected for {ionic_model}: "
        f"{len(early_files)} file(s)")

        early_list_path = output_dir / f"{ionic_model}_early_upstroke_files.txt"
        with early_list_path.open("w") as f:
            for fname in sorted(early_files):
                f.write(fname + "\n")

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



if __name__ == "__main__":
    raise SystemExit(
        "This module must be imported and called from mainS1-S2protocol.py"
    )