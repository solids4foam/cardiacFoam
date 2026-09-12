#!/usr/bin/env python3
"""
plot_ecg_video.py
=================
Renders the pseudoECG.dat traces as a frame-by-frame animation and
assembles it into an MP4 with ffmpeg.

Each frame shows all 6 leads (V1–V6) with the ECG trace built up to
time t — matching the same frame count as a 3D simulation video so
both can be played side-by-side.

Usage
-----
    python3 plot_ecg_video.py --case-dir PATHOS/LBBB
    python3 plot_ecg_video.py --case-dir PATHOS/LBBB --fps 30 --output lbbb_ecg.mp4
    python3 plot_ecg_video.py --case-dir PATHOS/LBBB --fps 30 --tmax 0.5

Requirements
------------
    pip install matplotlib numpy
    ffmpeg must be on PATH
"""

import argparse
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec


LEADS = ["V1", "V2", "V3", "V4", "V5", "V6"]
LAYOUT = [(0, 0), (0, 1), (1, 0), (1, 1), (2, 0), (2, 1)]   # row, col for each lead


def load_ecg(path: Path):
    data = np.loadtxt(path, comments="#")
    time = data[:, 0]
    signals = {lead: data[:, i + 1] for i, lead in enumerate(LEADS)}
    return time, signals


def render_frame(fig, axes, time, signals, frame_idx, ylims):
    t_now = time[frame_idx]
    t_arr = time[:frame_idx + 1]

    for ax, lead in zip(axes, LEADS):
        ax.cla()
        y = signals[lead][:frame_idx + 1]
        ax.plot(t_arr * 1000, y, lw=1.2, color="#1f77b4")   # time in ms
        ax.set_xlim(time[0] * 1000, time[-1] * 1000)
        ax.set_ylim(ylims[lead])
        ax.axvline(t_now * 1000, color="red", lw=0.8, alpha=0.6)
        ax.set_title(lead, fontsize=9, fontweight="bold", pad=2)
        ax.set_ylabel("pseudoECG", fontsize=7)
        ax.tick_params(labelsize=7)
        ax.grid(True, alpha=0.2)

    axes[-1].set_xlabel("time (ms)", fontsize=8)
    axes[-2].set_xlabel("time (ms)", fontsize=8)
    fig.suptitle(f"Pseudo-ECG   t = {t_now * 1000:.1f} ms", fontsize=10)


def main():
    parser = argparse.ArgumentParser(description="Render pseudoECG animation")
    parser.add_argument("--case-dir", default=".", help="Case directory")
    parser.add_argument("--fps", type=int, default=30, help="Frames per second (default 30)")
    parser.add_argument("--output", default=None, help="Output MP4 path")
    parser.add_argument("--tmin", type=float, default=None)
    parser.add_argument("--tmax", type=float, default=None)
    parser.add_argument("--dpi", type=int, default=120)
    args = parser.parse_args()

    case_dir = Path(args.case_dir).resolve()
    ecg_path = case_dir / "postProcessing" / "pseudoECG.dat"
    if not ecg_path.exists():
        sys.exit(f"Cannot find {ecg_path}")

    out_path = Path(args.output) if args.output else case_dir / "postProcessing" / "ecg_video.mp4"

    print(f"Loading {ecg_path} ...")
    time, signals = load_ecg(ecg_path)

    # Filter time range
    mask = np.ones(len(time), dtype=bool)
    if args.tmin is not None:
        mask &= time >= args.tmin
    if args.tmax is not None:
        mask &= time <= args.tmax
    time = time[mask]
    signals = {k: v[mask] for k, v in signals.items()}
    n_frames = len(time)

    print(f"  {n_frames} frames  →  {n_frames / args.fps:.1f} s at {args.fps} fps")

    # Pre-compute y limits with 5% padding
    ylims = {}
    for lead in LEADS:
        ymin, ymax = signals[lead].min(), signals[lead].max()
        pad = max(abs(ymax - ymin) * 0.05, 1e-9)
        ylims[lead] = (ymin - pad, ymax + pad)

    # Set up figure  — 3 rows × 2 cols
    fig = plt.figure(figsize=(10, 7))
    gs = gridspec.GridSpec(3, 2, figure=fig, hspace=0.55, wspace=0.35)
    axes = [fig.add_subplot(gs[r, c]) for r, c in LAYOUT]

    tmp_dir = tempfile.mkdtemp(prefix="ecg_frames_")
    print(f"Rendering frames to {tmp_dir} ...")

    try:
        for i in range(n_frames):
            render_frame(fig, axes, time, signals, i, ylims)
            frame_path = os.path.join(tmp_dir, f"frame_{i:05d}.png")
            fig.savefig(frame_path, dpi=args.dpi, bbox_inches="tight")
            if i % 50 == 0 or i == n_frames - 1:
                print(f"  frame {i+1}/{n_frames}", end="\r")

        print(f"\nAssembling video → {out_path} ...")
        out_path.parent.mkdir(parents=True, exist_ok=True)
        cmd = [
            "ffmpeg", "-y",
            "-framerate", str(args.fps),
            "-i", os.path.join(tmp_dir, "frame_%05d.png"),
            "-c:v", "libx264",
            "-pix_fmt", "yuv420p",
            "-crf", "20",
            str(out_path),
        ]
        subprocess.run(cmd, check=True)
        print(f"Done: {out_path}")

    finally:
        shutil.rmtree(tmp_dir, ignore_errors=True)
        plt.close(fig)


if __name__ == "__main__":
    main()
