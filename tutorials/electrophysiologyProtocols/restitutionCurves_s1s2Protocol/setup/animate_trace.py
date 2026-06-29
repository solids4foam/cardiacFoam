"""
Animate a single-cell voltage trace as an MP4 video.

Each frame draws the Vm trace progressively from left to right, giving a
"live recording" effect of the action potential waveform.

Dependencies: matplotlib, numpy, ffmpeg (for MP4) or Pillow (for GIF fallback).

Standalone usage:
    python animate_trace.py <trace.txt> <output.mp4>
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")  # headless — no display required
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np


# ── tuneable defaults ────────────────────────────────────────────────────────

N_FRAMES   = 300    # number of animation frames
FPS        = 30     # frames per second in the output video
LINE_COLOR = "steelblue"
LINE_WIDTH = 1.5


# ── core function ────────────────────────────────────────────────────────────

def create_voltage_animation(
    data_file: Path | str,
    output_path: Path | str,
    *,
    title: str = "",
    n_frames: int = N_FRAMES,
    fps: int = FPS,
) -> None:
    """Generate an animated voltage-trace video from a simulation output file.

    The animation shows Vm being drawn progressively across the full time axis.
    The file must have a whitespace-separated header with at least ``time`` and
    ``Vm`` columns (standard cardiacFoamEP output format).

    Parameters
    ----------
    data_file:
        Path to the ``.txt`` output produced by cardiacFoam.
    output_path:
        Destination video file.  Extension determines format:
        ``.mp4`` (requires ffmpeg) or ``.gif`` (requires Pillow).
    title:
        Optional plot title (e.g. the case_id).
    n_frames:
        Number of frames in the animation.
    fps:
        Playback speed in frames per second.
    """
    data_file = Path(data_file)
    output_path = Path(output_path)

    data = np.genfromtxt(data_file, names=True)
    time = data["time"]
    vm   = data["Vm"]

    # ── subsample to n_frames evenly-spaced indices ─────────────────────────
    indices = np.unique(
        np.linspace(1, len(time) - 1, n_frames, dtype=int)
    )

    # ── figure setup ─────────────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(10, 4))
    ax.set_xlim(time[0], time[-1])
    v_range = vm.max() - vm.min()
    ax.set_ylim(
        vm.min() - 0.05 * v_range,
        vm.max() + 0.05 * v_range,
    )
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Vm (mV)")
    if title:
        ax.set_title(title)
    (line,) = ax.plot([], [], color=LINE_COLOR, lw=LINE_WIDTH)
    fig.tight_layout()

    # ── animation callbacks ───────────────────────────────────────────────────
    def _init():
        line.set_data([], [])
        return (line,)

    def _update(frame_idx: int):
        n = indices[frame_idx]
        line.set_data(time[:n], vm[:n])
        return (line,)

    ani = animation.FuncAnimation(
        fig,
        _update,
        init_func=_init,
        frames=len(indices),
        interval=1000 // fps,
        blit=True,
    )

    # ── save ─────────────────────────────────────────────────────────────────
    suffix = output_path.suffix.lower()
    try:
        if suffix == ".gif":
            ani.save(output_path, writer="pillow", fps=fps)
        else:
            writer = animation.FFMpegWriter(fps=fps, bitrate=1800)
            ani.save(output_path, writer=writer)
        print(f"Video saved: {output_path}")
    except Exception as exc:
        print(f"Warning: could not save video ({exc}). Is ffmpeg installed?")
    finally:
        plt.close(fig)


# ── CLI entry point ───────────────────────────────────────────────────────────

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python animate_trace.py <trace.txt> <output.mp4>")
        sys.exit(1)
    create_voltage_animation(sys.argv[1], sys.argv[2])
