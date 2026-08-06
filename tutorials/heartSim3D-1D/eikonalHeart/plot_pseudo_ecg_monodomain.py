#!/usr/bin/env python3
"""Plot each precordial electrode (V1-V6) from postProcessing/pseudoECG.dat.

Usage:
    python3 plot_pseudo_ecg_monodomain.py [path/to/pseudoECG.dat]

Defaults to postProcessing/pseudoECG.dat next to this script.
"""

import sys
from pathlib import Path

import matplotlib.pyplot as plt


def load(path: Path):
    with open(path) as f:
        lines = [line.strip() for line in f if line.strip()]

    header = lines[0].lstrip("#").split()
    rows = [list(map(float, line.split())) for line in lines[1:]]
    columns = list(zip(*rows))
    return header, columns


def main() -> None:
    default_path = Path(__file__).parent / "postProcessing" / "pseudoECG.dat"
    path = Path(sys.argv[1]) if len(sys.argv) > 1 else default_path

    header, columns = load(path)
    time = columns[0]
    electrode_names = header[1:]
    electrode_values = columns[1:]

    fig, axes = plt.subplots(len(electrode_names), 1, figsize=(9, 12), sharex=True)

    for ax, name, values in zip(axes, electrode_names, electrode_values):
        ax.plot(time, values, color="black", linewidth=1.2)
        ax.axhline(0, color="gray", linewidth=0.6, linestyle="--")
        ax.set_ylabel(name)
        ax.grid(True, alpha=0.3)

    axes[-1].set_xlabel("time (s)")
    fig.suptitle(f"Pseudo-ECG: {path}")
    fig.tight_layout()

    out_path = path.with_suffix(".png")
    fig.savefig(out_path, dpi=150)
    print(f"Saved plot to {out_path}")

    plt.show()


if __name__ == "__main__":
    main()
