#!/usr/bin/env python3
from __future__ import annotations

import argparse
import math
from pathlib import Path
from typing import List, Sequence, Tuple

def ensure_matplotlib():
    try:
        import matplotlib.pyplot as plt  # type: ignore
    except ModuleNotFoundError as exc:
        raise ModuleNotFoundError(
            "matplotlib is required for plotting. Install it with: "
            "python3 -m pip install matplotlib"
        ) from exc
    return plt


def discover_input(case_dir: Path) -> Path:
    candidates: List[Path] = []

    root_file = case_dir / "postProcessing" / "pseudoECG.dat"
    if root_file.exists():
        candidates.append(root_file)

    proc_files = sorted(case_dir.glob("processor*/postProcessing/pseudoECG.dat"))
    candidates.extend(proc_files)

    if not candidates:
        raise FileNotFoundError(
            f"Could not find pseudoECG.dat in {case_dir} (root or processor*/postProcessing)."
        )

    non_empty = [p for p in candidates if p.stat().st_size > 0]
    if non_empty:
        return non_empty[0]

    return candidates[0]


def parse_pseudo_ecg(path: Path) -> Tuple[List[str], List[List[float]]]:
    header: List[str] | None = None
    rows: List[List[float]] = []

    with path.open("r", encoding="utf-8", errors="ignore") as f:
        for raw in f:
            line = raw.strip()
            if not line:
                continue

            if line.startswith("#"):
                tokens = line.lstrip("#").strip().split()
                if len(tokens) >= 2 and tokens[0].lower() == "time":
                    header = tokens
                continue

            tokens = line.split()
            if header is None:
                header = ["time"] + [f"signal_{i}" for i in range(1, len(tokens))]

            if len(tokens) != len(header):
                # Ignore incomplete/truncated lines while a run is still writing.
                continue

            try:
                rows.append([float(t) for t in tokens])
            except ValueError:
                continue

    if header is None:
        raise ValueError(f"No valid header found in {path}.")

    if not rows:
        raise ValueError(f"No numeric data rows found in {path}.")

    return header, rows


def select_signals(all_names: Sequence[str], requested: str | None) -> List[str]:
    if requested is None or not requested.strip():
        return list(all_names)

    names = [n.strip() for n in requested.split(",") if n.strip()]
    missing = [n for n in names if n not in all_names]
    if missing:
        raise ValueError(
            f"Requested electrodes not found: {missing}. Available: {list(all_names)}"
        )
    return names


def chunked(items: Sequence[str], group_size: int) -> List[List[str]]:
    if group_size <= 0 or group_size >= len(items):
        return [list(items)]
    return [list(items[i : i + group_size]) for i in range(0, len(items), group_size)]


def plot_group(
    time: Sequence[float],
    names: Sequence[str],
    values: Sequence[Sequence[float]],
    ncols: int,
    title: str,
    output: Path,
    y_limits: Tuple[float, float] | None = None,
) -> None:
    plt = ensure_matplotlib()
    n = len(names)
    if ncols <= 0:
        ncols = n  # Default: horizontal layout (single row)
    ncols = max(1, min(ncols, n))
    nrows = int(math.ceil(n / ncols))

    fig, axes = plt.subplots(
        nrows=nrows,
        ncols=ncols,
        figsize=(1.9 * ncols, 3.2 * nrows),
        sharex=True,
        sharey=True,
        squeeze=False,
    )

    time_ms = [t * 1000.0 for t in time]
    flat_axes = axes.flatten()
    for i, name in enumerate(names):
        ax = flat_axes[i]
        y = [row[i] * 1000.0 for row in values]
        ax.plot(time_ms, y, lw=1.1)
        ax.set_title(name, fontsize=9)
        ax.grid(True, alpha=0.25)
        ax.tick_params(labelsize=8)
        ax.set_ylabel("mV", fontsize=8)
        if y_limits is not None:
            ax.set_ylim(y_limits[0] * 1000.0, y_limits[1] * 1000.0)

    for j in range(n, len(flat_axes)):
        flat_axes[j].axis("off")

    for ax in flat_axes[:n]:
        ax.set_xlabel("time (ms)", fontsize=8)

    fig.suptitle(title, fontsize=10)
    fig.tight_layout()
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=170)
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Read pseudoECG.dat and plot one small time-series graph per electrode. "
            "Default layout is horizontal (one row)."
        )
    )
    parser.add_argument(
        "--case-dir",
        default=".",
        help="Case directory (default: current directory).",
    )
    parser.add_argument(
        "--input",
        default=None,
        help="Path to pseudoECG.dat. If omitted, auto-detect in case-dir.",
    )
    parser.add_argument(
        "--electrodes",
        default=None,
        help="Comma-separated electrode names to plot (subset), e.g. V1,V2,V3.",
    )
    parser.add_argument(
        "--tmin",
        type=float,
        default=None,
        help="Minimum time to plot (seconds).",
    )
    parser.add_argument(
        "--tmax",
        type=float,
        default=None,
        help="Maximum time to plot (seconds).",
    )
    parser.add_argument(
        "--ncols",
        type=int,
        default=0,
        help=(
            "Number of subplot columns. "
            "Default: 0 (auto horizontal, one row)."
        ),
    )
    parser.add_argument(
        "--group-size",
        type=int,
        default=0,
        help=(
            "Split selected electrodes into subsets of this size and write one image per subset. "
            "0 means all selected electrodes in one figure."
        ),
    )
    parser.add_argument(
        "--output",
        default="postProcessing/pseudoECG_plots.png",
        help="Output image path (or base name when --group-size > 0).",
    )
    parser.add_argument(
        "--show",
        action="store_true",
        help="Display the figure interactively in addition to saving.",
    )
    args = parser.parse_args()

    case_dir = Path(args.case_dir).resolve()
    input_path = Path(args.input).resolve() if args.input else discover_input(case_dir)

    header, rows = parse_pseudo_ecg(input_path)
    signal_names = header[1:]
    filtered_rows: List[List[float]] = []
    for row in rows:
        t = row[0]
        if args.tmin is not None and t < args.tmin:
            continue
        if args.tmax is not None and t > args.tmax:
            continue
        filtered_rows.append(row)

    if not filtered_rows:
        raise ValueError("No samples remain after applying tmin/tmax.")

    chosen = select_signals(signal_names, args.electrodes)
    chosen_idx = [signal_names.index(name) for name in chosen]
    time = [r[0] for r in filtered_rows]
    chosen_values = [[r[i + 1] for i in chosen_idx] for r in filtered_rows]

    # Use one common y-scale for all plots.
    all_vals = [v for row in chosen_values for v in row]
    ymin = min(all_vals)
    ymax = max(all_vals)
    if math.isclose(ymin, ymax, rel_tol=0.0, abs_tol=1e-14):
        pad = max(1e-9, abs(ymin) * 0.05 + 1e-9)
        y_limits = (ymin - pad, ymax + pad)
    else:
        pad = 0.05 * (ymax - ymin)
        y_limits = (ymin - pad, ymax + pad)

    output_path = (case_dir / args.output).resolve()
    groups = chunked(chosen, args.group_size)

    for gi, group_names in enumerate(groups, start=1):
        cols = [chosen.index(name) for name in group_names]
        group_values = [[row[c] for c in cols] for row in chosen_values]

        if len(groups) == 1:
            out = output_path
        else:
            out = output_path.with_name(
                f"{output_path.stem}_subset{gi:02d}{output_path.suffix}"
            )

        title = f"pseudoECG ({input_path.name})"
        plot_group(
            time,
            group_names,
            group_values,
            args.ncols,
            title,
            out,
            y_limits=y_limits,
        )
        print(f"Wrote {out}")

    if args.show:
        plt = ensure_matplotlib()
        plt.show()

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
