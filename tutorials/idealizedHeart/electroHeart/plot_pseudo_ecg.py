#!/usr/bin/env python3
from __future__ import annotations

import argparse
import math
from pathlib import Path
from typing import List, Optional, Sequence, Tuple

_FILENAMES = {
    "monodomain": "pseudoECG.dat",
    "eikonal":    "eikonalECG.dat",
}
_DEFAULT_SCALE = {
    "monodomain": 1000.0,
    "eikonal":    1000.0,
}


def ensure_matplotlib():
    try:
        import matplotlib.pyplot as plt  # type: ignore
    except ModuleNotFoundError as exc:
        raise ModuleNotFoundError(
            "matplotlib is required for plotting. Install it with: "
            "python3 -m pip install matplotlib"
        ) from exc
    return plt


def discover_input(case_dir: Path, filename: str) -> Path:
    candidates: List[Path] = []

    root_file = case_dir / "postProcessing" / filename
    if root_file.exists():
        candidates.append(root_file)

    proc_files = sorted(case_dir.glob(f"processor*/postProcessing/{filename}"))
    candidates.extend(proc_files)

    if not candidates:
        raise FileNotFoundError(
            f"Could not find {filename} in {case_dir} (root or processor*/postProcessing)."
        )

    non_empty = [p for p in candidates if p.stat().st_size > 0]
    if non_empty:
        return non_empty[0]

    return candidates[0]


def parse_pseudo_ecg(path: Path) -> Tuple[List[str], List[List[float]]]:
    header: Optional[List[str]] = None
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


def select_signals(all_names: Sequence[str], requested: Optional[str]) -> List[str]:
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
    scale: float = 1000.0,
    y_limits: Optional[Tuple[float, float]] = None,
) -> None:
    plt = ensure_matplotlib()
    n = len(names)
    if ncols <= 0:
        ncols = n
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
        y = [row[i] * scale for row in values]
        ax.plot(time_ms, y, lw=1.1)
        ax.set_title(name, fontsize=9)
        ax.grid(True, alpha=0.25)
        ax.tick_params(labelsize=8)
        ax.set_ylabel("mV", fontsize=8)
        if y_limits is not None:
            ax.set_ylim(y_limits[0] * scale, y_limits[1] * scale)

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
            "Plot ECG time series from pseudoECG.dat or eikonalECG.dat. "
            "Use --monodomain or --eikonal to select the source type."
        )
    )
    parser.add_argument(
        "--case-dir",
        default=".",
        help="Case directory (default: current directory).",
    )

    mode_group = parser.add_mutually_exclusive_group()
    mode_group.add_argument(
        "--monodomain",
        action="store_true",
        help="Read pseudoECG.dat (monodomain solver output, data in V, scaled ×1000 to mV).",
    )
    mode_group.add_argument(
        "--eikonal",
        action="store_true",
        help="Read eikonalECG.dat (eikonal surrogate output, data in V post-fix, scaled ×1000 to mV).",
    )

    parser.add_argument(
        "--input",
        default=None,
        help="Explicit path to ECG .dat file; overrides --monodomain/--eikonal auto-discovery.",
    )
    parser.add_argument(
        "--scale",
        type=float,
        default=None,
        help=(
            "Multiply raw signal values by this factor before plotting. "
            "Defaults to 1000.0 (V → mV) for both modes. "
            "Use --scale 1.0 for pre-fix eikonalECG.dat whose values are already in mV."
        ),
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
        help="Number of subplot columns. Default: 0 (one row).",
    )
    parser.add_argument(
        "--group-size",
        type=int,
        default=0,
        help=(
            "Split selected electrodes into subsets of this size, one image per subset. "
            "0 means all in one figure."
        ),
    )
    parser.add_argument(
        "--output",
        default=None,
        help=(
            "Output image path. Defaults to postProcessing/pseudoECG_plots.png "
            "or postProcessing/eikonalECG_plots.png depending on mode."
        ),
    )
    parser.add_argument(
        "--show",
        action="store_true",
        help="Display the figure interactively in addition to saving.",
    )
    args = parser.parse_args()

    case_dir = Path(args.case_dir).resolve()

    # Resolve mode
    if args.eikonal:
        mode = "eikonal"
    else:
        mode = "monodomain"

    # Resolve input file
    if args.input:
        input_path = Path(args.input).resolve()
    else:
        input_path = discover_input(case_dir, _FILENAMES[mode])

    # Resolve scale
    scale = args.scale if args.scale is not None else _DEFAULT_SCALE[mode]

    # Resolve output path
    if args.output:
        output_path = (case_dir / args.output).resolve()
    else:
        stem = input_path.stem
        output_path = (case_dir / "postProcessing" / f"{stem}_plots.png").resolve()

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

    all_vals = [v for row in chosen_values for v in row]
    ymin = min(all_vals)
    ymax = max(all_vals)
    if math.isclose(ymin, ymax, rel_tol=0.0, abs_tol=1e-14):
        pad = max(1e-9, abs(ymin) * 0.05 + 1e-9)
        y_limits = (ymin - pad, ymax + pad)
    else:
        pad = 0.05 * (ymax - ymin)
        y_limits = (ymin - pad, ymax + pad)

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

        title = f"{mode} ECG ({input_path.name})"
        plot_group(
            time,
            group_names,
            group_values,
            args.ncols,
            title,
            out,
            scale=scale,
            y_limits=y_limits,
        )
        print(f"Wrote {out}")

    if args.show:
        plt = ensure_matplotlib()
        plt.show()

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
