#!/usr/bin/env python3
import argparse
import sys
from pathlib import Path
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# Ensure plot_pseudo_ecg can be imported
curr_dir = Path(__file__).resolve().parent
if str(curr_dir) not in sys.path:
    sys.path.insert(0, str(curr_dir))

try:
    from plot_pseudo_ecg import discover_input, parse_pseudo_ecg, select_signals, _FILENAMES, _DEFAULT_SCALE
except ImportError as e:
    print(f"Error importing from plot_pseudo_ecg: {e}")
    print("Ensure plot_pseudo_ecg.py is in the same directory.")
    sys.exit(1)


def main():
    parser = argparse.ArgumentParser(description="Animate ECG time series to a video.")
    parser.add_argument("--case-dir", default=".", help="Case directory.")
    
    mode_group = parser.add_mutually_exclusive_group()
    mode_group.add_argument("--monodomain", action="store_true", help="Read pseudoECG.dat")
    mode_group.add_argument("--eikonal", action="store_true", help="Read eikonalECG.dat")
    
    parser.add_argument("--input", default=None, help="Explicit path to ECG .dat file.")
    parser.add_argument("--scale", type=float, default=None, help="Scale factor (default 1000.0 for V->mV).")
    parser.add_argument("--electrodes", default=None, help="Comma-separated electrode names to plot.")
    parser.add_argument("--tmin", type=float, default=None, help="Minimum time to plot.")
    parser.add_argument("--tmax", type=float, default=None, help="Maximum time to plot.")
    
    # New options for animation
    parser.add_argument("--fps", type=int, default=30, help="Frames per second for the output video.")
    parser.add_argument("--step", type=int, default=1, help="Plot every N-th time step to speed up animation generation.")
    parser.add_argument("--output", default=None, help="Output video file (e.g., .mp4 or .gif).")
    
    args = parser.parse_args()
    
    mode = "eikonal" if args.eikonal else "monodomain"
    case_dir = Path(args.case_dir).resolve()
    
    if args.input:
        input_path = Path(args.input).resolve()
    else:
        input_path = discover_input(case_dir, _FILENAMES[mode])
        
    scale = args.scale if args.scale is not None else _DEFAULT_SCALE[mode]
    
    if args.output:
        output_path = (case_dir / args.output).resolve()
    else:
        stem = input_path.stem
        output_path = (case_dir / "postProcessing" / f"{stem}_video.mp4").resolve()
        
    header, rows = parse_pseudo_ecg(input_path)
    signal_names = header[1:]
    
    filtered_rows = []
    for row in rows:
        t = row[0]
        if args.tmin is not None and t < args.tmin:
            continue
        if args.tmax is not None and t > args.tmax:
            continue
        filtered_rows.append(row)
        
    if not filtered_rows:
        raise ValueError("No samples remain after filtering time.")
        
    chosen = select_signals(signal_names, args.electrodes)
    chosen_idx = [signal_names.index(name) for name in chosen]
    
    time = [r[0] for r in filtered_rows]
    chosen_values = [[r[i + 1] for i in chosen_idx] for r in filtered_rows]
    
    # Calculate limits
    all_vals = [v * scale for row in chosen_values for v in row]
    ymin, ymax = min(all_vals), max(all_vals)
    if ymin == ymax:
        pad = 0.1
    else:
        pad = 0.05 * (ymax - ymin)
    y_limits = (ymin - pad, ymax + pad)
    
    time_ms = [t * 1000.0 for t in time]
    
    n = len(chosen)
    # Configure subplots (side-by-side)
    fig, axes = plt.subplots(1, n, figsize=(2.5 * n, 4), sharey=True, squeeze=False)
    fig.suptitle(f"{mode} ECG Animation", fontsize=12)
    
    lines = []
    for i, name in enumerate(chosen):
        ax = axes[0, i]
        ax.set_xlim(min(time_ms), max(time_ms))
        ax.set_ylim(y_limits)
        ax.set_title(name, fontsize=10)
        if i == 0:
            ax.set_ylabel("mV", fontsize=8)
        ax.set_xlabel("time (ms)", fontsize=10)
        ax.grid(True, alpha=0.3)
        # Empty line to update during animation
        line, = ax.plot([], [], lw=1.5)
        # Store full y values
        lines.append((line, [row[i] * scale for row in chosen_values]))
        
    fig.tight_layout()
    
    # Filter by step to reduce frames if requested
    frame_indices = list(range(0, len(time_ms), args.step))
    
    def init():
        for line, _ in lines:
            line.set_data([], [])
        return [line for line, _ in lines]
        
    def update(frame_idx):
        for line, vals in lines:
            # Plot from 0 up to current frame_idx
            line.set_data(time_ms[:frame_idx+1], vals[:frame_idx+1])
        return [line for line, _ in lines]

    ani = animation.FuncAnimation(
        fig, update, frames=frame_indices, init_func=init, blit=True
    )
    
    output_path.parent.mkdir(parents=True, exist_ok=True)
    print(f"Generating video with {len(frame_indices)} frames at {args.fps} FPS...")
    print(f"Saving to {output_path}...")
    
    if output_path.suffix.lower() == '.gif':
        writer = animation.PillowWriter(fps=args.fps)
    else:
        # Default to ffmpeg for mp4 or other formats
        try:
            writer = animation.FFMpegWriter(fps=args.fps)
        except Exception as e:
            print(f"Failed to use FFMpegWriter: {e}. Falling back to Pillow (will save as GIF).")
            output_path = output_path.with_suffix('.gif')
            writer = animation.PillowWriter(fps=args.fps)
            
    try:
        ani.save(str(output_path), writer=writer)
        print(f"Done! Saved {output_path}")
    except Exception as e:
        print(f"Failed to save video: {e}")
        print("Note: If saving as .mp4 fails, you might need to install ffmpeg (e.g. `brew install ffmpeg` or `apt-get install ffmpeg`).")
        print("Alternatively, try changing the --output extension to .gif")

if __name__ == "__main__":
    main()
