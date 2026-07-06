from __future__ import annotations

from pathlib import Path


def collect_plots(case_dir: Path) -> list[Path]:
    """Return existing PNG figures under the case's postProcessing tree."""
    pp = Path(case_dir) / "postProcessing"
    if not pp.is_dir():
        return []
    return sorted(p for p in pp.rglob("*.png") if p.is_file())


def generate_plot(dat_path: Path, out_png: Path) -> Path | None:
    """Draw a simple line plot from a .dat file. Returns None if matplotlib
    is unavailable or the data cannot be plotted."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception:
        return None
    from .metrics import load_dat
    try:
        cols, arr = load_dat(dat_path)
    except (ValueError, OSError):
        return None
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots()
    for j in range(1, arr.shape[1]):
        ax.plot(arr[:, 0], arr[:, j], label=cols[j])
    ax.set_xlabel(cols[0])
    ax.legend(fontsize="small")
    fig.savefig(out_png, dpi=100, bbox_inches="tight")
    plt.close(fig)
    return out_png
