from __future__ import annotations

from collections.abc import Mapping
from pathlib import Path
from typing import Any


def apply_plotly_layout(
    fig: Any,
    *,
    title: str,
    xaxis_title: str,
    yaxis_title: str,
    template: str = "plotly_white",
    legend_title: str | None = None,
    showlegend: bool = True,
    legend: Mapping[str, Any] | None = None,
    updatemenus: list[dict[str, Any]] | None = None,
    width: int | None = None,
    height: int | None = None,
    extra_layout: Mapping[str, Any] | None = None,
) -> None:
    layout: dict[str, Any] = {
        "title": title,
        "xaxis_title": xaxis_title,
        "yaxis_title": yaxis_title,
        "template": template,
        "showlegend": showlegend,
    }
    if legend_title is not None:
        layout["legend_title"] = legend_title
    if legend is not None:
        layout["legend"] = legend
    if updatemenus is not None:
        layout["updatemenus"] = updatemenus
    if width is not None:
        layout["width"] = width
    if height is not None:
        layout["height"] = height
    if extra_layout:
        layout.update(dict(extra_layout))
    fig.update_layout(**layout)


def write_plotly_html(fig: Any, output_path: str | Path) -> Path:
    destination = Path(output_path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    fig.write_html(str(destination))
    return destination


def configure_matplotlib_defaults() -> None:
    import matplotlib as mpl

    mpl.rcParams.update(
        {
            "axes.grid": True,
            "grid.linestyle": "--",
            "grid.alpha": 0.6,
            "axes.titlesize": 12,
            "axes.labelsize": 11,
            "legend.frameon": False,
        }
    )


def style_matplotlib_axes(
    ax: Any,
    *,
    title: str,
    xlabel: str,
    ylabel: str,
    legend: bool = True,
    grid: bool = True,
    grid_kwargs: Mapping[str, Any] | None = None,
) -> None:
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    if grid:
        kwargs = {"which": "both"}
        if grid_kwargs:
            kwargs.update(dict(grid_kwargs))
        ax.grid(True, **kwargs)
    if legend:
        ax.legend()


def finalize_matplotlib_figure(
    fig: Any,
    *,
    save_path: str | Path | None = None,
    show: bool = True,
    close: bool = False,
    dpi: int = 300,
) -> None:
    fig.tight_layout()
    if save_path is not None:
        destination = Path(save_path)
        destination.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(destination, dpi=dpi)
    if show:
        import matplotlib.pyplot as plt

        plt.show()
    if close:
        import matplotlib.pyplot as plt

        plt.close(fig)
