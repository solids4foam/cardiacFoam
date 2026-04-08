"""Shared utilities for tutorial post-processing scripts."""
from __future__ import annotations

from typing import Protocol, runtime_checkable

from .driver import PostprocessTask, run_postprocess_tasks
from .plot_builder import (
    DEFAULT_PALETTE,
    GroupShadedColors,
    PlotSpec,
    TraceSpec,
    build_line_traces,
    load_csv_folder,
    make_toggle_button,
    write_plot,
)
from .style import (
    apply_plotly_layout,
    configure_matplotlib_defaults,
    finalize_matplotlib_figure,
    style_matplotlib_axes,
    write_plotly_html,
)


@runtime_checkable
class PostprocessingProtocol(Protocol):
    """Typing stub for tutorial post-processing entry points.

    Every tutorial postprocessing script must expose a function named
    ``run_postprocessing`` that satisfies this signature.  Scripts are not
    required to subclass this Protocol — IDEs and mypy will flag mismatches
    when type-checking is enabled.
    """

    def __call__(
        self,
        *,
        output_dir: str,
        setup_root: str | None = None,
        **kwargs: object,
    ) -> list[dict]: ...


__all__ = [
    # contract
    "PostprocessingProtocol",
    # driver
    "PostprocessTask",
    "run_postprocess_tasks",
    # plot_builder — declarative Plotly helpers
    "DEFAULT_PALETTE",
    "GroupShadedColors",
    "PlotSpec",
    "TraceSpec",
    "build_line_traces",
    "load_csv_folder",
    "make_toggle_button",
    "write_plot",
    # style
    "apply_plotly_layout",
    "configure_matplotlib_defaults",
    "finalize_matplotlib_figure",
    "style_matplotlib_axes",
    "write_plotly_html",
]
