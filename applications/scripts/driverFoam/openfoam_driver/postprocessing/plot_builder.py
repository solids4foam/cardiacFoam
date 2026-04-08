"""
plot_builder — declarative Plotly figure construction for post-processing scripts.

Motivation
----------
Tutorial post-processing scripts share a common pattern:

1. Load a set of CSV (or similar) files from the output directory.
2. Extract metadata (e.g. mesh size DX, time step DT) from each filename.
3. Build Plotly traces with a group-based color strategy (one base color per
   group, lightened by variant position within the group).
4. Apply a shared layout and write an interactive HTML file.
5. Return an artifact dict for inclusion in ``plots.json``.

``plot_builder`` provides the building blocks that make this pattern
re-usable across tutorials without duplicating code.

Public API
----------
GroupShadedColors
    Assigns a distinct base color to each group (e.g. DX value) and lightens
    it proportionally for each shade variant within that group (e.g. DT value).

TraceSpec
    Describes which DataFrame columns map to X and Y, including optional unit
    scaling.

PlotSpec
    Full figure description: title, axis labels, output filename, and artifact
    metadata — everything needed to apply layout and register the output.

load_csv_folder(folder, glob_pattern, *, drop_columns)
    Load all CSV files matching a glob pattern from a directory, optionally
    dropping irrelevant columns.

build_line_traces(data, *, trace, name_fn, colors, group_fn, shade_fn)
    Turn a list of (filename, DataFrame) pairs into Plotly Scatter traces,
    with optional group-shaded coloring.

write_plot(fig, spec, output_dir, *, updatemenus, show)
    Apply PlotSpec layout to a figure, write HTML, and return the artifact dict.

Typical script pattern
----------------------
::

    from openfoam_driver.postprocessing.plot_builder import (
        GroupShadedColors, TraceSpec, PlotSpec,
        load_csv_folder, build_line_traces, write_plot,
    )
    from openfoam_driver.postprocessing.plotting_common import extract_dx_dt

    def run_postprocessing(*, output_dir, setup_root=None, **_):
        folder = Path(output_dir)
        data = load_csv_folder(folder, "*line*.csv", drop_columns=["Point ID"])

        group_fn = lambda f: extract_dx_dt(f)[0]  # DX
        shade_fn = lambda f: extract_dx_dt(f)[1]  # DT
        colors = GroupShadedColors.from_groups_and_traces(
            [f for f, _ in data], group_fn, shade_fn
        )

        trace = TraceSpec(x_col=1, y_col=0, x_scale=1000.0, y_scale=1000.0)
        name_fn = lambda f: f"DX={extract_dx_dt(f)[0]:.1f} mm, DT={extract_dx_dt(f)[1]:.3f} ms"
        traces = build_line_traces(
            data, trace=trace, name_fn=name_fn,
            colors=colors, group_fn=group_fn, shade_fn=shade_fn,
        )

        spec = PlotSpec(
            title="Activation time along diagonal",
            xaxis_title="Distance (mm)",
            yaxis_title="Activation Time (ms)",
            output_filename="allSimulations.html",
            label="All simulations",
        )

        import plotly.graph_objects as go
        fig = go.Figure(traces)
        artifact = write_plot(fig, spec, folder)
        return [artifact]
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable, Iterable

import pandas as pd
import plotly.graph_objects as go

from .plotting_common import build_visibility_mask, lighten_hex_color, ordered_unique
from .style import apply_plotly_layout, write_plotly_html


# ---------------------------------------------------------------------------
# Default color palette
# ---------------------------------------------------------------------------

#: Default base-color palette used by :class:`GroupShadedColors` when no
#: explicit palette is supplied.  Colors are chosen to be visually distinct
#: on a white background.
DEFAULT_PALETTE: list[str] = [
    "#1f77b4",  # blue
    "#2ca02c",  # green
    "#eb1616",  # red
    "#9467bd",  # purple
    "#8c564b",  # brown
    "#e377c2",  # pink
    "#17becf",  # cyan
    "#bcbd22",  # yellow-green
]


# ---------------------------------------------------------------------------
# GroupShadedColors
# ---------------------------------------------------------------------------

class GroupShadedColors:
    """Assign colors to traces using a two-level group + shade strategy.

    Each unique *group* value (e.g. a DX mesh size) receives a distinct base
    color from the palette.  Within a group, each *shade* value (e.g. a DT
    time step) is mapped to a lightened variant of that base color, so that
    coarser/finer variants are visually distinguishable while remaining
    obviously part of the same group.

    Parameters
    ----------
    palette:
        List of hex color strings (``"#rrggbb"``).  Defaults to
        :data:`DEFAULT_PALETTE`.
    max_shade_amount:
        Maximum lightening fraction applied to the darkest shade variant.
        ``0.0`` means no lightening; ``1.0`` would produce white.  The default
        ``0.4`` produces a noticeable but still readable lightening.

    Example
    -------
    ::

        from openfoam_driver.postprocessing.plotting_common import extract_dx_dt

        filenames = ["case_DX1_DT005.csv", "case_DX1_DT01.csv", "case_DX2_DT005.csv"]
        group_fn = lambda f: extract_dx_dt(f)[0]   # DX value
        shade_fn = lambda f: extract_dx_dt(f)[1]   # DT value

        colors = GroupShadedColors.from_groups_and_traces(filenames, group_fn, shade_fn)
        c = colors.color_for(0.1, 0.005)  # darkest blue for smallest DX, smallest DT
    """

    def __init__(
        self,
        palette: list[str] | None = None,
        max_shade_amount: float = 0.4,
    ) -> None:
        self._palette = palette or DEFAULT_PALETTE
        self._max_shade = max_shade_amount
        self._group_order: list[Any] = []
        self._shades: dict[Any, list[Any]] = {}

    # ------------------------------------------------------------------
    # Registration
    # ------------------------------------------------------------------

    def register_groups(self, groups: Iterable[Any]) -> None:
        """Register group keys in display order (first seen = first color)."""
        for g in ordered_unique(list(groups)):
            if g not in self._group_order:
                self._group_order.append(g)

    def register_shades(self, group: Any, shades: Iterable[Any]) -> None:
        """Register shade variants for one group, sorted ascending."""
        self._shades[group] = sorted(set(shades))

    # ------------------------------------------------------------------
    # Color query
    # ------------------------------------------------------------------

    def color_for(self, group: Any, shade: Any) -> str:
        """Return the ``rgb(r,g,b)`` color string for a given (group, shade) pair."""
        if group in self._group_order:
            group_index = self._group_order.index(group)
        else:
            group_index = 0
        base_color = self._palette[group_index % len(self._palette)]

        shade_list = self._shades.get(group, [shade])
        try:
            shade_index = shade_list.index(shade)
        except ValueError:
            shade_index = 0
        shade_count = max(len(shade_list) - 1, 1)
        amount = (shade_index / shade_count) * self._max_shade

        return lighten_hex_color(base_color, amount)

    def base_color_for(self, group: Any) -> str:
        """Return the unmodified base color for a group (shade position 0)."""
        if group in self._group_order:
            group_index = self._group_order.index(group)
        else:
            group_index = 0
        return self._palette[group_index % len(self._palette)]

    # ------------------------------------------------------------------
    # Factory
    # ------------------------------------------------------------------

    @classmethod
    def from_groups_and_traces(
        cls,
        filenames: Iterable[str],
        group_fn: Callable[[str], Any],
        shade_fn: Callable[[str], Any],
        palette: list[str] | None = None,
        max_shade_amount: float = 0.4,
    ) -> "GroupShadedColors":
        """Build and fully register a :class:`GroupShadedColors` from filenames.

        Parameters
        ----------
        filenames:
            Iterable of filename strings (basenames or full paths — only the
            filename part is passed to the key functions).
        group_fn:
            Callable ``filename -> group_key`` (e.g. returns the DX float).
        shade_fn:
            Callable ``filename -> shade_key`` (e.g. returns the DT float).
        palette:
            Optional explicit palette; defaults to :data:`DEFAULT_PALETTE`.
        max_shade_amount:
            Forwarded to the constructor.

        Returns
        -------
        Fully populated :class:`GroupShadedColors` ready for :meth:`color_for`.
        """
        colors = cls(palette=palette, max_shade_amount=max_shade_amount)
        names = list(filenames)
        groups = [group_fn(f) for f in names]
        colors.register_groups(groups)
        for g in ordered_unique(groups):
            shades = [shade_fn(f) for f, gg in zip(names, groups) if gg == g]
            colors.register_shades(g, shades)
        return colors


# ---------------------------------------------------------------------------
# CSV loading
# ---------------------------------------------------------------------------

def load_csv_folder(
    folder: str | Path,
    glob_pattern: str,
    *,
    drop_columns: list[str] | None = None,
) -> list[tuple[str, pd.DataFrame]]:
    """Load CSV files matching a glob pattern from *folder*.

    Parameters
    ----------
    folder:
        Directory to search.
    glob_pattern:
        Glob pattern relative to *folder* (e.g. ``"*line*.csv"`` or
        ``"*points_DT*_DX*.csv"``).
    drop_columns:
        Column names to discard from each DataFrame.  Columns that are absent
        in a particular file are silently skipped.

    Returns
    -------
    Sorted list of ``(basename, DataFrame)`` pairs.  Returns an empty list (with
    a warning printed) if no files match the pattern.

    Notes
    -----
    Common OpenFOAM/ParaView export columns that carry no simulation meaning
    (``"Point ID"``, ``"vtkValidPointMask"``, ``"Points_Magnitude"``) can be
    passed as *drop_columns* to keep DataFrames clean.
    """
    folder = Path(folder)
    paths = sorted(folder.glob(glob_pattern))
    if not paths:
        print(f"[plot_builder] No files match '{glob_pattern}' in {folder}")
        return []

    drop = set(drop_columns or [])
    result: list[tuple[str, pd.DataFrame]] = []
    for p in paths:
        df = pd.read_csv(p)
        to_drop = [c for c in drop if c in df.columns]
        if to_drop:
            df = df.drop(columns=to_drop)
        result.append((p.name, df))

    print(f"[plot_builder] Loaded {len(result)} file(s) matching '{glob_pattern}'")
    return result


# ---------------------------------------------------------------------------
# TraceSpec
# ---------------------------------------------------------------------------

@dataclass
class TraceSpec:
    """Describes how to extract X and Y series from a DataFrame for one trace.

    Attributes
    ----------
    x_col:
        Column name (``str``) or integer index selecting the X column.
    y_col:
        Column name (``str``) or integer index selecting the Y column.
    x_scale:
        Multiplicative scale applied to all X values after loading.
        Use ``1000.0`` to convert metres to millimetres, etc.
    y_scale:
        Multiplicative scale applied to all Y values.
    mode:
        Plotly trace mode string: ``"lines"``, ``"markers"``, or
        ``"lines+markers"`` (default).
    """

    x_col: str | int = 0
    y_col: str | int = 1
    x_scale: float = 1.0
    y_scale: float = 1.0
    mode: str = "lines+markers"

    def x(self, df: pd.DataFrame) -> pd.Series:
        """Return the scaled X series from *df*."""
        col = df.columns[self.x_col] if isinstance(self.x_col, int) else self.x_col
        return df[col] * self.x_scale

    def y(self, df: pd.DataFrame) -> pd.Series:
        """Return the scaled Y series from *df*."""
        col = df.columns[self.y_col] if isinstance(self.y_col, int) else self.y_col
        return df[col] * self.y_scale


# ---------------------------------------------------------------------------
# Trace building
# ---------------------------------------------------------------------------

def build_line_traces(
    data: list[tuple[str, pd.DataFrame]],
    *,
    trace: TraceSpec,
    name_fn: Callable[[str], str],
    colors: GroupShadedColors | None = None,
    group_fn: Callable[[str], Any] | None = None,
    shade_fn: Callable[[str], Any] | None = None,
) -> list[go.Scatter]:
    """Build Plotly ``Scatter`` traces from a list of ``(filename, DataFrame)`` pairs.

    Parameters
    ----------
    data:
        Output of :func:`load_csv_folder` — sorted ``(basename, df)`` pairs.
    trace:
        :class:`TraceSpec` describing which columns to use and any scaling.
    name_fn:
        Callable ``filename -> trace_name`` used to populate the legend.
    colors:
        Optional :class:`GroupShadedColors` instance.  When provided along
        with *group_fn* and *shade_fn*, each trace receives a color determined
        by its group and shade position.  When omitted, Plotly's default color
        cycle is used.
    group_fn:
        Callable ``filename -> group_key`` (e.g. returns the DX float).
        Required when *colors* is provided.
    shade_fn:
        Callable ``filename -> shade_key`` (e.g. returns the DT float).
        Required when *colors* is provided.

    Returns
    -------
    List of :class:`plotly.graph_objects.Scatter` traces, one per input file.
    """
    use_colors = colors is not None and group_fn is not None and shade_fn is not None
    result: list[go.Scatter] = []

    for filename, df in data:
        x = trace.x(df)
        y = trace.y(df)
        name = name_fn(filename)
        kw: dict[str, Any] = {"mode": trace.mode}

        if use_colors:
            color = colors.color_for(group_fn(filename), shade_fn(filename))
            kw["line"] = dict(color=color)
            kw["marker"] = dict(color=color)

        result.append(go.Scatter(x=x, y=y, name=name, **kw))

    return result


# ---------------------------------------------------------------------------
# PlotSpec
# ---------------------------------------------------------------------------

@dataclass
class PlotSpec:
    """Complete description of a single Plotly HTML figure.

    :class:`PlotSpec` captures everything needed to apply a layout and register
    the output as a post-processing artifact.  It does *not* hold trace data —
    traces are built separately and added to a ``go.Figure`` before calling
    :func:`write_plot`.

    Attributes
    ----------
    title:
        Figure title displayed at the top.
    xaxis_title:
        Label for the X axis.
    yaxis_title:
        Label for the Y axis.
    output_filename:
        Basename of the HTML file written to the output directory
        (e.g. ``"cardiacFoam_allSimulations.html"``).
    label:
        Human-readable label used in the ``plots.json`` artifact entry.
    legend_title:
        Optional text shown above the legend.
    kind:
        Artifact kind for ``plots.json``: ``"plot"`` (default), ``"table"``,
        ``"data"``, or ``"report"``.
    height:
        Optional figure height in pixels.
    width:
        Optional figure width in pixels.
    template:
        Plotly layout template.  Defaults to ``"plotly_white"``.
    """

    title: str
    xaxis_title: str
    yaxis_title: str
    output_filename: str
    label: str
    legend_title: str | None = None
    kind: str = "plot"
    height: int | None = None
    width: int | None = None
    template: str = "plotly_white"

    def artifact(self, output_dir: Path) -> dict[str, Any]:
        """Return the artifact dict for this figure (for inclusion in ``plots.json``)."""
        return {
            "path": str(output_dir / self.output_filename),
            "label": self.label,
            "kind": self.kind,
            "format": "html",
        }


# ---------------------------------------------------------------------------
# write_plot
# ---------------------------------------------------------------------------

def write_plot(
    fig: go.Figure,
    spec: PlotSpec,
    output_dir: Path,
    *,
    updatemenus: list[dict[str, Any]] | None = None,
    show: bool = False,
) -> dict[str, Any]:
    """Apply layout from *spec*, write HTML, and return the artifact dict.

    Parameters
    ----------
    fig:
        Plotly figure with all traces already added.
    spec:
        :class:`PlotSpec` that provides title, axis labels, and output path.
    output_dir:
        Directory into which the HTML file is written.
    updatemenus:
        Optional list of Plotly ``updatemenus`` entries (toggle buttons).
        Pass the output of :func:`make_toggle_button` here.
    show:
        If ``True``, also call ``fig.show()`` after writing.

    Returns
    -------
    Artifact dict with keys ``path``, ``label``, ``kind``, ``format`` — ready
    to be returned from a ``run_postprocessing`` function.
    """
    apply_plotly_layout(
        fig,
        title=spec.title,
        xaxis_title=spec.xaxis_title,
        yaxis_title=spec.yaxis_title,
        legend_title=spec.legend_title,
        template=spec.template,
        showlegend=True,
        updatemenus=updatemenus,
        width=spec.width,
        height=spec.height,
    )
    output_path = output_dir / spec.output_filename
    write_plotly_html(fig, output_path)
    if show:
        fig.show()
    return spec.artifact(output_dir)


# ---------------------------------------------------------------------------
# Toggle button helper
# ---------------------------------------------------------------------------

def make_toggle_button(
    fig: go.Figure,
    *,
    label_a: str,
    label_b: str,
    indices_a: Iterable[int],
    indices_b: Iterable[int],
    names_a: list[str] | None = None,
    names_b: list[str] | None = None,
    x: float = 0.2,
    y: float = 1.0,
) -> list[dict[str, Any]]:
    """Build a two-state Plotly toggle button (``updatemenus`` entry).

    The button switches between two sets of visible traces.  This is the
    general form of the "All simulations / Niederer vs cardiacFoam" toggle
    in the Niederer benchmark scripts.

    Parameters
    ----------
    fig:
        The figure whose traces the button controls.  Used to read the
        current total trace count.
    label_a:
        Button label for the first view (typically the "all" view).
    label_b:
        Button label for the second view (typically the "comparison" view).
    indices_a:
        Trace indices that should be visible in the *label_a* view.
    indices_b:
        Trace indices that should be visible in the *label_b* view.
    names_a:
        Optional list of trace name strings to apply when switching to view A.
        Must have the same length as ``len(fig.data)`` if provided.
    names_b:
        Optional list of trace name strings to apply when switching to view B.
    x, y:
        Button placement in figure coordinates (Plotly ``paper`` units).

    Returns
    -------
    A list with one dict, suitable for passing as ``updatemenus`` to
    :func:`write_plot` or ``fig.update_layout``.
    """
    n = len(fig.data)
    set_a = set(indices_a)
    set_b = set(indices_b)
    vis_a = build_visibility_mask(set_a, n)
    vis_b = build_visibility_mask(set_b, n)

    args_a: dict[str, Any] = {"visible": vis_a}
    args_b: dict[str, Any] = {"visible": vis_b}
    if names_a is not None:
        args_a["name"] = names_a
    if names_b is not None:
        args_b["name"] = names_b

    return [
        {
            "type": "buttons",
            "direction": "right",
            "x": x,
            "y": y,
            "xanchor": "left",
            "yanchor": "top",
            "showactive": False,
            "buttons": [
                {"label": label_a, "method": "update", "args": [args_a]},
                {"label": label_b, "method": "update", "args": [args_b]},
            ],
        }
    ]
