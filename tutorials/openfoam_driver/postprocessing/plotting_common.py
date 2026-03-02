from __future__ import annotations

import re
from collections.abc import Iterable, Mapping

import plotly.colors as plotly_colors


_DX_PATTERN = re.compile(r"DX(\d+)")
_DT_PATTERN = re.compile(r"DT(\d+)")


def extract_dx_dt(filename: str) -> tuple[float, float]:
    """Extract DX and DT values from filenames like ...DX5...DT005..."""
    dx_match = _DX_PATTERN.search(filename)
    dt_match = _DT_PATTERN.search(filename)
    if not dx_match or not dt_match:
        return float("inf"), float("inf")

    dx = int(dx_match.group(1)) / 10.0
    dt_token = dt_match.group(1)
    dt = int(dt_token) / (10 ** (len(dt_token) - 1))
    return dx, dt


def lighten_hex_color(color: str, amount: float) -> str:
    """Lighten a hex color by amount (0=no change, 1=white)."""
    bounded_amount = min(max(float(amount), 0.0), 1.0)
    red, green, blue = plotly_colors.hex_to_rgb(color)
    light_red = int(red + (255 - red) * bounded_amount)
    light_green = int(green + (255 - green) * bounded_amount)
    light_blue = int(blue + (255 - blue) * bounded_amount)
    return f"rgb({light_red},{light_green},{light_blue})"


def ordered_unique(values: Iterable[str]) -> list[str]:
    seen: set[str] = set()
    result: list[str] = []
    for value in values:
        if value in seen:
            continue
        seen.add(value)
        result.append(value)
    return result


def parse_model_and_cell(
    filename: str,
    *,
    model_map: Mapping[str, str],
    cell_map: Mapping[str, str],
) -> tuple[str, str]:
    file_base = filename.rsplit(".", maxsplit=1)[0]
    parts = file_base.split("_")

    raw_model = parts[0] if parts else "UnknownModel"
    raw_cell = parts[1] if len(parts) > 1 else "UnknownCell"

    model = model_map.get(raw_model, raw_model)
    cell = cell_map.get(raw_cell, raw_cell)
    return model, cell


def build_visibility_mask(trace_indices: Iterable[int], total_traces: int) -> list[bool]:
    selected = set(trace_indices)
    return [index in selected for index in range(total_traces)]


def rename_cardiacfoam_trace(name: str) -> str:
    if ", ΔT=" in name:
        return name.split(", ΔT=")[0] + " cardiacFoam"
    return name
