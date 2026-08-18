#----------------------------------------------------------------------------#
# License
#     This file is part of cardiacFoam.
#
#     cardiacFoam is free software: you can redistribute it and/or modify it
#     under the terms of the GNU General Public License as published by the
#     Free Software Foundation, either version 3 of the License, or (at your
#     option) any later version.
#
#     cardiacFoam is distributed in the hope that it will be useful, but
#     WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#     General Public License for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with cardiacFoam.  If not, see <http://www.gnu.org/licenses/>.
#
# Module
#     plotting_common
#
# Description
#     Provides shared utilities for trace plotting and visualization.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

from __future__ import annotations

import re
from collections.abc import Iterable, Mapping

import plotly.colors as plotly_colors





_DX_PATTERN = re.compile(r"DX(\d+)")
_DT_PATTERN = re.compile(r"DT(\d+)")


def extract_dx_dt(filename: str) -> tuple[float, float]:
    """Extract DX and DT values from filenames like ...DX5...DT005..."""
    from pathlib import Path
    basename = Path(filename).name
    dx_match = _DX_PATTERN.search(basename)
    dt_match = _DT_PATTERN.search(basename)
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


def parse_two_part_stem(
    filename: str,
    *,
    first_map: Mapping[str, str],
    second_map: Mapping[str, str],
) -> tuple[str, str]:
    """Split a `<first>_<second>_...` filename stem and map both tokens.

    Tokens absent from their map pass through unchanged, so callers only supply
    the renames they care about.
    """
    from pathlib import Path
    file_base = Path(filename).stem
    parts = file_base.split("_")

    raw_first = parts[0] if parts else "UnknownFirst"
    raw_second = parts[1] if len(parts) > 1 else "UnknownSecond"

    return first_map.get(raw_first, raw_first), second_map.get(raw_second, raw_second)


def build_visibility_mask(trace_indices: Iterable[int], total_traces: int) -> list[bool]:
    selected = set(trace_indices)
    return [index in selected for index in range(total_traces)]
