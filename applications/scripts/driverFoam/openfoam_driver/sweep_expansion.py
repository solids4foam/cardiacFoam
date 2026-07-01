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
#     sweep_expansion
#
# Description
#     Pure axis expansion for parameter sweeps (cross product / zip).
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

from __future__ import annotations

from dataclasses import dataclass
from itertools import product
from math import prod
from typing import Any


class SweepValidationError(ValueError):
    """Raised for any sweep spec problem caught before materialization runs."""


@dataclass(frozen=True)
class ResolvedCase:
    case_id: str
    resolved_axis_values: dict[str, Any]


def _sweep_block(sweep_spec: dict[str, Any]) -> dict[str, Any]:
    try:
        sweep = sweep_spec["sweep"]
    except KeyError:
        raise SweepValidationError("sweep spec is missing required 'sweep' object") from None
    if not isinstance(sweep, dict):
        raise SweepValidationError("sweep must be a JSON object")
    return sweep


def _independent_axes(sweep: dict[str, Any]) -> dict[str, list[Any]]:
    independent = sweep.get("independent")
    if not isinstance(independent, dict):
        raise SweepValidationError("sweep.independent must be a JSON object")
    for name, values in independent.items():
        if not isinstance(name, str) or not name:
            raise SweepValidationError("every independent axis name must be a non-empty string")
        if not isinstance(values, list):
            raise SweepValidationError(f"sweep.independent.{name} must be a list")
        if not values:
            raise SweepValidationError(f"sweep.independent.{name} is empty")
    return independent


def _mode(sweep: dict[str, Any]) -> str:
    mode = sweep.get("mode")
    if mode not in {"cross_product", "zip"}:
        raise SweepValidationError(
            f"sweep.mode must be 'cross_product' or 'zip', got {mode!r}"
        )
    return str(mode)


def compute_case_count(sweep_spec: dict[str, Any]) -> int:
    sweep = _sweep_block(sweep_spec)
    mode = _mode(sweep)
    independent = _independent_axes(sweep)
    if not independent:
        return 1
    lengths = [len(values) for values in independent.values()]
    if mode == "zip":
        _validate_zip_lengths(independent)
        return lengths[0]
    return prod(lengths)


def _validate_zip_lengths(independent: dict[str, list[Any]]) -> None:
    lengths = {name: len(values) for name, values in independent.items()}
    distinct = set(lengths.values())
    if len(distinct) > 1:
        detail = ", ".join(f"{name}={n}" for name, n in lengths.items())
        raise SweepValidationError(
            f"zip mode requires all independent axes to have equal length, got: {detail}"
        )


def _combinations(mode: str, independent: dict[str, list[Any]]) -> list[dict[str, Any]]:
    names = list(independent.keys())
    if not names:
        return [{}]
    if mode == "zip":
        _validate_zip_lengths(independent)
        rows = zip(*(independent[name] for name in names))
        return [dict(zip(names, row)) for row in rows]
    rows = product(*(independent[name] for name in names))
    return [dict(zip(names, row)) for row in rows]


def expand_sweep(sweep_spec: dict[str, Any]) -> list[ResolvedCase]:
    sweep = _sweep_block(sweep_spec)
    mode = _mode(sweep)
    independent = _independent_axes(sweep)

    combinations = _combinations(mode, independent)

    cases: list[ResolvedCase] = []
    for index, combo in enumerate(combinations, start=1):
        cases.append(ResolvedCase(case_id=f"case_{index:04d}", resolved_axis_values=dict(combo)))
    return cases
