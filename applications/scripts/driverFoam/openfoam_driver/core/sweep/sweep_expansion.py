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
from typing import Any, Callable


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


def _validate_dependent_declarations(
    dependent: list[dict[str, Any]],
    independent_names: set[str],
) -> None:
    known_names = set(independent_names)
    for entry in dependent:
        if not isinstance(entry, dict):
            raise SweepValidationError("each dependent entry must be a JSON object")
        name = entry["name"]
        if not isinstance(name, str) or not name:
            raise SweepValidationError("dependent entry names must be non-empty strings")
        if "derive" not in entry:
            raise SweepValidationError(f"dependent entry '{name}' is missing required 'derive'")
        if "of" not in entry:
            raise SweepValidationError(f"dependent entry '{name}' is missing required 'of'")
        of = entry["of"] if isinstance(entry["of"], list) else [entry["of"]]
        for ref in of:
            if ref not in known_names:
                raise SweepValidationError(
                    f"dependent entry '{name}' references '{ref}', which is not "
                    "an independent variable or an earlier dependent entry "
                    "(forward references are not allowed)"
                )
        if name in known_names:
            raise SweepValidationError(
                f"dependent entry '{name}' collision: an independent variable "
                "or earlier dependent entry already has this name"
            )
        known_names.add(name)


def _validate_path_safe_case_id(label: str) -> str:
    if (
        not label
        or label in {".", ".."}
        or label.strip() != label
        or "/" in label
        or "\x00" in label
    ):
        raise SweepValidationError(
            f"caseId {label!r} is not path-safe and cannot be used as a case directory"
        )
    return label


def expand_sweep(
    sweep_spec: dict[str, Any],
    *,
    get_derivation: Callable[[str], Callable[[dict[str, Any]], dict[str, Any]]] | None = None,
) -> list[ResolvedCase]:
    sweep = _sweep_block(sweep_spec)
    mode = _mode(sweep)
    independent = _independent_axes(sweep)
    dependent: list[dict[str, Any]] = sweep.get("dependent", [])
    if dependent and get_derivation is None:
        raise SweepValidationError(
            "expand_sweep requires a get_derivation lookup function when dependent entries are present"
        )

    _validate_dependent_declarations(dependent, set(independent.keys()))

    combinations = _combinations(mode, independent)

    resolved_values: list[dict[str, Any]] = []
    for combo in combinations:
        values: dict[str, Any] = dict(combo)
        for entry in dependent:
            of = entry["of"] if isinstance(entry["of"], list) else [entry["of"]]
            derivation_inputs = {ref: values[ref] for ref in of}
            fn = get_derivation(entry["derive"])
            result = fn(derivation_inputs)
            for key in result:
                if key in values:
                    raise SweepValidationError(
                        f"derivation '{entry['derive']}' returned key '{key}', which "
                        "collides with an existing value for this case"
                    )
            values.update(result)
        resolved_values.append(values)

    labels: list[str] = []
    for index, values in enumerate(resolved_values, start=1):
        label = str(values["caseId"]) if "caseId" in values else f"case_{index:04d}"
        labels.append(_validate_path_safe_case_id(label))

    if len(set(labels)) != len(labels):
        seen: dict[str, int] = {}
        duplicates = set()
        for label in labels:
            seen[label] = seen.get(label, 0) + 1
            if seen[label] > 1:
                duplicates.add(label)
        raise SweepValidationError(
            f"caseId values are not unique across the sweep (duplicates: "
            f"{sorted(duplicates)}); the caseId derivation must reference enough "
            "axes to disambiguate every combination"
        )

    return [
        ResolvedCase(case_id=label, resolved_axis_values=dict(values))
        for label, values in zip(labels, resolved_values)
    ]


DEFAULT_MAX_CASES = 200


def check_case_count_cap(sweep_spec: dict[str, Any], *, max_cases: int = DEFAULT_MAX_CASES) -> None:
    """Raise before any expansion/materialization if N exceeds the cap.

    Both sweep-plan and sweep-run must call this before expanding or
    materializing even one case, since materialization writes real files to
    disk (see design doc's "Per-case filesystem layout").
    """
    n = compute_case_count(sweep_spec)
    if n > max_cases:
        raise SweepValidationError(
            f"sweep expands to {n} cases, exceeding max_cases={max_cases}. "
            "Pass an explicit --max-cases to proceed."
        )
