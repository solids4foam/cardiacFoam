from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from typing import Iterable, Mapping


def _canonical_label(value: str) -> str:
    return value.strip().lower()


@dataclass(frozen=True)
class TreeSample:
    """
    One sampled element from a Purkinje tree representation.

    - `tissue`: expected labels like "endocardium" / "myocardium".
    - `is_terminal`: True for terminal/leaf elements.
    - `section_tag`: anatomical region id/name for metric #4.
    - `weight`: contribution weight (defaults to 1.0).
      For metric #1 this can be segment length if available.
    """

    tissue: str
    is_terminal: bool = False
    section_tag: str | None = None
    weight: float = 1.0


@dataclass(frozen=True)
class TreeDensityMetrics:
    """
    Output bundle for the 4 requested quantities.
    """

    any_part_probability_by_tissue: dict[str, float]
    terminal_probability_by_tissue: dict[str, float]
    terminal_density_by_tissue_mm2: dict[str, float | None]
    terminal_density_by_section_mm2: dict[str, float | None]
    terminal_count_by_section: dict[str, float]
    total_any_part_weight: float
    total_terminal_weight: float


@dataclass(frozen=True)
class TreePresenceMetrics:
    """
    Minimal output bundle for early-stage modeling:
    1) Probability of any tree part in each tissue.
    2) Probability of terminal tree parts in each tissue.
    """

    any_part_probability_by_tissue: dict[str, float]
    terminal_probability_by_tissue: dict[str, float]
    total_any_part_weight: float
    total_terminal_weight: float


def _probabilities_from_counter(counter: Counter[str]) -> dict[str, float]:
    total = float(sum(counter.values()))
    if total <= 0.0:
        return {}
    return {k: float(v) / total for k, v in counter.items()}


def _densities(
    counts: Mapping[str, float],
    areas_mm2: Mapping[str, float] | None,
) -> dict[str, float | None]:
    if areas_mm2 is None:
        return {k: None for k in counts}

    out: dict[str, float | None] = {}
    for key, count in counts.items():
        area = areas_mm2.get(key)
        if area is None or area <= 0.0:
            out[key] = None
            continue
        out[key] = float(count) / float(area)
    return out


def compute_tree_presence_metrics(
    samples: Iterable[TreeSample],
    expected_tissues: tuple[str, ...] = ("endocardium", "myocardium"),
) -> TreePresenceMetrics:
    """
    Compute only the two base probabilities:

    1) Probability of any tree part in each tissue.
    2) Probability of terminal tree parts in each tissue.
    """
    samples_list = list(samples)
    if not samples_list:
        raise ValueError("No tree samples were provided.")

    any_part_counter: Counter[str] = Counter()
    terminal_counter: Counter[str] = Counter()

    for sample in samples_list:
        label = _canonical_label(sample.tissue)
        weight = float(sample.weight)
        if weight < 0.0:
            raise ValueError("Sample weight must be non-negative.")

        any_part_counter[label] += weight
        if sample.is_terminal:
            terminal_counter[label] += weight

    # Ensure expected tissues are present in outputs.
    for tissue in expected_tissues:
        any_part_counter[_canonical_label(tissue)] += 0.0
        terminal_counter[_canonical_label(tissue)] += 0.0

    return TreePresenceMetrics(
        any_part_probability_by_tissue=_probabilities_from_counter(any_part_counter),
        terminal_probability_by_tissue=_probabilities_from_counter(terminal_counter),
        total_any_part_weight=float(sum(any_part_counter.values())),
        total_terminal_weight=float(sum(terminal_counter.values())),
    )


def compute_tree_density_metrics(
    samples: Iterable[TreeSample],
    tissue_surface_areas_mm2: Mapping[str, float] | None = None,
    section_surface_areas_mm2: Mapping[str, float] | None = None,
    expected_tissues: tuple[str, ...] = ("endocardium", "myocardium"),
) -> TreeDensityMetrics:
    """
    Compute the 4 requested metrics:

    1) Probability of any tree part in each tissue.
    2) Probability of terminal tree parts in each tissue.
    3) Terminal density (ends/mm^2) in each tissue.
    4) Terminal density (ends/mm^2) per tagged section.
    """
    samples_list = list(samples)
    presence = compute_tree_presence_metrics(
        samples=samples_list,
        expected_tissues=expected_tissues,
    )

    terminal_counter: Counter[str] = Counter()
    terminal_by_section_counter: Counter[str] = Counter()

    for sample in samples_list:
        weight = float(sample.weight)
        label = _canonical_label(sample.tissue)
        if sample.is_terminal:
            terminal_counter[label] += weight
            if sample.section_tag is not None:
                terminal_by_section_counter[str(sample.section_tag)] += weight

    for tissue in expected_tissues:
        terminal_counter[_canonical_label(tissue)] += 0.0

    if tissue_surface_areas_mm2 is not None:
        tissue_area_map = {
            _canonical_label(k): float(v)
            for k, v in tissue_surface_areas_mm2.items()
        }
    else:
        tissue_area_map = None

    section_area_map = (
        {str(k): float(v) for k, v in section_surface_areas_mm2.items()}
        if section_surface_areas_mm2 is not None
        else None
    )

    terminal_density_tissue = _densities(terminal_counter, tissue_area_map)
    terminal_density_section = _densities(terminal_by_section_counter, section_area_map)

    return TreeDensityMetrics(
        any_part_probability_by_tissue=presence.any_part_probability_by_tissue,
        terminal_probability_by_tissue=presence.terminal_probability_by_tissue,
        terminal_density_by_tissue_mm2=terminal_density_tissue,
        terminal_density_by_section_mm2=terminal_density_section,
        terminal_count_by_section=dict(terminal_by_section_counter),
        total_any_part_weight=presence.total_any_part_weight,
        total_terminal_weight=presence.total_terminal_weight,
    )
