from __future__ import annotations

from typing import Callable

from .types import FeatureVector, LocationSample, ProbabilityContext

FeatureExtractor = Callable[[LocationSample, ProbabilityContext | None], FeatureVector]


def _bin_flag(value: str | None, positive_values: set[str]) -> float:
    if value is None:
        return 0.0
    return 1.0 if value.strip().lower() in positive_values else 0.0


def default_feature_extractor(
    sample: LocationSample,
    context: ProbabilityContext | None = None,
) -> FeatureVector:
    """
    Baseline extractor:
    - copies numeric scalar fields
    - copies numeric tags
    - adds binary anatomy flags from `chamber` and `wall`
    """
    features: FeatureVector = {}

    for name, value in sample.scalar_fields.items():
        try:
            features[name] = float(value)
        except (TypeError, ValueError):
            continue

    for name, value in sample.tags.items():
        try:
            features[name] = float(value)
        except (TypeError, ValueError):
            continue

    features["is_lv"] = _bin_flag(sample.chamber, {"lv", "left", "left_ventricle"})
    features["is_rv"] = _bin_flag(sample.chamber, {"rv", "right", "right_ventricle"})
    features["is_septal"] = _bin_flag(sample.wall, {"septal", "septum"})
    features["is_freewall"] = _bin_flag(sample.wall, {"freewall", "free_wall", "lateral"})

    if context is not None and context.metadata:
        # Optional pass-through numeric metadata for global controls/calibration.
        for name, value in context.metadata.items():
            if isinstance(value, (int, float)):
                features[f"meta_{name}"] = float(value)

    return features


def compose_extractors(*extractors: FeatureExtractor) -> FeatureExtractor:
    """
    Merge features from multiple extractors (later extractors override earlier keys).
    """
    if not extractors:
        return default_feature_extractor

    def _composed(
        sample: LocationSample,
        context: ProbabilityContext | None = None,
    ) -> FeatureVector:
        merged: FeatureVector = {}
        for extractor in extractors:
            merged.update(extractor(sample, context))
        return merged

    return _composed
