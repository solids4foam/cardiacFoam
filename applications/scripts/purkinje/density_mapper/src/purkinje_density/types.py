from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Mapping


FeatureVector = dict[str, float]


@dataclass(frozen=True)
class LocationSample:
    """
    Minimal location container for probability inference.

    Notes:
    - `scalar_fields` is where most normalized geometric descriptors should go
      (for example: `uvc_longitudinal`, `transmural`, `theta`).
    - `chamber` and `wall` become binary convenience features in the default
      extractor.
    - `tags` can carry numeric values from your anatomical tagging pipeline.
    """

    xyz: tuple[float, float, float] | None = None
    chamber: str | None = None
    wall: str | None = None
    zone: str | None = None
    scalar_fields: Mapping[str, float] = field(default_factory=dict)
    categorical_fields: Mapping[str, str] = field(default_factory=dict)
    tags: Mapping[str, float] = field(default_factory=dict)


@dataclass(frozen=True)
class ProbabilityContext:
    """
    Extra runtime information passed to feature extraction/model evaluation.

    Keep this generic so the engine remains standalone.
    """

    geometry: Mapping[str, Any] = field(default_factory=dict)
    metadata: Mapping[str, Any] = field(default_factory=dict)
