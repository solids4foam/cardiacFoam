from __future__ import annotations

from dataclasses import dataclass, field

from .features import FeatureExtractor, default_feature_extractor
from .registry import TagModelRegistry
from .types import FeatureVector, LocationSample, ProbabilityContext


@dataclass
class ProbabilityEngine:
    """
    Standalone probability API:
    - extracts a shared feature vector once
    - routes to model registry per anatomical tag
    """

    extractor: FeatureExtractor = default_feature_extractor
    registry: TagModelRegistry = field(default_factory=TagModelRegistry)

    def features(
        self,
        sample: LocationSample,
        context: ProbabilityContext | None = None,
    ) -> FeatureVector:
        return self.extractor(sample, context)

    def probability(
        self,
        tag: str,
        sample: LocationSample,
        context: ProbabilityContext | None = None,
    ) -> float:
        feats = self.features(sample, context)
        return self.registry.predict(tag, feats)

    def probabilities(
        self,
        sample: LocationSample,
        context: ProbabilityContext | None = None,
    ) -> dict[str, float]:
        feats = self.features(sample, context)
        return self.registry.predict_all(feats)
