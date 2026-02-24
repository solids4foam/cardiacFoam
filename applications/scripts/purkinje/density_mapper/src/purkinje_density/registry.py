from __future__ import annotations

from dataclasses import dataclass, field

from .models import ProbabilityModel
from .types import FeatureVector


@dataclass
class TagModelRegistry:
    """
    Maps anatomical tag names to probability models.
    """

    _models: dict[str, ProbabilityModel] = field(default_factory=dict)
    _fallback: ProbabilityModel | None = None

    def register(self, tag: str, model: ProbabilityModel) -> None:
        self._models[tag] = model

    def set_fallback(self, model: ProbabilityModel | None) -> None:
        self._fallback = model

    def has(self, tag: str) -> bool:
        return tag in self._models

    def predict(self, tag: str, features: FeatureVector) -> float:
        model = self._models.get(tag, self._fallback)
        if model is None:
            raise KeyError(f"No probability model registered for tag '{tag}'.")
        return model.predict(features)

    def predict_all(self, features: FeatureVector) -> dict[str, float]:
        return {tag: model.predict(features) for tag, model in self._models.items()}

    def tags(self) -> tuple[str, ...]:
        return tuple(self._models.keys())
