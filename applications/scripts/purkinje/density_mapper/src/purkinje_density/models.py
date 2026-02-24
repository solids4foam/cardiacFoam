from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Mapping, Protocol

from .types import FeatureVector


def _sigmoid(x: float) -> float:
    if x >= 0:
        z = math.exp(-x)
        return 1.0 / (1.0 + z)
    z = math.exp(x)
    return z / (1.0 + z)


class ProbabilityModel(Protocol):
    def predict(self, features: FeatureVector) -> float:
        ...


@dataclass(frozen=True)
class LinearLogisticModel:
    """
    p = sigmoid(sum_i(w_i * x_i) + bias)
    """

    weights: Mapping[str, float]
    bias: float = 0.0

    def predict(self, features: FeatureVector) -> float:
        score = float(self.bias)
        for key, weight in self.weights.items():
            score += float(weight) * float(features.get(key, 0.0))
        return _sigmoid(score)


@dataclass(frozen=True)
class RangeRule:
    feature: str
    min_value: float | None = None
    max_value: float | None = None
    weight: float = 1.0
    missing_penalty: float = 0.5


@dataclass(frozen=True)
class RangeRuleModel:
    """
    Rule scorer:
    - Adds `weight` when a feature is inside [min, max].
    - Subtracts `weight` when outside.
    - Subtracts `missing_penalty` when feature is missing.
    - Final probability is sigmoid(total_score + bias).
    """

    rules: tuple[RangeRule, ...]
    bias: float = 0.0

    def predict(self, features: FeatureVector) -> float:
        score = float(self.bias)
        for rule in self.rules:
            value = features.get(rule.feature)
            if value is None:
                score -= float(rule.missing_penalty)
                continue

            is_inside = True
            if rule.min_value is not None and value < rule.min_value:
                is_inside = False
            if rule.max_value is not None and value > rule.max_value:
                is_inside = False

            score += float(rule.weight) if is_inside else -float(rule.weight)

        return _sigmoid(score)
