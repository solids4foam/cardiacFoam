from __future__ import annotations

from typing import Any, Mapping

from .engine import ProbabilityEngine
from .features import FeatureExtractor, default_feature_extractor
from .models import LinearLogisticModel, RangeRule, RangeRuleModel


def build_engine_from_dict(
    config: Mapping[str, Any],
    extractor: FeatureExtractor | None = None,
) -> ProbabilityEngine:
    """
    Build engine from dict configuration.

    Supported model types:
    - `linear_logistic`
    - `range_rule`
    """
    engine = ProbabilityEngine(extractor=extractor or default_feature_extractor)
    models = config.get("models", {})

    for tag_name, model_cfg in models.items():
        model_type = str(model_cfg.get("type", "linear_logistic")).strip().lower()
        bias = float(model_cfg.get("bias", 0.0))

        if model_type == "linear_logistic":
            weights = {
                str(k): float(v)
                for k, v in dict(model_cfg.get("weights", {})).items()
            }
            model = LinearLogisticModel(weights=weights, bias=bias)
        elif model_type == "range_rule":
            rules_in = model_cfg.get("rules", [])
            rules = []
            for r in rules_in:
                rules.append(
                    RangeRule(
                        feature=str(r["feature"]),
                        min_value=float(r["min"]) if "min" in r and r["min"] is not None else None,
                        max_value=float(r["max"]) if "max" in r and r["max"] is not None else None,
                        weight=float(r.get("weight", 1.0)),
                        missing_penalty=float(r.get("missing_penalty", 0.5)),
                    )
                )
            model = RangeRuleModel(rules=tuple(rules), bias=bias)
        else:
            raise ValueError(
                f"Unsupported model type '{model_type}' for tag '{tag_name}'."
            )

        engine.registry.register(str(tag_name), model)

    return engine
