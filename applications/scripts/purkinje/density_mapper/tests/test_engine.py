import unittest

from purkinje_density.defaults import build_engine_from_dict
from purkinje_density.features import default_feature_extractor
from purkinje_density.models import RangeRule, RangeRuleModel
from purkinje_density.types import LocationSample


class FeatureExtractorTests(unittest.TestCase):
    def test_default_extractor_includes_binary_anatomy_flags(self) -> None:
        sample = LocationSample(
            chamber="LV",
            wall="septal",
            scalar_fields={"uvc_longitudinal": 0.2},
            tags={"existing_tag_score": 0.75},
        )
        feats = default_feature_extractor(sample)

        self.assertAlmostEqual(feats["is_lv"], 1.0)
        self.assertAlmostEqual(feats["is_rv"], 0.0)
        self.assertAlmostEqual(feats["is_septal"], 1.0)
        self.assertAlmostEqual(feats["is_freewall"], 0.0)
        self.assertAlmostEqual(feats["uvc_longitudinal"], 0.2)
        self.assertAlmostEqual(feats["existing_tag_score"], 0.75)


class ModelTests(unittest.TestCase):
    def test_range_rule_model_prefers_inside_range(self) -> None:
        model = RangeRuleModel(
            rules=(RangeRule(feature="uvc_longitudinal", min_value=0.1, max_value=0.3, weight=1.0),),
            bias=0.0,
        )
        p_inside = model.predict({"uvc_longitudinal": 0.2})
        p_outside = model.predict({"uvc_longitudinal": 0.7})

        self.assertGreater(p_inside, p_outside)


class EngineTests(unittest.TestCase):
    def test_build_and_predict(self) -> None:
        config = {
            "models": {
                "lv_apical_septal": {
                    "type": "linear_logistic",
                    "weights": {
                        "uvc_longitudinal": -2.0,
                        "is_lv": 1.0,
                        "is_septal": 1.0,
                    },
                    "bias": 0.5,
                },
                "rv_mid_freewall": {
                    "type": "range_rule",
                    "bias": 0.0,
                    "rules": [
                        {"feature": "uvc_longitudinal", "min": 0.3, "max": 0.7, "weight": 1.0},
                        {"feature": "is_rv", "min": 1.0, "max": 1.0, "weight": 1.0},
                    ],
                },
            }
        }

        engine = build_engine_from_dict(config)

        lv_sample = LocationSample(
            chamber="LV",
            wall="septal",
            scalar_fields={"uvc_longitudinal": 0.15},
        )
        rv_sample = LocationSample(
            chamber="RV",
            wall="freewall",
            scalar_fields={"uvc_longitudinal": 0.5},
        )

        p_lv = engine.probability("lv_apical_septal", lv_sample)
        p_rv_on_rv_model = engine.probability("rv_mid_freewall", rv_sample)
        all_lv = engine.probabilities(lv_sample)

        self.assertGreaterEqual(p_lv, 0.0)
        self.assertLessEqual(p_lv, 1.0)
        self.assertGreaterEqual(p_rv_on_rv_model, 0.0)
        self.assertLessEqual(p_rv_on_rv_model, 1.0)
        self.assertIn("lv_apical_septal", all_lv)
        self.assertIn("rv_mid_freewall", all_lv)


if __name__ == "__main__":
    unittest.main()
