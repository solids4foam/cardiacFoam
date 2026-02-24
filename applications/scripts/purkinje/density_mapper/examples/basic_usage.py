from purkinje_density import LocationSample
from purkinje_density.defaults import build_engine_from_dict


def main() -> None:
    config = {
        "models": {
            "lv_apical_septal": {
                "type": "linear_logistic",
                "weights": {
                    "uvc_longitudinal": -3.0,
                    "is_lv": 1.5,
                    "is_septal": 1.2,
                },
                "bias": 1.0,
            },
            "rv_mid_freewall": {
                "type": "range_rule",
                "bias": -0.2,
                "rules": [
                    {"feature": "uvc_longitudinal", "min": 0.3, "max": 0.7, "weight": 1.0},
                    {"feature": "is_rv", "min": 1.0, "max": 1.0, "weight": 1.3},
                    {"feature": "is_freewall", "min": 1.0, "max": 1.0, "weight": 1.1},
                ],
            },
        }
    }

    engine = build_engine_from_dict(config)

    sample = LocationSample(
        chamber="LV",
        wall="septal",
        scalar_fields={"uvc_longitudinal": 0.15, "transmural": 0.2},
    )

    p_one = engine.probability("lv_apical_septal", sample)
    p_all = engine.probabilities(sample)

    print(f"lv_apical_septal = {p_one:.4f}")
    print("all tags:", {k: round(v, 4) for k, v in p_all.items()})


if __name__ == "__main__":
    main()
