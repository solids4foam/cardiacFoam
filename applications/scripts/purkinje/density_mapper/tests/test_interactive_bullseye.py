import unittest

from purkinje_density.interactive_bullseye import (
    DEFAULT_AHA17_REFERENCE_CFG,
    build_bullseye_segment_geometry,
    segment_ids_from_reference,
    validate_segment_percentages,
)


class InteractiveBullseyeTests(unittest.TestCase):
    def test_segment_ids_include_apex_for_default_aha17(self) -> None:
        ids = segment_ids_from_reference(DEFAULT_AHA17_REFERENCE_CFG)
        self.assertEqual(ids[0], 1)
        self.assertEqual(ids[-1], 17)
        self.assertEqual(len(ids), 17)

    def test_geometry_has_expected_count(self) -> None:
        geometry = build_bullseye_segment_geometry(DEFAULT_AHA17_REFERENCE_CFG)
        self.assertEqual(len(geometry), 17)
        self.assertTrue(any(g.is_apex_cap for g in geometry))

    def test_apical_ring_starts_at_45_deg_for_visualization(self) -> None:
        geometry = build_bullseye_segment_geometry(DEFAULT_AHA17_REFERENCE_CFG)
        apical_wedges = [g for g in geometry if g.ring_name == "apical" and not g.is_apex_cap]
        self.assertTrue(apical_wedges)
        first_apical = apical_wedges[0]
        self.assertAlmostEqual(first_apical.theta1_deg, 45.0)

    def test_validate_allows_partial(self) -> None:
        validate_segment_percentages({1: 20.0, 5: 35.0}, DEFAULT_AHA17_REFERENCE_CFG)

    def test_validate_rejects_unknown_segment(self) -> None:
        with self.assertRaises(ValueError):
            validate_segment_percentages({99: 10.0}, DEFAULT_AHA17_REFERENCE_CFG)

    def test_validate_rejects_out_of_range(self) -> None:
        with self.assertRaises(ValueError):
            validate_segment_percentages({1: 120.0}, DEFAULT_AHA17_REFERENCE_CFG)

    def test_validate_respects_custom_max_value(self) -> None:
        with self.assertRaises(ValueError):
            validate_segment_percentages(
                {1: 12.0},
                DEFAULT_AHA17_REFERENCE_CFG,
                max_value=10.0,
            )

    def test_validate_requires_all_when_requested(self) -> None:
        with self.assertRaises(ValueError):
            validate_segment_percentages(
                {1: 20.0, 2: 30.0},
                DEFAULT_AHA17_REFERENCE_CFG,
                allow_partial=False,
            )


if __name__ == "__main__":
    unittest.main()
