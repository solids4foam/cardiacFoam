import unittest

from purkinje_density.tree_metrics import (
    TreeSample,
    compute_tree_density_metrics,
)


class TreeMetricsTests(unittest.TestCase):
    def test_requested_four_metrics(self) -> None:
        samples = [
            TreeSample(tissue="endocardium", is_terminal=False, section_tag="s1"),
            TreeSample(tissue="myocardium", is_terminal=False, section_tag="s2"),
            TreeSample(tissue="endocardium", is_terminal=False, section_tag="s1"),
            TreeSample(tissue="endocardium", is_terminal=True, section_tag="s1"),
            TreeSample(tissue="endocardium", is_terminal=True, section_tag="s1"),
            TreeSample(tissue="myocardium", is_terminal=True, section_tag="s2"),
        ]

        tissue_areas = {"endocardium": 100.0, "myocardium": 200.0}
        section_areas = {"s1": 50.0, "s2": 80.0}

        m = compute_tree_density_metrics(
            samples=samples,
            tissue_surface_areas_mm2=tissue_areas,
            section_surface_areas_mm2=section_areas,
        )

        # 1) Any-part probability
        self.assertAlmostEqual(m.any_part_probability_by_tissue["endocardium"], 4.0 / 6.0)
        self.assertAlmostEqual(m.any_part_probability_by_tissue["myocardium"], 2.0 / 6.0)

        # 2) Terminal probability
        self.assertAlmostEqual(m.terminal_probability_by_tissue["endocardium"], 2.0 / 3.0)
        self.assertAlmostEqual(m.terminal_probability_by_tissue["myocardium"], 1.0 / 3.0)

        # 3) Terminal density by tissue (ends/mm^2)
        self.assertAlmostEqual(m.terminal_density_by_tissue_mm2["endocardium"], 2.0 / 100.0)
        self.assertAlmostEqual(m.terminal_density_by_tissue_mm2["myocardium"], 1.0 / 200.0)

        # 4) Terminal density by section (ends/mm^2)
        self.assertAlmostEqual(m.terminal_density_by_section_mm2["s1"], 2.0 / 50.0)
        self.assertAlmostEqual(m.terminal_density_by_section_mm2["s2"], 1.0 / 80.0)

    def test_missing_area_returns_none_density(self) -> None:
        samples = [
            TreeSample(tissue="endocardium", is_terminal=True, section_tag="s1"),
            TreeSample(tissue="myocardium", is_terminal=True, section_tag="s2"),
        ]
        m = compute_tree_density_metrics(samples=samples)

        self.assertIsNone(m.terminal_density_by_tissue_mm2["endocardium"])
        self.assertIsNone(m.terminal_density_by_tissue_mm2["myocardium"])
        self.assertIsNone(m.terminal_density_by_section_mm2["s1"])
        self.assertIsNone(m.terminal_density_by_section_mm2["s2"])


if __name__ == "__main__":
    unittest.main()
