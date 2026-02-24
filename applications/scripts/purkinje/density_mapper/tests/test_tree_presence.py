import unittest

from purkinje_density.tree_metrics import (
    TreeSample,
    compute_tree_presence_metrics,
)


class TreePresenceTests(unittest.TestCase):
    def test_any_and_terminal_probabilities(self) -> None:
        samples = [
            TreeSample(tissue="endocardium", is_terminal=False),
            TreeSample(tissue="myocardium", is_terminal=False),
            TreeSample(tissue="endocardium", is_terminal=False),
            TreeSample(tissue="endocardium", is_terminal=True),
            TreeSample(tissue="myocardium", is_terminal=True),
            TreeSample(tissue="endocardium", is_terminal=True),
        ]

        m = compute_tree_presence_metrics(samples)
        self.assertAlmostEqual(m.any_part_probability_by_tissue["endocardium"], 4.0 / 6.0)
        self.assertAlmostEqual(m.any_part_probability_by_tissue["myocardium"], 2.0 / 6.0)
        self.assertAlmostEqual(m.terminal_probability_by_tissue["endocardium"], 2.0 / 3.0)
        self.assertAlmostEqual(m.terminal_probability_by_tissue["myocardium"], 1.0 / 3.0)

    def test_negative_weight_raises(self) -> None:
        with self.assertRaises(ValueError):
            compute_tree_presence_metrics(
                [TreeSample(tissue="endocardium", is_terminal=True, weight=-1.0)]
            )


if __name__ == "__main__":
    unittest.main()
