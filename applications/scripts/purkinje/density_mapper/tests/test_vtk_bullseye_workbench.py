import unittest

from purkinje_density.vtk_bullseye_workbench import (
    load_division_reference,
    map_percentages_to_segments,
)


class VtkWorkbenchLogicTests(unittest.TestCase):
    def test_load_local_division_reference(self) -> None:
        division = load_division_reference(config_source="local", reference_name="aha17")
        self.assertEqual(division.reference_name, "aha17")
        self.assertIn("ring_layout", division.reference_cfg)

    def test_map_percentages_to_segments(self) -> None:
        seg = [1, 1, 2, 3, 0, -1]
        mapping = {1: 20.0, 2: 60.0}
        out = map_percentages_to_segments(seg, mapping, default_value=0.0)
        self.assertEqual([float(x) for x in out], [20.0, 20.0, 60.0, 0.0, 0.0, 0.0])

    def test_map_percentages_without_ignoring_non_positive(self) -> None:
        seg = [0, -1, 2]
        mapping = {0: 5.0, -1: 3.0, 2: 10.0}
        out = map_percentages_to_segments(
            seg,
            mapping,
            default_value=0.0,
            ignore_non_positive=False,
        )
        self.assertEqual([float(x) for x in out], [5.0, 3.0, 10.0])


if __name__ == "__main__":
    unittest.main()
