import unittest

from purkinje_density.vtk_bullseye_workbench import aggregate_surface_area_by_segment


class SurfaceAreaByTagTests(unittest.TestCase):
    def test_aggregate_surface_area_by_segment(self) -> None:
        seg = [1, 1, 2, 2, 3]
        areas = [2.0, 3.0, 10.0, 1.0, 4.0]
        out = aggregate_surface_area_by_segment(seg, areas)
        self.assertEqual(out, {1: 5.0, 2: 11.0, 3: 4.0})

    def test_aggregate_surface_area_with_mask(self) -> None:
        seg = [1, 1, 2, 2, 3]
        areas = [2.0, 3.0, 10.0, 1.0, 4.0]
        mask = [True, False, True, False, True]
        out = aggregate_surface_area_by_segment(seg, areas, cell_mask=mask)
        self.assertEqual(out, {1: 2.0, 2: 10.0, 3: 4.0})

    def test_aggregate_ignores_non_positive_segment_ids(self) -> None:
        seg = [-1, 0, 1, 2]
        areas = [2.0, 3.0, 4.0, 5.0]
        out = aggregate_surface_area_by_segment(seg, areas)
        self.assertEqual(out, {1: 4.0, 2: 5.0})


if __name__ == "__main__":
    unittest.main()
