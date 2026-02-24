import json
import tempfile
import unittest
from pathlib import Path

from purkinje_density.apply_solution_json import load_solution_json


class ApplySolutionJsonTests(unittest.TestCase):
    def test_load_solution_json_parses_segment_maps(self) -> None:
        payload = {
            "purkinje": {"1": 0.2, "2": 0.5},
            "pmj": {"1": 0.1, "2": 0.9},
            "thickness_bundles": {"1": 0.3, "2": 0.7},
        }
        with tempfile.TemporaryDirectory() as td:
            path = Path(td) / "solution.json"
            path.write_text(json.dumps(payload), encoding="utf-8")
            purkinje, pmj, thickness_bundles = load_solution_json(path)

        self.assertEqual(purkinje, {1: 0.2, 2: 0.5})
        self.assertEqual(pmj, {1: 0.1, 2: 0.9})
        self.assertEqual(thickness_bundles, {1: 0.3, 2: 0.7})

    def test_load_solution_json_accepts_missing_tag_map(self) -> None:
        payload = {"purkinje": {"1": 0.25}}
        with tempfile.TemporaryDirectory() as td:
            path = Path(td) / "solution.json"
            path.write_text(json.dumps(payload), encoding="utf-8")
            purkinje, pmj, thickness_bundles = load_solution_json(path)

        self.assertEqual(purkinje, {1: 0.25})
        self.assertEqual(pmj, {})
        self.assertEqual(thickness_bundles, {})


if __name__ == "__main__":
    unittest.main()
