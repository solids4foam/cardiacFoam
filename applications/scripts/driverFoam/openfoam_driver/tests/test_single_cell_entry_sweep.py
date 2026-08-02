#----------------------------------------------------------------------------#
# License
#     This file is part of cardiacFoam.
#
#     cardiacFoam is free software: you can redistribute it and/or modify it
#     under the terms of the GNU General Public License as published by the
#     Free Software Foundation, either version 3 of the License, or (at your
#     option) any later version.
#
#     cardiacFoam is distributed in the hope that it will be useful, but
#     WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#     General Public License for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with cardiacFoam.  If not, see <http://www.gnu.org/licenses/>.
#
# Module
#     test_single_cell_entry_sweep
#
# Description
#     Tests the ionic_model/tissue single-case override make_spec() gains so
#     an entry-mode sweep can select exactly one (model, tissue) case per
#     resolved axis combination -- the invariant _materialize_entry_case
#     requires (build_cases() must collapse to exactly 1 case).
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

from __future__ import annotations

import unittest
from pathlib import Path

from openfoam_driver.plugins.cardiacfoam.tutorials.single_cell import make_spec


def _repo_root() -> Path:
    here = Path(__file__).resolve()
    for parent in here.parents:
        if (parent / "tutorials").exists() and (parent / "applications").exists():
            return parent
    raise AssertionError("Could not locate repository root from test path")


class TestSingleCellIonicModelTissueOverride(unittest.TestCase):
    def test_single_pair_override_collapses_to_exactly_one_case(self) -> None:
        spec = make_spec(
            tutorials_root=_repo_root() / "tutorials",
            ionic_model="TNNP",
            tissue="epicardialCells",
        )
        cases = spec.build_cases()
        self.assertEqual(len(cases), 1)
        self.assertEqual(cases[0].case_id, "TNNP_epicardialCells")
        self.assertEqual(cases[0].params, {"ionicModel": "TNNP", "tissue": "epicardialCells"})

    def test_default_ionic_models_still_produces_the_full_catalog_sweep(self) -> None:
        spec = make_spec(tutorials_root=_repo_root() / "tutorials")
        cases = spec.build_cases()
        self.assertGreater(len(cases), 1)

    def test_ionic_model_without_tissue_is_rejected(self) -> None:
        with self.assertRaises(ValueError):
            make_spec(tutorials_root=_repo_root() / "tutorials", ionic_model="TNNP")

    def test_tissue_without_ionic_model_is_rejected(self) -> None:
        with self.assertRaises(ValueError):
            make_spec(tutorials_root=_repo_root() / "tutorials", tissue="epicardialCells")


if __name__ == "__main__":
    unittest.main()
