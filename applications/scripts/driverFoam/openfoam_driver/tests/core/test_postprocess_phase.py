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
#     test_postprocess_phase
#
# Description
#     Tests the post-DAG hand-off stub and its JSON serialization.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

"""Tests for the post-DAG hand-off placeholder.

`run_postprocess_phase` currently does no real work -- it proves the hand-off
point exists and returns a stable, serializable shape. These tests pin that
shape so the eventual real implementation has a contract to keep.
"""
from __future__ import annotations

import unittest
from pathlib import Path

from openfoam_driver.core.runtime.postprocess_phase import (
    PostprocessOutcome,
    run_postprocess_phase,
)


class PostprocessPhaseTests(unittest.TestCase):
    def test_returns_stub_outcome(self) -> None:
        outcome = run_postprocess_phase(entry="singleCell", output_dir=Path("/tmp/out"))
        self.assertIsInstance(outcome, PostprocessOutcome)
        self.assertEqual(outcome.status, "stub")
        self.assertIn("singleCell", outcome.message)
        self.assertIn("/tmp/out", outcome.message)

    def test_to_json_round_trips_status_and_message(self) -> None:
        outcome = run_postprocess_phase(entry="niederer2012", output_dir=Path("/tmp/x"))
        payload = outcome.to_json()
        self.assertEqual(payload, {"status": outcome.status, "message": outcome.message})


if __name__ == "__main__":
    unittest.main()
