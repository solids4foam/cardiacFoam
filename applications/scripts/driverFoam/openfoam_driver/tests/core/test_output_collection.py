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
#     test_output_collection
#
# Description
#     Tests the generic snapshot/diff output-collection mechanism used by
#     entry-mode sweeps to archive each case's postProcessing/ output, one
#     subfolder per case_id.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

from __future__ import annotations

import tempfile
import time
import unittest
from pathlib import Path

from openfoam_driver.core.runtime.output_collection import (
    OutputCollisionError,
    collect_new_outputs,
    snapshot_postprocessing,
)


class TestSnapshotAndCollect(unittest.TestCase):
    def test_new_file_after_snapshot_is_collected_under_its_case_id(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            case_root = Path(tmp) / "case"
            (case_root / "postProcessing").mkdir(parents=True)
            archive_dir = Path(tmp) / "sweepCases"

            before = snapshot_postprocessing(case_root)
            (case_root / "postProcessing" / "1D_10_cells_implicit.dat").write_text("case A data")

            collected = collect_new_outputs(case_root, before, archive_dir, case_id="10")

            self.assertEqual(
                [p.name for p in collected], ["1D_10_cells_implicit.dat"]
            )
            self.assertEqual(
                (archive_dir / "10" / "1D_10_cells_implicit.dat").read_text(), "case A data"
            )

    def test_preexisting_unchanged_file_is_not_recollected(self) -> None:
        """A leftover file from an earlier, unrelated run that this case's
        solve doesn't touch is not "this case's output" and must not be
        copied again -- only new-or-changed-since-snapshot counts."""
        with tempfile.TemporaryDirectory() as tmp:
            case_root = Path(tmp) / "case"
            (case_root / "postProcessing").mkdir(parents=True)
            (case_root / "postProcessing" / "stale.dat").write_text("old leftover")
            archive_dir = Path(tmp) / "sweepCases"

            before = snapshot_postprocessing(case_root)
            # Nothing changes this "run" -- stale.dat is untouched.
            collected = collect_new_outputs(case_root, before, archive_dir, case_id="10")

            self.assertEqual(collected, [])
            self.assertFalse((archive_dir / "10" / "stale.dat").exists())

    def test_two_sequential_cases_get_their_own_subfolder(self) -> None:
        """The real scenario: case_root is shared and reused sequentially
        (apply_case mutates it in place between cases). Each case's
        postProcessing/ output lands under its own case_id subfolder --
        organizing by ran case, not a flat merge, so which case produced what
        is always unambiguous regardless of whether its filenames happen to
        be case-parameter-qualified."""
        with tempfile.TemporaryDirectory() as tmp:
            case_root = Path(tmp) / "case"
            (case_root / "postProcessing").mkdir(parents=True)
            archive_dir = Path(tmp) / "sweepCases"

            before_1 = snapshot_postprocessing(case_root)
            (case_root / "postProcessing" / "1D_10_cells_implicit.dat").write_text("N=10 result")
            collect_new_outputs(case_root, before_1, archive_dir, case_id="10")

            before_2 = snapshot_postprocessing(case_root)
            (case_root / "postProcessing" / "1D_20_cells_implicit.dat").write_text("N=20 result")
            collect_new_outputs(case_root, before_2, archive_dir, case_id="20")

            self.assertEqual(
                (archive_dir / "10" / "1D_10_cells_implicit.dat").read_text(), "N=10 result"
            )
            self.assertEqual(
                (archive_dir / "20" / "1D_20_cells_implicit.dat").read_text(), "N=20 result"
            )
            # Per-case subfolders make a cross-case name collision structurally
            # impossible even for a non-case-qualified functionObject name.
            self.assertEqual(sorted(p.name for p in archive_dir.iterdir()), ["10", "20"])

    def test_rerunning_the_same_case_id_with_different_content_raises(self) -> None:
        """A retry that produces genuinely different content at the same
        relative path under the SAME case_id's own subfolder is suspicious --
        e.g. a concurrent/partial retry -- and must not be silently
        overwritten."""
        with tempfile.TemporaryDirectory() as tmp:
            case_root = Path(tmp) / "case"
            (case_root / "postProcessing").mkdir(parents=True)
            archive_dir = Path(tmp) / "sweepCases"

            before_1 = snapshot_postprocessing(case_root)
            (case_root / "postProcessing" / "result.dat").write_text("first attempt")
            collect_new_outputs(case_root, before_1, archive_dir, case_id="10")

            time.sleep(0.01)
            before_2 = snapshot_postprocessing(case_root)
            (case_root / "postProcessing" / "result.dat").write_text("second attempt")
            with self.assertRaises(OutputCollisionError):
                collect_new_outputs(case_root, before_2, archive_dir, case_id="10")

    def test_identical_rerun_content_is_not_a_collision(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            case_root = Path(tmp) / "case"
            (case_root / "postProcessing").mkdir(parents=True)
            archive_dir = Path(tmp) / "sweepCases"

            before_1 = snapshot_postprocessing(case_root)
            (case_root / "postProcessing" / "result.dat").write_text("same content")
            collect_new_outputs(case_root, before_1, archive_dir, case_id="10")

            time.sleep(0.01)
            before_2 = snapshot_postprocessing(case_root)
            (case_root / "postProcessing" / "result.dat").write_text("same content")
            collected = collect_new_outputs(case_root, before_2, archive_dir, case_id="10")
            self.assertEqual(collected, [])

    def test_nested_functionobject_directories_preserve_relative_structure(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            case_root = Path(tmp) / "case"
            (case_root / "postProcessing" / "Niedererpoints" / "10").mkdir(parents=True)
            archive_dir = Path(tmp) / "sweepCases"

            before = snapshot_postprocessing(case_root)
            (case_root / "postProcessing" / "Niedererpoints" / "10" / "data.xy").write_text("xy")
            collected = collect_new_outputs(case_root, before, archive_dir, case_id="dx0.5")

            self.assertEqual(len(collected), 1)
            self.assertTrue(
                (archive_dir / "dx0.5" / "Niedererpoints" / "10" / "data.xy").exists()
            )

    def test_missing_postprocessing_directory_snapshots_as_empty(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            case_root = Path(tmp) / "case"
            case_root.mkdir()
            before = snapshot_postprocessing(case_root)
            self.assertEqual(before, {})


if __name__ == "__main__":
    unittest.main()
