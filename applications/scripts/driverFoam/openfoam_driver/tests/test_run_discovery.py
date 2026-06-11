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
#     test_run_discovery
#
# Description
#     Tests run discovery logic and specification contracts.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

"""Tests for run discovery (autonomous-agent completion plan Task 5).

`list_runs(root)` walks a directory and returns one parsed manifest dict
per `run_manifest.json` found. Used by agents to inspect past runs
without manual filesystem traversal.
"""
from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path


class TestListRunsModule(unittest.TestCase):
    def test_module_exports_list_runs(self) -> None:
        from openfoam_driver.core.runtime.run_discovery import list_runs
        self.assertTrue(callable(list_runs))


class TestListRunsBehaviour(unittest.TestCase):
    def test_empty_root_returns_empty_list(self) -> None:
        from openfoam_driver.core.runtime.run_discovery import list_runs
        with tempfile.TemporaryDirectory() as temp:
            self.assertEqual(list(list_runs(Path(temp))), [])

    def test_returns_one_entry_per_manifest_found(self) -> None:
        from openfoam_driver.core.runtime.run_discovery import list_runs
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            # Two synthesized runs under different sub-paths.
            for run_id, sub in [("r1", "alpha/output"), ("r2", "beta/output")]:
                (root / sub).mkdir(parents=True)
                (root / sub / "run_manifest.json").write_text(
                    json.dumps({"run_id": run_id, "status": "completed"})
                )
            # And one decoy non-manifest file that must be ignored.
            (root / "alpha" / "other.json").write_text('{"random": true}')

            runs = list(list_runs(root))
            self.assertEqual(len(runs), 2)
            run_ids = sorted(r["run_id"] for r in runs)
            self.assertEqual(run_ids, ["r1", "r2"])

    def test_returns_manifest_path_alongside_payload(self) -> None:
        """Each yielded entry carries the absolute path to the manifest
        file so agents can locate sibling sidecars."""
        from openfoam_driver.core.runtime.run_discovery import list_runs
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            (root / "out").mkdir()
            manifest = root / "out" / "run_manifest.json"
            manifest.write_text(json.dumps({"run_id": "r"}))
            runs = list(list_runs(root))
            self.assertEqual(len(runs), 1)
            self.assertEqual(runs[0]["_manifest_path"], str(manifest))

    def test_malformed_manifest_is_skipped_silently(self) -> None:
        """Malformed JSON should not crash the iterator; agents should
        still see the well-formed manifests."""
        from openfoam_driver.core.runtime.run_discovery import list_runs
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            (root / "broken").mkdir()
            (root / "broken" / "run_manifest.json").write_text("{ not json")
            (root / "good").mkdir()
            (root / "good" / "run_manifest.json").write_text('{"run_id": "ok"}')
            runs = list(list_runs(root))
            ids = [r["run_id"] for r in runs]
            self.assertEqual(ids, ["ok"])


if __name__ == "__main__":
    unittest.main()
