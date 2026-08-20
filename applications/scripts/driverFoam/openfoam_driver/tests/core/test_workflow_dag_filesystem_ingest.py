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
#     test_workflow_dag_filesystem_ingest
#
# Description
#     Tests workflow dag filesystem ingest logic and specification contracts.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

"""Tests for filesystem case workflow ownership.

Plain case folders are owned by their on-disk Allrun. Registry discovery may
still find a non-runnable cardiac-marked folder, but it must not invent a
workflow_dag unless Allrun exists.
"""
from __future__ import annotations

import unittest
from pathlib import Path
import tempfile

from openfoam_driver.core.runtime.registry import load_tutorial_spec, resolve_entry


class TestFilesystemCaseWorkflowOwnership(unittest.TestCase):
    """Filesystem case folders own their run definition through Allrun."""

    def _write_case_files(self, case_root: Path) -> None:
        (case_root / "constant").mkdir(parents=True, exist_ok=True)
        (case_root / "constant" / "electroProperties").write_text(
            "myocardiumSolver singleCellSolver;\n"
        )
        (case_root / "constant" / "physicsProperties").write_text(
            "type electroModel;\n"
        )

    def test_filesystem_case_with_allrun_uses_allrun_fallback(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            tutorials_root = Path(temp_dir)
            case_root = tutorials_root / "myCase"
            self._write_case_files(case_root)
            (case_root / "Allrun").write_text("#!/bin/sh\n")

            spec = load_tutorial_spec(
                "myCase",
                overrides={"tutorials_root": tutorials_root},
            )

            dag = spec.metadata.get("workflow_dag")
            self.assertEqual(
                dag,
                {"steps": [{"id": "run", "command": "Allrun", "depends_on": []}]},
            )

    def test_filesystem_case_without_allrun_has_no_workflow_dag(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            tutorials_root = Path(temp_dir)
            case_root = tutorials_root / "bareCase"
            self._write_case_files(case_root)

            spec = load_tutorial_spec(
                "bareCase",
                overrides={"tutorials_root": tutorials_root},
            )

            dag = spec.metadata.get("workflow_dag")
            self.assertIsNone(dag, "workflow_dag must be None when Allrun is absent")

    def test_variant_electro_properties_case_is_discoverable_and_runnable(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            tutorials_root = Path(temp_dir)
            case_root = tutorials_root / "variantCase"
            (case_root / "constant").mkdir(parents=True, exist_ok=True)
            (case_root / "system").mkdir(parents=True, exist_ok=True)
            (case_root / "constant" / "electroProperties.monodomain").write_text(
                "myocardiumSolver monodomainSolver;\n"
            )
            (case_root / "constant" / "physicsProperties").write_text("type electroModel;\n")
            for relpath in ("controlDict", "fvSchemes", "fvSolution"):
                (case_root / "system" / relpath).write_text("\n")

            resolution = resolve_entry(
                "variantCase",
                overrides={"tutorials_root": tutorials_root},
            )

            self.assertEqual(resolution["resolution"], "case_folder")
            self.assertTrue(resolution["is_runnable"])


if __name__ == "__main__":
    unittest.main()
