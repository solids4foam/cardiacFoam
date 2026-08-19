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
#     test_execution_context
#
# Description
#     Tests resolve_execution_context(), the neutral case_root/setup_root/
#     output_dir/workflow_state_path resolver that replaces strict_plan's reuse of
#     describe_launch("sim", ...) for path calculation.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

from __future__ import annotations

import unittest
from pathlib import Path

from openfoam_driver.core.runtime.execution_context import resolve_execution_context
from openfoam_driver.core.runtime.registry import load_entry_spec

from openfoam_driver.tests.conftest import monorepo_root, skip_without_monorepo



@skip_without_monorepo
class TestResolveExecutionContext(unittest.TestCase):
    def test_reports_case_setup_output_and_workflow_state_paths(self) -> None:
        tutorials_root = monorepo_root / "tutorials"  # type: ignore[operator]
        spec = load_entry_spec("singleCell", overrides={"tutorials_root": str(tutorials_root)})

        context = resolve_execution_context(spec)

        self.assertEqual(context.case_root, Path(spec.case_root))
        self.assertEqual(context.setup_root, Path(spec.setup_root))
        self.assertEqual(context.output_dir, Path(spec.output_dir))
        self.assertEqual(
            context.workflow_state_path,
            Path(spec.output_dir) / "workflow_state.json",
        )

    def test_never_re_resolves_the_entry(self) -> None:
        """Takes an already-built spec directly -- no resolve_entry/factory call of
        its own, unlike describe_launch (which strict_plan used to call a second
        time on the same entry purely to get these four paths)."""
        tutorials_root = monorepo_root / "tutorials"  # type: ignore[operator]
        spec = load_entry_spec("singleCell", overrides={"tutorials_root": str(tutorials_root)})

        context = resolve_execution_context(spec)

        # No entry/entry_kind/overrides/config_path parameter exists to pass --
        # the only input is the spec itself.
        self.assertIsInstance(context.case_root, Path)


if __name__ == "__main__":
    unittest.main()
