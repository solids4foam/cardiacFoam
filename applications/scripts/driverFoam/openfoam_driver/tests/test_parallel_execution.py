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
#     test_parallel_execution
#
# Description
#     Tests the shared decomposePar/mpirun/reconstructPar workflow_dag step
#     builder used by every manufactured-solution tutorial's _workflow_dag_for.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

from __future__ import annotations

import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

from openfoam_driver.core.runtime.parallel_execution import (
    read_number_of_subdomains,
    solve_steps,
)

_DECOMPOSE_PAR_DICT = """FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      decomposeParDict;
}
numberOfSubdomains  6;
method          scotch;
"""


class TestReadNumberOfSubdomains(unittest.TestCase):
    def test_reads_committed_value(self) -> None:
        with TemporaryDirectory() as tmp:
            case_root = Path(tmp)
            (case_root / "system").mkdir()
            (case_root / "system" / "decomposeParDict").write_text(_DECOMPOSE_PAR_DICT)
            self.assertEqual(read_number_of_subdomains(case_root), 6)


class TestSolveSteps(unittest.TestCase):
    def test_serial_is_a_single_solve_step(self) -> None:
        steps, last_id = solve_steps(
            solve_id="solve", solve_command="cardiacFoam",
            depends_on=["setConductivity"], run_in_parallel=False, case_root=None,
        )
        self.assertEqual(steps, [
            {"id": "solve", "command": "cardiacFoam", "depends_on": ["setConductivity"]},
        ])
        self.assertEqual(last_id, "solve")

    def test_parallel_wraps_decompose_mpirun_reconstruct(self) -> None:
        with TemporaryDirectory() as tmp:
            case_root = Path(tmp)
            (case_root / "system").mkdir()
            (case_root / "system" / "decomposeParDict").write_text(_DECOMPOSE_PAR_DICT)

            steps, last_id = solve_steps(
                solve_id="solve", solve_command="cardiacFoam",
                depends_on=["setConductivity"], run_in_parallel=True, case_root=case_root,
            )
            self.assertEqual([s["id"] for s in steps], ["decomposePar", "solve", "reconstructPar"])
            self.assertEqual(steps[0]["command"], "decomposePar")
            self.assertEqual(steps[0]["depends_on"], ["setConductivity"])
            self.assertEqual(steps[1]["command"], "mpirun")
            self.assertEqual(steps[1]["args"], ["-np", "6", "cardiacFoam", "-parallel"])
            self.assertEqual(steps[1]["depends_on"], ["decomposePar"])
            self.assertEqual(steps[2]["command"], "reconstructPar")
            self.assertEqual(steps[2]["depends_on"], ["solve"])
            self.assertEqual(last_id, "reconstructPar")

    def test_parallel_requires_case_root(self) -> None:
        with self.assertRaises(ValueError):
            solve_steps(
                solve_id="solve", solve_command="cardiacFoam",
                depends_on=[], run_in_parallel=True, case_root=None,
            )


if __name__ == "__main__":
    unittest.main()
