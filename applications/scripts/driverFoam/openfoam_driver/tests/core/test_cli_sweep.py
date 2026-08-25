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
#     test_cli_sweep
#
# Description
#     Tests CLI wiring for sweep-plan/sweep-run actions.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

import json
import unittest
from unittest import mock

from openfoam_driver.cli import main


class TestCliSweepActions(unittest.TestCase):
    def test_sweep_plan_dispatches_to_sweep_runner(self):
        captured = []

        def fake_print(*args, **kwargs):
            captured.append(" ".join(str(a) for a in args))

        with mock.patch("openfoam_driver.cli.sweep_plan", return_value={"case_count": 2, "cases": []}) as mock_fn, \
             mock.patch("builtins.print", side_effect=fake_print):
            code = main(["sweep-plan", "--spec", "sweep.json", "--output-dir", "/tmp/out"])

        assert code == 0
        mock_fn.assert_called_once()
        payload = json.loads(captured[0])
        assert payload["case_count"] == 2

    def test_sweep_plan_defaults_to_repo_scratch_output_dir(self):
        with mock.patch(
            "openfoam_driver.cli.sweep_plan",
            return_value={"case_count": 0, "cases": []},
        ) as mock_fn, mock.patch("builtins.print"):
            assert main(["sweep-plan", "--spec", "paperI_methods.json"]) == 0

        assert ".tmp/driverfoam/sweeps/paperI_methods" in str(
            mock_fn.call_args.kwargs["output_dir"]
        )

    def test_sweep_run_defaults_to_repo_scratch_output_dir(self):
        with mock.patch(
            "openfoam_driver.cli.sweep_run",
            return_value={"case_count": 0, "completed_count": 0, "failed_count": 0},
        ) as mock_fn, mock.patch("builtins.print"):
            assert main(["sweep-run", "--spec", "paperI_methods.json"]) == 0

        assert ".tmp/driverfoam/sweeps/paperI_methods" in str(
            mock_fn.call_args.kwargs["output_dir"]
        )

    def test_sweep_run_passes_max_cases_and_retry_flag(self):
        captured = []

        def fake_print(*args, **kwargs):
            captured.append(" ".join(str(a) for a in args))

        with mock.patch(
            "openfoam_driver.cli.sweep_run",
            return_value={"case_count": 1, "completed_count": 1, "failed_count": 0, "skipped_count": 0},
        ) as mock_fn, mock.patch("builtins.print", side_effect=fake_print):
            code = main([
                "sweep-run", "--spec", "sweep.json", "--output-dir", "/tmp/out",
                "--max-cases", "500", "--retry-failed",
            ])

        assert code == 0
        mock_fn.assert_called_once()
        kwargs = mock_fn.call_args.kwargs
        assert kwargs["max_cases"] == 500
        assert kwargs["retry_failed"] is True

    def test_sweep_run_passes_fresh_flag(self):
        captured = []

        def fake_print(*args, **kwargs):
            captured.append(" ".join(str(a) for a in args))

        with mock.patch(
            "openfoam_driver.cli.sweep_run",
            return_value={"case_count": 1, "completed_count": 1, "failed_count": 0, "skipped_count": 0},
        ) as mock_fn, mock.patch("builtins.print", side_effect=fake_print):
            code = main([
                "sweep-run", "--spec", "sweep.json", "--output-dir", "/tmp/out", "--fresh",
            ])

        assert code == 0
        mock_fn.assert_called_once()
        assert mock_fn.call_args.kwargs["fresh"] is True

    def test_fresh_and_retry_failed_are_mutually_exclusive(self):
        with mock.patch("builtins.print"):
            with self.assertRaises(SystemExit):
                main([
                    "sweep-run", "--spec", "sweep.json", "--output-dir", "/tmp/out",
                    "--fresh", "--retry-failed",
                ])

    def test_fresh_rejected_for_describe_action(self):
        with mock.patch("builtins.print"):
            with self.assertRaises(SystemExit):
                main(["describe", "--entry", "singleCell", "--fresh"])

    def test_sweep_plan_exits_nonzero_when_a_case_failed(self):
        with mock.patch(
            "openfoam_driver.cli.sweep_plan",
            return_value={
                "case_count": 2,
                "cases": [
                    {"case_id": "TNNP", "status": "ok"},
                    {"case_id": "BadModel", "status": "failed", "materialization_error": "bad ionicModel"},
                ],
            },
        ), mock.patch("builtins.print"):
            code = main(["sweep-plan", "--spec", "sweep.json", "--output-dir", "/tmp/out"])

        assert code == 1

    def test_sweep_plan_rejects_entry_flag(self):
        with mock.patch("builtins.print"):
            with self.assertRaises(SystemExit):
                main(["sweep-plan", "--spec", "sweep.json", "--output-dir", "/tmp/out", "--entry", "singleCell"])

    def test_non_sweep_action_rejects_spec_flag(self):
        with mock.patch("builtins.print"):
            with self.assertRaises(SystemExit):
                main(["plan", "--strict", "--entry", "singleCell", "--spec", "sweep.json"])


if __name__ == "__main__":
    unittest.main()
