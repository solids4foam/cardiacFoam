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
#     test_cli_legacy_preflight
#
# Description
#     Preflight gate for legacy sim/all/post actions in cli.py.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

import json
import unittest
from unittest import mock

from openfoam_driver.cli import main
from openfoam_driver.planning_types import diagnostic


class TestLegacyPreflight(unittest.TestCase):
    def _run_legacy(self, action: str, *, env_error: bool) -> tuple[int, str]:
        """Run cli main() for a legacy action; return (exit_code, stdout_capture)."""
        captured = []

        def fake_print(*args, **kwargs):
            captured.append(" ".join(str(a) for a in args))

        with mock.patch("openfoam_driver.cli._environment_diagnostics") as mock_diag, \
             mock.patch("builtins.print", side_effect=fake_print):

            if env_error:
                mock_diag.return_value = (
                    diagnostic(
                        "error",
                        "missing_openfoam_env",
                        "WM_PROJECT_DIR is not set. OpenFOAM environment not sourced.",
                        source="environment",
                    ),
                )
            else:
                mock_diag.return_value = ()

            with mock.patch("openfoam_driver.core.runtime.engine.DriverEngine.run_simulations") as mock_run, \
                 mock.patch("openfoam_driver.core.runtime.engine.DriverEngine.run_all") as mock_all, \
                 mock.patch("openfoam_driver.core.runtime.engine.DriverEngine.run_postprocess") as mock_post:
                mock_run.return_value = []
                mock_all.return_value = []
                try:
                    code = main([action, "--entry", "singleCell"])
                except SystemExit as exc:
                    code = int(exc.code) if exc.code is not None else 0

        return code, "\n".join(captured)

    def test_legacy_sim_blocked_by_env_error(self):
        code, output = self._run_legacy("sim", env_error=True)
        assert code == 1, f"expected exit code 1, got {code}"
        payload = json.loads(output)
        assert payload["status"] == "failed"
        assert "environment_diagnostics" in payload
        assert payload["action"] == "sim"

    def test_legacy_all_blocked_by_env_error(self):
        code, output = self._run_legacy("all", env_error=True)
        assert code == 1
        payload = json.loads(output)
        assert payload["status"] == "failed"
        assert payload["action"] == "all"

    def test_legacy_sim_passes_through_on_clean_env(self):
        code, _ = self._run_legacy("sim", env_error=False)
        assert code == 0, "clean env must not be blocked"

    def test_legacy_preflight_output_has_entry_field(self):
        code, output = self._run_legacy("sim", env_error=True)
        payload = json.loads(output)
        assert payload["entry"] == "singleCell"


if __name__ == "__main__":
    unittest.main()
