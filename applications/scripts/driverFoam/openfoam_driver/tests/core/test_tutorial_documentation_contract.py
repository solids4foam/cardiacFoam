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
#     test_tutorial_documentation_contract
#
# Description
#     Keeps the tutorial index aligned with case paths and executable registries.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

from pathlib import Path
import re
import subprocess
import pytest
from openfoam_driver.tests.conftest import skip_without_monorepo
pytestmark = skip_without_monorepo

from openfoam_driver.core.runtime.registry import list_tutorials, _normalized_registry


REPO_ROOT = Path(__file__).resolve().parents[6]
TUTORIALS_ROOT = REPO_ROOT / "tutorials"
README = TUTORIALS_ROOT / "README.md"
REGRESSION_RUNNER = TUTORIALS_ROOT / "Alltest-regression"


def _documented_rows() -> dict[str, dict[str, str]]:
    rows = {}
    for line in README.read_text().splitlines():
        if not line.startswith("| `"):
            continue
        cells = [cell.strip() for cell in line.strip("|").split("|")]
        if len(cells) != 7:
            continue
        path, purpose, solver, mode, driver, regression, outputs = cells
        path = path.strip("`")
        rows[path] = {
            "purpose": purpose,
            "solver": solver,
            "mode": mode,
            "driver": driver.strip("`"),
            "regression": regression,
            "outputs": outputs,
        }
    return rows


def _runner_cases() -> set[str]:
    text = REGRESSION_RUNNER.read_text()
    match = re.search(r"REGRESSION_TESTS=\(\n(?P<body>.*?)\n\)", text, re.DOTALL)
    assert match, "Alltest-regression has no parseable REGRESSION_TESTS array"
    return set(re.findall(r'^\s+"([^"]+)"$', match.group("body"), re.MULTILINE))


def _registered_case_paths() -> dict[str, str]:
    paths = {}
    for entry in list_tutorials():
        spec = _normalized_registry()[entry.casefold()](tutorials_root=TUTORIALS_ROOT)
        paths[entry] = str(Path(spec.case_root).relative_to(TUTORIALS_ROOT))
    return paths


def test_documented_canonical_paths_exist():
    rows = _documented_rows()
    assert rows, "tutorial index contains no canonical case rows"
    for relative_path in rows:
        case_root = TUTORIALS_ROOT / relative_path
        assert case_root.is_dir(), f"documented tutorial path does not exist: {relative_path}"
        assert (case_root / "Allrun").is_file(), f"documented tutorial has no Allrun: {relative_path}"
        tracked = subprocess.run(
            ["git", "ls-files", "--error-unmatch", f"tutorials/{relative_path}/Allrun"],
            cwd=REPO_ROOT,
            capture_output=True,
            text=True,
            check=False,
        )
        assert tracked.returncode == 0, f"documented tutorial is not committed: {relative_path}"


def test_documented_driver_registration_is_exact():
    registered_paths = _registered_case_paths()
    for relative_path, row in _documented_rows().items():
        driver = row["driver"]
        if driver == "—":
            assert relative_path not in registered_paths.values(), (
                f"{relative_path} is registered but documented without its driver entry"
            )
        else:
            assert driver in registered_paths, f"unknown documented driver entry: {driver}"
            assert registered_paths[driver] == relative_path


def test_documented_regression_coverage_matches_runner():
    rows = _documented_rows()
    documented_covered = {
        path for path, row in rows.items()
        if row["regression"].startswith("`Alltest-regression`")
    }
    runner_covered = _runner_cases()
    assert documented_covered == runner_covered
