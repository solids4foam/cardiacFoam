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
#     conftest
#
# Description
#     Configures shared pytest fixtures and runtime dependencies.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

import os
import pytest
from pathlib import Path

os.environ["SKIP_ENV_DIAGNOSTICS"] = "1"


def _find_monorepo_root() -> Path | None:
    """Walk parent directories looking for the cardiacFoam monorepo root.

    Returns the first ancestor that has both ``tutorials/`` and
    ``applications/`` siblings, or ``None`` when running in a standalone
    (temp-folder / CI) checkout that does not include the full monorepo tree.
    """
    current = Path(__file__).resolve()
    for parent in current.parents:
        if (parent / "tutorials").exists() and (parent / "applications").exists():
            return parent
    return None


#: The monorepo root resolved once at collection time.  ``None`` in standalone.
monorepo_root: Path | None = _find_monorepo_root()

#: Apply this decorator to any test class/function that reads real tutorial
#: case directories from the monorepo ``tutorials/`` tree.  The test is
#: automatically skipped in standalone clones and CI environments.
skip_without_monorepo = pytest.mark.skipif(
    monorepo_root is None,
    reason=(
        "Requires the full cardiacFoam monorepo tree (tutorials/ + applications/). "
        "Clone the full repository or run with --tutorials-root to enable this test."
    ),
)
