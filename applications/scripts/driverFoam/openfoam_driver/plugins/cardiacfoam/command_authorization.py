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
#     command_authorization
#
# Description
#     Declares the workflow commands this plugin authorizes: the cardiacFoam
#     solver binary, its solver-specific utilities, and the utility manifests
#     discovered under this plugin's utilities root.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

from __future__ import annotations

from functools import lru_cache
from pathlib import Path
from typing import Any

# The cardiacFoam solver binary and this plugin's solver-specific commands.
# These are NOT core knowledge: a different plugin authorizes different
# commands. bathBidomainInterfaceMetrics is listed here rather than coming
# through utility_manifests() because it ships no utility.manifest.toml, so
# UTILITY_CATALOG does not contain it -- yet it is a live workflow step in the
# manufacturedFDABathBidomain tutorial.
CARDIAC_SOLVER_COMMANDS = frozenset(
    {"cardiacFoam", "bathBidomainInterfaceMetrics"}
)

# Repo-relative root holding this plugin's `utility.manifest.toml` sidecars.
_UTILITIES_ROOT = (
    Path(__file__).resolve().parents[6] / "applications" / "utilities"
)


def solver_commands() -> frozenset[str]:
    return CARDIAC_SOLVER_COMMANDS


def utility_roots() -> tuple[Path, ...]:
    return (_UTILITIES_ROOT,) if _UTILITIES_ROOT.is_dir() else ()


@lru_cache(maxsize=1)
def utility_manifests() -> dict[str, Any]:
    """Parse this plugin's utility manifests. Cached: parsing walks the tree."""
    from openfoam_driver.utility_catalog import load_utility_manifests

    manifests: dict[str, Any] = {}
    for root in utility_roots():
        manifests.update(load_utility_manifests(root))
    return manifests
