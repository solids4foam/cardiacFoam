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

from collections.abc import Mapping
from functools import lru_cache
from pathlib import Path
from types import MappingProxyType
from typing import Any

# Authorization is split into two accessors because the two facts are
# different, and core consumes them differently:
#
#   solver_commands()    -- this plugin's artifact-producing solver binaries.
#                           Authorized to run AND eligible to be credited with
#                           a run's unclaimed expected artifacts.
#   auxiliary_commands() -- other plugin-specific commands authorized to run
#                           but which do not produce the run's artifacts.
#
# Their union is the plugin's accept-surface contribution; only
# solver_commands() feeds the artifact-producer heuristic in
# normalize_workflow_dag. Conflating them would credit a post-processing
# utility with the solver's outputs, and `produces` is enforced post-step
# (workflow_runner's missing_artifacts check), so a silent solver would blame
# the utility.
#
# These are NOT core knowledge: a different plugin authorizes different
# commands.
CARDIAC_SOLVER_COMMANDS = frozenset({"cardiacFoam"})

# bathBidomainInterfaceMetrics is listed here rather than coming through
# utility_manifests() because it ships no utility.manifest.toml, so
# UTILITY_CATALOG does not contain it -- yet it is a live workflow step in the
# manufacturedFDABathBidomain tutorial. It is post-processing: it reads the
# reconstructed final-time solution the solver already wrote.
#
# gradientReconstructionOrder (applications/test/gradientReconstructionOrder)
# is the same shape of thing: no utility.manifest.toml, but a live workflow
# step in manufactured_eikonal_ecg.py's gradient_reconstruction=True path,
# appended after the solve step (see that module's _workflow_dag_for).
#
# The error_localisation_analysis=True path's other two steps need no entry
# here: `postProcess` is already core-authorized generically
# (CORE_NEUTRAL_COMMANDS in core/runtime/workflow.py), and the analysis
# script itself is invoked by its own relative path (containing "/"), which
# command_authorization does not gate at all -- see
# workflow_runner._resolve_command: a "/" in the command is used verbatim as
# an explicit opt-in, the same as any case's own Allrun/Allclean script.
CARDIAC_AUXILIARY_COMMANDS = frozenset({
    "bathBidomainInterfaceMetrics",
    "gradientReconstructionOrder",
})


def solver_commands() -> frozenset[str]:
    """This plugin's artifact-producing solver commands."""
    return CARDIAC_SOLVER_COMMANDS


def auxiliary_commands() -> frozenset[str]:
    """Authorized plugin commands that do not produce the run's artifacts."""
    return CARDIAC_AUXILIARY_COMMANDS


def utility_roots() -> tuple[Path, ...]:
    """Roots holding this plugin's ``utility.manifest.toml`` sidecars.

    Derived from ``utility_catalog``'s own constant rather than recomputed
    from ``__file__``, so the two cannot drift apart (``utility_roots()``
    returns ``()`` on a missing directory, which would degrade silently).
    """
    from openfoam_driver.core.utility_catalog import UTILITIES_ROOT

    return (UTILITIES_ROOT,) if UTILITIES_ROOT.is_dir() else ()


@lru_cache(maxsize=1)
def utility_manifests() -> Mapping[str, Any]:
    """Parse this plugin's utility manifests. Cached: parsing walks the tree.

    Returns a read-only view: the cache hands the same object to every caller,
    so a mutable dict would let one consumer corrupt the authorization input of
    all the others.
    """
    from openfoam_driver.core.utility_catalog import load_utility_manifests

    manifests: dict[str, Any] = {}
    for root in utility_roots():
        manifests.update(load_utility_manifests(root))
    return MappingProxyType(manifests)
