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
#     sweep_materialize
#
# Description
#     Materializes a resolved sweep case via build_and_launch and Allrun.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

from __future__ import annotations

from pathlib import Path
from typing import Any


def _materialize_case_legacy(*, case_dir: Path, routed: dict[str, Any]) -> None:
    """Write a resolved+routed case's dict files and Allrun script.

    Raises ValueError (propagated from build_and_launch/
    build_electro_properties) if the routed selectors are structurally invalid
    — the caller treats that as this case's failure, not a crash of the whole
    sweep. Per-run workflow intent/state is persisted later by strict planning
    and execution as ``run_document.json`` and ``workflow_state.json``.
    """
    from .plugins.cardiacfoam.sweep import materialize_case

    materialize_case(case_dir=case_dir, routed=routed)


def materialize_case(
    *, case_dir: Path, routed: dict[str, Any], driver_context=None,
) -> None:
    """Materialize via the selected plugin while preserving the public API."""

    from .core.compatibility import resolve_public_driver_context
    from .core.plugin_capabilities import SweepMaterializationRequest

    driver_context = resolve_public_driver_context(driver_context)
    driver_context.capabilities.sweep_materializer.materialize(
        SweepMaterializationRequest(case_dir=case_dir, routed=routed),
    )
