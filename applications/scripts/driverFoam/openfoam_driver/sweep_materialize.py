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
#     Materializes a resolved sweep case via build_and_launch, Allrun, and workflow_contract.json.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

from __future__ import annotations

import json
import stat
from pathlib import Path
from typing import Any

from .specs.dict_builder import build_and_launch


def materialize_case(*, case_dir: Path, routed: dict[str, Any]) -> None:
    """Write a resolved+routed case's dict files, Allrun script, and
    workflow_contract.json. Raises ValueError (propagated from
    build_and_launch/build_electro_properties) if the routed selectors are
    structurally invalid — the caller treats that as this case's failure,
    not a crash of the whole sweep.
    """
    result = build_and_launch(
        electro_selectors=routed["electro_selectors"],
        physics_selectors=routed["physics_selectors"],
        case_dir=case_dir,
        electro_overrides=routed["electro_overrides"] or None,
        physics_overrides=routed["physics_overrides"] or None,
        delta_t=routed["delta_t"],
        end_time=routed["end_time"],
        dx=routed.get("dx"),
        dry_run=True,
        overwrite=True,
    )

    allrun_body = "blockMesh\ncardiacFoam\n" if result.get("needs_block_mesh") else "cardiacFoam\n"
    allrun_path = case_dir / "Allrun"
    allrun_path.write_text("#!/bin/sh\n" + allrun_body)
    allrun_path.chmod(allrun_path.stat().st_mode | stat.S_IEXEC | stat.S_IXGRP | stat.S_IXOTH)

    contract_path = case_dir / "workflow_contract.json"
    contract_path.write_text(json.dumps({"steps": [{"id": "run", "command": "Allrun", "depends_on": []}]}))
