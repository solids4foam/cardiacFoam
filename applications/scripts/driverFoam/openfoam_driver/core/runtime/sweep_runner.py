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
#     sweep_runner
#
# Description
#     Orchestrates sweep-plan: expand, materialize, and audit each resolved case.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from ...strict_planning import strict_plan
from ...sweep_derivation_catalog import get_derivation
from ...sweep_expansion import check_case_count_cap, expand_sweep
from ...sweep_materialize import materialize_case
from ...sweep_routing import route_case_values


def _load_spec(spec_path: str | Path) -> dict[str, Any]:
    return json.loads(Path(spec_path).read_text())


def sweep_plan(
    spec_path: str | Path,
    *,
    output_dir: str | Path,
    max_cases: int = 200,
) -> dict[str, Any]:
    sweep_spec = _load_spec(spec_path)
    check_case_count_cap(sweep_spec, max_cases=max_cases)

    output_dir = Path(output_dir)
    resolved_cases = expand_sweep(sweep_spec, get_derivation=get_derivation)
    base = sweep_spec.get("base", {})

    case_reports = []
    for case in resolved_cases:
        routed = route_case_values(base=base, resolved_axis_values=case.resolved_axis_values)
        case_dir = output_dir / case.case_id
        try:
            materialize_case(case_dir=case_dir, routed=routed)
        except (OSError, ValueError) as exc:
            case_reports.append(
                {
                    "case_id": case.case_id,
                    "resolved_axis_values": case.resolved_axis_values,
                    "status": "failed",
                    "materialization_error": str(exc),
                }
            )
            continue

        report = strict_plan(case.case_id, entry_kind="case_folder", overrides={"tutorials_root": str(output_dir)})
        report_payload = report.to_json()
        case_reports.append(
            {
                "case_id": case.case_id,
                "resolved_axis_values": case.resolved_axis_values,
                "status": report.status,
                "plan": report_payload,
            }
        )

    return {"case_count": len(resolved_cases), "cases": case_reports}
