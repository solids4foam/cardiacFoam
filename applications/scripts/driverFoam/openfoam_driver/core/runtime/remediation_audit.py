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
#     remediation_audit
#
# Description
#     Append-only audit trail for agent-chosen remediations applied via
#     `step --strict --apply`. One JSONL record per apply, under the strict-plan
#     output directory. Best-effort: never raises, so an audit-append problem can
#     never crash the rerun the agent is waiting on.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


def append_remediation_record(
    output_dir: Path,
    *,
    step_id: str,
    attempt: int,
    applied_overrides: list[dict[str, Any]],
    resulting_status: str,
) -> None:
    """Append one audit line. Best-effort: never raises (must not crash a rerun)."""
    record = {
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "step_id": step_id,
        "attempt": attempt,
        "applied_overrides": applied_overrides,
        "resulting_status": resulting_status,
    }
    try:
        path = Path(output_dir) / "remediation_history.jsonl"
        with path.open("a", encoding="utf-8") as handle:
            handle.write(json.dumps(record) + "\n")
    except OSError:
        return
