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
#     planning_types
#
# Description
#     Shared strict-planning diagnostic and audit data types.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import Any

from .runtime.models import DataArtifact


@dataclass(frozen=True)
class StrictDiagnostic:
    level: str
    code: str
    message: str
    source: str = ""
    field: str = ""


@dataclass(frozen=True)
class SimulationAuditItem:
    stage: str
    status: str
    points: int
    max_points: int
    summary: str
    evidence: dict[str, Any] = field(default_factory=dict)


def diagnostic(
    level: str,
    code: str,
    message: str,
    *,
    source: str = "",
    field: str = "",
) -> StrictDiagnostic:
    return StrictDiagnostic(
        level=level,
        code=code,
        message=message,
        source=source,
        field=field,
    )


def has_error(diagnostics: tuple[StrictDiagnostic, ...]) -> bool:
    return any(diagnostic.level == "error" for diagnostic in diagnostics)


def has_warning(diagnostics: tuple[StrictDiagnostic, ...]) -> bool:
    return any(diagnostic.level == "warning" for diagnostic in diagnostics)


def artifact_to_json(artifact: DataArtifact) -> dict[str, Any]:
    payload = asdict(artifact)
    payload["variables"] = list(artifact.variables)
    return payload

