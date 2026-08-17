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
#     run_model
#
# Description
#     Defines programmatic models for cardiacFoam execution payload documents.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

"""Python model for the cardiacFoam Run document.

Its shape is defined in ``schemas/run-document.json`` (the single source
of truth); this module provides a Python dataclass for code that wants
to construct, validate, or round-trip a Run programmatically.
"""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass, field
from importlib import resources
from typing import Any, Literal

import jsonschema

Phase = Literal["anatomy", "physics", "stimulus", "solver"]
Status = Literal["draft", "queued", "planning", "planned", "running", "completed", "failed"]

_SCHEMA = json.loads(
    resources.files("openfoam_driver.schemas")
    .joinpath("run-document.json")
    .read_text()
)


@dataclass
class RunDocument:
    """Run document as defined by ``schemas/run-document.json``.

    Construction does not validate; call :meth:`to_json` to produce a
    schema-validated dict, or :meth:`from_json` to parse with validation.
    """

    id: str
    name: str
    status: Status
    # Plugin-defined: the core schema constrains ``config`` to an object but
    # imposes no shape on the per-phase values (P2.2). Annotating the values
    # as ``dict`` would assert a guarantee the schema no longer makes;
    # ``specs.validation.validate_run`` enforces the mapping shape and
    # reports violations as diagnostics.
    config: dict[str, Any]
    version: str = "3"
    createdAt: str = ""
    lastModified: str = ""
    intent: dict[str, Any] = field(default_factory=dict)
    plugin: dict[str, str] | None = None
    resolvedEntry: dict[str, Any] | None = None
    workflowDag: dict[str, Any] | None = None
    workflowState: dict[str, Any] | None = None
    launch: dict[str, Any] | None = None
    expectedArtifacts: list[dict[str, Any]] = field(default_factory=list)
    validation: dict[str, Any] = field(default_factory=dict)
    results: dict[str, Any] | None = None
    reports: dict[str, Any] | None = None
    terminalStatusValues: list[str] = field(
        default_factory=lambda: ["completed", "failed"]
    )

    def to_json(self) -> dict[str, Any]:
        data = asdict(self)
        if data.get("version") != "3":
            raise ValueError(
                "RunDocument.to_json emits only version '3'; use "
                "RunDocument.migrate_v1(...) or RunDocument.migrate_v2(...) "
                "before serializing old documents"
            )
        if data.get("reports") is None:
            data.pop("reports", None)
        if data.get("plugin") is None:
            data.pop("plugin", None)
        jsonschema.validate(data, _SCHEMA)
        return data

    @classmethod
    def from_json(cls, data: dict[str, Any]) -> "RunDocument":
        if data.get("version") in ("1", "2"):
            raise ValueError(
                f"RunDocument.from_json expects version '3'. Use "
                f"RunDocument.migrate_v1(data) or RunDocument.migrate_v2(data) "
                f"to migrate version {data.get('version')!r} input explicitly."
            )
        jsonschema.validate(data, _SCHEMA)
        return cls(
            id=data["id"],
            name=data["name"],
            status=data["status"],
            config=data["config"],
            version=data.get("version", "3"),
            createdAt=data.get("createdAt", ""),
            lastModified=data.get("lastModified", ""),
            intent=data.get("intent", {}),
            plugin=data.get("plugin"),
            resolvedEntry=data.get("resolvedEntry"),
            workflowDag=data.get("workflowDag"),
            workflowState=data.get("workflowState"),
            launch=data.get("launch"),
            expectedArtifacts=data.get("expectedArtifacts", []),
            validation=data.get("validation", {}),
            results=data.get("results"),
            reports=data.get("reports"),
            terminalStatusValues=data.get(
                "terminalStatusValues", ["completed", "failed"]
            ),
        )

    @classmethod
    def migrate_v1(cls, data: dict[str, Any]) -> "RunDocument":
        """Return a v3 RunDocument from the legacy v1 shape.

        The migration is deliberately conservative: it preserves the config,
        validation, results, reports, and timestamps, then adds empty v3 planning
        fields. Callers must still run the strict planner to populate
        resolvedEntry, workflowDag, launch, and expectedArtifacts.
        """
        if data.get("version") != "1":
            raise ValueError("migrate_v1 expects a RunDocument with version '1'")
        migrated = {
            "version": "3",
            "id": data["id"],
            "name": data["name"],
            "createdAt": data.get("createdAt", ""),
            "lastModified": data.get("lastModified", ""),
            "status": data.get("status", "draft"),
            "intent": {},
            "plugin": None,
            "config": data["config"],
            "resolvedEntry": None,
            "workflowDag": None,
            "workflowState": None,
            "launch": None,
            "expectedArtifacts": [],
            "validation": data.get("validation", {}),
            "results": data.get("results"),
            "terminalStatusValues": ["completed", "failed"],
        }
        if "reports" in data:
            migrated["reports"] = data["reports"]
        return cls.from_json(migrated)

    @classmethod
    def migrate_v2(cls, data: dict[str, Any]) -> "RunDocument":
        """Return a v3 RunDocument from the v2 shape.

        v2's config was schema-constrained by core to the cardiac phase
        envelope; v3 makes config an open object validated by the plugin's
        own schema instead. The migration is a pure version-field bump --
        config, validation, results, reports, and timestamps carry over
        unchanged, since v2's config already satisfies any v3 plugin schema
        shaped like the (now plugin-owned) v2 constraint.
        """
        if data.get("version") != "2":
            raise ValueError("migrate_v2 expects a RunDocument with version '2'")
        migrated = dict(data)
        migrated["version"] = "3"
        return cls.from_json(migrated)
