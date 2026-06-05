"""Python model for the cardiacFoam Run document.

Its shape is defined in ``schemas/run-document.json`` (the single source
of truth); this module provides a Python dataclass for code that wants
to construct, validate, or round-trip a Run programmatically.
"""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Literal

import jsonschema

Phase = Literal["anatomy", "physics", "stimulus", "solver"]
Status = Literal["draft", "queued", "planning", "planned", "running", "completed", "failed"]

_SCHEMA_PATH = (
    Path(__file__).resolve().parents[3] / "schemas" / "run-document.json"
)
_SCHEMA = json.loads(_SCHEMA_PATH.read_text())


@dataclass
class RunDocument:
    """Run document as defined by ``schemas/run-document.json``.

    Construction does not validate; call :meth:`to_json` to produce a
    schema-validated dict, or :meth:`from_json` to parse with validation.
    """

    id: str
    name: str
    status: Status
    config: dict[str, dict[str, Any]]
    version: str = "2"
    createdAt: str = ""
    lastModified: str = ""
    intent: dict[str, Any] = field(default_factory=dict)
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
        if data.get("version") != "2":
            raise ValueError(
                "RunDocument.to_json emits only version '2'; use "
                "RunDocument.migrate_v1(...) before serializing old documents"
            )
        if data.get("reports") is None:
            data.pop("reports", None)
        jsonschema.validate(data, _SCHEMA)
        return data

    @classmethod
    def from_json(cls, data: dict[str, Any]) -> "RunDocument":
        if data.get("version") == "1":
            raise ValueError(
                "RunDocument.from_json expects version '2'. Use "
                "RunDocument.migrate_v1(data) to migrate v1 input explicitly."
            )
        jsonschema.validate(data, _SCHEMA)
        return cls(
            id=data["id"],
            name=data["name"],
            status=data["status"],
            config=data["config"],
            version=data.get("version", "2"),
            createdAt=data.get("createdAt", ""),
            lastModified=data.get("lastModified", ""),
            intent=data.get("intent", {}),
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
        """Return a v2 RunDocument from the legacy v1 shape.

        The migration is deliberately conservative: it preserves the config,
        validation, results, reports, and timestamps, then adds empty v2 planning
        fields. Callers must still run the strict planner to populate
        resolvedEntry, workflowDag, launch, and expectedArtifacts.
        """
        if data.get("version") != "1":
            raise ValueError("migrate_v1 expects a RunDocument with version '1'")
        migrated = {
            "version": "2",
            "id": data["id"],
            "name": data["name"],
            "createdAt": data.get("createdAt", ""),
            "lastModified": data.get("lastModified", ""),
            "status": data.get("status", "draft"),
            "intent": {},
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
