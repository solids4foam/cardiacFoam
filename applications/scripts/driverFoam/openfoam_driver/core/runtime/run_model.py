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
Status = Literal["draft", "queued", "running", "completed", "failed"]

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
    version: str = "1"
    createdAt: str = ""
    lastModified: str = ""
    validation: dict[str, Any] = field(default_factory=dict)
    results: dict[str, Any] | None = None

    def to_json(self) -> dict[str, Any]:
        data = asdict(self)
        jsonschema.validate(data, _SCHEMA)
        return data

    @classmethod
    def from_json(cls, data: dict[str, Any]) -> "RunDocument":
        jsonschema.validate(data, _SCHEMA)
        return cls(
            id=data["id"],
            name=data["name"],
            status=data["status"],
            config=data["config"],
            version=data.get("version", "1"),
            createdAt=data.get("createdAt", ""),
            lastModified=data.get("lastModified", ""),
            validation=data.get("validation", {}),
            results=data.get("results"),
        )
