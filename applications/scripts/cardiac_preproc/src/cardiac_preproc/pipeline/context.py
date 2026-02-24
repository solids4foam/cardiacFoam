"""Shared context/result objects for pipeline steps."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any


@dataclass
class StepContext:
    """Generic step context passed to step runners."""

    project_root: Path
    current_mesh: str | None = None
    output_dir: Path | None = None
    options: dict[str, Any] = field(default_factory=dict)


@dataclass
class StepResult:
    """Standard step return type."""

    output_mesh: str | None = None
    metadata: dict[str, Any] = field(default_factory=dict)
