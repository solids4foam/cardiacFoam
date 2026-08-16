"""Generic dictionary-catalog value objects.

The core owns this shape and its generic constraint vocabulary. Individual
plugins own the entries and document-specific catalogues built from it.
"""

from __future__ import annotations

from dataclasses import dataclass, field, replace
from typing import Any

from ..runtime.run_model import Phase


@dataclass(frozen=True)
class DictEntry:
    driver_path: str
    description: str
    source_refs: tuple[str, ...] = ()
    notes: str = ""
    value_kind: str = "openfoam_literal"
    enum_values: tuple[str, ...] = ()
    examples: tuple[str, ...] = ()
    dynamic_path: bool = False
    required: bool = False
    constraints: tuple[str, ...] = ()
    unit: str = ""
    typical_value: str = ""
    phases: frozenset[str] = frozenset()
    applicable_when: dict[str, str | tuple[str, ...]] = field(default_factory=dict)
    forbidden_when: dict[str, str | tuple[str, ...]] = field(default_factory=dict)
    required_when: dict[str, str | tuple[str, ...]] = field(default_factory=dict)
    mutually_exclusive_with: tuple[str, ...] = ()


def build_group(
    defaults: dict[str, Any],
    entries: tuple[DictEntry, ...],
) -> tuple[DictEntry, ...]:
    """Apply group defaults without mutating plugin-owned entry objects."""
    out = []
    for entry in entries:
        changes = {}
        for key, value in defaults.items():
            current = getattr(entry, key)
            if not current:
                changes[key] = value
            elif isinstance(current, dict) and isinstance(value, dict):
                changes[key] = {**value, **current}
        out.append(replace(entry, **changes) if changes else entry)
    return tuple(out)
