"""Immutable, plugin-owned catalogue partitioned by OpenFOAM document."""

from __future__ import annotations

from dataclasses import dataclass, field
from types import MappingProxyType
from typing import Mapping

from .dictionary import DictEntry


@dataclass(frozen=True)
class DictionaryCatalog:
    """Dictionary entries grouped by their concrete OpenFOAM document.

    Document names are plugin vocabulary (for example ``electroProperties``),
    not a core requirement. The core only relies on deterministic grouping and
    globally unique driver paths.
    """

    documents: Mapping[str, tuple[DictEntry, ...]]
    _entries: tuple[DictEntry, ...] = field(init=False, repr=False)

    def __post_init__(self) -> None:
        normalised: dict[str, tuple[DictEntry, ...]] = {}
        paths: set[str] = set()
        all_entries: list[DictEntry] = []
        for name, entries in self.documents.items():
            if not isinstance(name, str) or not name:
                raise ValueError("DictionaryCatalog document names must be non-empty strings")
            values = tuple(entries)
            for entry in values:
                if not isinstance(entry, DictEntry) or not entry.driver_path:
                    raise ValueError("DictionaryCatalog entries must be DictEntry values with paths")
                if entry.driver_path in paths:
                    raise ValueError(f"DictionaryCatalog has duplicate path {entry.driver_path!r}")
                paths.add(entry.driver_path)
            normalised[name] = values
            all_entries.extend(values)
        object.__setattr__(self, "documents", MappingProxyType(normalised))
        object.__setattr__(self, "_entries", tuple(all_entries))

    @property
    def entries(self) -> tuple[DictEntry, ...]:
        return self._entries

    def entries_for(self, document: str) -> tuple[DictEntry, ...]:
        return self.documents.get(document, ())
