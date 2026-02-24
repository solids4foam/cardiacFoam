"""cardiac_preproc package."""

from __future__ import annotations

from typing import Any

__all__ = ["CardiacPreprocEngine"]


def __getattr__(name: str) -> Any:
    if name == "CardiacPreprocEngine":
        from .engine import CardiacPreprocEngine
        return CardiacPreprocEngine
    raise AttributeError(name)
