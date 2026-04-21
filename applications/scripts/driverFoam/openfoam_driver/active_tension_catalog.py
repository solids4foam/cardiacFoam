"""
Active Tension Model Catalog

A static, build-time-generated catalog of active tension models.
This module exposes the exact variables each model supports so an autonomous
agent can plan active-tension output without running the solver.

All variable names are extracted from C++ source files and are guaranteed to be
exact.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Final


@dataclass(frozen=True)
class ActiveTensionModelEntry:
    """Metadata for a single active tension model."""

    states: tuple[str, ...]
    """State variables."""

    algebraic: tuple[str, ...]
    """Algebraic variables."""

    constants: tuple[str, ...]
    """Constant parameters."""

    rates: tuple[str, ...]
    """Rate variables (derivatives of states)."""

    recommended_exports: tuple[str, ...]
    """Minimum useful export set."""

    description: str
    """Human-readable description."""

    notes: str = ""
    """Additional notes or warnings."""

    aliases: tuple[str, ...] = ()
    """Alternative names for this model."""


ACTIVE_TENSION_MODEL_CATALOG: Final[dict[str, ActiveTensionModelEntry]] = {
    "GoktepeKuhl": ActiveTensionModelEntry(
        states=("Ta",),
        algebraic=("AV_e", "AV_Vm", "AV_u"),
        constants=("AC_Vr", "AC_eInfty", "AC_e0", "AC_eXi", "AC_Vshift", "AC_kTa"),
        rates=("Ta",),
        recommended_exports=("Ta",),
        description="Goktepe-Kuhl active tension model (2004).",
        aliases=("Goktepe-Kuhl", "active stress model"),
    ),
    "NashPanfilov": ActiveTensionModelEntry(
        states=("Ta",),
        algebraic=("AV_u", "AV_e"),
        constants=("AC_Vp", "AC_Vr", "AC_Vth", "AC_e0", "AC_kTa", "AC_a"),
        rates=("Ta",),
        recommended_exports=("Ta",),
        description="Nash-Panfilov active tension model (2004).",
        aliases=("Nash-Panfilov", "phenomenological active stress"),
    ),
}


def get_active_tension_entry(name: str) -> ActiveTensionModelEntry:
    """
    Return the catalog entry for the named active tension model.

    Args:
        name: The active tension model name (e.g. 'GoktepeKuhl').

    Returns:
        The ActiveTensionModelEntry for that model.

    Raises:
        KeyError: If the model is not in the catalog.
    """
    if name not in ACTIVE_TENSION_MODEL_CATALOG:
        raise KeyError(
            f"Unknown active tension model '{name}'. "
            f"Available models: {', '.join(ACTIVE_TENSION_MODEL_CATALOG.keys())}"
        )
    return ACTIVE_TENSION_MODEL_CATALOG[name]
