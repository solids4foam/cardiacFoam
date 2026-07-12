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
#     active_tension_catalog
#
# Description
#     Manages metadata and parameter constraints for active tension models.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

"""
Active Tension Model Catalog

A static, build-time-generated catalog of active tension models.
This module exposes the exact variables each model supports so an autonomous
agent can plan active-tension output without running the solver.

All variable names are extracted from C++ source files and are guaranteed to be
exact.

The ``constants`` field lists only the *user-facing*, dict-overridable parameters.
Constants that a model *derives* from these at ``initConsts`` (e.g. the
Land-Niederer transition rates ``AC_k_uw``/``AC_k_ws``/``AC_k_wu``/``AC_k_su``,
``AC_cds``/``AC_cdw``, ``AC_ktm_block``, ``AC_A``, ``AC_XSSS``, ``AC_XWSS``,
``AC_fPKA_TnI``, ``AC_PKAForceMultiplier``) are intentionally omitted: overriding
them has no effect because the solver recomputes them from their inputs. An agent
may still *reason about* such derived values, but to change one it must override
the user-facing constant(s) it derives from.
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
    "GoktepeKuhlBatched": ActiveTensionModelEntry(
        states=("Ta",),
        algebraic=("AV_e", "AV_Vm", "AV_u"),
        constants=("AC_Vr", "AC_eInfty", "AC_e0", "AC_eXi", "AC_Vshift", "AC_kTa"),
        rates=("Ta",),
        recommended_exports=("Ta",),
        description="Goktepe-Kuhl active tension model (2004) - GPU batched implementation.",
        aliases=("Goktepe-Kuhl GPU",),
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
    "NashPanfilovBatched": ActiveTensionModelEntry(
        states=("Ta",),
        algebraic=("AV_u", "AV_e"),
        constants=("AC_Vp", "AC_Vr", "AC_Vth", "AC_e0", "AC_kTa", "AC_a"),
        rates=("Ta",),
        recommended_exports=("Ta",),
        description="Nash-Panfilov active tension model (2004) - GPU batched implementation.",
        aliases=("Nash-Panfilov GPU",),
    ),
    "LandNiederer": ActiveTensionModelEntry(
        states=("Ca_TRPN", "TmBlocked", "XW", "XS", "ZETAS", "ZETAW"),
        algebraic=(
            "AV_XU", "AV_gamma_rate", "AV_gamma_rate_w", "AV_xb_uw", "AV_xb_ws",
            "AV_xb_su", "AV_xb_wu", "AV_xb_su_gamma", "AV_xb_wu_gamma", "AV_Ta",
            "AV_d_Ca_TRPN", "AV_Cai", "AV_lambda", "AV_lambda_rate", "AV_Lfac",
            "AV_ca50",
        ),
        # 21 user-facing (dict-overridable) parameters; 12 further derived
        # constants (AC_A, AC_XSSS, AC_XWSS, AC_fPKA_TnI, AC_k_uw, AC_k_ws,
        # AC_k_wu, AC_k_su, AC_cds, AC_cdw, AC_ktm_block, AC_PKAForceMultiplier)
        # are computed in initConsts and are not set via the case dict.
        constants=(
            "AC_TOT_A", "AC_TRPN_n", "AC_Tref", "AC_beta_0", "AC_beta_1", "AC_dr",
            "AC_fracTnIpo", "AC_contraction_gamma", "AC_gamma_wu", "AC_koff",
            "AC_ktm_unblock", "AC_lambda_max", "AC_lambda_min", "AC_mu", "AC_nperm",
            "AC_nu", "AC_perm50", "AC_phi", "AC_wfrac", "AC_fMyBPC_PKA", "AC_fTnI_PKA",
        ),
        rates=("Ca_TRPN", "TmBlocked", "XW", "XS", "ZETAS", "ZETAW"),
        recommended_exports=("Ta",),
        description=(
            "Land-Niederer myofilament contraction model (2017): 6-state "
            "crossbridge/tropomyosin system driven by intracellular Ca and "
            "sarcomere stretch."
        ),
        notes=(
            "Length- and velocity-dependent active tension (Frank-Starling); "
            "requires AV_Cai, AV_lambda, AV_lambda_rate to be supplied by the "
            "caller. Dict also accepts 'preconditioningTime' (default 1000 ms)."
        ),
        aliases=("Land-Niederer", "Land2017"),
    ),
    "LandNiedererBatched": ActiveTensionModelEntry(
        states=("Ca_TRPN", "TmBlocked", "XW", "XS", "ZETAS", "ZETAW"),
        algebraic=(
            "AV_XU", "AV_gamma_rate", "AV_gamma_rate_w", "AV_xb_uw", "AV_xb_ws",
            "AV_xb_su", "AV_xb_wu", "AV_xb_su_gamma", "AV_xb_wu_gamma", "AV_Ta",
            "AV_d_Ca_TRPN", "AV_Cai", "AV_lambda", "AV_lambda_rate", "AV_Lfac",
            "AV_ca50",
        ),
        constants=(
            "AC_TOT_A", "AC_TRPN_n", "AC_Tref", "AC_beta_0", "AC_beta_1", "AC_dr",
            "AC_fracTnIpo", "AC_contraction_gamma", "AC_gamma_wu", "AC_koff",
            "AC_ktm_unblock", "AC_lambda_max", "AC_lambda_min", "AC_mu", "AC_nperm",
            "AC_nu", "AC_perm50", "AC_phi", "AC_wfrac", "AC_fMyBPC_PKA", "AC_fTnI_PKA",
        ),
        rates=("Ca_TRPN", "TmBlocked", "XW", "XS", "ZETAS", "ZETAW"),
        recommended_exports=("Ta",),
        description=(
            "Land-Niederer myofilament contraction model (2017) - GPU batched "
            "implementation."
        ),
        notes=(
            "GPU batched variant of LandNiederer; same states/constants. Dict "
            "also accepts 'preconditioningTime' (default 1000 ms)."
        ),
        aliases=("Land-Niederer GPU", "Land2017 GPU"),
    ),
    "ManufacturedElectromechanics": ActiveTensionModelEntry(
        states=("Ta",),
        algebraic=("AV_Vm", "AV_lambda", "AV_voltageActivation", "AV_lengthFactor"),
        constants=("AC_Tmax", "AC_V0", "AC_gamma"),
        rates=("Ta",),
        recommended_exports=("Ta",),
        description="Manufactured active-tension model for electromechanics verification.",
        notes="Verification-only model used by the total-Lagrangian electromechanics MMS.",
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
