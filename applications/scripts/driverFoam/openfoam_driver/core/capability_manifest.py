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
#     capability_manifest
#
# Description
#     Assembles the machine-readable surface of what the driver will accept:
#     the allowed workflow-command set and the field names each solver can
#     sample. Pure assembly from explicit inputs supplied by the calling
#     plugin (commands, utility manifests, samplable fields) plus the
#     core-neutral workflow allowlist, so the manifest cannot drift from the
#     enforcers that actually gate execution and this module holds no solver
#     knowledge of its own.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

from __future__ import annotations

from pathlib import Path
from typing import Any, Iterable

from .runtime.workflow import (
    CASE_SCRIPT_COMMANDS,
    CORE_NEUTRAL_COMMANDS,
)


def resolve_case_models(
    case_root: str | Path,
) -> tuple[str | None, str | None, str | None]:
    """Deprecated 3-tuple shim over the cardiac plugin's
    ``CaseIntrospectionCapability.resolve_case_models``. ``(solver,
    ionic_model, active_tension)``, any of which may be ``None``; never
    raises. Kept for callers that imported this function directly before it
    became plugin-owned; new code should go through
    ``driver_context.capabilities.case_introspection`` instead."""
    from .plugin_interface import default_driver_context

    introspection = default_driver_context().capabilities.case_introspection
    resolved = introspection.resolve_case_models(Path(case_root))
    return resolved.get("solver"), resolved.get("ionic_model"), resolved.get("active_tension")


def utility_produces(
    utility_manifests: dict[str, Any],
) -> dict[str, tuple[str, ...]]:
    """The plugin utility commands the allowlist accepts, keyed to what they
    produce. Mirrors the acceptance rule in ``validate_workflow_commands``
    (only utilities that declare ``produces`` are accepted).

    Single owner: the advertised accept-surface (this module's manifest) and
    the enforced one (strict planning's artifact coverage) must not drift, so
    both derive from this one function rather than each keeping a copy.
    """

    return {
        command: tuple(produce.artifact_id for produce in manifest.produces)
        for command, manifest in utility_manifests.items()
        if manifest.produces
    }


def _utility_commands(utility_manifests: dict[str, Any]) -> dict[str, list[str]]:
    """JSON-shaped view of :func:`utility_produces` (lists, not tuples)."""

    return {
        command: list(artifacts)
        for command, artifacts in utility_produces(utility_manifests).items()
    }


def build_capability_manifest(
    *,
    plugin_commands: Iterable[str] = (),
    utility_manifests: dict[str, Any] | None = None,
    samplable_fields: dict[str, tuple[str, ...]] | None = None,
) -> dict[str, Any]:
    """Return the driver's accept-surface as a plain JSON-able dict.

    ``plugin_commands``, ``utility_manifests``, and ``samplable_fields`` are
    all supplied by the calling plugin so that core names no solver here;
    together with :data:`CORE_NEUTRAL_COMMANDS` the commands reproduce
    exactly what ``validate_workflow_commands`` accepts for that plugin.

    ``allowed_commands`` names exactly what a workflow DAG step may invoke;
    ``samplable_fields`` names the fields a function object may sample for the
    resolved model, split by region -- resolving the model and naming its
    fields is entirely the plugin's ``CaseIntrospectionCapability``, not this
    module's concern. A caller with nothing resolved passes ``None`` and gets
    an empty field set rather than this function raising.
    """

    fields = samplable_fields or {}

    return {
        "allowed_commands": {
            "core": sorted(set(CORE_NEUTRAL_COMMANDS) | set(plugin_commands)),
            "case_scripts": sorted(CASE_SCRIPT_COMMANDS),
            "utilities": _utility_commands(utility_manifests or {}),
            "installed_openfoam_apps_note": (
                "When OpenFOAM is sourced, any executable under $FOAM_APPBIN or "
                "$FOAM_USER_APPBIN is also accepted (core apps + your compiled "
                "utilities). Unsourced, only core + case_scripts + utilities apply."
            ),
        },
        "samplable_fields": {
            **{region: sorted(names) for region, names in fields.items()},
            "note": (
                "Function objects are OpenFOAM's; these are the field NAMES this "
                "solver exposes. Sampling a name not listed here is silently "
                "dropped by the solver (strict planning emits unknown_sampled_field "
                "warnings for such names)."
            ),
        },
    }
