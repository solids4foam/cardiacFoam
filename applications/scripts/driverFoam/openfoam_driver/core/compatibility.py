"""Named compatibility boundaries: what a capability returns with no hook.

These adapters preserve observable behavior. They produce no warnings and make
no policy changes. Keeping them named and documented stops legacy decisions
from being rediscovered deep inside solver-neutral code, and gives the
capability seams an explicit place at which behaviour may later change.

**Every ``legacy_*`` below is now plugin-neutral.** Nineteen of them used to
branch on ``plugin_id == "org.cardiacfoam"`` and import from
``plugins/cardiacfoam/`` -- a string comparison that let core reach into one
named plugin's package. Those branches were removed on 2026-08-26 once the
cardiac plugin implemented all fifteen optional hooks itself, and only after
instrumenting all nineteen and confirming across the full 1703-test suite that
not one was reachable. Deleted on measured evidence, not on the argument that
the "v1 plugin" population is empty.

Two functions still name the cardiac plugin, both deliberately and both
documented at their definitions: :func:`legacy_default_driver_context` encodes
the product decision that cardiacFoam is the plugin you get when you name none,
and :func:`legacy_generic_case_mutation` serves direct callers of core
``make_spec``.
"""

from __future__ import annotations

import contextvars
import functools
from contextlib import contextmanager
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from .plugin_interface import DriverContext

_fallback_call_log: contextvars.ContextVar[list[str] | None] = contextvars.ContextVar(
    "_fallback_call_log", default=None,
)


@contextmanager
def track_fallback_calls():
    """Yield a list that fills with the name of every legacy_* fallback
    invoked inside the ``with`` block, in call order. Empty means none fired
    -- the P2.4 assertion an explicit non-cardiac v2 context should satisfy."""
    token = _fallback_call_log.set([])
    try:
        yield _fallback_call_log.get()
    finally:
        _fallback_call_log.reset(token)


def _instrumented(func):
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        log = _fallback_call_log.get()
        if log is not None:
            log.append(func.__name__)
        return func(*args, **kwargs)
    return wrapper


@_instrumented
def legacy_default_driver_context() -> "DriverContext":
    """Return the historical built-in cardiacFoam context.

    Why: public CLI and Python callers have always selected cardiacFoam when no
    plugin/context is supplied.
    Activation: the public boundary receives no explicit plugin context.
    Preserved by: core plugin-context, CLI matrix, validation, and strict-plan
    tests.
    Plan 2 seam: default selection or deprecation policy may change there.
    """

    from .plugin_interface import driver_context
    from ..plugins.cardiacfoam_plugin import CardiacFoamPlugin

    target = "openfoam_driver.plugins.cardiacfoam_plugin:CardiacFoamPlugin"
    return driver_context(
        CardiacFoamPlugin(), source=f"trusted-import:{target}",
    )


def resolve_public_driver_context(
    driver_context: "DriverContext | None",
) -> "DriverContext":
    """Resolve the unchanged optional-context public API convention once."""

    return driver_context if driver_context is not None else legacy_default_driver_context()


@_instrumented
def legacy_generic_case_mutation(*args, **kwargs) -> None:
    """Preserve direct callers of the formerly cardiac-owned generic factory.

    Why: the historical public ``make_spec`` mutated the cardiac electro and
    physics dictionaries.  Activation: a caller imports core ``make_spec``
    directly rather than the solver-neutral registry alias.  Tests in the
    cardiac generic-case and template suites preserve it.  Plan 2 may replace
    this with explicit generic dictionary mutations.
    """

    from ..plugins.cardiacfoam.generic_case_mutation import apply_case_mutation

    apply_case_mutation(*args, **kwargs)


# The dictionary files the historical generic-case factory addressed, keyed by
# the generic ``dict_file_relpaths`` names.  Insertion order matters: the first
# entry is the "primary" file whose presence marks a folder as non-generic, and
# ``electroProperties`` has always been that marker.
_LEGACY_GENERIC_CASE_DICT_FILES = (
    ("electro", "constant/electroProperties"),
    ("physics", "constant/physicsProperties"),
)

# Historical cardiac-named ``make_spec`` keyword arguments, mapped onto the
# generic ``(bucket, dict-file name)`` they now address.
_LEGACY_GENERIC_CASE_ALIASES = {
    "electro_properties_relpath": ("relpaths", "electro"),
    "physics_properties_relpath": ("relpaths", "physics"),
    "electro_property_overrides": ("overrides", "electro"),
    "physics_property_overrides": ("overrides", "physics"),
}


@_instrumented
def legacy_generic_case_dict_file_relpaths() -> dict[str, str]:
    """Return the dictionary files core ``make_spec`` has always defaulted to.

    Why: the historical signature defaulted ``electro_properties_relpath`` and
    ``physics_properties_relpath`` to fixed cardiac paths, and the generic-case
    detection keyed off the first of them.  Activation: a caller of core
    ``make_spec`` declares no ``dict_file_relpaths``.  Preserved by the core
    generic-case and strict-plan suites.  Plan 2 seam: a plugin declaring its
    own dictionary files makes this default unnecessary.
    """

    return dict(_LEGACY_GENERIC_CASE_DICT_FILES)


def legacy_generic_case_alias_names() -> frozenset[str]:
    """Names :func:`legacy_generic_case_dict_file_aliases` recognises.

    Uninstrumented on purpose: callers use it to *decide* whether a legacy
    alias is present at all, so consulting it is not itself a fallback.
    """

    return frozenset(_LEGACY_GENERIC_CASE_ALIASES)


@_instrumented
def legacy_generic_case_dict_file_aliases(
    payload,
) -> tuple[dict, dict, list[str]]:
    """Translate deprecated cardiac-named generic-case keywords.

    Why: ``electro_property_overrides`` and friends are advertised as common
    override keys and reach ``make_spec`` verbatim from ``--config``/``--set``.
    Activation: any such key appears in a ``make_spec`` call or in a ``cases``
    entry.  Returns ``(relpaths, overrides, unknown_keys)`` so the caller keeps
    ownership of rejecting genuinely unknown keywords.  Plan 2 seam: the
    aliases may be dropped once callers migrate to ``dict_file_relpaths`` and
    ``dict_file_overrides``.
    """

    relpaths: dict = {}
    overrides: dict = {}
    unknown: list[str] = []
    for key, value in dict(payload).items():
        target = _LEGACY_GENERIC_CASE_ALIASES.get(key)
        if target is None:
            unknown.append(key)
            continue
        bucket, name = target
        (relpaths if bucket == "relpaths" else overrides)[name] = value
    return relpaths, overrides, unknown


@_instrumented
def legacy_case_marker(plugin, case_root) -> bool:
    """Whether a case folder belongs to this plugin, for a plugin that does not
    implement ``has_case_marker()``.

    ``False``: absent filesystem evidence of its own, a plugin claims nothing.
    Core then relies on an executable ``Allrun``, which is plugin-neutral."""

    del plugin, case_root

    return False


@_instrumented
def legacy_case_runnable_without_workflow(plugin, case_root) -> bool:
    """Whether an uncontracted case is runnable, for a plugin that does not
    implement ``is_case_runnable_without_workflow()``.

    ``False``: core falls back to an executable ``Allrun``, which is
    plugin-neutral filesystem evidence rather than a guess."""

    del plugin, case_root

    return False


@_instrumented
def legacy_run_document_config(plugin, spec):
    """The RunDocument ``config`` object for a plugin that does not implement
    ``build_run_document_config()``.

    An empty config and no diagnostics -- the plugin constrains nothing,
    matching the fully open schema :func:`legacy_run_document_config_schema`
    hands it. RunDocument v3 removed any fixed phase vocabulary from core, so
    there is no shape to supply on a plugin's behalf."""

    del plugin, spec

    return {}, ()


@_instrumented
def legacy_run_document_config_schema(plugin) -> dict:
    """The RunDocument ``config`` JSON Schema for a plugin that does not
    implement ``get_run_document_config_schema()``.

    A fully open object: this plugin constrains nothing."""

    del plugin

    return {"type": "object", "additionalProperties": True}


@_instrumented
def legacy_nondimensional_case(plugin, spec) -> bool:
    """Whether SI mesh-scale diagnostics should be skipped, for a plugin that
    does not implement ``is_nondimensional_case()``.

    ``False``, which keeps the diagnostics on. Silence is not the safe
    default; a plugin that wants them skipped says so."""

    del plugin, spec

    return False


@_instrumented
def legacy_route_sweep_case(plugin, *, base, resolved_axis_values, driver_context):
    """Plugins predating route_sweep_case_values().

    Unlike every other fallback here, a neutral empty return is not available:
    routing produces the values a case is then materialized from, so an empty
    routing silently yields a case that is not the one the sweep asked for.
    The honest neutral is to refuse, naming the hook the plugin must
    implement.

    Historically this delegated to the cardiac router, which validates axes
    against ``electroProperties``/``physicsProperties`` vocabulary -- so
    ungated it rejected a non-cardiac plugin's axes in cardiac terms, or,
    worse, accepted them."""

    from openfoam_driver.core.sweep.sweep_expansion import SweepValidationError

    raise SweepValidationError(
        f"plugin {getattr(plugin, 'plugin_id', '<unknown>')!r} does not implement "
        "route_sweep_case_values(); driverFOAM cannot route sweep axes for it. "
        "Implement route_sweep_case_values(base, resolved_axis_values, "
        "driver_context) on the plugin to support sweeps."
    )


@_instrumented
def legacy_materialize_sweep_case(plugin, *, case_dir, routed) -> None:
    """Plugins predating materialize_sweep_case(). Refuses for the same
    reason as :func:`legacy_route_sweep_case`.

    This is the fallback with real teeth. It used to delegate to the cardiac
    materializer, which writes an ``Allrun`` containing a hardcoded
    ``cardiacFoam`` command -- so ungated it generated a case invoking the
    cardiacFoam binary under whichever plugin was loaded (reproduced against
    GenericOpenFOAMPlugin, 2026-08-19)."""

    from openfoam_driver.core.sweep.sweep_expansion import SweepValidationError

    raise SweepValidationError(
        f"plugin {getattr(plugin, 'plugin_id', '<unknown>')!r} does not implement "
        "materialize_sweep_case(); driverFOAM cannot materialize sweep cases "
        "for it. Implement materialize_sweep_case(case_dir, routed) on the "
        "plugin to support sweeps."
    )


@_instrumented
def legacy_solver_commands(plugin) -> frozenset[str]:
    """Solver binaries for a plugin that does not implement
    ``get_solver_commands()``.

    Empty: core authorizes no binary it was not told about."""

    del plugin

    return frozenset()


@_instrumented
def legacy_auxiliary_commands(plugin) -> frozenset[str]:
    """Auxiliary binaries for a plugin that does not implement
    ``get_auxiliary_commands()``.

    Empty, for the same reason as :func:`legacy_solver_commands`."""

    del plugin

    return frozenset()


@_instrumented
def legacy_utility_manifests(plugin) -> dict:
    """Per-utility declarations for a plugin that does not implement
    ``get_utility_manifests()``. Empty: no utility is pre-authorized."""

    del plugin

    return {}


@_instrumented
def legacy_utility_roots(plugin) -> tuple:
    """Utility source roots for a plugin that does not implement
    ``get_utility_roots()``. Empty: nothing extra to fingerprint."""

    del plugin

    return ()


@_instrumented
def legacy_resolve_case_models(plugin, case_root) -> dict:
    """On-disk model selections for a plugin that does not implement
    ``resolve_case_models()``. Empty: core reads no dictionary on its own."""

    del plugin, case_root

    return {}


@_instrumented
def legacy_samplable_fields(plugin, resolved) -> dict:
    """Samplable fields for a plugin that does not implement
    ``get_samplable_fields()``. Empty: core invents no field names."""

    del plugin, resolved

    return {}


@_instrumented
def legacy_override_schema(plugin, tutorial_name: str, make_spec_info: dict) -> dict:
    """The ``--config`` override schema for a plugin that does not implement
    ``get_override_schema()``. Empty: nothing is advertised as overridable."""

    del plugin, tutorial_name, make_spec_info

    return {}


@_instrumented
def legacy_dict_entry_catalog(plugin) -> dict:
    """Dictionary entries by document name, for a plugin that does not
    implement ``get_dict_entry_catalog()``. Empty: core owns no vocabulary."""

    del plugin

    return {}


@_instrumented
def legacy_describe_config_resolution(plugin) -> str:
    """Plugins that do not implement get_config_resolution_description() get a
    plugin-neutral sentence. The built-in cardiac plugin now implements the hook
    itself, so no plugin is named here."""

    del plugin
    return "The plugin's configuration files resolve into a valid RunDocument config."


@_instrumented
def legacy_report_catalog(plugin) -> tuple:
    """Plugins that do not implement get_report_catalog() have no post-run
    reports. The built-in cardiac plugin now implements the hook itself, so core
    no longer imports any plugin's report module."""

    del plugin
    return ()


@_instrumented
def legacy_named_catalogs(plugin) -> dict:
    """Plugin-chosen catalogs for a plugin that does not implement
    ``get_named_catalogs()``. Empty: ``describe`` namespaces nothing extra."""

    del plugin

    return {}


@_instrumented
def legacy_override_scopes(plugin) -> tuple:
    """``$TOKEN.`` override scopes for a plugin that does not implement
    ``get_override_scopes()``. Empty: no scoped patch targets."""

    del plugin

    return ()


@_instrumented
def legacy_dict_regeneration_scopes(plugin) -> tuple:
    """Regeneration scopes for a plugin that does not implement
    ``get_regeneration_scopes()``. Empty: no override regenerates a dict."""

    del plugin

    return ()


@_instrumented
def legacy_phases(plugin) -> tuple[str, ...]:
    """The dictionary phases for a plugin that does not implement
    ``get_phases()``: those its own ``DictEntry`` values declare, sorted for
    determinism.

    Sorted, not ordered -- and the order is the semantics, since
    ``primary_phase()`` returns the first phase in it that an entry claims. A
    plugin with multi-phase entries should implement ``get_phases()`` rather
    than accept an alphabetical guess. What this must never do is hand back
    cardiacFoam's four to a plugin that never declared them: that was the
    silent defect this replaces."""

    declared: set[str] = set()
    for entry in plugin.get_dict_entries():
        declared.update(entry.phases)
    return tuple(sorted(declared))
