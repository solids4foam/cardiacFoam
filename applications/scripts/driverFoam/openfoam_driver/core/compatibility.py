"""Named Plan-1 compatibility boundaries.

These adapters intentionally preserve observable behavior.  They produce no
warnings and make no policy changes.  Keeping them named and documented stops
legacy decisions from being rediscovered deep inside solver-neutral code and
gives Plan 2 explicit seams at which behaviour may later change.
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
    """Plugins predating has_case_marker(). Only the built-in cardiac plugin
    has authored filesystem evidence for its own cases; every other plugin
    gets ``False`` and must declare its own marker.

    Gating matters here even though ``False`` is also what the cardiac rule
    returns for a case with no ``constant/electroProperties``: ungated, the
    *reason* a non-cardiac case was rejected was that it failed a cardiac
    test, so a non-cardiac case that happened to carry an
    ``electroProperties`` file was claimed by whichever plugin was loaded."""

    if getattr(plugin, "plugin_id", "") == "org.cardiacfoam":
        from ..plugins.cardiacfoam.case_compatibility import has_case_marker

        return has_case_marker(case_root)
    return False


@_instrumented
def legacy_case_runnable_without_workflow(plugin, case_root) -> bool:
    """Plugins predating is_case_runnable_without_workflow(). Same rule as
    :func:`legacy_case_marker`: only the built-in cardiac plugin can judge an
    uncontracted case runnable, because the judgement reads cardiac
    dictionaries. Others get ``False`` -- core then falls back to an
    executable ``Allrun``, which is plugin-neutral filesystem evidence."""

    if getattr(plugin, "plugin_id", "") == "org.cardiacfoam":
        from ..plugins.cardiacfoam.case_compatibility import is_runnable_without_workflow

        return is_runnable_without_workflow(case_root)
    return False


@_instrumented
def legacy_run_document_config(plugin, spec):
    """Plugins predating build_run_document_config(). Only the built-in
    cardiac plugin has an authored RunDocument config builder; others get an
    empty config and no diagnostics -- they constrain nothing, exactly as
    :func:`legacy_run_document_config_schema` hands them a fully open schema.

    The pre-gate return for a non-cardiac plugin was the cardiac *phase*
    vocabulary (``anatomy``/``physics``/``stimulus``/``solver``). That
    vocabulary is precisely what RunDocument v3 removed from core, where
    ``config`` is an open object with no fixed phases, so returning it for a
    plugin that never declared those phases contradicts the schema."""

    if getattr(plugin, "plugin_id", "") == "org.cardiacfoam":
        from ..plugins.cardiacfoam.run_document_config import build_config

        return build_config(spec)
    return {}, ()


@_instrumented
def legacy_run_document_config_schema(plugin) -> dict:
    """v1 plugins predate get_run_document_config_schema(). Only the built-in
    cardiac plugin has an authored config schema; other v1 plugins get a fully
    open schema (no constraint) and must declare their own by migrating to v2."""

    if getattr(plugin, "plugin_id", "") == "org.cardiacfoam":
        from ..plugins.cardiacfoam.config_schema import get_run_document_config_schema

        return get_run_document_config_schema()
    return {"type": "object", "additionalProperties": True}


@_instrumented
def legacy_nondimensional_case(plugin, spec) -> bool:
    """Plugins predating is_nondimensional_case(). The cardiac exemption is
    read out of ``constant/electroProperties`` (a singleCell or verification
    model), so only the built-in cardiac plugin can answer it. Others get
    ``False``: their meshes are dimensional until they say otherwise, which
    is the conservative answer -- it keeps mesh-scale diagnostics ON rather
    than silently exempting a case from them."""

    if getattr(plugin, "plugin_id", "") == "org.cardiacfoam":
        from ..plugins.cardiacfoam.planning_policy import is_nondimensional_case

        return is_nondimensional_case(spec)
    return False


@_instrumented
def legacy_route_sweep_case(plugin, *, base, resolved_axis_values, driver_context):
    """Plugins predating route_sweep_case_values().

    Unlike every other fallback here, a neutral empty return is not available:
    routing produces the values a case is then materialized from, so an empty
    routing silently yields a case that is not the one the sweep asked for.
    The honest neutral is to refuse, naming the hook the plugin must
    implement.

    The cardiac router validates axes against ``electroProperties``/
    ``physicsProperties`` vocabulary, so ungated it rejected a non-cardiac
    plugin's axes in cardiac terms -- or, worse, accepted them."""

    if getattr(plugin, "plugin_id", "") == "org.cardiacfoam":
        from ..plugins.cardiacfoam.sweep import route_case_values

        return route_case_values(
            base=base,
            resolved_axis_values=resolved_axis_values,
            driver_context=driver_context,
        )
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

    This is the fallback with real teeth. The cardiac materializer writes an
    ``Allrun`` containing a hardcoded ``cardiacFoam`` command, so ungated it
    generated a case invoking the cardiacFoam binary under whichever plugin
    was loaded (reproduced against GenericOpenFOAMPlugin, 2026-08-19)."""

    if getattr(plugin, "plugin_id", "") == "org.cardiacfoam":
        from ..plugins.cardiacfoam.sweep import materialize_case

        materialize_case(case_dir=case_dir, routed=routed)
        return
    from openfoam_driver.core.sweep.sweep_expansion import SweepValidationError

    raise SweepValidationError(
        f"plugin {getattr(plugin, 'plugin_id', '<unknown>')!r} does not implement "
        "materialize_sweep_case(); driverFOAM cannot materialize sweep cases "
        "for it. Implement materialize_sweep_case(case_dir, routed) on the "
        "plugin to support sweeps."
    )


@_instrumented
def legacy_solver_commands(plugin) -> frozenset[str]:
    """v1 plugins predate get_solver_commands(). Only the built-in cardiac
    plugin can be given a solver name; a third-party v1 plugin gets none and
    must declare its commands by migrating to v2."""

    if getattr(plugin, "plugin_id", "") == "org.cardiacfoam":
        from ..plugins.cardiacfoam.command_authorization import solver_commands

        return solver_commands()
    return frozenset()


@_instrumented
def legacy_auxiliary_commands(plugin) -> frozenset[str]:
    """v1 plugins predate get_auxiliary_commands(). Same rule as
    :func:`legacy_solver_commands`: only the built-in cardiac plugin gets its
    non-solver commands authorized."""

    if getattr(plugin, "plugin_id", "") == "org.cardiacfoam":
        from ..plugins.cardiacfoam.command_authorization import auxiliary_commands

        return auxiliary_commands()
    return frozenset()


@_instrumented
def legacy_utility_manifests(plugin) -> dict:
    """Preserve the cardiac utility catalog for plugins without the new hook."""

    if getattr(plugin, "plugin_id", "") == "org.cardiacfoam":
        from ..plugins.cardiacfoam.command_authorization import utility_manifests

        # Cached read-only view; copy so a caller cannot reach the shared cache.
        return dict(utility_manifests())
    return {}


@_instrumented
def legacy_utility_roots(plugin) -> tuple:
    """Preserve the cardiac utilities root for plugins without the new hook."""

    if getattr(plugin, "plugin_id", "") == "org.cardiacfoam":
        from ..plugins.cardiacfoam.command_authorization import utility_roots

        return utility_roots()
    return ()


@_instrumented
def legacy_resolve_case_models(plugin, case_root) -> dict:
    """v1 plugins predate resolve_case_models(). Only the built-in cardiac
    plugin can resolve a case's models; other v1 plugins get nothing and must
    declare their own resolution by migrating to v2."""

    if getattr(plugin, "plugin_id", "") == "org.cardiacfoam":
        from ..plugins.cardiacfoam.case_introspection import resolve_case_models

        return resolve_case_models(case_root)
    return {}


@_instrumented
def legacy_samplable_fields(plugin, resolved) -> dict:
    """v1 plugins predate get_samplable_fields(). Same rule as
    :func:`legacy_resolve_case_models`: only the built-in cardiac plugin
    names any fields."""

    if getattr(plugin, "plugin_id", "") == "org.cardiacfoam":
        from ..plugins.cardiacfoam.case_introspection import samplable_fields

        return samplable_fields(resolved)
    return {}


@_instrumented
def legacy_override_schema(plugin, tutorial_name: str, make_spec_info: dict) -> dict:
    """v1 plugins predate get_override_schema(). Only the built-in cardiac
    plugin has an authored configuration vocabulary; other v1 plugins get an
    empty schema and must declare their own by migrating to v2."""

    if getattr(plugin, "plugin_id", "") == "org.cardiacfoam":
        from ..plugins.cardiacfoam.override_schema import config_schema

        return config_schema(tutorial_name, make_spec_info)
    return {}


@_instrumented
def legacy_dict_entry_catalog(plugin) -> dict:
    """v1 plugins predate get_dict_entry_catalog(). Same rule as
    :func:`legacy_override_schema`: only the built-in cardiac plugin knows the
    electro/physics document shape."""

    if getattr(plugin, "plugin_id", "") == "org.cardiacfoam":
        from ..plugins.cardiacfoam.override_schema import dict_entry_catalog

        return dict_entry_catalog(
            plugin.get_dictionary_catalog(), plugin.get_dict_groups(),
        )
    return {}


@_instrumented
def legacy_describe_config_resolution(plugin) -> str:
    """v1 plugins predate describe_config_resolution(). Only the built-in
    cardiac plugin has an authored description; other v1 plugins get a
    plugin-neutral sentence."""

    if getattr(plugin, "plugin_id", "") == "org.cardiacfoam":
        return "physicsProperties and electroProperties resolve into a valid RunDocument config."
    return "The plugin's configuration files resolve into a valid RunDocument config."


@_instrumented
def legacy_report_catalog(plugin) -> tuple:
    """v1 plugins predate get_report_catalog(). Same rule as
    :func:`legacy_override_schema`: only the built-in cardiac plugin has an
    authored post-run report catalog; other v1 plugins get no reports and
    must declare their own by migrating to v2."""

    if getattr(plugin, "plugin_id", "") == "org.cardiacfoam":
        from ..plugins.cardiacfoam.reports import CARDIAC_REPORTS

        return CARDIAC_REPORTS
    return ()


@_instrumented
def legacy_named_catalogs(plugin) -> dict:
    """v1 plugins predate get_named_catalogs(). Same rule as
    :func:`legacy_override_schema`: only the built-in cardiac plugin has
    ionic-model/active-tension catalogs; other v1 plugins get none and must
    declare their own by migrating to v2."""

    if getattr(plugin, "plugin_id", "") == "org.cardiacfoam":
        from ..plugins.cardiacfoam.named_catalogs import named_catalogs

        return named_catalogs(plugin.get_capabilities())
    return {}


@_instrumented
def legacy_override_scopes(plugin) -> tuple:
    """v1/v2 plugins predate get_override_scopes(). Only the built-in
    cardiac plugin has an authored override scope ($ELECTRO_MODEL_COEFFS);
    other plugins get none and must declare their own by implementing
    get_override_scopes()."""

    if getattr(plugin, "plugin_id", "") == "org.cardiacfoam":
        from ..plugins.cardiacfoam.overrides import electro_model_coeffs_scope

        return (electro_model_coeffs_scope(),)
    return ()


@_instrumented
def legacy_dict_regeneration_scopes(plugin) -> tuple:
    """v1/v2 plugins predate get_regeneration_scopes(). Only the built-in
    cardiac plugin has an authored regeneration scope (myocardiumSolver ->
    constant/electroProperties); other plugins get none and must declare
    their own by implementing get_regeneration_scopes()."""

    if getattr(plugin, "plugin_id", "") == "org.cardiacfoam":
        from ..plugins.cardiacfoam.overrides import electro_properties_regeneration_scope

        return (electro_properties_regeneration_scope(),)
    return ()
