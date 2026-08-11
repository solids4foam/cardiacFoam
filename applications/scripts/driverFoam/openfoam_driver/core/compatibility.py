"""Named Plan-1 compatibility boundaries.

These adapters intentionally preserve observable behavior.  They produce no
warnings and make no policy changes.  Keeping them named and documented stops
legacy decisions from being rediscovered deep inside solver-neutral code and
gives Plan 2 explicit seams at which behaviour may later change.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from .plugin_interface import DriverContext


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


def legacy_generic_case_mutation(*args, **kwargs) -> None:
    """Preserve direct callers of the formerly cardiac-owned generic factory.

    Why: the historical public ``make_spec`` accepted electro/physics override
    arguments.  Activation: a caller imports core ``make_spec`` directly rather
    than the solver-neutral registry alias.  Tests in the cardiac generic-case
    and template suites preserve it.  Plan 2 may replace this with explicit
    generic dictionary mutations.
    """

    from ..plugins.cardiacfoam.generic_case_mutation import apply_case_mutation

    apply_case_mutation(*args, **kwargs)


def legacy_case_marker(case_root) -> bool:
    """Preserve cardiac filesystem evidence for plugins without the new hook."""

    from ..plugins.cardiacfoam.case_compatibility import has_case_marker

    return has_case_marker(case_root)


def legacy_case_runnable_without_workflow(case_root) -> bool:
    """Preserve historical uncontracted-case runnability for legacy plugins."""

    from ..plugins.cardiacfoam.case_compatibility import is_runnable_without_workflow

    return is_runnable_without_workflow(case_root)


def legacy_run_document_config(spec):
    """Preserve the cardiac-shaped RunDocument-v2 parser for legacy plugins."""

    from ..plugins.cardiacfoam.run_document_config import build_config

    return build_config(spec)


def legacy_nondimensional_case(spec) -> bool:
    """Preserve cardiac mesh-diagnostic exemptions for legacy plugins."""

    from ..plugins.cardiacfoam.planning_policy import is_nondimensional_case

    return is_nondimensional_case(spec)


def legacy_route_sweep_case(*, base, resolved_axis_values, driver_context):
    """Preserve the cardiac-shaped generic sweep router for legacy plugins."""

    from ..plugins.cardiacfoam.sweep import route_case_values

    return route_case_values(
        base=base,
        resolved_axis_values=resolved_axis_values,
        driver_context=driver_context,
    )


def legacy_materialize_sweep_case(*, case_dir, routed) -> None:
    """Preserve build_and_launch-based sweep materialization."""

    from ..plugins.cardiacfoam.sweep import materialize_case

    materialize_case(case_dir=case_dir, routed=routed)


def legacy_solver_commands(plugin) -> frozenset[str]:
    """v1 plugins predate get_solver_commands(). Only the built-in cardiac
    plugin can be given a solver name; a third-party v1 plugin gets none and
    must declare its commands by migrating to v2."""

    if getattr(plugin, "plugin_id", "") == "org.cardiacfoam":
        from ..plugins.cardiacfoam.command_authorization import solver_commands

        return solver_commands()
    return frozenset()


def legacy_utility_manifests(plugin) -> dict:
    """Preserve the cardiac utility catalog for plugins without the new hook."""

    if getattr(plugin, "plugin_id", "") == "org.cardiacfoam":
        from ..plugins.cardiacfoam.command_authorization import utility_manifests

        return utility_manifests()
    return {}


def legacy_utility_roots(plugin) -> tuple:
    """Preserve the cardiac utilities root for plugins without the new hook."""

    if getattr(plugin, "plugin_id", "") == "org.cardiacfoam":
        from ..plugins.cardiacfoam.command_authorization import utility_roots

        return utility_roots()
    return ()
