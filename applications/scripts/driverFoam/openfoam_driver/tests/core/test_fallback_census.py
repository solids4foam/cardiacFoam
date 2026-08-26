"""An explicitly-contexted operation must never fall back to the cardiac default.

Before Phase 1, `legacy_default_driver_context` fired 51,540 times across one
full test run: core resolved its own plugin context whenever a caller omitted
one, and what it resolved to was always cardiacFoam. With `driver_context` now a
required parameter throughout `core/`, an operation driven by a named plugin must
fire it zero times.

The static guard in `test_core_context_is_explicit.py` proves core contains no
implicit resolution *syntactically*. This proves the runtime consequence.
"""
from __future__ import annotations

import collections
import os
import traceback

from openfoam_driver.core import compatibility
from openfoam_driver.core.generic_plugin import GenericOpenFOAMPlugin
from openfoam_driver.core.plugin_interface import (
    driver_context,
    generic_openfoam_context,
)
from openfoam_driver.plugins.cardiacfoam_plugin import CardiacFoamPlugin


def assert_no_default_context_fallback(operation) -> None:
    """Run `operation()` and assert it never resolved an implicit context."""
    with compatibility.track_fallback_calls() as calls:
        operation()
        fired = [name for name in calls if name == "legacy_default_driver_context"]
    assert fired == [], (
        f"operation resolved the built-in cardiac context {len(fired)} time(s); "
        "it should use the DriverContext it was given"
    )


def test_capability_reads_under_an_explicit_cardiac_context_use_no_default() -> None:
    ctx = driver_context(CardiacFoamPlugin(), source="test")

    def op() -> None:
        caps = ctx.capabilities
        caps.dictionaries.entries()
        caps.dictionaries.groups()
        caps.manifest.manifest()
        caps.tutorials.displays()

    assert_no_default_context_fallback(op)


def test_capability_reads_under_an_explicit_generic_context_use_no_default() -> None:
    ctx = driver_context(GenericOpenFOAMPlugin(), source="test")

    def op() -> None:
        caps = ctx.capabilities
        caps.dictionaries.entries()
        caps.dictionaries.groups()
        caps.manifest.manifest()
        caps.tutorials.displays()

    assert_no_default_context_fallback(op)


def _default_context_callers(operation) -> collections.Counter:
    """Run `operation()`, returning a count of who resolved the cardiac default.

    Keyed by the first frame outside ``compatibility.py`` -- i.e. the module that
    actually asked for a context it did not have.
    """
    callers: collections.Counter = collections.Counter()
    original = compatibility.legacy_default_driver_context

    def traced(*args, **kwargs):
        for frame in reversed(traceback.extract_stack()[:-1]):
            if "compatibility.py" not in frame.filename:
                callers[frame.filename] += 1
                break
        return original(*args, **kwargs)

    globals_ = compatibility.resolve_public_driver_context.__globals__
    globals_["legacy_default_driver_context"] = traced
    try:
        operation()
    finally:
        globals_["legacy_default_driver_context"] = original
    return callers


def test_a_generic_plan_never_resolves_the_cardiac_default(tmp_path) -> None:
    """The guarantee Phase 1 actually makes, and the regression that matters.

    A plan driven by a non-cardiac plugin must never reach for the built-in
    cardiac context -- not once. This is the end-to-end version of the static
    guard in test_core_context_is_explicit.py.
    """
    from openfoam_driver.core.strict_planning import strict_plan

    case_root = tmp_path / "plainOpenFoamCase"
    case_root.mkdir()
    (case_root / "Allrun").write_text("#!/bin/sh\nexit 0\n")

    callers = _default_context_callers(
        lambda: strict_plan(
            "plainOpenFoamCase",
            overrides={"tutorials_root": str(tmp_path)},
            driver_context=generic_openfoam_context(),
        )
    )
    assert callers == collections.Counter(), (
        "a generic plan resolved the built-in cardiac context from: "
        f"{sorted(callers)}"
    )


def test_no_core_module_resolves_the_cardiac_default_during_a_cardiac_plan() -> None:
    """A cardiac plan still resolves the default, but only from cardiac's own code.

    `plugins/cardiacfoam/dict_builder.py` and `validation.py` reach their own
    catalogs through the solver-neutral helpers in `openfoam_driver/dict_entries.py`,
    which default to cardiac. That yields the right answer by an ugly route, and
    it is confined to the plugin -- it is NOT core guessing its plugin.

    This test pins that confinement: the moment a module under `core/` appears in
    this list, Phase 1 has regressed.
    """
    from openfoam_driver.core.strict_planning import strict_plan

    ctx = driver_context(CardiacFoamPlugin(), source="test")
    callers = _default_context_callers(
        lambda: strict_plan("singleCell", driver_context=ctx)
    )

    core_offenders = sorted(
        path for path in callers
        if f"openfoam_driver{os.sep}core{os.sep}" in path
    )
    assert core_offenders == [], (
        "core modules resolved an implicit cardiac context during a plan: "
        f"{core_offenders}"
    )
