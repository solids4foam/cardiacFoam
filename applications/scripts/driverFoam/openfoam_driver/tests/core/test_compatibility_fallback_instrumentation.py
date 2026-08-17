"""P2.4: compatibility fallbacks must be individually observable by tests,
and an explicit non-cardiac v2 context must call none of them."""
from __future__ import annotations

from openfoam_driver.core import compatibility


def test_recorder_captures_a_fallback_call() -> None:
    with compatibility.track_fallback_calls() as calls:
        compatibility.legacy_resolve_case_models(object(), case_root=None)
    assert calls == ["legacy_resolve_case_models"]


def test_recorder_is_empty_when_no_fallback_is_invoked() -> None:
    with compatibility.track_fallback_calls() as calls:
        pass
    assert calls == []


def test_explicit_v2_plugin_calls_no_legacy_fallback() -> None:
    from openfoam_driver.tests.plugins.minimal_plugin import MinimalOpenFOAMPlugin
    from openfoam_driver.core.plugin_interface import driver_context

    context = driver_context(MinimalOpenFOAMPlugin(), source="test:minimal")
    with compatibility.track_fallback_calls() as calls:
        context.capabilities.run_document_configuration.schema()
    assert calls == []
