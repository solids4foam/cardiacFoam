from __future__ import annotations

from dataclasses import fields
from pathlib import Path

import pytest

from openfoam_driver.core.plugin_capabilities import (
    ArtifactPredictionRequest,
    CaseCompatibilityRequest,
    ConfigurationValidationRequest,
    RunDocumentConfigurationRequest,
    RunSemanticValidationRequest,
)
from openfoam_driver.core.plugin_interface import DriverContext, driver_context
from openfoam_driver.core.runtime.models import TutorialSpec
from openfoam_driver.tests.plugins.minimal_plugin import MinimalOpenFOAMPlugin


def _spec(tmp_path: Path) -> TutorialSpec:
    return TutorialSpec(
        name="minimal",
        case_root=tmp_path,
        setup_root=tmp_path,
        output_dir=tmp_path / "postProcessing",
        build_cases=lambda: [],
        apply_case=lambda *_args: None,
        metadata={"generic_case": True},
    )


def test_context_exposes_focused_adapters_without_replacing_public_plugin(
    tmp_path: Path,
) -> None:
    plugin = MinimalOpenFOAMPlugin()
    context = driver_context(plugin, source="test")
    spec = _spec(tmp_path)

    assert context.plugin is plugin
    assert context.capabilities.tutorials.catalog() == plugin.get_tutorial_catalog()
    assert context.capabilities.dictionaries.entries() == plugin.get_dict_entries()
    assert context.capabilities.manifest.manifest() == plugin.get_capabilities()
    assert context.capabilities.configuration_validator.validate(
        ConfigurationValidationRequest(spec),
    ) == ()
    assert context.capabilities.run_semantic_validator.validate(
        RunSemanticValidationRequest({}),
    ) == ()
    assert context.capabilities.artifacts.predict(
        ArtifactPredictionRequest(tmp_path, spec),
    ) == ()
    config, diagnostics = context.capabilities.run_document_configuration.build(
        RunDocumentConfigurationRequest(spec),
    )
    # A non-cardiac plugin with no build_run_document_config() hook now gets
    # an empty config rather than the cardiac phase vocabulary. Those four
    # phase names are exactly what RunDocument v3 removed from core, where
    # `config` is an open object with no fixed phases (schemas/run-document.json),
    # so handing them to a plugin that never declared them contradicted the
    # schema. Matches legacy_run_document_config_schema, which already handed
    # non-cardiac plugins a fully open schema.
    assert config == {}
    assert diagnostics == ()

    # Existing callers that constructed DriverContext(plugin, identity)
    # directly retain the same constructor shape.
    reconstructed = DriverContext(plugin, context.identity)
    assert reconstructed.plugin is plugin
    assert reconstructed.capabilities.tutorials.catalog() == plugin.get_tutorial_catalog()
    assert [item.name for item in fields(reconstructed)] == ["plugin", "identity"]


def test_non_cardiac_plugin_does_not_inherit_cardiac_case_evidence(
    tmp_path: Path,
) -> None:
    """Same rule as :func:`test_report_catalog_is_empty_for_non_cardiac_plugin`,
    applied to case compatibility.

    This test previously asserted the opposite, under the name
    ``test_legacy_plugin_case_evidence_preserves_pre_capability_behavior``: a
    non-cardiac plugin DID claim a case carrying ``electroProperties*``,
    because legacy_case_marker/legacy_case_runnable_without_workflow called
    the cardiac implementation without checking plugin_id -- unlike the
    thirteen sibling fallbacks, which all gate on ``org.cardiacfoam``.

    Preserving that behaviour was never the intent; it was the Plan-1
    fallback's unexamined default, and it meant a third-party plugin was
    silently judged by cardiac filesystem evidence."""
    plugin = MinimalOpenFOAMPlugin()
    context = driver_context(plugin, source="test")
    case_root = tmp_path / "case"
    for relative in (
        "constant/electroProperties.variant",
        "constant/physicsProperties",
        "system/controlDict",
        "system/fvSchemes",
        "system/fvSolution",
    ):
        path = case_root / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("")

    request = CaseCompatibilityRequest(case_root)
    assert not context.capabilities.case_compatibility.has_case_marker(request)
    assert not context.capabilities.case_compatibility.is_runnable_without_workflow(
        request
    )


def test_report_catalog_is_empty_for_non_cardiac_plugin() -> None:
    """P2.7: report_catalog.py's former REPORTS tuple was cardiac-specific
    data consumed unconditionally. A non-cardiac v1 plugin must get an empty
    report catalog, not the built-in "Vm field"/"activation map" reports."""
    plugin = MinimalOpenFOAMPlugin()
    context = driver_context(plugin, source="test")

    reports = context.capabilities.report_catalog.reports()

    assert reports == ()
    assert not any("Vm field" in r.title for r in reports)
    assert not any("activation map" in r.title.lower() for r in reports)


def test_capability_adapter_preserves_plugin_exceptions(tmp_path: Path) -> None:
    class ThrowingPlugin(MinimalOpenFOAMPlugin):
        def validate_configuration(self, spec):
            del spec
            raise RuntimeError("same failure")

    context = driver_context(ThrowingPlugin(), source="test")
    with pytest.raises(RuntimeError, match="same failure"):
        context.capabilities.configuration_validator.validate(
            ConfigurationValidationRequest(_spec(tmp_path)),
        )
