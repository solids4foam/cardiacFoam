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
#     test_run_document_config_schema
#
# Description
#     Tests plugin-declared RunDocument.config schema validation (P2.2).
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

"""Plugin-declared RunDocument.config schema validation (P2.2).

Covers both directions of the contract: the *emission* path (a plugin-built
config, checked in ``run_document_adapter``) and the *ingestion* path (an
agent-authored document read off disk, checked in ``run_document_exec``).
The two must stay symmetric -- a config the planner would refuse to emit is
a config the executor must refuse to ingest.
"""
from __future__ import annotations

import json
import tempfile
from pathlib import Path

from openfoam_driver.core.plugin_interface import default_driver_context
from openfoam_driver.plugins.cardiacfoam.run_document_config import _read_physics_type
from openfoam_driver.core.strict_planning import strict_plan


def test_cardiac_plugin_declares_a_config_schema() -> None:
    context = default_driver_context()
    schema = context.plugin.get_run_document_config_schema()
    assert schema["required"] == ["anatomy", "physics", "stimulus", "solver"]


def test_strict_plan_reports_a_structured_diagnostic_for_schema_violation(monkeypatch) -> None:
    """A plugin that builds a config violating its own declared schema must
    surface a StrictDiagnostic an agent can read and act on -- not a raw
    jsonschema traceback and not a silent pass."""
    from openfoam_driver.core.plugin_capabilities import RunDocumentConfigurationRequest
    from openfoam_driver.plugins import cardiacfoam_plugin

    context = default_driver_context()

    def _broken_build(spec):
        # Deliberately omit the required "solver" phase key.
        return {"anatomy": {}, "physics": {}, "stimulus": {}}, ()

    monkeypatch.setattr(
        context.plugin, "build_run_document_config", _broken_build, raising=False,
    )
    report = strict_plan("singleCell", driver_context=context)
    codes = {d.code for d in report.validation_diagnostics}
    assert "plugin_config_schema_violation" in codes
    messages = [d.message for d in report.validation_diagnostics if d.code == "plugin_config_schema_violation"]
    assert any("solver" in message for message in messages)


def _document_json(config: dict) -> dict:
    return {
        "version": "3",
        "id": "ingested",
        "name": "ingested",
        "status": "planned",
        "config": config,
        "launch": {"caseRoot": "/nonexistent/case", "outputDir": "/nonexistent/case/out"},
        "workflowDag": {
            "schema_version": "1",
            "step_status_values": [
                "pending", "running", "completed", "failed", "skipped",
            ],
            "steps": [{
                "id": "solve", "command": "Allrun", "args": [], "cwd": ".",
                "depends_on": [], "produces": [], "consumes": [],
                "retry_policy": {}, "command_display": "Allrun",
            }],
        },
    }


def _ingest(config: dict) -> tuple[dict, ...]:
    """Load a hand-authored document off disk and adapt it for execution."""
    from openfoam_driver.core.runtime.run_document_exec import (
        build_execution_inputs,
        load_run_document,
    )

    with tempfile.TemporaryDirectory() as temp:
        path = Path(temp) / "run.json"
        path.write_text(json.dumps(_document_json(config)))
        run_doc = load_run_document(path)
    _inputs, diagnostics = build_execution_inputs(run_doc)
    return diagnostics


def test_ingested_document_is_checked_against_the_plugin_config_schema() -> None:
    """An agent-authored config that violates the plugin's own schema must be
    rejected at *ingestion*, with the same diagnostic code the emission path
    uses. The ingestion path is the untrusted one: before this gate it saw
    only ``validate_run`` and never consulted the plugin schema at all.

    ``tissue`` carries a closed enum in the cardiac plugin's config schema
    but no catalog enum that ``validate_run`` would independently reject, so
    an out-of-enum value here isolates the plugin-schema gate.
    """
    diagnostics = _ingest({
        "anatomy": {},
        "physics": {"tissue": "notATissue"},
        "stimulus": {},
        "solver": {},
    })
    violations = [
        d for d in diagnostics if d["code"] == "plugin_config_schema_violation"
    ]
    assert violations, diagnostics
    assert violations[0]["field"] == "physics.tissue"
    assert "notATissue" in violations[0]["message"]


def test_ingested_document_with_a_schema_valid_config_raises_no_violation() -> None:
    """The gate must not fire on a config the plugin schema accepts (the
    all-empty phase config a core generic case emits)."""
    diagnostics = _ingest({
        "anatomy": {}, "physics": {}, "stimulus": {}, "solver": {},
    })
    assert not [
        d for d in diagnostics if d["code"] == "plugin_config_schema_violation"
    ], diagnostics


def test_read_physics_type_ignores_a_nested_type_key(tmp_path):
    """A nested block's own 'type' key must not shadow the real top-level one.

    The pre-migration scanner matches the first line starting with 'type'
    anywhere in the file, with no nesting awareness -- a hypothetical
    nested block declared before the real entry would silently win.
    """
    path = tmp_path / "physicsProperties"
    path.write_text(
        "FoamFile{ version 2.0; format ascii; class dictionary; object physicsProperties; }\n"
        "someSubBlock\n{\n    type notTheRealAnswer;\n}\n"
        "type monodomain;\n"
    )
    assert _read_physics_type(path) == "monodomain"


def test_read_physics_type_returns_none_when_file_missing(tmp_path):
    assert _read_physics_type(tmp_path / "nope") is None
