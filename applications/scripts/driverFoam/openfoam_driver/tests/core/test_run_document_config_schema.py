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

"""Plugin-declared RunDocument.config schema validation (P2.2)."""
from __future__ import annotations

from openfoam_driver.core.plugin_interface import default_driver_context
from openfoam_driver.strict_planning import strict_plan


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
