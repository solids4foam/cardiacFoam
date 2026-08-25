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
#     test_case_provenance_capability
#
# Description
#     CaseProvenanceCapability is routed through the adapter exactly like
#     every Phase 1 capability, with an empty fallback -- which under I1's
#     precedence means "everything unknown is a required input", the safe
#     default for a plugin (or plugin version) that declares nothing.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

from __future__ import annotations

from pathlib import Path

from openfoam_driver.core.plugin_interface import (
    default_driver_context,
    driver_context,
    generic_openfoam_context,
)
from openfoam_driver.tests.plugins.minimal_plugin import MinimalOpenFOAMPlugin


def test_generic_plugin_declares_no_required_inputs(tmp_path: Path) -> None:
    generic = generic_openfoam_context().capabilities.case_provenance
    assert generic.required_inputs(tmp_path, {}, "0") == ()


def test_generic_plugin_declares_no_generated_outputs(tmp_path: Path) -> None:
    generic = generic_openfoam_context().capabilities.case_provenance
    assert generic.generated_output_globs(tmp_path, {}, "0") == ()


def test_a_v1_plugin_with_no_hooks_gets_the_empty_fallback(tmp_path: Path) -> None:
    """A plugin that predates this capability -- v1 or a v2 third-party
    plugin that never implemented it -- must still load and adapt cleanly.
    CaseProvenanceCapability is not a mandatory SolverPlugin member."""
    context = driver_context(MinimalOpenFOAMPlugin(), source="test")
    assert context.capabilities.case_provenance.required_inputs(tmp_path, {}, "0") == ()
    assert context.capabilities.case_provenance.generated_output_globs(tmp_path, {}, "0") == ()


def test_cardiac_declares_the_mesh_diagnostic_fields_as_generated_outputs(
    tmp_path: Path,
) -> None:
    """constant/C, Cx, Cy, Cz, skewness are consumed by nothing -- confirmed
    by exhaustive grep across src/ and applications/: zero hits."""
    cardiac = default_driver_context().capabilities.case_provenance
    globs = cardiac.generated_output_globs(tmp_path, {}, "0")
    assert set(globs) == {
        "constant/C",
        "constant/Cx",
        "constant/Cy",
        "constant/Cz",
        "constant/skewness",
    }


def test_cardiac_required_inputs_defers_to_the_safe_default(tmp_path: Path) -> None:
    """Model-dependent per-field required-input resolution is Task 2b's job
    (input enumeration). Returning () here is safe under I1's precedence:
    an unclassified file still defaults to required_input upstream."""
    cardiac = default_driver_context().capabilities.case_provenance
    assert cardiac.required_inputs(tmp_path, {}, "0") == ()
