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
#     test_runtime_evidence
#
# Description
#     Plugins declare where runtime evidence lives; core consumes it later.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

"""Declaration surface consumed by Phase 2 (provenance), Phase 4 (telemetry),
and Phase 5 (observables). Nothing in Phase 1 reads these yet -- they exist so
those phases need not reopen the plugin contract mid-flight."""

from __future__ import annotations

from pathlib import Path

from openfoam_driver.core.plugin_interface import (
    default_driver_context,
    generic_openfoam_context,
)


def test_cardiac_declares_its_solver_as_a_solve_step() -> None:
    evidence = default_driver_context().capabilities.runtime_evidence
    assert "cardiacFoam" in evidence.solve_step_commands()


def test_a_post_processing_utility_is_not_a_solve_step() -> None:
    """bathBidomainInterfaceMetrics is authorized to run but does not solve;
    Phase 4 must not expect solver telemetry from it."""
    evidence = default_driver_context().capabilities.runtime_evidence
    assert "bathBidomainInterfaceMetrics" not in evidence.solve_step_commands()


def test_generic_declares_no_solve_steps() -> None:
    evidence = generic_openfoam_context().capabilities.runtime_evidence
    assert evidence.solve_step_commands() == frozenset()


def test_allrun_declares_a_log_glob_so_redirected_output_is_findable() -> None:
    """OpenFOAM's runApplication redirects solver output to log.<app>, so an
    Allrun step produces no parseable driver-captured stdout."""
    evidence = default_driver_context().capabilities.runtime_evidence
    globs = evidence.telemetry_source_globs("Allrun")
    assert any("log." in glob for glob in globs)


def test_a_command_with_no_declared_globs_returns_empty() -> None:
    evidence = default_driver_context().capabilities.runtime_evidence
    assert evidence.telemetry_source_globs("blockMesh") == ()


def test_extra_provenance_paths_default_to_empty(tmp_path: Path) -> None:
    generic = generic_openfoam_context().capabilities.runtime_evidence
    assert generic.extra_provenance_paths(tmp_path) == ()


def test_cardiac_extra_provenance_paths_declares_the_solver_as_a_dependency(
    tmp_path: Path,
) -> None:
    """End-to-end through the adapter: the whole reason extra_provenance_paths
    was replaced with a typed RuntimeDependency form is so the cardiacFoam
    binary itself -- never named by an Allrun-driven step's command -- is
    still declared as something Phase 2 must fingerprint."""
    from openfoam_driver.core.plugin_capabilities import RuntimeDependency

    evidence = default_driver_context().capabilities.runtime_evidence
    dependencies = evidence.extra_provenance_paths(tmp_path)
    assert dependencies
    assert all(isinstance(dependency, RuntimeDependency) for dependency in dependencies)
    by_name = {dependency.name: dependency for dependency in dependencies}
    assert "cardiacFoam" in by_name
    assert by_name["cardiacFoam"].required is True


def test_an_unknown_artifact_format_has_no_reader() -> None:
    evidence = default_driver_context().capabilities.runtime_evidence
    assert evidence.artifact_value_reader("not_a_real_format") is None


def test_generic_plugin_provides_no_artifact_readers() -> None:
    evidence = generic_openfoam_context().capabilities.runtime_evidence
    assert evidence.artifact_value_reader("openfoam_log") is None
