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
#     test_provenance_inputs
#
# Description
#     Canonical input enumeration (Task 2b): which on-disk files, case
#     scripts, and runtime dependencies a workflow run actually consumes.
#     Classification is by consumption, not authorship (I1) -- an
#     unclassified file defaults to required_input, since a spurious refusal
#     is recoverable and a silent stale replay is the incident this phase
#     exists to prevent.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

from __future__ import annotations

import os
import stat
from pathlib import Path

import pytest

from openfoam_driver.core.plugin_capabilities import ResolvedInput, RuntimeDependency
from openfoam_driver.core.plugin_interface import default_driver_context, driver_context, generic_openfoam_context
from openfoam_driver.core.runtime.provenance_inputs import enumerate_case_inputs
from openfoam_driver.tests.plugins.minimal_plugin import MinimalOpenFOAMPlugin


def _paths(components, *, kind: str | None = None) -> set[str]:
    return {c.path for c in components if kind is None or c.kind == kind}


def _by_path(components, path: str):
    matches = [c for c in components if c.path == path]
    assert len(matches) == 1, f"expected exactly one component for {path!r}, got {matches}"
    return matches[0]


def _write_control_dict(case_root: Path, *, start_from: str, start_time: str = "0") -> None:
    system = case_root / "system"
    system.mkdir(parents=True, exist_ok=True)
    (system / "controlDict").write_text(
        f"startFrom       {start_from};\nstartTime       {start_time};\n"
        "endTime         1;\ndeltaT          0.01;\n"
    )


def _make_executable(path: Path, content: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(content)
    path.chmod(path.stat().st_mode | stat.S_IEXEC | stat.S_IXGRP | stat.S_IXOTH)


class _FakePlugin(MinimalOpenFOAMPlugin):
    """A v1 plugin that declares CaseProvenanceCapability / RuntimeEvidence
    hooks inline, so precedence can be exercised without a tutorial."""

    def __init__(self, *, required_inputs=(), generated_output_globs=(), extra_provenance_paths=()):
        self._required_inputs = required_inputs
        self._generated_output_globs = generated_output_globs
        self._extra_provenance_paths = extra_provenance_paths

    def get_required_inputs(self, case_root, resolved_case, selected_start_time):
        return self._required_inputs

    def get_generated_output_globs(self, case_root, resolved_case, selected_start_time):
        return self._generated_output_globs

    def get_extra_provenance_paths(self, case_root):
        return self._extra_provenance_paths


def test_system_and_constant_are_required_and_diagnostic_outputs_are_excluded(tmp_path: Path) -> None:
    """cardiacFoam's own CaseProvenanceCapability excludes the mesh-diagnostic
    byproducts nothing reads (I1's worked example)."""
    _write_control_dict(tmp_path, start_from="startTime", start_time="0")
    (tmp_path / "constant").mkdir()
    (tmp_path / "constant" / "electroProperties").write_text("solver monodomain;\n")
    (tmp_path / "constant" / "C").write_bytes(b"mesh-diagnostic-byproduct")
    (tmp_path / "constant" / "skewness").write_bytes(b"mesh-diagnostic-byproduct")

    components = enumerate_case_inputs(
        tmp_path, workflow_dag={"steps": []}, driver_context=default_driver_context(),
    )
    included = _paths(components, kind="case_file")

    assert "system/controlDict" in included
    assert "constant/electroProperties" in included
    assert "constant/C" not in included
    assert "constant/skewness" not in included


def test_selected_start_time_directory_is_included_others_excluded(tmp_path: Path) -> None:
    _write_control_dict(tmp_path, start_from="startTime", start_time="0")
    for time_name in ("0", "0.5", "1"):
        time_dir = tmp_path / time_name
        time_dir.mkdir()
        (time_dir / "Vm").write_text(f"field-at-{time_name}")

    components = enumerate_case_inputs(
        tmp_path, workflow_dag={"steps": []}, driver_context=generic_openfoam_context(),
    )
    included = _paths(components, kind="case_file")

    assert "0/Vm" in included
    assert "0.5/Vm" not in included
    assert "1/Vm" not in included


def test_latest_time_selects_the_latest_written_time_directory(tmp_path: Path) -> None:
    """startFrom latestTime selects the latest written time, not 0 -- 0/
    still holds initial/boundary conditions but is not itself the input."""
    _write_control_dict(tmp_path, start_from="latestTime", start_time="0")
    for time_name in ("0", "0.5"):
        time_dir = tmp_path / time_name
        time_dir.mkdir()
        (time_dir / "Vm").write_text(f"field-at-{time_name}")

    components = enumerate_case_inputs(
        tmp_path, workflow_dag={"steps": []}, driver_context=generic_openfoam_context(),
    )
    included = _paths(components, kind="case_file")

    assert "0.5/Vm" in included
    assert "0/Vm" not in included


def test_postprocessing_and_workflow_logs_are_excluded(tmp_path: Path) -> None:
    _write_control_dict(tmp_path, start_from="startTime", start_time="0")
    (tmp_path / "postProcessing" / "workflow_logs").mkdir(parents=True)
    (tmp_path / "postProcessing" / "workflow_logs" / "solve.attempt1.stdout.log").write_text("log")
    (tmp_path / "postProcessing" / "ecgProbes" / "0").mkdir(parents=True)
    (tmp_path / "postProcessing" / "ecgProbes" / "0" / "data").write_text("probe data")

    components = enumerate_case_inputs(
        tmp_path, workflow_dag={"steps": []}, driver_context=generic_openfoam_context(),
    )
    included = _paths(components, kind="case_file")

    assert not any(path.startswith("postProcessing/") for path in included)


def test_a_case_with_no_constant_does_not_raise(tmp_path: Path) -> None:
    _write_control_dict(tmp_path, start_from="startTime", start_time="0")

    components = enumerate_case_inputs(
        tmp_path, workflow_dag={"steps": []}, driver_context=default_driver_context(),
    )

    assert "system/controlDict" in _paths(components, kind="case_file")


def test_generic_plugin_still_requires_unknown_files(tmp_path: Path) -> None:
    """Under the generic plugin (declares nothing), the same file cardiacFoam
    would exclude stays required -- an unclassified file always defaults to
    required_input (I1)."""
    _write_control_dict(tmp_path, start_from="startTime", start_time="0")
    (tmp_path / "constant").mkdir()
    (tmp_path / "constant" / "C").write_bytes(b"mesh-diagnostic-byproduct")

    components = enumerate_case_inputs(
        tmp_path, workflow_dag={"steps": []}, driver_context=generic_openfoam_context(),
    )

    assert "constant/C" in _paths(components, kind="case_file")


def test_an_allrun_named_by_the_dag_is_included(tmp_path: Path) -> None:
    _write_control_dict(tmp_path, start_from="startTime", start_time="0")
    _make_executable(tmp_path / "Allrun", b"#!/bin/sh\ncardiacFoam\n")

    workflow_dag = {"steps": [{"id": "solve", "command": "Allrun", "depends_on": []}]}
    components = enumerate_case_inputs(
        tmp_path, workflow_dag=workflow_dag, driver_context=generic_openfoam_context(),
    )

    allrun = _by_path(components, "Allrun")
    assert allrun.kind == "case_file"
    assert allrun.strength == "content"


def test_a_parallel_step_includes_both_mpirun_and_its_payload(tmp_path: Path) -> None:
    _write_control_dict(tmp_path, start_from="startTime", start_time="0")
    bin_dir = tmp_path / "fakebin"
    _make_executable(bin_dir / "mpirun", b"#!/bin/sh\necho mpirun\n")
    _make_executable(bin_dir / "cardiacFoam", b"#!/bin/sh\necho solve\n")

    workflow_dag = {
        "steps": [
            {
                "id": "solve",
                "command": "mpirun",
                "args": ["-np", "4", "cardiacFoam", "-parallel"],
                "depends_on": [],
            }
        ]
    }
    env = {"PATH": str(bin_dir)}
    components = enumerate_case_inputs(
        tmp_path, workflow_dag=workflow_dag, driver_context=generic_openfoam_context(), env=env,
    )

    mpirun = _by_path(components, "mpirun")
    payload = _by_path(components, "cardiacFoam")
    assert mpirun.kind == "runtime_dependency"
    assert mpirun.strength == "content"
    assert payload.kind == "runtime_dependency"
    assert payload.strength == "content"


def test_a_required_but_unresolved_step_executable_is_unavailable_not_omitted(tmp_path: Path) -> None:
    _write_control_dict(tmp_path, start_from="startTime", start_time="0")
    workflow_dag = {"steps": [{"id": "mesh", "command": "blockMesh", "depends_on": []}]}
    env = {"PATH": str(tmp_path / "nowhere")}

    components = enumerate_case_inputs(
        tmp_path, workflow_dag=workflow_dag, driver_context=generic_openfoam_context(), env=env,
    )

    block_mesh = _by_path(components, "blockMesh")
    assert block_mesh.kind == "runtime_dependency"
    assert block_mesh.strength == "unavailable"


def test_plugin_runtime_dependency_entries_appear_and_missing_required_is_unavailable(
    tmp_path: Path,
) -> None:
    _write_control_dict(tmp_path, start_from="startTime", start_time="0")
    present_lib = tmp_path / "libbin" / "libBar.dylib"
    present_lib.parent.mkdir(parents=True)
    present_lib.write_bytes(b"library contents")

    plugin = _FakePlugin(
        extra_provenance_paths=(
            RuntimeDependency(name="libFoo", path=None, required=True),
            RuntimeDependency(name="libBar", path=present_lib, required=False),
        )
    )
    context = driver_context(plugin, source="test")

    components = enumerate_case_inputs(
        tmp_path, workflow_dag={"steps": []}, driver_context=context,
    )

    lib_foo = _by_path(components, "libFoo")
    lib_bar = _by_path(components, "libBar")
    assert lib_foo.kind == "runtime_dependency"
    assert lib_foo.strength == "unavailable"
    assert lib_bar.strength == "content"


def test_dag_consumes_declaration_wins_over_a_generated_output_glob(tmp_path: Path) -> None:
    """I1's resolution precedence: a DAG step's consumes declaration beats a
    plugin's generated_output_globs exclusion."""
    _write_control_dict(tmp_path, start_from="startTime", start_time="0")
    (tmp_path / "constant").mkdir()
    (tmp_path / "constant" / "C").write_bytes(b"mesh-diagnostic-byproduct")

    workflow_dag = {
        "steps": [
            {
                "id": "verify",
                "command": "postProcess",
                "depends_on": [],
                "consumes": ["constant/C"],
            }
        ]
    }
    components = enumerate_case_inputs(
        tmp_path, workflow_dag=workflow_dag, driver_context=default_driver_context(),
    )

    assert "constant/C" in _paths(components, kind="case_file")


def test_plugin_required_inputs_entry_wins_over_a_generated_output_glob(tmp_path: Path) -> None:
    """I1's resolution precedence: a plugin required_inputs() entry beats its
    own generated_output_globs exclusion."""
    _write_control_dict(tmp_path, start_from="startTime", start_time="0")
    (tmp_path / "constant").mkdir()
    target = tmp_path / "constant" / "C"
    target.write_bytes(b"actually consumed this time")

    plugin = _FakePlugin(
        required_inputs=(ResolvedInput(name="C", path=target, required=True, consumer="test"),),
        generated_output_globs=("constant/C",),
    )
    context = driver_context(plugin, source="test")

    components = enumerate_case_inputs(
        tmp_path, workflow_dag={"steps": []}, driver_context=context,
    )

    assert "constant/C" in _paths(components, kind="case_file")


def test_optional_required_input_that_is_absent_is_not_added(tmp_path: Path) -> None:
    """READ_IF_PRESENT and absent is not the same as MUST_READ and missing --
    nothing was going to be consumed, so nothing is fingerprinted."""
    _write_control_dict(tmp_path, start_from="startTime", start_time="0")
    plugin = _FakePlugin(
        required_inputs=(ResolvedInput(name="Optional", path=None, required=False, consumer="test"),),
    )
    context = driver_context(plugin, source="test")

    components = enumerate_case_inputs(
        tmp_path, workflow_dag={"steps": []}, driver_context=context,
    )

    assert "Optional" not in _paths(components)


def test_processor_selected_time_is_included_other_processor_times_excluded(tmp_path: Path) -> None:
    """I9: processor*/<selected-time>/** is a required input on the same
    footing as the serial case; other times under processor*/ are outputs."""
    _write_control_dict(tmp_path, start_from="startTime", start_time="0")
    proc0 = tmp_path / "processor0"
    (proc0 / "0").mkdir(parents=True)
    (proc0 / "0" / "Vm").write_text("decomposed restart field")
    (proc0 / "0.5").mkdir(parents=True)
    (proc0 / "0.5" / "Vm").write_text("later time, not an input")

    components = enumerate_case_inputs(
        tmp_path, workflow_dag={"steps": []}, driver_context=generic_openfoam_context(),
    )
    included = _paths(components, kind="case_file")

    assert "processor0/0/Vm" in included
    assert "processor0/0.5/Vm" not in included
