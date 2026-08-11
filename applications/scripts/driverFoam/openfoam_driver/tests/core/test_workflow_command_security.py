"""Command-boundary security tests for RunDocument execution."""
from __future__ import annotations

import os
import stat
import tempfile
import unittest
from pathlib import Path

from openfoam_driver.core.plugin_interface import default_driver_context
from openfoam_driver.core.runtime.workflow import validate_workflow_commands
from openfoam_driver.core.runtime.workflow_runner import (
    _resolve_case_cwd,
    _resolve_command,
)


def _make_executable(path: Path) -> None:
    path.write_text("#!/bin/sh\necho shadow\n")
    path.chmod(path.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)


class TestValidateWorkflowCommands(unittest.TestCase):
    def setUp(self) -> None:
        # The allowlist is now sourced from the active plugin context; the
        # cardiac context is what these cases have always exercised.
        self.context = default_driver_context()

    def test_known_openfoam_command_is_allowed(self) -> None:
        dag = {"steps": [{"id": "s", "command": "cardiacFoam"}]}
        self.assertEqual(validate_workflow_commands(dag, driver_context=self.context), ())

    def test_case_script_command_is_allowed(self) -> None:
        dag = {"steps": [{"id": "s", "command": "Allrun"}]}
        self.assertEqual(validate_workflow_commands(dag, driver_context=self.context), ())

    def test_gmsh_is_allowed(self) -> None:
        # Tet-mesh sweep workflows (mesh_family="tet") run gmsh/gmshToFoam/
        # checkMesh as explicit workflow steps; the allowlist must be static
        # (not gated on whether the caller's shell happens to have OpenFOAM
        # sourced), so these are named directly, same as blockMesh.
        dag = {"steps": [{"id": "s", "command": "gmsh"}]}
        self.assertEqual(validate_workflow_commands(dag, driver_context=self.context), ())

    def test_gmsh_to_foam_is_allowed(self) -> None:
        dag = {"steps": [{"id": "s", "command": "gmshToFoam"}]}
        self.assertEqual(validate_workflow_commands(dag, driver_context=self.context), ())

    def test_check_mesh_is_allowed(self) -> None:
        dag = {"steps": [{"id": "s", "command": "checkMesh"}]}
        self.assertEqual(validate_workflow_commands(dag, driver_context=self.context), ())

    def test_bath_bidomain_interface_metrics_is_allowed(self) -> None:
        # bath_tet's canonical reported metrics come from this utility (a
        # post-hoc pass over the reconstructed mesh, since the live verifier
        # can't do heart/bath mesh-subsetting during a parallel-decomposed
        # solve) run as its own workflow step after solve -- authorized by the
        # cardiac plugin (it ships no utility.manifest.toml), not by core.
        dag = {"steps": [{"id": "s", "command": "bathBidomainInterfaceMetrics"}]}
        self.assertEqual(validate_workflow_commands(dag, driver_context=self.context), ())

    def test_mpirun_is_allowed(self) -> None:
        # run_in_parallel=True wraps the solve step as
        # `mpirun -np <N> cardiacFoam -parallel`; only the bare command is
        # allowlist-checked (args are not re-validated as commands), same
        # reasoning as the other explicit workflow-step entries above.
        dag = {"steps": [{"id": "s", "command": "mpirun", "args": ["-np", "6", "cardiacFoam", "-parallel"]}]}
        self.assertEqual(validate_workflow_commands(dag, driver_context=self.context), ())

    def test_unknown_command_is_rejected(self) -> None:
        dag = {"steps": [{"id": "s", "command": "rm"}]}
        codes = {d.code for d in validate_workflow_commands(dag, driver_context=self.context)}
        self.assertIn("unknown_workflow_command", codes)

    def test_empty_command_is_rejected(self) -> None:
        dag = {"steps": [{"id": "s", "command": ""}]}
        codes = {d.code for d in validate_workflow_commands(dag, driver_context=self.context)}
        self.assertIn("workflow_step_without_command", codes)

    def test_case_script_allclean_is_allowed(self) -> None:
        dag = {"steps": [{"id": "s", "command": "Allclean"}]}
        self.assertEqual(validate_workflow_commands(dag, driver_context=self.context), ())

    def test_explicit_relative_case_script_is_allowed(self) -> None:
        dag = {"steps": [{"id": "s", "command": "./Allrun"}]}
        self.assertEqual(validate_workflow_commands(dag, driver_context=self.context), ())

    def test_dotted_case_script_is_allowed(self) -> None:
        dag = {"steps": [{"id": "s", "command": "Allrun.pre"}]}
        self.assertEqual(validate_workflow_commands(dag, driver_context=self.context), ())

    def test_explicit_relative_dotted_case_script_is_allowed(self) -> None:
        dag = {"steps": [{"id": "s", "command": "./Allrun.post"}]}
        self.assertEqual(validate_workflow_commands(dag, driver_context=self.context), ())

    def test_explicit_relative_non_case_script_is_rejected(self) -> None:
        dag = {"steps": [{"id": "s", "command": "./notAllrun"}]}
        codes = {d.code for d in validate_workflow_commands(dag, driver_context=self.context)}
        self.assertIn("unknown_workflow_command", codes)

    def test_absolute_path_is_rejected(self) -> None:
        dag = {"steps": [{"id": "s", "command": "/usr/bin/checkMesh"}]}
        codes = {d.code for d in validate_workflow_commands(dag, driver_context=self.context)}
        self.assertIn("unknown_workflow_command", codes)

    def test_none_dag_is_empty(self) -> None:
        self.assertEqual(validate_workflow_commands(None, driver_context=self.context), ())


class TestValidateWorkflowCommandsFoamApp(unittest.TestCase):
    def setUp(self) -> None:
        self.context = default_driver_context()

    def test_installed_openfoam_app_is_allowed(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            appbin = Path(temp) / "bin"
            appbin.mkdir()
            _make_executable(appbin / "checkMeshFake")
            old_path = os.environ.get("PATH", "")
            old_appbin = os.environ.get("FOAM_APPBIN")
            os.environ["PATH"] = f"{appbin}{os.pathsep}{old_path}"
            os.environ["FOAM_APPBIN"] = str(appbin)
            try:
                dag = {"steps": [{"id": "s", "command": "checkMeshFake"}]}
                self.assertEqual(validate_workflow_commands(dag, driver_context=self.context), ())
            finally:
                os.environ["PATH"] = old_path
                if old_appbin is None:
                    os.environ.pop("FOAM_APPBIN", None)
                else:
                    os.environ["FOAM_APPBIN"] = old_appbin

    def test_command_outside_foam_bins_is_rejected(self) -> None:
        old_appbin = os.environ.pop("FOAM_APPBIN", None)
        old_userbin = os.environ.pop("FOAM_USER_APPBIN", None)
        try:
            dag = {"steps": [{"id": "s", "command": "definitelyNotAFoamApp"}]}
            codes = {d.code for d in validate_workflow_commands(dag, driver_context=self.context)}
            self.assertIn("unknown_workflow_command", codes)
        finally:
            if old_appbin is not None:
                os.environ["FOAM_APPBIN"] = old_appbin
            if old_userbin is not None:
                os.environ["FOAM_USER_APPBIN"] = old_userbin


class TestResolveCommandShadowing(unittest.TestCase):
    def test_bare_binary_name_is_never_resolved_to_case_local_file(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            cwd = Path(temp)
            _make_executable(cwd / "cardiacFoam")
            self.assertEqual(_resolve_command("cardiacFoam", cwd), "cardiacFoam")

    def test_recognized_case_script_resolves_locally_when_present(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            cwd = Path(temp)
            _make_executable(cwd / "Allrun")
            self.assertEqual(_resolve_command("Allrun", cwd), str(cwd / "Allrun"))

    def test_case_script_falls_through_when_absent(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            self.assertEqual(_resolve_command("Allrun", Path(temp)), "Allrun")

    def test_explicit_relative_path_passes_through(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            self.assertEqual(_resolve_command("./Allrun", Path(temp)), "./Allrun")

    def test_absolute_path_passes_through(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            self.assertEqual(_resolve_command("/usr/bin/env", Path(temp)), "/usr/bin/env")


class TestResolveCaseCwdContainment(unittest.TestCase):
    def test_symlinked_cwd_escape_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp) / "case"
            root.mkdir()
            outside = Path(temp) / "outside"
            outside.mkdir()
            (root / "evil").symlink_to(outside, target_is_directory=True)
            with self.assertRaises(ValueError):
                _resolve_case_cwd(root, "evil")

    def test_in_tree_cwd_is_allowed(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp) / "case"
            (root / "sub").mkdir(parents=True)
            resolved = _resolve_case_cwd(root, "sub")
            self.assertEqual(resolved.name, "sub")


if __name__ == "__main__":
    unittest.main()
