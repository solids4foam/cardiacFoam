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
#     test_runtime_dependencies
#
# Description
#     Resolves cardiacFoam's runtime dependencies -- the solver binary, its
#     linked/runtime-loaded libraries, and anything the case's own
#     controlDict pulls in via libs (...) -- entirely from an injected
#     environment, so these tests need no sourced OpenFOAM. The live
#     verification module cross-checks the same resolver against a real
#     sourced install.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

from __future__ import annotations

import re
import stat
from pathlib import Path

import pytest

from openfoam_driver.plugins.cardiacfoam.runtime_evidence import (
    _LIBRARY_CATALOG,
    _parse_control_dict_libs,
    resolve_runtime_dependencies,
)
from openfoam_driver.tests.conftest import monorepo_root, skip_without_monorepo


def _write_control_dict(tmp_path: Path, body: str) -> Path:
    p = tmp_path / "controlDict"
    p.write_text(
        "FoamFile\n{\n    version 2.0;\n    format ascii;\n"
        "    class dictionary;\n    object controlDict;\n}\n\n" + body
    )
    return p


class TestParseControlDictLibs:
    """Unit coverage for _parse_control_dict_libs, which reads a
    controlDict's libs ( ... ) list via foamlib's structural, read-only
    parsing rather than hand-rolled paren-counting. The higher-level
    resolve_runtime_dependencies tests above exercise this indirectly for
    the realistic multi-line layout; these pin its edge cases directly."""

    def test_returns_empty_tuple_when_libs_is_absent(self, tmp_path: Path) -> None:
        p = _write_control_dict(tmp_path, "application cardiacFoam;\n")
        assert _parse_control_dict_libs(p) == ()

    def test_returns_empty_tuple_when_libs_is_an_empty_list(self, tmp_path: Path) -> None:
        p = _write_control_dict(tmp_path, "libs\n(\n);\n")
        assert _parse_control_dict_libs(p) == ()

    def test_strips_quotes_from_a_single_entry(self, tmp_path: Path) -> None:
        p = _write_control_dict(tmp_path, 'libs\n(\n    "libverificationModels.so"\n);\n')
        assert _parse_control_dict_libs(p) == ("libverificationModels.so",)

    def test_parses_multiple_entries_in_declaration_order(self, tmp_path: Path) -> None:
        p = _write_control_dict(
            tmp_path,
            'libs\n(\n    "libfoo.so"\n    "libbar.so"\n);\n',
        )
        assert _parse_control_dict_libs(p) == ("libfoo.so", "libbar.so")

    def test_returns_empty_tuple_when_file_does_not_exist(self, tmp_path: Path) -> None:
        assert _parse_control_dict_libs(tmp_path / "missing" / "controlDict") == ()


def _make_lib(directory: Path, name: str, extension: str, content: bytes = b"x") -> Path:
    directory.mkdir(parents=True, exist_ok=True)
    lib = directory / f"lib{name}{extension}"
    lib.write_bytes(content)
    return lib


def _make_executable(path: Path, content: str) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content)
    path.chmod(path.stat().st_mode | stat.S_IEXEC | stat.S_IXGRP | stat.S_IXOTH)
    return path


def test_finds_a_library_in_whichever_defined_directory_holds_it(tmp_path: Path) -> None:
    """libphysicsModel is declared under FOAM_MODULE_LIBBIN in Make/files but
    verified live to live in FOAM_USER_LIBBIN -- so the resolver must search
    every defined lib directory for every name, never one fixed mapping."""
    user_libbin = tmp_path / "user_libbin"
    _make_lib(user_libbin, "physicsModel", ".dylib")
    env = {"FOAM_USER_LIBBIN": str(user_libbin), "PATH": ""}

    deps = {d.name: d for d in resolve_runtime_dependencies(tmp_path / "case", env=env)}
    assert deps["physicsModel"].path == user_libbin / "libphysicsModel.dylib"


def test_missing_foam_module_libbin_does_not_crash_the_search(tmp_path: Path) -> None:
    """FOAM_MODULE_LIBBIN is unset in a plain sourced shell -- the resolver
    must simply skip an undefined variable, not fail."""
    env = {"FOAM_USER_LIBBIN": str(tmp_path / "nonexistent"), "PATH": ""}
    deps = resolve_runtime_dependencies(tmp_path / "case", env=env)
    assert deps  # did not raise


def test_tries_both_platform_extensions(tmp_path: Path) -> None:
    lib_dir = tmp_path / "libbin"
    _make_lib(lib_dir, "electroModels", ".so")
    env = {"FOAM_LIBBIN": str(lib_dir), "PATH": ""}
    deps = {d.name: d for d in resolve_runtime_dependencies(tmp_path / "case", env=env)}
    assert deps["electroModels"].path == lib_dir / "libelectroModels.so"


def test_a_required_library_absent_everywhere_is_reported_not_omitted(tmp_path: Path) -> None:
    env = {"PATH": ""}
    deps = {d.name: d for d in resolve_runtime_dependencies(tmp_path / "case", env=env)}
    assert deps["electroModels"].path is None
    assert deps["electroModels"].required is True


def test_an_optional_library_absent_in_lightweight_mode_is_still_reported(
    tmp_path: Path,
) -> None:
    """libelectroMechanicalModels is absent entirely in the maintainer's
    lightweight (no-solids4foam) default. Its absence is a normal state, not
    an error -- but it must still appear, not be silently dropped, per I3b."""
    env = {"PATH": ""}
    deps = {d.name: d for d in resolve_runtime_dependencies(tmp_path / "case", env=env)}
    assert "electroMechanicalModels" in deps
    assert deps["electroMechanicalModels"].path is None
    assert deps["electroMechanicalModels"].required is False


def test_the_solver_binary_resolves_on_path(tmp_path: Path) -> None:
    bin_dir = tmp_path / "bin"
    fake_solver = _make_executable(bin_dir / "cardiacFoam", "#!/bin/sh\necho v1\n")
    env = {"PATH": str(bin_dir)}
    deps = {d.name: d for d in resolve_runtime_dependencies(tmp_path / "case", env=env)}
    assert deps["cardiacFoam"].path == fake_solver
    assert deps["cardiacFoam"].required is True


def test_the_solver_binary_absent_from_path_is_unavailable(tmp_path: Path) -> None:
    env = {"PATH": ""}
    deps = {d.name: d for d in resolve_runtime_dependencies(tmp_path / "case", env=env)}
    assert deps["cardiacFoam"].path is None
    assert deps["cardiacFoam"].required is True


def test_a_bare_library_named_in_controldict_libs_becomes_required(tmp_path: Path) -> None:
    """libverificationModels is not linked into cardiacFoam at all -- every
    manufactured-solution tutorial depends on its own controlDict libs (...)
    entry to run at all (I2c)."""
    case_root = tmp_path / "case"
    lib_dir = tmp_path / "user_libbin"
    _make_lib(lib_dir, "verificationModels", ".dylib")
    (case_root / "system").mkdir(parents=True)
    (case_root / "system" / "controlDict").write_text(
        'application cardiacFoam;\n\nlibs\n(\n    "libverificationModels.so"\n);\n'
    )
    env = {"FOAM_USER_LIBBIN": str(lib_dir), "PATH": ""}
    deps = {d.name: d for d in resolve_runtime_dependencies(case_root, env=env)}
    assert deps["verificationModels"].required is True
    assert deps["verificationModels"].path == lib_dir / "libverificationModels.dylib"


def test_a_case_local_library_path_expands_foam_case_and_wm_options(tmp_path: Path) -> None:
    """monodomainTotalLagrangianEM loads
    $FOAM_CASE/platforms/$WM_OPTIONS/lib/lib<name>.so -- a case-local library
    built from sources inside the case. $FOAM_CASE is the case being
    resolved; $WM_OPTIONS must come from the executor's own environment, not
    a guess."""
    case_root = tmp_path / "case"
    built_lib_dir = case_root / "platforms" / "myOptions" / "lib"
    built_lib_dir.mkdir(parents=True)
    built_lib = built_lib_dir / "libmanufacturedMonodomainTotalLagrangianEM.so"
    built_lib.write_bytes(b"compiled")
    (case_root / "system").mkdir(parents=True, exist_ok=True)
    (case_root / "system" / "controlDict").write_text(
        "libs\n(\n"
        '    "$FOAM_CASE/platforms/$WM_OPTIONS/lib/libmanufacturedMonodomainTotalLagrangianEM.so"\n'
        ");\n"
    )
    env = {"WM_OPTIONS": "myOptions", "PATH": ""}
    deps = {d.name: d for d in resolve_runtime_dependencies(case_root, env=env)}
    assert deps["manufacturedMonodomainTotalLagrangianEM"].path == built_lib
    assert deps["manufacturedMonodomainTotalLagrangianEM"].required is True


def test_a_case_local_library_path_with_unset_wm_options_is_unavailable_not_guessed(
    tmp_path: Path,
) -> None:
    """If the executor's own environment cannot resolve $WM_OPTIONS, the
    dependency must come back unavailable -- never a guessed path."""
    case_root = tmp_path / "case"
    (case_root / "system").mkdir(parents=True, exist_ok=True)
    (case_root / "system" / "controlDict").write_text(
        "libs\n(\n"
        '    "$FOAM_CASE/platforms/$WM_OPTIONS/lib/libmanufacturedMonodomainTotalLagrangianEM.so"\n'
        ");\n"
    )
    env = {"PATH": ""}  # WM_OPTIONS deliberately absent
    deps = {d.name: d for d in resolve_runtime_dependencies(case_root, env=env)}
    assert deps["manufacturedMonodomainTotalLagrangianEM"].path is None
    assert deps["manufacturedMonodomainTotalLagrangianEM"].required is True


def test_a_case_with_no_controldict_still_resolves_the_fixed_catalog(tmp_path: Path) -> None:
    env = {"PATH": ""}
    deps = resolve_runtime_dependencies(tmp_path / "case_without_controldict", env=env)
    assert {d.name for d in deps} >= set(_LIBRARY_CATALOG)


@skip_without_monorepo
def test_declared_library_catalog_covers_every_make_files_library() -> None:
    """Ratchet: a newly added library under src/*/Make/files must not be
    able to silently escape fingerprinting."""
    declared: set[str] = set()
    for make_file in (monorepo_root / "src").glob("*/Make/files"):
        text = make_file.read_text()
        match = re.search(r"^LIB\s*=\s*\S+/lib(\w+)\s*$", text, re.MULTILINE)
        if match:
            declared.add(match.group(1))

    assert declared, "expected at least one src/*/Make/files LIB declaration"
    missing = declared - set(_LIBRARY_CATALOG)
    assert not missing, (
        f"src/*/Make/files declares libraries the runtime-dependency catalog "
        f"does not cover: {sorted(missing)}. Add them to _LIBRARY_CATALOG in "
        "runtime_evidence.py so a rebuilt copy is fingerprinted."
    )
