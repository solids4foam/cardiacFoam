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
#     environment_preflight
#
# Description
#     Runtime environment diagnostics for strict workflow execution.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

from __future__ import annotations

import os
import shlex
import shutil
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from openfoam_driver.core.planning_types import StrictDiagnostic, diagnostic
from .openfoam_environment import load_openfoam_environment
from .workflow import CASE_SCRIPT_COMMANDS


_MPI_LAUNCHERS = frozenset({"mpirun", "mpiexec", "orterun"})
_INTERPRETER_SKIP = frozenset({"python", "python3"})
_MPI_VALUE_FLAGS = frozenset({"-np", "-n", "--np"})


def _unwrap_mpi_program(args: tuple[str, ...]) -> str | None:
    """Return the wrapped program from an MPI launcher's args, or None."""
    index = 0
    while index < len(args):
        token = args[index]
        if token in _MPI_VALUE_FLAGS:
            index += 2
            continue
        if token.startswith("-"):
            index += 1
            continue
        return token
    return None


@dataclass(frozen=True)
class _ExecutableRequirements:
    executables: tuple[str, ...]
    is_parallel: bool
    mpi_launcher_in_dag: bool


def _required_executables(workflow_dag: dict[str, Any] | None) -> _ExecutableRequirements:
    """Derive the executables a plan will invoke from its workflow DAG."""
    executables: list[str] = []
    is_parallel = False
    mpi_launcher_in_dag = False

    def _add(name: str) -> None:
        if (
            name
            and name not in _INTERPRETER_SKIP
            and name not in CASE_SCRIPT_COMMANDS
            and name not in executables
        ):
            executables.append(name)

    for step in (workflow_dag or {}).get("steps", ()):
        raw_command = str(step.get("command", "")).strip()
        args = tuple(str(arg) for arg in step.get("args", ()))
        if not raw_command:
            continue
        # A step may carry its entire invocation in "command", for example
        # "postProcess -func Niedererpoints -latestTime". Only the leading
        # token names an executable; the remainder are arguments and must not
        # be passed to a PATH lookup.
        try:
            command, *inline_args = shlex.split(raw_command)
        except ValueError:
            command, *inline_args = raw_command.split()
        if not command:
            continue
        args = tuple(inline_args) + args
        if command in _MPI_LAUNCHERS:
            is_parallel = True
            mpi_launcher_in_dag = True
            _add(command)
            wrapped = _unwrap_mpi_program(args)
            if wrapped is not None:
                _add(wrapped)
            continue
        if command == "decomposePar" or "-parallel" in args:
            is_parallel = True
        _add(command)

    return _ExecutableRequirements(
        executables=tuple(executables),
        is_parallel=is_parallel,
        mpi_launcher_in_dag=mpi_launcher_in_dag,
    )


_SOURCE_SUFFIXES = frozenset({".C", ".H", ".cu", ".cuh"})


def _discover_src_root() -> Path | None:
    """Locate the repository ``src/`` tree from this module's location, or
    ``None`` if it cannot be found (e.g. an installed-package layout)."""
    for parent in Path(__file__).resolve().parents:
        if (parent / "src").is_dir() and (parent / "tutorials").is_dir():
            return parent / "src"
    return None


def _newest_source_mtime(src_root: Path) -> float | None:
    """Return the mtime of the most recently modified C++/CUDA source under
    ``src_root``, or ``None`` if there are no source files."""
    newest: float | None = None
    for dirpath, _dirnames, filenames in os.walk(src_root):
        for name in filenames:
            if os.path.splitext(name)[1] in _SOURCE_SUFFIXES:
                mtime = os.path.getmtime(os.path.join(dirpath, name))
                if newest is None or mtime > newest:
                    newest = mtime
    return newest


def _build_staleness_diagnostics(
    workflow_dag: dict[str, Any] | None,
    checked_env: dict[str, str],
    *,
    src_root: Path | str | None,
) -> tuple[StrictDiagnostic, ...]:
    """Warn (never block) when a user-compiled utility the plan invokes is older
    than the newest C++/CUDA source under ``src_root`` -- i.e. the binary was
    not rebuilt after the source changed (the classic stale-``libso`` footgun).

    Only executables that resolve under ``$FOAM_USER_APPBIN`` are policed; core
    OpenFOAM apps and system binaries are never flagged.
    """
    if src_root is None:
        return ()
    src_root = Path(src_root)
    if not src_root.exists():
        return ()

    user_appbin = checked_env.get("FOAM_USER_APPBIN")
    if not user_appbin:
        return ()
    user_appbin_resolved = Path(user_appbin).resolve()

    # Find plan binaries that actually live under the user appbin before doing
    # the (potentially large) source-tree walk.
    candidates: list[tuple[str, str]] = []
    path = checked_env.get("PATH")
    for executable in _required_executables(workflow_dag).executables:
        resolved = shutil.which(executable, path=path)
        if not resolved:
            continue
        try:
            Path(resolved).resolve().relative_to(user_appbin_resolved)
        except ValueError:
            continue
        candidates.append((executable, resolved))

    if not candidates:
        return ()

    newest_source = _newest_source_mtime(src_root)
    if newest_source is None:
        return ()

    diagnostics: list[StrictDiagnostic] = []
    for executable, resolved in candidates:
        if os.path.getmtime(resolved) < newest_source:
            diagnostics.append(diagnostic(
                "warning",
                "stale_build",
                f"{executable} is older than the newest source under {src_root}; "
                "rebuild (e.g. wmake / wmake libso) before running.",
                source="environment",
                field=executable,
            ))
    return tuple(diagnostics)


def _environment_diagnostics(
    workflow_dag: dict[str, Any] | None,
    *,
    env: dict[str, str] | None = None,
    openfoam_bashrc: str | None = None,
    driver_context: Any | None = None,
) -> tuple[StrictDiagnostic, ...]:
    """Preflight the runtime environment against the plan's actual commands."""
    if "SKIP_ENV_DIAGNOSTICS" in os.environ:
        return ()
    diagnostics: list[StrictDiagnostic] = []
    checked_env = env
    loaded_environment = None
    if checked_env is None:
        loaded_environment = load_openfoam_environment(
            explicit_bashrc=openfoam_bashrc,
            driver_context=driver_context,
        )
        checked_env = loaded_environment.env

    if loaded_environment is not None and loaded_environment.error:
        diagnostics.append(diagnostic(
            "error",
            "openfoam_env_source_failed",
            loaded_environment.error,
            source="environment",
            field=loaded_environment.bashrc or openfoam_bashrc or "",
        ))

    if "WM_PROJECT_DIR" not in checked_env:
        diagnostics.append(diagnostic(
            "error",
            "missing_openfoam_env",
            "WM_PROJECT_DIR is not set. OpenFOAM environment not sourced.",
            source="environment",
        ))
    else:
        for var in ("WM_PROJECT_VERSION", "FOAM_USER_LIBBIN"):
            if var not in checked_env:
                diagnostics.append(diagnostic(
                    "warning",
                    "partial_openfoam_env",
                    f"{var} is not set. OpenFOAM environment may be partially sourced.",
                    source="environment",
                    field=var,
                ))

    requirements = _required_executables(workflow_dag)
    for executable in requirements.executables:
        if not shutil.which(executable, path=checked_env.get("PATH")):
            diagnostics.append(diagnostic(
                "error",
                "missing_executable",
                f"{executable} not found on PATH.",
                source="environment",
                field=executable,
            ))

    if (
        requirements.is_parallel
        and not requirements.mpi_launcher_in_dag
        and not (
            shutil.which("mpirun", path=checked_env.get("PATH"))
            or shutil.which("mpiexec", path=checked_env.get("PATH"))
        )
    ):
        diagnostics.append(diagnostic(
            "error",
            "missing_mpi",
            "Plan is parallel but no MPI launcher (mpirun/mpiexec) found on PATH.",
            source="environment",
            field="mpirun",
        ))

    diagnostics.extend(
        _build_staleness_diagnostics(
            workflow_dag, checked_env, src_root=_discover_src_root()
        )
    )

    return tuple(diagnostics)
