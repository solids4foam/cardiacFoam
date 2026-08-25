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
#     runtime_evidence
#
# Description
#     Where cardiacFoam's runtime evidence lives: which steps actually solve,
#     where their logs land, what extra inputs belong in a provenance
#     snapshot, and how to read values out of solver-specific artifacts.
#     ``solve_step_commands``/``telemetry_source_globs``/
#     ``artifact_value_reader`` stay declaration-only for later phases;
#     ``extra_provenance_paths`` is consumed by Phase 2 now.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

from __future__ import annotations

import os
import re
import shutil
from collections.abc import Mapping
from pathlib import Path

from openfoam_driver.core.plugin_capabilities import RuntimeDependency

# Commands that run the solver. Deliberately excludes
# bathBidomainInterfaceMetrics, which is authorized to run but post-processes
# rather than solves -- expecting solver telemetry from it would be wrong.
_SOLVE_STEP_COMMANDS = frozenset({"cardiacFoam"})

# OpenFOAM's runApplication redirects solver output to log.<app>, so a step
# that runs an Allrun script produces no parseable driver-captured stdout.
# Phase 4's telemetry collector uses these globs to find the real log.
_TELEMETRY_GLOBS: dict[str, tuple[str, ...]] = {
    "Allrun": ("log.cardiacFoam", "log.*"),
    "Allrun.pre": ("log.*",),
    "Allrun.post": ("log.*",),
    "cardiacFoam": (),
}


def solve_step_commands() -> frozenset[str]:
    return _SOLVE_STEP_COMMANDS


def telemetry_source_globs(command: str) -> tuple[str, ...]:
    """Case-relative globs where ``command``'s solver log may land, beyond
    driver-captured stdout. Empty for commands that write no solver log."""
    return _TELEMETRY_GLOBS.get(command, ())


# Libraries cardiacFoam links or loads, and whether their absence is
# abnormal.
#
# required=True  -- linked unconditionally into cardiacFoam's EXE_LIBS
#                    (applications/solvers/cardiacFoam/Make/options); the
#                    solver cannot start without it.
# required=False -- either mode-dependent (exactly one of physicsModel /
#                    electroMechanicalModels is linked, selected by
#                    USE_LIGHTWEIGHT_PHYSICSMODEL at build time, and which
#                    one is baked into *this* installed binary cannot be
#                    determined from the outside) or not linked into
#                    cardiacFoam at all (verificationModels -- six
#                    manufactured-solution tutorials pull it in at runtime
#                    through their own controlDict libs (...) instead, which
#                    overrides this default to required=True for those
#                    cases; see _resolve_control_dict_libs_entry below).
#
# Verified live against a sourced OpenFOAM v2412 install (I3b), not assumed
# from Make/files: FOAM_MODULE_LIBBIN is unset in a plain sourced shell, and
# libphysicsModel is actually found in FOAM_USER_LIBBIN despite Make/files
# declaring FOAM_MODULE_LIBBIN. See _defined_lib_dirs -- the search below
# never maps a name to one fixed variable.
_LIBRARY_CATALOG: dict[str, bool] = {
    "electroModels": True,
    "ionicModels": True,
    "genericWriter": True,
    "activeTensionModels": True,
    "physicsModel": False,
    "electroMechanicalModels": False,
    "verificationModels": False,
}

# Search every *defined* one of these, never one fixed variable per name
# (I3b findings 1 and 2).
_LIB_DIR_ENV_VARS = ("FOAM_USER_LIBBIN", "FOAM_MODULE_LIBBIN", "FOAM_LIBBIN")

# The extension is platform-dependent -- .dylib here, .so on Linux (I3b
# finding 3). Try both rather than branching on sys.platform.
_LIB_EXTENSIONS = (".dylib", ".so")

_BARE_LIB_RE = re.compile(r"^lib(.+)\.(?:so|dylib)$")
_ENV_VAR_RE = re.compile(r"\$\{(\w+)\}|\$(\w+)")


def _defined_lib_dirs(env: Mapping[str, str]) -> tuple[Path, ...]:
    return tuple(Path(env[var]) for var in _LIB_DIR_ENV_VARS if env.get(var))


def _resolve_named_library(name: str, *, lib_dirs: tuple[Path, ...]) -> Path | None:
    for lib_dir in lib_dirs:
        for extension in _LIB_EXTENSIONS:
            candidate = lib_dir / f"lib{name}{extension}"
            if candidate.is_file():
                return candidate
    return None


def _parse_control_dict_libs(control_dict_path: Path) -> tuple[str, ...]:
    """Extract the entries of a controlDict's ``libs ( ... )`` list.

    Uses foamlib for structural, read-only parsing instead of hand-rolled
    paren-counting: this is a pure read of a list's shape with no write-
    back or provenance-digest concern (contrast core/runtime/mutators.py's
    read_foam_entry, which deliberately avoids foamlib because ITS values
    feed dict builders and must stay byte-verbatim), so none of that
    module's determinism caveats apply here -- foamlib parses in-process
    without evaluating ``#calc``/``#codeStream`` either way (see
    core/runtime/foam_backend.py's header comment), so it is safe on an
    arbitrary case's controlDict regardless.

    Unparseable or absent ``libs`` blocks yield an empty tuple rather than
    raising -- a case with no libs entry simply has none to resolve.
    """
    from foamlib import FoamFile

    try:
        libs = FoamFile(control_dict_path).get("libs")
    except (OSError, ValueError):
        return ()
    if not libs:
        return ()
    return tuple(_strip_quotes(str(entry)) for entry in libs)


def _strip_quotes(value: str) -> str:
    if len(value) >= 2 and value[0] == '"' and value[-1] == '"':
        return value[1:-1]
    return value


def _library_name_from_entry(entry: str) -> str:
    basename = entry.rsplit("/", 1)[-1]
    match = _BARE_LIB_RE.match(basename)
    return match.group(1) if match else basename


def _expand_openfoam_path(raw: str, *, case_root: Path, env: Mapping[str, str]) -> str | None:
    """Expand ``$FOAM_CASE`` and any other ``$VAR`` in a controlDict-style
    path using the executor's own environment.

    ``$FOAM_CASE`` is not a shell variable set by ``etc/bashrc`` -- OpenFOAM
    substitutes it at solver startup from the case being run (``-case``), so
    it is resolved here as the ``case_root`` this resolver was already given.
    Every other variable (``$WM_OPTIONS`` in practice) must come from
    ``env``; if it is undefined this returns ``None`` rather than guessing a
    value, so the caller can report the dependency unavailable instead of
    fabricating a path.
    """
    missing: list[str] = []

    def _substitute(match: "re.Match[str]") -> str:
        name = match.group(1) or match.group(2)
        if name == "FOAM_CASE":
            return str(case_root)
        value = env.get(name)
        if value is None:
            missing.append(name)
            return ""
        return value

    expanded = _ENV_VAR_RE.sub(_substitute, raw)
    return None if missing else expanded


def _resolve_control_dict_libs_entry(
    entry: str,
    *,
    case_root: Path,
    lib_dirs: tuple[Path, ...],
    env: Mapping[str, str],
) -> tuple[str, Path | None]:
    name = _library_name_from_entry(entry)
    if "/" not in entry:
        # A bare entry such as "libverificationModels.so" -- search the
        # defined lib directories exactly like the fixed catalog.
        return name, _resolve_named_library(name, lib_dirs=lib_dirs)

    # A case-local path such as
    # "$FOAM_CASE/platforms/$WM_OPTIONS/lib/lib<name>.so" (I2c).
    expanded = _expand_openfoam_path(entry, case_root=case_root, env=env)
    if expanded is None:
        return name, None
    candidate = Path(expanded)
    return name, candidate if candidate.is_file() else None


def resolve_runtime_dependencies(
    case_root: Path,
    *,
    env: Mapping[str, str] | None = None,
) -> tuple[RuntimeDependency, ...]:
    """The cardiacFoam binary, the libraries it links or loads, and anything
    this case's own controlDict pulls in through ``libs ( ... )``.

    This is the fix for the incident Phase 2 exists to prevent: most cases
    run through an ``Allrun`` script, so a workflow step's command
    fingerprints the script, never the solver binary the script invokes.
    Declaring the binary and its libraries here, independent of how a step
    happens to be launched, is what makes a rebuilt solver visible to
    provenance.

    ``env`` defaults to ``os.environ`` (the executor's own environment); a
    test may inject a synthetic mapping so these dependencies resolve
    without a sourced OpenFOAM install.
    """
    environment = os.environ if env is None else env
    case_root = Path(case_root)
    lib_dirs = _defined_lib_dirs(environment)

    dependencies: dict[str, RuntimeDependency] = {}

    solver_path = shutil.which("cardiacFoam", path=environment.get("PATH"))
    dependencies["cardiacFoam"] = RuntimeDependency(
        name="cardiacFoam",
        path=Path(solver_path) if solver_path else None,
        required=True,
    )

    for name, required in _LIBRARY_CATALOG.items():
        dependencies[name] = RuntimeDependency(
            name=name,
            path=_resolve_named_library(name, lib_dirs=lib_dirs),
            required=required,
        )

    control_dict = case_root / "system" / "controlDict"
    if control_dict.is_file():
        for entry in _parse_control_dict_libs(control_dict):
            name, resolved = _resolve_control_dict_libs_entry(
                entry, case_root=case_root, lib_dirs=lib_dirs, env=environment,
            )
            # A case that names a library in its own controlDict cannot run
            # without it -- required regardless of the fixed catalog's
            # default above (I2c: six manufactured-solution tutorials depend
            # on exactly this for libverificationModels).
            dependencies[name] = RuntimeDependency(name=name, path=resolved, required=True)

    return tuple(dependencies.values())


def extra_provenance_paths(
    case_root: Path,
    *,
    env: Mapping[str, str] | None = None,
) -> tuple[RuntimeDependency, ...]:
    """Runtime dependencies Phase 2 must fingerprint beyond ``system/`` and
    ``constant/``.

    Replaces the earlier ``tuple[Path, ...]`` stub, which could not express
    "this was required and I could not find it" -- a missing library would
    have been silently omitted rather than reported (I3)."""
    return resolve_runtime_dependencies(case_root, env=env)


def artifact_value_reader(artifact_format: str):
    """Reader for a solver-specific artifact format, or ``None``.

    Empty today. Phase 5 registers readers here for cardiac formats such as
    ECG traces and Purkinje time series. Returning ``None`` must make Phase 5
    report ``not_evaluated`` with a reason -- never an implicit pass."""
    del artifact_format
    return None
