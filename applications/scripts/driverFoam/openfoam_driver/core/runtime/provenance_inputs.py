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
#     provenance_inputs
#
# Description
#     Canonical input enumeration: which files, case scripts, and runtime
#     dependencies a workflow run actually consumes, so Task 3 can
#     checkpoint the answer. Reuses Task 1's snapshot model
#     (component_for_path) and Task 2a's runtime-dependency composition
#     (component_for_runtime_dependency) rather than building a parallel
#     one, and resolves step executables through the executor's own
#     _resolve_command rather than a second PATH rule.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

"""Enumerate a case's canonical provenance inputs.

Classification is by **consumption, not authorship** (I1). ``constant/polyMesh``
is written by ``blockMesh`` and then read by the solver -- generated, and
still a mandatory input. The resolution precedence, first match wins:

1. A DAG step's ``consumes`` declaration.
2. A plugin ``CaseProvenanceCapability.required_inputs()`` entry.
3. A plugin ``CaseProvenanceCapability.generated_output_globs()`` match
   (excludes).
4. Otherwise: **required input**. This is the safety property, not laziness
   -- an unclassified file that turns out to matter causes a spurious
   refusal (recoverable); the reverse is the silent stale replay this whole
   phase exists to prevent.

Concretely: ``system/**``, ``constant/**``, the **selected** start-time
directory (read from ``system/controlDict``'s ``startFrom``/``startTime`` --
not always ``0/``), and ``processor*/<selected-time>/**`` during a decomposed
restart (I9) are walked and classified by that precedence. The exact
Allrun-family scripts named by the DAG are added directly. Other time
directories, ``postProcessing/`` (which holds ``workflow_logs/``), and
state/manifest files are never walked, so they are excluded by construction
rather than by an exclusion rule that could be gotten wrong.

Step executables -- including an MPI launcher's payload
(``mpirun -np 4 cardiacFoam`` -> both) -- are resolved through
``workflow_runner._resolve_command``, the executor's own resolution, so a
provenance digest is never computed against a different binary than the one
that actually runs. A bare command not found locally is looked up on
``PATH`` (via the same MPI-payload-unwrapping rule
``environment_preflight`` already uses) and folded into a
:class:`RuntimeDependency`, then fingerprinted through Task 2a's
``component_for_runtime_dependency`` -- so a required-but-unresolved
executable surfaces as ``unavailable``, never silently omitted, exactly like
a plugin-declared library.
"""

from __future__ import annotations

import fnmatch
import os
import shlex
import shutil
from dataclasses import replace
from pathlib import Path
from typing import Any, Iterable, Mapping, TYPE_CHECKING

from ..plugin_capabilities import ResolvedInput, RuntimeDependency
from .environment_preflight import _MPI_LAUNCHERS, _unwrap_mpi_program
from .mutators import read_foam_entry
from .provenance import ProvenanceComponent, component_for_path
from .provenance_dependencies import component_for_runtime_dependency
from .workflow import CASE_SCRIPT_COMMANDS
from .workflow_runner import _resolve_case_cwd, _resolve_command

if TYPE_CHECKING:
    from ..plugin_interface import DriverContext

_DEFAULT_START_TIME = "0"


def _select_start_time(case_root: Path) -> str:
    """The **selected** start-time directory, per ``system/controlDict``.

    ``startFrom latestTime`` selects the latest written time, not ``0``;
    ``firstTime`` selects the earliest. Anything else (including an absent
    ``controlDict``) falls back to the literal ``startTime`` value, or
    ``"0"`` if that too is absent -- the common case, but never assumed
    without checking.
    """
    control_dict = case_root / "system" / "controlDict"
    start_from = (read_foam_entry(control_dict, "startFrom") or "startTime").strip()

    if start_from in ("latestTime", "firstTime"):
        candidates = _list_time_dir_names(case_root)
        if not candidates:
            return _DEFAULT_START_TIME
        selector = max if start_from == "latestTime" else min
        return selector(candidates, key=float)

    start_time = read_foam_entry(control_dict, "startTime")
    return start_time.strip() if start_time is not None else _DEFAULT_START_TIME


def _list_time_dir_names(case_root: Path) -> list[str]:
    names: list[str] = []
    try:
        children = list(case_root.iterdir())
    except OSError:
        return names
    for child in children:
        if not child.is_dir():
            continue
        try:
            float(child.name)
        except ValueError:
            continue
        names.append(child.name)
    return names


def _walk_files(root: Path) -> list[Path]:
    """Every file (or symlink, valid or broken) under ``root``, recursively.

    Plain sub-directories are not emitted -- only leaves. A symlink is always
    emitted regardless of what it resolves to, so a required-but-broken link
    still surfaces (via ``component_for_path``'s own unavailable handling)
    instead of silently vanishing from the walk.
    """
    if not root.is_dir():
        return []
    found: list[Path] = []
    for path in root.rglob("*"):
        if path.is_dir() and not path.is_symlink():
            continue
        found.append(path)
    return sorted(found)


def _matches_any_glob(rel_path: str, patterns: Iterable[str]) -> bool:
    return any(fnmatch.fnmatch(rel_path, pattern) for pattern in patterns)


def _collect_consumed_relpaths(workflow_dag: Mapping[str, Any] | None) -> set[str]:
    relpaths: set[str] = set()
    for step in (workflow_dag or {}).get("steps", ()):
        for entry in step.get("consumes", ()) or ():
            relpaths.add(str(entry))
    return relpaths


def _component_for_resolved_input(
    resolved_input: ResolvedInput, case_root: Path
) -> ProvenanceComponent | None:
    """A plugin-resolved required input, fingerprinted -- or, for a
    ``required`` input that failed to resolve, an explicit ``unavailable``
    component (I2/I5). A ``required=False`` input with no resolved path was
    genuinely absent and optional: nothing was going to be consumed, so
    nothing is fingerprinted."""
    path = resolved_input.path
    if path is None:
        if not resolved_input.required:
            return None
        return ProvenanceComponent(
            kind="case_file",
            path=resolved_input.name,
            role="required_input",
            method="unavailable",
            strength="unavailable",
        )
    try:
        path.relative_to(case_root)
    except ValueError:
        # Resolved outside the case tree entirely (not even via a symlink
        # component_for_path would classify as external_link) -- identify it
        # by the plugin's own declared name, the same pattern
        # provenance_dependencies.py uses for runtime dependencies.
        component = component_for_path(path, kind="case_file", relative_to=path.parent)
        return replace(component, path=resolved_input.name)
    return component_for_path(path, kind="case_file", relative_to=case_root)


def _is_case_local_script(
    command: str, executable: str, *, resolved_cwd: Path, case_root: Path
) -> Path | None:
    """If ``_resolve_command`` would run a script from inside the case tree
    for this exact ``command``, return its on-disk path; else ``None``.

    Mirrors ``_resolve_command``'s own precondition (an explicit path, or a
    recognized case-script name) rather than re-deciding independently --
    checking existence for every bare command name here would find files the
    real executor's PATH-only bare-name rule would never run.
    """
    if "/" not in command and command not in CASE_SCRIPT_COMMANDS:
        return None
    candidate = Path(executable)
    if not candidate.is_absolute():
        candidate = resolved_cwd / candidate
    if not candidate.is_file():
        return None
    try:
        within_case = candidate.resolve().is_relative_to(case_root.resolve())
    except OSError:
        return None
    return candidate if within_case else None


def _resolve_dependency_path(
    name: str, *, resolved_cwd: Path, env: Mapping[str, str]
) -> Path | None:
    """Where a step executable that is *not* a case-local script lives:
    an explicit path resolved against the step's cwd, or a PATH lookup for a
    bare name -- the same two cases ``_resolve_command`` distinguishes."""
    if "/" in name:
        candidate = Path(name)
        if not candidate.is_absolute():
            candidate = resolved_cwd / candidate
        return candidate if candidate.is_file() else None
    # An explicit empty/absent PATH must mean "nothing is found", matching
    # what the real subprocess call would see -- shutil.which's own fallback
    # to the real os.environ on a bare ``None`` would silently defeat a
    # caller's isolated ``env``.
    found = shutil.which(name, path=env.get("PATH", ""))
    return Path(found) if found else None


def _step_command_and_args(step: Mapping[str, Any]) -> tuple[str, tuple[str, ...]] | None:
    raw_command = str(step.get("command", "")).strip()
    if not raw_command:
        return None
    try:
        command, *inline_args = shlex.split(raw_command)
    except ValueError:
        command, *inline_args = raw_command.split()
    if not command:
        return None
    args = tuple(inline_args) + tuple(str(arg) for arg in step.get("args", ()))
    return command, args


def _register_step_executable(
    name: str,
    *,
    resolved_cwd: Path,
    case_root: Path,
    env: Mapping[str, str],
    add_case_file: "_ComponentAdder",
    dependencies: dict[str, RuntimeDependency],
) -> None:
    executable = _resolve_command(name, resolved_cwd)
    local_script = _is_case_local_script(
        name, executable, resolved_cwd=resolved_cwd, case_root=case_root
    )
    if local_script is not None:
        add_case_file(component_for_path(local_script, kind="case_file", relative_to=case_root))
        return
    dependencies[name] = RuntimeDependency(
        name=name,
        path=_resolve_dependency_path(name, resolved_cwd=resolved_cwd, env=env),
        required=True,
    )


class _ComponentAdder:
    """Callable wrapper so a bound helper can carry a type hint above."""

    def __init__(self, sink: dict[tuple[str, str], ProvenanceComponent]) -> None:
        self._sink = sink

    def __call__(self, component: ProvenanceComponent | None) -> None:
        if component is not None:
            self._sink[(component.kind, component.path)] = component


def enumerate_case_inputs(
    case_root: Path,
    *,
    workflow_dag: dict[str, Any],
    driver_context: "DriverContext",
    env: Mapping[str, str] | None = None,
) -> tuple[ProvenanceComponent, ...]:
    """Enumerate every provenance input a run of ``workflow_dag`` against
    ``case_root`` consumes: case files, case-local scripts, step executables,
    and plugin-declared runtime dependencies. Excludes generated outputs.

    Never raises on an incomplete or minimal case (e.g. no ``constant/``) --
    every filesystem read here degrades to an ``unavailable`` component via
    ``component_for_path``/``component_for_runtime_dependency`` rather than
    propagating an exception, consistent with those functions' own contract.
    """
    case_root = Path(case_root)
    environment = dict(os.environ) if env is None else dict(env)
    capabilities = driver_context.capabilities

    resolved_case = capabilities.case_introspection.resolve_case_models(case_root)
    selected_start_time = _select_start_time(case_root)

    consumed_relpaths = _collect_consumed_relpaths(workflow_dag)
    required_inputs = capabilities.case_provenance.required_inputs(
        case_root, resolved_case, selected_start_time
    )
    generated_globs = capabilities.case_provenance.generated_output_globs(
        case_root, resolved_case, selected_start_time
    )

    components: dict[tuple[str, str], ProvenanceComponent] = {}
    add = _ComponentAdder(components)

    # -- system/**, constant/**, the selected start-time dir, and
    # processor*/<selected-time>/** (I9) -- classified by precedence steps
    # 1 (consumes), 3 (generated_output_globs), 4 (fallback required). Step 2
    # (plugin required_inputs) is applied uniformly below instead, since a
    # resolved input's path need not fall under any of these directories.
    walk_roots = [case_root / "system", case_root / "constant", case_root / selected_start_time]
    for processor_dir in sorted(case_root.glob("processor*")):
        if processor_dir.is_dir():
            walk_roots.append(processor_dir / selected_start_time)

    for root in walk_roots:
        for path in _walk_files(root):
            rel = path.relative_to(case_root).as_posix()
            if rel in consumed_relpaths:
                add(component_for_path(path, kind="case_file", relative_to=case_root))
                continue
            if _matches_any_glob(rel, generated_globs):
                continue
            add(component_for_path(path, kind="case_file", relative_to=case_root))

    # -- Precedence step 2: plugin-resolved required inputs always win,
    # regardless of whether the walk above already included or excluded
    # them (e.g. a field resolved outside the selected time via backward
    # findInstance, or one explicitly re-included over a generated glob).
    for resolved_input in required_inputs:
        add(_component_for_resolved_input(resolved_input, case_root))

    # -- Precedence step 1, applied last so it always wins even for a path
    # the walk above excluded or never visited at all.
    for rel in consumed_relpaths:
        add(component_for_path(case_root / rel, kind="case_file", relative_to=case_root))

    # -- step executables, resolved the same way the executor resolves them
    # (workflow_runner._resolve_command), including an MPI launcher's
    # payload. Case-local scripts land as case_file components above; every
    # other command becomes a RuntimeDependency, deduplicated by name so a
    # binary named by two steps (or already declared by the plugin below)
    # fingerprints once.
    dependencies: dict[str, RuntimeDependency] = {}
    for step in (workflow_dag or {}).get("steps", ()):
        parsed = _step_command_and_args(step)
        if parsed is None:
            continue
        command, args = parsed
        resolved_cwd = _resolve_case_cwd(case_root, str(step.get("cwd", ".")))
        _register_step_executable(
            command,
            resolved_cwd=resolved_cwd,
            case_root=case_root,
            env=environment,
            add_case_file=add,
            dependencies=dependencies,
        )
        if command in _MPI_LAUNCHERS:
            payload = _unwrap_mpi_program(args)
            if payload:
                _register_step_executable(
                    payload,
                    resolved_cwd=resolved_cwd,
                    case_root=case_root,
                    env=environment,
                    add_case_file=add,
                    dependencies=dependencies,
                )

    # -- plugin-declared runtime dependencies (Task 2a): the solver binary,
    # its libraries, and anything the case's own controlDict libs (...)
    # pulls in. Authoritative over the generic PATH-only resolution above --
    # it knows the library search directories and required/optional split a
    # bare PATH lookup cannot.
    for dependency in capabilities.runtime_evidence.extra_provenance_paths(case_root):
        dependencies[dependency.name] = dependency

    for dependency in dependencies.values():
        add(component_for_runtime_dependency(dependency))

    return tuple(sorted(components.values(), key=lambda component: (component.kind, component.path)))
