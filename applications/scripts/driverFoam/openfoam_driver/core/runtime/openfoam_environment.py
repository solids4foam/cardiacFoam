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
#     openfoam_environment
#
# Description
#     Discovers and sources an OpenFOAM bashrc for driver-managed execution.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

from __future__ import annotations

import glob
import os
import shlex
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping


def _discover_openfoam_bashrcs() -> tuple[Path, ...]:
    """Discover OpenFOAM installations across common locations.

    Searches (in order):
      - /opt/openfoam* (Linux)
      - /usr/local/openfoam* (Linux/macOS)
      - /Volumes/OpenFOAM-v* (macOS)
      - Result of `which foamVersion` if available
    """
    candidates: set[Path] = set()

    for search_pattern in [
        "/opt/openfoam*/etc/bashrc",
        "/usr/local/openfoam*/etc/bashrc",
        "/Volumes/OpenFOAM-v*/etc/bashrc",
    ]:
        for path in glob.glob(search_pattern):
            candidates.add(Path(path))

    try:
        result = subprocess.run(
            ("which", "foamVersion"),
            capture_output=True,
            text=True,
            timeout=2,
            check=False,
        )
        if result.returncode == 0:
            foam_exe = Path(result.stdout.strip())
            if foam_exe.exists():
                wm_project_dir = foam_exe.parent.parent.parent
                bashrc = wm_project_dir / "etc" / "bashrc"
                if bashrc.exists():
                    candidates.add(bashrc)
    except (subprocess.TimeoutExpired, FileNotFoundError):
        pass

    return tuple(sorted(candidates, reverse=True))


@dataclass(frozen=True)
class OpenFOAMEnvironment:
    env: dict[str, str]
    bashrc: str | None = None
    error: str | None = None


def _parse_exported_environment(payload: str) -> dict[str, str]:
    env: dict[str, str] = {}
    for line in payload.splitlines():
        if not line.startswith("declare -x "):
            continue
        try:
            parts = shlex.split(line)
        except ValueError:
            continue
        if len(parts) < 3:
            continue
        assignment = parts[2]
        if "=" in assignment:
            key, value = assignment.split("=", 1)
        else:
            key, value = assignment, ""
        env[key] = value
    return env


def _candidate_bashrcs(
    *,
    explicit_bashrc: str | Path | None = None,
    base_env: Mapping[str, str] | None = None,
) -> tuple[Path, ...]:
    env = base_env or os.environ
    candidates: list[Path] = []

    if explicit_bashrc:
        return (Path(explicit_bashrc).expanduser(),)

    openfoam_bashrc = env.get("OPENFOAM_BASHRC")
    if openfoam_bashrc:
        candidates.append(Path(openfoam_bashrc).expanduser())

    wm_project_dir = env.get("WM_PROJECT_DIR")
    if wm_project_dir:
        candidates.append(Path(wm_project_dir).expanduser() / "etc" / "bashrc")

    candidates.extend(_discover_openfoam_bashrcs())

    seen: set[str] = set()
    unique: list[Path] = []
    for candidate in candidates:
        key = str(candidate)
        if key in seen:
            continue
        seen.add(key)
        unique.append(candidate)
    return tuple(unique)


def discover_openfoam_bashrc(
    *,
    explicit_bashrc: str | Path | None = None,
    base_env: Mapping[str, str] | None = None,
) -> Path | None:
    for candidate in _candidate_bashrcs(
        explicit_bashrc=explicit_bashrc,
        base_env=base_env,
    ):
        if candidate.is_file():
            return candidate
    return None


def load_openfoam_environment(
    *,
    explicit_bashrc: str | Path | None = None,
    base_env: Mapping[str, str] | None = None,
    driver_context: Any | None = None,
    timeout_s: float = 20.0,
) -> OpenFOAMEnvironment:
    """Return an environment suitable for strict OpenFOAM execution.

    If no bashrc can be found, return the current environment unchanged so the
    normal preflight diagnostics can report the missing OpenFOAM variables and
    executables.
    """
    env = dict(base_env or os.environ)
    if explicit_bashrc is None and driver_context is not None:
        bashrc_hook = getattr(driver_context.plugin, "get_openfoam_bashrc", None)
        if callable(bashrc_hook):
            try:
                explicit_bashrc = bashrc_hook(env)
            except Exception as exc:  # plugin diagnostics must not crash planning
                return OpenFOAMEnvironment(
                    env=env,
                    error=f"Plugin OpenFOAM runtime configuration failed: {exc}",
                )
    bashrc = discover_openfoam_bashrc(
        explicit_bashrc=explicit_bashrc,
        base_env=env,
    )
    if bashrc is None:
        if explicit_bashrc:
            return OpenFOAMEnvironment(
                env=env,
                error=f"OpenFOAM bashrc not found: {explicit_bashrc}",
            )
        sourced = OpenFOAMEnvironment(env=env)
        return _configure_plugin_environment(sourced, driver_context)

    script = (
        'set +e +u\n'
        'source "$_DRIVER_OPENFOAM_BASHRC" >/dev/null\n'
        'export -p > "$_DRIVER_ENV_FILE"'
    )
    env_file = tempfile.NamedTemporaryFile(delete=False)
    err_file = tempfile.NamedTemporaryFile(delete=False)
    env_file.close()
    err_file.close()
    temp_paths = (Path(env_file.name), Path(err_file.name))
    try:
        with open(err_file.name, "wb") as stderr_handle:
            source_env = {
                **env,
                "_DRIVER_OPENFOAM_BASHRC": str(bashrc),
                "_DRIVER_ENV_FILE": env_file.name,
            }
            completed = subprocess.run(
                ("bash", "-c", script),
                env=source_env,
                stdin=subprocess.DEVNULL,
                stdout=subprocess.DEVNULL,
                stderr=stderr_handle,
                timeout=timeout_s,
                check=False,
            )
    except subprocess.TimeoutExpired:
        for temp_path in temp_paths:
            temp_path.unlink(missing_ok=True)
        return OpenFOAMEnvironment(
            env=env,
            bashrc=str(bashrc),
            error=f"source {bashrc} timed out after {timeout_s:g} seconds",
        )
    except OSError as exc:
        for temp_path in temp_paths:
            temp_path.unlink(missing_ok=True)
        return OpenFOAMEnvironment(env=env, bashrc=str(bashrc), error=str(exc))

    if completed.returncode != 0:
        stderr = Path(err_file.name).read_text(errors="replace").strip()
        for temp_path in temp_paths:
            temp_path.unlink(missing_ok=True)
        return OpenFOAMEnvironment(
            env=env,
            bashrc=str(bashrc),
            error=stderr or f"source {bashrc} exited with {completed.returncode}",
        )

    raw_env = Path(env_file.name).read_text(errors="replace")
    for temp_path in temp_paths:
        temp_path.unlink(missing_ok=True)
    sourced_env = _parse_exported_environment(raw_env)

    sourced = OpenFOAMEnvironment(env=sourced_env, bashrc=str(bashrc))
    return _configure_plugin_environment(sourced, driver_context)


def _configure_plugin_environment(
    environment: OpenFOAMEnvironment,
    driver_context: Any | None,
) -> OpenFOAMEnvironment:
    """Apply an optional plugin-owned runtime environment contract.

    OpenFOAM sourcing is generic. Project-specific library selection belongs
    to the selected plugin, so the core only invokes the optional hook and
    carries its returned environment/error forward.
    """
    if environment.error or driver_context is None:
        return environment

    hook = getattr(driver_context.plugin, "configure_execution_environment", None)
    if hook is None:
        return environment

    try:
        configured_env, error = hook(dict(environment.env))
    except Exception as exc:  # plugin diagnostics must not crash planning
        return OpenFOAMEnvironment(
            env=environment.env,
            bashrc=environment.bashrc,
            error=f"Plugin runtime environment configuration failed: {exc}",
        )

    return OpenFOAMEnvironment(
        env=dict(configured_env),
        bashrc=environment.bashrc,
        error=error,
    )


def configure_plugin_environment(
    env: Mapping[str, str],
    driver_context: Any | None,
) -> OpenFOAMEnvironment:
    """Apply a plugin environment contract without sourcing OpenFOAM again."""
    return _configure_plugin_environment(
        OpenFOAMEnvironment(env=dict(env)),
        driver_context,
    )
