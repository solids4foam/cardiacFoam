"""Machine-readable cardiacFoam/OpenFOAM runtime backend contract."""

from __future__ import annotations

import hashlib
import json
import os
import shutil
import subprocess
from pathlib import Path
from typing import Any, Mapping

import yaml


_PLUGIN_ID = "org.cardiacfoam"
_PROFILE_PATH = Path(__file__).with_name("plugin.yaml")
_RUNTIME_CONFIG_ENV = "DRIVERFOAM_RUNTIME_CONFIG"
_BACKEND_ENV = "DRIVERFOAM_CARDIACFOAM_BACKEND"
_SOLIDS_ROOT_ENV = "DRIVERFOAM_CARDIACFOAM_SOLIDS4FOAM_ROOT"
_MANIFEST_ENV = "DRIVERFOAM_CARDIACFOAM_BUILD_MANIFEST"
_LIBRARY_EXTENSIONS = ("dylib", "so")


def _profile_contract() -> dict[str, Any]:
    payload = yaml.safe_load(_PROFILE_PATH.read_text(encoding="utf-8"))
    try:
        contract = payload["runtime"]["backend"]
    except (KeyError, TypeError) as exc:
        raise ValueError(f"CardiacFoam plugin profile has no runtime backend contract: {exc}") from exc
    if not isinstance(contract, dict) or not isinstance(contract.get("options"), dict):
        raise ValueError("CardiacFoam runtime.backend.options must be a mapping")
    return contract


def _user_selection(env: Mapping[str, str]) -> dict[str, Any]:
    config_name = env.get(_RUNTIME_CONFIG_ENV)
    if not config_name:
        return {}
    config_path = Path(os.path.expandvars(config_name)).expanduser().resolve()
    if not config_path.is_file():
        raise ValueError(f"Configured driverFOAM runtime file does not exist: {config_path}")
    try:
        payload = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    except yaml.YAMLError as exc:
        raise ValueError(f"Invalid driverFOAM runtime YAML {config_path}: {exc}") from exc
    if not isinstance(payload, dict):
        raise ValueError(f"DriverFOAM runtime file must contain a mapping: {config_path}")
    plugins = payload.get("plugins", {})
    if not isinstance(plugins, dict):
        raise ValueError("DriverFOAM runtime 'plugins' must be a mapping")
    selection = plugins.get(_PLUGIN_ID, {})
    if not isinstance(selection, dict):
        raise ValueError(f"DriverFOAM runtime selection for {_PLUGIN_ID} must be a mapping")
    return selection


def configured_openfoam_bashrc(env: Mapping[str, str]) -> str | None:
    """Return the OpenFOAM bashrc from the selected host runtime file."""
    config_name = env.get(_RUNTIME_CONFIG_ENV)
    if not config_name:
        return None
    config_path = Path(os.path.expandvars(config_name)).expanduser().resolve()
    if not config_path.is_file():
        return None
    payload = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        return None
    openfoam = payload.get("openfoam", {})
    value = openfoam.get("bashrc") if isinstance(openfoam, dict) else None
    return str(value) if value else None


def configure_runtime_environment(env: Mapping[str, str]) -> tuple[dict[str, str], str | None]:
    """Validate and export the selected cardiacFoam runtime backend."""
    configured_env = dict(env)
    try:
        contract = _profile_contract()
        selection = _user_selection(configured_env)
    except (OSError, ValueError, yaml.YAMLError) as exc:
        return configured_env, str(exc)

    backend = selection.get("backend") or configured_env.get(_BACKEND_ENV)
    options = contract["options"]
    if not isinstance(backend, str) or backend not in options:
        return configured_env, (
            f"{_BACKEND_ENV} must select one of: {', '.join(sorted(options))}. "
            "Set it directly or provide DRIVERFOAM_RUNTIME_CONFIG."
        )

    option = options[backend]
    if not isinstance(option, dict):
        return configured_env, f"Invalid cardiacFoam backend declaration: {backend}"

    solids_root_value = selection.get("solids4foam_root") or configured_env.get(_SOLIDS_ROOT_ENV)
    solids_root: Path | None = None
    if backend == "full":
        if not solids_root_value:
            return configured_env, (
                f"Full cardiacFoam backend requires {_SOLIDS_ROOT_ENV} or "
                "plugins.org.cardiacfoam.solids4foam_root."
            )
        solids_root = Path(os.path.expandvars(str(solids_root_value))).expanduser().resolve()
        source_header = solids_root / "src/solids4FoamModels/physicsModel/physicsModel.H"
        ln_include = solids_root / "src/solids4FoamModels/lnInclude/physicsModel.H"
        if not solids_root.is_dir() or not source_header.is_file():
            return configured_env, f"Invalid solids4foam source root: {solids_root}"
        if not ln_include.is_file():
            return configured_env, f"solids4foam is not built (missing {ln_include})"

    manifest_value = (
        selection.get("build_manifest")
        or configured_env.get(_MANIFEST_ENV)
        or (
            str(Path(configured_env["FOAM_USER_LIBBIN"]) / "cardiacFoam.build.json")
            if configured_env.get("FOAM_USER_LIBBIN")
            else None
        )
    )
    if not manifest_value:
        return configured_env, (
            "No cardiacFoam build manifest configured. Set "
            f"{_MANIFEST_ENV} or plugins.org.cardiacfoam.build_manifest."
        )

    manifest_path = Path(os.path.expandvars(str(manifest_value))).expanduser().resolve()

    regeneration_error = _ensure_build_manifest(manifest_path, configured_env, contract, solids_root)
    if regeneration_error:
        return configured_env, regeneration_error

    error = _validate_build_manifest(
        manifest_path,
        backend=backend,
        openfoam_root=Path(configured_env["WM_PROJECT_DIR"]).resolve()
        if configured_env.get("WM_PROJECT_DIR") else None,
        solids_root=solids_root,
        required_libraries=tuple(option.get("required_libraries", ())),
        forbidden_libraries=tuple(option.get("forbidden_libraries", ())),
    )
    if error:
        return configured_env, error

    configured_env[_BACKEND_ENV] = backend
    configured_env[_MANIFEST_ENV] = str(manifest_path)
    if solids_root is not None:
        configured_env[_SOLIDS_ROOT_ENV] = str(solids_root)
        configured_env["SOLIDS4FOAM_INST_DIR"] = str(solids_root)
    return configured_env, None


def _library_search_dirs(env: Mapping[str, str]) -> tuple[Path, ...]:
    dirs = []
    for key in ("FOAM_USER_LIBBIN", "FOAM_MODULE_LIBBIN", "FOAM_LIBBIN"):
        value = env.get(key)
        if value:
            dirs.append(Path(value))
    return tuple(dirs)


def _find_library(bare_name: str, search_dirs: tuple[Path, ...]) -> Path | None:
    for directory in search_dirs:
        for extension in _LIBRARY_EXTENSIONS:
            candidate = directory / f"lib{bare_name}.{extension}"
            if candidate.is_file():
                return candidate
    return None


def _linked_library_names(binary: Path) -> tuple[str, ...]:
    """Return the basenames of shared libraries linked into `binary`.

    Uses `otool -L` on macOS or `ldd` on Linux — whichever is on PATH.
    """
    if shutil.which("otool"):
        result = subprocess.run(["otool", "-L", str(binary)], capture_output=True, text=True)
        if result.returncode != 0:
            return ()
        lines = result.stdout.splitlines()[1:]
        return tuple(Path(line.split()[0]).name for line in lines if line.strip())
    if shutil.which("ldd"):
        result = subprocess.run(["ldd", str(binary)], capture_output=True, text=True)
        if result.returncode != 0:
            return ()
        names = []
        for line in result.stdout.splitlines():
            lhs = line.strip().split(" =>")[0].strip()
            if lhs.startswith("lib"):
                names.append(Path(lhs).name)
        return tuple(names)
    return ()


def _infer_backend(linked: tuple[str, ...], options: Mapping[str, Any]) -> str | None:
    """Infer which backend was compiled from the solver's linked libraries.

    Each backend option in the plugin contract declares the libraries that
    must (`required_libraries`) and must not (`forbidden_libraries`) appear
    linked into the solver. Exactly one backend matching both conditions is
    the compiled backend; anything else (zero or more than one match) is
    ambiguous.
    """
    matches = []
    for backend_name, option in options.items():
        if not isinstance(option, dict):
            continue
        required = [Path(lib).name for lib in option.get("required_libraries", ())]
        forbidden = [Path(lib).name for lib in option.get("forbidden_libraries", ())]
        if not required:
            continue
        has_all_required = all(any(name.startswith(req) for name in linked) for req in required)
        has_any_forbidden = any(any(name.startswith(forb) for name in linked) for forb in forbidden)
        if has_all_required and not has_any_forbidden:
            matches.append(backend_name)
    return matches[0] if len(matches) == 1 else None


def _ensure_build_manifest(
    manifest_path: Path,
    env: Mapping[str, str],
    contract: Mapping[str, Any],
    solids_root_hint: Path | None,
) -> str | None:
    """Regenerate `cardiacFoam.build.json` from the compiled artifacts.

    Runs only when the manifest is missing or older than the compiled
    solver, so a fresh (re)build is always reflected without any build-time
    hook. The backend is inferred from what the solver actually links
    against, not asserted by the caller, so the manifest can never disagree
    with reality — it IS the inspection of reality.

    Returns an error string if regeneration was attempted and failed.
    Returns None if regeneration succeeded, was unnecessary (already
    up to date), or could not be attempted (e.g. nothing compiled yet) —
    in the last case the existing missing/stale manifest is reported by the
    caller's own validation step.
    """
    user_appbin = env.get("FOAM_USER_APPBIN")
    if not user_appbin:
        return None
    solver = Path(user_appbin) / "cardiacFoam"
    if not solver.is_file():
        return None

    if manifest_path.is_file() and manifest_path.stat().st_mtime >= solver.stat().st_mtime:
        return None

    linked = _linked_library_names(solver)
    if not linked:
        return (
            f"Could not determine libraries linked into {solver} "
            "(otool/ldd unavailable or produced no output)."
        )

    options = contract["options"]
    backend = _infer_backend(linked, options)
    if backend is None:
        return (
            f"Could not infer the cardiacFoam build backend from {solver}'s "
            f"linked libraries: {', '.join(sorted(linked))}"
        )

    search_dirs = _library_search_dirs(env)
    library_names = tuple(contract.get("common_libraries", ())) + tuple(
        options[backend].get("required_libraries", ())
    )
    artifact_paths: list[tuple[str, Path]] = [("cardiacFoam", solver)]
    for library_name in library_names:
        path = _find_library(library_name.removeprefix("lib"), search_dirs)
        if path is None:
            return f"Required library {library_name} was not found for backend {backend!r}."
        artifact_paths.append((library_name, path))

    solids_root: Path | None = None
    source_revision: str | None = None
    if backend == "full":
        solids_root = solids_root_hint
        if solids_root is None and env.get("SOLIDS4FOAM_INST_DIR"):
            solids_root = Path(env["SOLIDS4FOAM_INST_DIR"])
        if solids_root is not None and (solids_root / ".git").exists():
            result = subprocess.run(
                ["git", "-C", str(solids_root), "rev-parse", "HEAD"],
                capture_output=True,
                text=True,
            )
            if result.returncode == 0:
                source_revision = result.stdout.strip() or None

    payload = {
        "schema_version": 1,
        "plugin": _PLUGIN_ID,
        "backend": backend,
        "openfoam": {
            "root": str(Path(env["WM_PROJECT_DIR"]).resolve()) if env.get("WM_PROJECT_DIR") else None,
            "version": env.get("WM_PROJECT_VERSION", ""),
            "options": env.get("WM_OPTIONS", ""),
        },
        "solids4foam": {
            "root": str(solids_root.resolve()) if solids_root else None,
            "revision": source_revision,
        },
        "linked_libraries": sorted(linked),
        "artifacts": [
            {
                "name": name,
                "path": str(path),
                "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            }
            for name, path in artifact_paths
        ],
    }
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    manifest_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    return None


def _validate_build_manifest(
    path: Path,
    *,
    backend: str,
    openfoam_root: Path | None,
    solids_root: Path | None,
    required_libraries: tuple[str, ...],
    forbidden_libraries: tuple[str, ...],
) -> str | None:
    if not path.is_file():
        return f"CardiacFoam build manifest does not exist: {path}"
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        return f"Invalid cardiacFoam build manifest {path}: {exc}"
    if payload.get("backend") != backend:
        return f"Build manifest backend is {payload.get('backend')!r}, requested {backend!r}"
    manifest_openfoam = payload.get("openfoam", {}).get("root")
    if openfoam_root is not None and manifest_openfoam and Path(manifest_openfoam).resolve() != openfoam_root:
        return f"Build manifest OpenFOAM root does not match {openfoam_root}"
    manifest_solids = payload.get("solids4foam", {}).get("root")
    if solids_root is not None and manifest_solids and Path(manifest_solids).resolve() != solids_root:
        return f"Build manifest solids4foam root does not match {solids_root}"

    linked = set(payload.get("linked_libraries", ()))
    for library in required_libraries:
        if not any(library in item for item in linked):
            return f"Build manifest is missing required linked library {library}"
    for library in forbidden_libraries:
        if any(library in item for item in linked):
            return f"Build manifest links forbidden library {library}"

    for artifact in payload.get("artifacts", ()):
        artifact_path = Path(artifact.get("path", ""))
        expected = artifact.get("sha256")
        if not artifact_path.is_file() or not expected:
            return f"Build manifest artifact is unavailable: {artifact_path}"
        digest = hashlib.sha256(artifact_path.read_bytes()).hexdigest()
        if digest != expected:
            return f"Build artifact changed since manifest creation: {artifact_path}"
    return None
