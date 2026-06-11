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
#     common
#
# Description
#     Provides shared definitions and defaults for specification templates.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

from __future__ import annotations

import importlib.util
import re
import shutil
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

from ..core.runtime.mutators import ensure_foam_dict, remove_foam_dict, update_foam_entry


def repo_root_default() -> Path:
    current = Path(__file__).resolve()
    for parent in current.parents:
        if (parent / "tutorials").exists() and (parent / "src").exists():
            return parent
    return current.parents[2]


def tutorials_root_default() -> Path:
    repo_root = repo_root_default()
    tutorials_root = repo_root / "tutorials"
    if tutorials_root.exists():
        return tutorials_root
    return repo_root


def default_setup_dir_name(case_dir_name: str) -> str:
    normalized_case_dir = case_dir_name.strip()
    if not normalized_case_dir:
        raise ValueError("case_dir_name cannot be empty")
    case_leaf = Path(normalized_case_dir).name
    return f"setup{case_leaf[:1].upper()}{case_leaf[1:]}"


def detect_myocardium_solver_name(electro_properties_path: Path) -> str:
    for line in electro_properties_path.read_text().splitlines():
        stripped = line.split("//", 1)[0].strip()
        if not stripped.startswith("myocardiumSolver"):
            continue
        tokens = stripped.rstrip(";").split()
        if len(tokens) < 2:
            break
        return tokens[1]
    raise KeyError(f"Could not determine myocardiumSolver from {electro_properties_path}")


def detect_electro_coeffs_scope(electro_properties_path: Path) -> str:
    return f"{detect_myocardium_solver_name(electro_properties_path)}Coeffs"


def detect_ionic_model_name(electro_properties_path: Path) -> str:
    """Return the ionicModel value from the active <solver>Coeffs block.

    Uses the same line-based scan as :func:`detect_myocardium_solver_name`
    rather than a full OpenFOAM dictionary parser. The scope is identified
    by entering the brace following the relevant ``<solver>Coeffs`` header
    and looking for the first ``ionicModel`` key inside it. Raises
    ``KeyError`` if ionicModel is not declared.
    """
    scope = detect_electro_coeffs_scope(electro_properties_path)
    in_scope = False
    depth = 0
    for line in electro_properties_path.read_text().splitlines():
        stripped = line.split("//", 1)[0].strip()
        if not in_scope:
            if stripped == scope or stripped.startswith(f"{scope} ") or stripped.startswith(f"{scope}{{"):
                in_scope = True
            continue
        if "{" in stripped:
            depth += stripped.count("{")
        if stripped.startswith("ionicModel") and depth == 1:
            tokens = stripped.rstrip(";").split()
            if len(tokens) >= 2:
                return tokens[1]
        if "}" in stripped:
            depth -= stripped.count("}")
            if depth <= 0:
                break
    raise KeyError(
        f"Could not determine ionicModel from {electro_properties_path} "
        f"(scope {scope!r})"
    )


_IONIC_EXPORT_RE = re.compile(r"\bionic\s*\{[^}]*\bexport\s*\(([^)]*)\)", re.DOTALL)


def detect_ionic_export_list(
    electro_properties_path: Path,
) -> tuple[str, ...] | None:
    """Return the names declared in ``outputVariables.ionic.export ( ... )``.

    Returns ``None`` when no export declaration is present. The parser is
    intentionally permissive: it strips ``//`` line comments and matches the
    first ``export ( ... )`` block in the file. The myocardium ionic export
    appears before any purkinje-coupling or sub-domain ``export`` by
    OpenFOAM convention (top-level ``<solver>Coeffs`` precedes
    ``domainCouplings`` and other nested blocks), so first-match is correct
    in practice. If this assumption ever breaks, scope this to the active
    ``<solver>Coeffs`` block using brace-depth tracking.
    """
    text = electro_properties_path.read_text()
    cleaned = "\n".join(line.split("//", 1)[0] for line in text.splitlines())
    match = _IONIC_EXPORT_RE.search(cleaned)
    if match is None:
        return None
    tokens = tuple(t for t in match.group(1).split() if t)
    return tokens if tokens else None


_BLOCK_DECL_RE = re.compile(
    r"^\s*(?P<name>[A-Za-z_][A-Za-z0-9_]*)\s*(?:\{|$)",
)


def electro_properties_has_block(
    electro_properties_path: Path,
    block_name: str,
) -> bool:
    """Return True if ``electro_properties_path`` declares a top-level
    OpenFOAM block with the given ``block_name``.

    Matches both inline (``foo { ... }``) and multi-line
    (``foo\\n{\\n...\\n}``) declarations. Substring matches on keys named
    similarly (e.g. ``ecgDomainsCount``) are rejected — the regex anchors
    to a whole identifier followed by an open brace or line end.
    """
    text = electro_properties_path.read_text()
    cleaned_lines = [line.split("//", 1)[0] for line in text.splitlines()]
    for i, line in enumerate(cleaned_lines):
        match = _BLOCK_DECL_RE.match(line)
        if not match or match.group("name") != block_name:
            continue
        # The brace may be consumed by the regex (inline `foo {`) or appear
        # after the match or on the next non-empty line.
        if "{" in match.group(0):
            return True
        rest = line[match.end():].lstrip()
        if rest.startswith("{"):
            return True
        for following in cleaned_lines[i + 1:]:
            stripped = following.strip()
            if not stripped:
                continue
            if stripped.startswith("{"):
                return True
            break
    return False


def detect_verification_model_type(
    electro_properties_path: Path,
) -> str | None:
    """Return the value of ``verificationModel.type`` inside the active
    ``<solver>Coeffs`` block, or None when no verificationModel is
    declared."""
    scope = detect_electro_coeffs_scope(electro_properties_path)
    in_scope = False
    in_verification = False
    depth = 0
    for line in electro_properties_path.read_text().splitlines():
        stripped = line.split("//", 1)[0].strip()
        if not in_scope:
            if (
                stripped == scope
                or stripped.startswith(f"{scope} ")
                or stripped.startswith(f"{scope}{{")
            ):
                in_scope = True
            continue
        if "{" in stripped:
            depth += stripped.count("{")
        if (
            not in_verification
            and stripped.startswith("verificationModel")
            and depth == 1
        ):
            in_verification = True
        if in_verification and stripped.startswith("type") and depth == 2:
            tokens = stripped.rstrip(";").split()
            if len(tokens) >= 2:
                return tokens[1]
        if "}" in stripped:
            depth -= stripped.count("}")
            if in_verification and depth <= 1:
                in_verification = False
            if depth <= 0:
                break
    return None


_AT_EXPORT_RE = re.compile(
    r"\bactiveTension\s*\{[^}]*\bexport\s*\(([^)]*)\)", re.DOTALL
)


def detect_active_tension_model_name(
    electro_properties_path: Path,
) -> str | None:
    """Return the ``activeTensionModel`` value from inside ``<solver>Coeffs``.

    Enters the ``<solver>Coeffs`` block via brace-depth tracking, then enters
    the ``activeTensionModel`` sub-block and reads the ``activeTensionModel``
    key inside it. Returns ``None`` when no ``activeTensionModel`` block is
    declared (coupling disabled).
    """
    scope = detect_electro_coeffs_scope(electro_properties_path)
    in_scope = False
    in_at_block = False
    depth = 0
    for line in electro_properties_path.read_text().splitlines():
        stripped = line.split("//", 1)[0].strip()
        if not in_scope:
            if (
                stripped == scope
                or stripped.startswith(f"{scope} ")
                or stripped.startswith(f"{scope}{{")
            ):
                in_scope = True
            continue
        if "{" in stripped:
            depth += stripped.count("{")
        if not in_at_block and stripped.startswith("activeTensionModel") and depth == 1:
            in_at_block = True
        if in_at_block and stripped.startswith("activeTensionModel") and depth == 2:
            tokens = stripped.rstrip(";").split()
            if len(tokens) >= 2:
                return tokens[1]
        if "}" in stripped:
            depth -= stripped.count("}")
            if in_at_block and depth <= 1:
                in_at_block = False
            if depth <= 0:
                break
    return None


def detect_active_tension_export_list(
    electro_properties_path: Path,
) -> tuple[str, ...] | None:
    """Return the names declared in ``outputVariables.activeTension.export ( ... )``.

    Returns ``None`` when no active-tension export declaration is present.
    Mirrors :func:`detect_ionic_export_list` but scoped to the
    ``activeTension`` sub-block.
    """
    text = electro_properties_path.read_text()
    cleaned = "\n".join(line.split("//", 1)[0] for line in text.splitlines())
    match = _AT_EXPORT_RE.search(cleaned)
    if match is None:
        return None
    tokens = tuple(t for t in match.group(1).split() if t)
    return tokens if tokens else None


def _resolve_scope_tokens(
    path: str,
    *,
    electro_properties_path: Path | None = None,
) -> tuple[str, ...]:
    resolved_parts: list[str] = []
    for token in path.split("."):
        if token == "$ELECTRO_MODEL_COEFFS":
            if electro_properties_path is None:
                raise ValueError(
                    "Scope token '$ELECTRO_MODEL_COEFFS' requires electro_properties_path"
                )
            resolved_parts.append(detect_electro_coeffs_scope(electro_properties_path))
            continue
        resolved_parts.append(token)
    return tuple(part for part in resolved_parts if part)


def normalize_entry_overrides(
    overrides: Mapping[str, Any] | Sequence[Mapping[str, Any]] | None,
    *,
    electro_properties_path: Path | None = None,
) -> list[dict[str, Any]]:
    if overrides is None:
        return []

    def normalize_from_key_value(key_path: str, value: Any) -> dict[str, Any]:
        parts = _resolve_scope_tokens(
            str(key_path),
            electro_properties_path=electro_properties_path,
        )
        if not parts:
            raise ValueError("Override path cannot be empty")
        if len(parts) == 1:
            return {"key": parts[0], "value": value, "scope": None}
        return {"key": parts[-1], "value": value, "scope": parts[:-1]}

    if isinstance(overrides, Mapping):
        return [normalize_from_key_value(key_path, value) for key_path, value in overrides.items()]

    normalized: list[dict[str, Any]] = []
    for item in overrides:
        if not isinstance(item, Mapping):
            raise TypeError("Entry overrides must be a mapping or sequence of mappings")
        if "key" not in item or "value" not in item:
            raise KeyError("Override items must define 'key' and 'value'")

        key = str(item["key"])
        value = item["value"]
        if "scope" in item:
            raw_scope = item["scope"]
            if raw_scope is None:
                scope = None
            elif isinstance(raw_scope, str):
                scope = _resolve_scope_tokens(
                    raw_scope,
                    electro_properties_path=electro_properties_path,
                )
            else:
                scope = tuple(
                    part
                    for token in raw_scope
                    for part in _resolve_scope_tokens(
                        str(token),
                        electro_properties_path=electro_properties_path,
                    )
                )
        else:
            normalized_item = normalize_from_key_value(key, value)
            normalized.append(normalized_item)
            continue

        normalized.append(
            {
                "key": key,
                "value": value,
                "scope": scope,
            }
        )

    return normalized


def apply_entry_overrides(
    file_path: Path,
    overrides: Mapping[str, Any] | Sequence[Mapping[str, Any]] | None,
    *,
    electro_properties_path: Path | None = None,
) -> None:
    for item in normalize_entry_overrides(
        overrides,
        electro_properties_path=electro_properties_path,
    ):
        update_foam_entry(
            file_path,
            item["key"],
            item["value"],
            scope=item["scope"],
        )


def apply_electro_property_overrides(
    electro_properties_path: Path,
    overrides: Mapping[str, Any] | Sequence[Mapping[str, Any]] | None,
) -> None:
    apply_entry_overrides(
        electro_properties_path,
        overrides,
        electro_properties_path=electro_properties_path,
    )


def apply_physics_property_overrides(
    physics_properties_path: Path,
    overrides: Mapping[str, Any] | Sequence[Mapping[str, Any]] | None,
) -> None:
    apply_entry_overrides(physics_properties_path, overrides)


def remove_electro_property_dict(
    electro_properties_path: Path,
    dict_name: str,
    *,
    scope: str | Sequence[str] | None = None,
    missing_ok: bool = False,
) -> None:
    resolved_scope = None
    if scope is not None:
        raw_scope = (scope,) if isinstance(scope, str) else tuple(scope)
        resolved_scope = tuple(
            part
            for token in raw_scope
            for part in _resolve_scope_tokens(
                str(token),
                electro_properties_path=electro_properties_path,
            )
        )

    remove_foam_dict(
        electro_properties_path,
        dict_name,
        scope=resolved_scope,
        missing_ok=missing_ok,
    )


def ensure_electro_property_dict(
    electro_properties_path: Path,
    dict_name: str,
    block_text: str,
    *,
    scope: str | Sequence[str] | None = None,
) -> bool:
    resolved_scope = None
    if scope is not None:
        raw_scope = (scope,) if isinstance(scope, str) else tuple(scope)
        resolved_scope = tuple(
            part
            for token in raw_scope
            for part in _resolve_scope_tokens(
                str(token),
                electro_properties_path=electro_properties_path,
            )
        )

    return ensure_foam_dict(
        electro_properties_path,
        dict_name,
        block_text,
        scope=resolved_scope,
    )


def resolve_spec_paths(
    *,
    tutorials_root: Path | None,
    case_dir_name: str,
    setup_dir_name: str | None = None,
    output_dir_name: str | Path | None = None,
    default_output_dir_name: str | Path | None = None,
) -> tuple[Path, Path, Path]:
    resolved_tutorials_root = (
        Path(tutorials_root) if tutorials_root is not None else tutorials_root_default()
    )
    resolved_case_dir = case_dir_name.strip()
    if not resolved_case_dir:
        raise ValueError("case_dir_name cannot be empty")
    resolved_setup_dir = setup_dir_name or default_setup_dir_name(resolved_case_dir)
    resolved_output_dir = output_dir_name or default_output_dir_name
    if resolved_output_dir is None:
        raise ValueError("output_dir_name or default_output_dir_name must be provided")

    case_root = resolved_tutorials_root / resolved_case_dir
    setup_root = case_root / resolved_setup_dir
    output_dir = case_root / Path(resolved_output_dir)
    return case_root, setup_root, output_dir


def resolve_run_script_path(
    *,
    tutorials_root: Path | None,
    run_script_relpath: Path,
) -> Path:
    if run_script_relpath.is_absolute():
        return run_script_relpath

    candidate_roots: list[Path] = []
    if tutorials_root is not None:
        candidate_roots.append(Path(tutorials_root))
    candidate_roots.append(repo_root_default())
    candidate_roots.append(tutorials_root_default())

    checked_paths: list[Path] = []
    seen: set[Path] = set()
    for root in candidate_roots:
        resolved_root = root.resolve()
        if resolved_root in seen:
            continue
        seen.add(resolved_root)
        candidate = resolved_root / run_script_relpath
        checked_paths.append(candidate)
        if candidate.exists():
            return candidate

    checked_str = ", ".join(str(path) for path in checked_paths)
    raise FileNotFoundError(
        f"Run script not found for '{run_script_relpath}'. Checked: {checked_str}. "
        "Use an absolute path or a path relative to repository root."
    )


def load_python_module(module_path: Path, *, module_name: str):
    spec = importlib.util.spec_from_file_location(module_name, module_path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Unable to load module from {module_path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def collect_outputs_by_pattern(case_root: Path, output_dir: Path, *, pattern: str) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    matching_files = sorted(case_root.glob(pattern))
    for source in matching_files:
        destination = output_dir / source.name
        if destination.exists():
            destination.unlink()
        shutil.move(str(source), str(destination))
        print(f"Moved output: {source.name} -> {destination}")


def set_delta_t(control_dict_path: Path, delta_t_seconds: float) -> None:
    update_foam_entry(control_dict_path, "deltaT", delta_t_seconds)


def set_end_time(control_dict_path: Path, t_s: float) -> None:
    update_foam_entry(control_dict_path, "endTime", t_s)
