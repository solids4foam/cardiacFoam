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
#     registry
#
# Description
#     Manages tutorial discovery, registration, and specification loading.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

from __future__ import annotations

import json
import os
from dataclasses import replace
from pathlib import Path
from typing import Callable

from .models import TutorialSpec
from ...specs.common import tutorials_root_default
from ...specs.tutorials.manufactured_fda_bidomain import (
    make_spec as make_manufactured_fda_bidomain_spec,
)
from ...specs.tutorials.manufactured_fda_bath_bidomain import (
    make_spec as make_manufactured_fda_bath_bidomain_spec,
)
from ...specs.tutorials.manufactured_eikonal_ecg import (
    make_spec as make_manufactured_eikonal_ecg_spec,
)
from ...specs.tutorials.manufactured_electromechanics_bc import (
    make_spec as make_manufactured_electromechanics_bc_spec,
)
from ...specs.tutorials.generic_case import make_spec as make_generic_case_spec
from ...specs.tutorials.manufactured_fda import make_spec as make_manufactured_fda_spec
from ...specs.tutorials.niederer_2012 import make_spec as make_niederer_2012_spec
from ...specs.tutorials.restitution_curves import make_spec as make_restitution_curves_spec
from ...specs.tutorials.single_cell import make_spec as make_single_cell_spec

SPEC_FACTORIES = {
    "singleCell": make_single_cell_spec,
    "singlecell": make_single_cell_spec,
    "niederer2012": make_niederer_2012_spec,
    "niedereretal2012": make_niederer_2012_spec,
    "manufacturedFDA": make_manufactured_fda_spec,
    "manufacturedfda": make_manufactured_fda_spec,
    "manufacturedFDABidomain": make_manufactured_fda_bidomain_spec,
    "manufacturedfdabidomain": make_manufactured_fda_bidomain_spec,
    "manufacturedFDABathBidomain": make_manufactured_fda_bath_bidomain_spec,
    "manufacturedfdabathbidomain": make_manufactured_fda_bath_bidomain_spec,
    "manufacturedEikonalECG": make_manufactured_eikonal_ecg_spec,
    "manufacturedeikonalecg": make_manufactured_eikonal_ecg_spec,
    "manufacturedElectromechanicsBC": make_manufactured_electromechanics_bc_spec,
    "manufacturedelectromechanicsbc": make_manufactured_electromechanics_bc_spec,
    "restitutionCurves": make_restitution_curves_spec,
    "restitutioncurves": make_restitution_curves_spec,
}

SpecFactory = Callable[..., TutorialSpec]

REGISTERED_TUTORIALS = (
    "singleCell",
    "niederer2012",
    "manufacturedFDA",
    "manufacturedFDABidomain",
    "manufacturedFDABathBidomain",
    "manufacturedEikonalECG",
    "manufacturedElectromechanicsBC",
    "restitutionCurves",
)

ENTRY_KIND_VALUES = (
    "registered_tutorial",
    "workflow_template",
    "workflow_case",
    "case_folder",
)

_ENTRY_HINTS: dict[str, dict[str, object]] = {}

_CORE_REQUIRED_FILES = (
    "constant/electroProperties",
    "constant/physicsProperties",
)

_SOLVER_REQUIRED_FILES = (
    "system/controlDict",
    "system/fvSchemes",
    "system/fvSolution",
)


def _is_case_directory(path: Path) -> bool:
    if not path.is_dir() or path.name.startswith(".") or path.name == "__pycache__":
        return False
    return (path / "constant" / "electroProperties").exists()


def _read_json_if_exists(path: Path) -> dict[str, object] | None:
    if not path.exists():
        return None
    return json.loads(path.read_text())


def _case_is_runnable(
    case_root: Path,
    authoring_contract: dict[str, object] | None = None,
) -> bool:
    if authoring_contract is not None:
        status = authoring_contract.get("status")
        if isinstance(status, dict):
            runnable_without_substitution = status.get("runnable_without_substitution")
            if runnable_without_substitution is False:
                return False

    required_paths = (*_CORE_REQUIRED_FILES, *_SOLVER_REQUIRED_FILES)
    return all((case_root / relpath).exists() for relpath in required_paths)


def _iter_case_directories_recursive(tutorials_root: Path) -> list[Path]:
    if not tutorials_root.exists():
        return []

    discovered: list[Path] = []
    for current_root, dirnames, _filenames in os.walk(tutorials_root):
        path = Path(current_root)
        dirnames[:] = [
            dirname
            for dirname in dirnames
            if not dirname.startswith(".")
            and dirname != "__pycache__"
            and not dirname.startswith("processor")
            and dirname not in {"postProcessing", "logs"}
        ]
        if _is_case_directory(path):
            discovered.append(path)
            dirnames[:] = []
    return discovered


def _registered_tutorial_entry(
    tutorial: str,
    tutorials_root: Path,
) -> dict[str, object]:
    factory = _normalized_registry()[tutorial.casefold()]
    spec = factory(tutorials_root=tutorials_root)
    case_root = Path(spec.case_root)
    try:
        entry_path = str(case_root.relative_to(tutorials_root))
    except ValueError:
        entry_path = case_root.name
    return {
        "entry_name": tutorial,
        "entry_kind": "registered_tutorial",
        "entry_path": entry_path,
        "is_runnable": True,
        "source_type": "spec_factory",
        "workflow_family": None,
    }


def _classify_case_entry(case_root: Path, tutorials_root: Path) -> dict[str, object]:
    relative_path = str(case_root.relative_to(tutorials_root))
    normalized_relative_path = relative_path.casefold()
    authoring_contract = _read_json_if_exists(case_root / "workflow_contract.json")

    hint = _ENTRY_HINTS.get(normalized_relative_path)
    if hint is None:
        hint = _ENTRY_HINTS.get(case_root.name.casefold())

    if hint is not None:
        entry_kind = str(hint["entry_kind"])
        source_type = str(hint["source_type"])
        workflow_family = hint.get("workflow_family")
        is_runnable = bool(hint["is_runnable"])
    elif authoring_contract is not None:
        status = authoring_contract.get("status", {})
        template_kind = status.get("template_kind") if isinstance(status, dict) else None
        runnable_without_substitution = (
            status.get("runnable_without_substitution") if isinstance(status, dict) else None
        )
        if template_kind == "symbolic_authoring_template" or runnable_without_substitution is False:
            entry_kind = "workflow_template"
            source_type = "workflow_contract"
            workflow_family = authoring_contract.get("tutorial_family")
            is_runnable = False
        else:
            entry_kind = "case_folder"
            source_type = "filesystem_case"
            workflow_family = authoring_contract.get("tutorial_family")
            is_runnable = _case_is_runnable(case_root, authoring_contract)
    else:
        entry_kind = "case_folder"
        source_type = "filesystem_case"
        workflow_family = None
        is_runnable = _case_is_runnable(case_root)

    # Extract workflow_dag from the on-disk contract if a steps array is present.
    workflow_dag: dict[str, object] | None = None
    if authoring_contract is not None:
        raw_steps = authoring_contract.get("steps")
        if isinstance(raw_steps, list) and raw_steps:
            workflow_dag = {"steps": raw_steps}

    return {
        "entry_name": case_root.name,
        "entry_kind": entry_kind,
        "entry_path": relative_path,
        "is_runnable": is_runnable,
        "source_type": source_type,
        "workflow_family": workflow_family,
        "workflow_dag": workflow_dag,
    }


def _entry_catalog_for_root(tutorials_root: Path) -> list[dict[str, object]]:
    entries: list[dict[str, object]] = [
        _registered_tutorial_entry(tutorial, tutorials_root)
        for tutorial in REGISTERED_TUTORIALS
    ]
    known_registered = {tutorial.casefold() for tutorial in REGISTERED_TUTORIALS}
    for case_root in _iter_case_directories_recursive(tutorials_root):
        classified = _classify_case_entry(case_root, tutorials_root)
        if classified["entry_name"].casefold() in known_registered:
            continue
        entries.append(classified)

    return sorted(
        entries,
        key=lambda entry: (
            ENTRY_KIND_VALUES.index(str(entry["entry_kind"])),
            str(entry["entry_name"]).casefold(),
            str(entry["entry_path"]).casefold(),
        ),
    )


def list_tutorials() -> list[str]:
    return list(REGISTERED_TUTORIALS)


def list_case_directories(tutorials_root: Path | None = None) -> list[str]:
    resolved_root = Path(tutorials_root) if tutorials_root is not None else tutorials_root_default()
    if not resolved_root.exists():
        return []
    return sorted(child.name for child in resolved_root.iterdir() if _is_case_directory(child))


def list_available_tutorials(tutorials_root: Path | None = None) -> list[str]:
    available = list_tutorials()
    known = {name.casefold() for name in available}
    for case_dir in list_case_directories(tutorials_root):
        if case_dir.casefold() in known:
            continue
        available.append(case_dir)
        known.add(case_dir.casefold())
    return available


def list_entries(tutorials_root: Path | None = None) -> list[dict[str, object]]:
    resolved_root = Path(tutorials_root) if tutorials_root is not None else tutorials_root_default()
    return _entry_catalog_for_root(resolved_root)


def _normalized_registry() -> dict[str, object]:
    return {name.casefold(): factory for name, factory in SPEC_FACTORIES.items()}


def _find_case_dir_match(name: str, tutorials_root: Path) -> str | None:
    requested = name.strip()
    if not requested:
        return None

    direct = tutorials_root / requested
    if _is_case_directory(direct):
        return requested

    normalized = requested.casefold()
    for child in tutorials_root.iterdir():
        if not _is_case_directory(child):
            continue
        if child.name.casefold() == normalized:
            return child.name
    return None


def load_tutorial_spec(name: str, overrides: dict | None = None) -> TutorialSpec:
    resolution = resolve_tutorial(name, overrides=overrides)
    spec = resolution["factory"](**resolution["factory_overrides"])
    return _with_entry_metadata(spec, resolution)


def load_entry_spec(
    name: str,
    *,
    entry_kind: str | None = None,
    overrides: dict | None = None,
) -> TutorialSpec:
    resolution = resolve_entry(name, entry_kind=entry_kind, overrides=overrides)
    spec = resolution["factory"](**resolution["factory_overrides"])
    return _with_entry_metadata(spec, resolution)


def _with_entry_metadata(
    spec: TutorialSpec,
    resolution: dict[str, object],
) -> TutorialSpec:
    metadata = dict(spec.metadata)
    metadata.update(
        {
            "entry_name": resolution["entry_name"],
            "entry_kind": resolution["entry_kind"],
            "entry_path": resolution["entry_path"],
            "source_type": resolution["source_type"],
            "workflow_family": resolution["workflow_family"],
            "resolution": resolution["resolution"],
        }
    )
    # For filesystem cases, the on-disk workflow_contract.json is authoritative.
    # When the registry found a 'steps' array there, set it unconditionally so
    # it overrides any generic-spec fallback. When there is no on-disk DAG
    # (resolution key absent or explicitly None from a contract without steps),
    # leave the spec's own metadata untouched — spec-factory DAGs are preserved.
    if "workflow_dag" in resolution:
        on_disk_dag = resolution["workflow_dag"]
        if on_disk_dag is not None:
            # On-disk steps win; overwrite spec-factory default.
            metadata["workflow_dag"] = on_disk_dag
        else:
            # Contract present but no steps array (or absent contract) — clear
            # any generic-spec placeholder so callers see None.
            metadata["workflow_dag"] = None
    return replace(spec, metadata=metadata)


def _match_entry(
    name: str,
    entry_kind: str | None,
    tutorials_root: Path,
) -> dict[str, object] | None:
    normalized_name = name.strip().casefold()
    matches = [
        entry
        for entry in list_entries(tutorials_root)
        if (
            normalized_name in {
                str(entry["entry_name"]).casefold(),
                str(entry["entry_path"]).casefold(),
            }
            and (entry_kind is None or str(entry["entry_kind"]) == entry_kind)
        )
    ]

    if not matches:
        return None
    if len(matches) == 1:
        return matches[0]

    exact_path_matches = [
        entry
        for entry in matches
        if str(entry["entry_path"]).casefold() == normalized_name
    ]
    if len(exact_path_matches) == 1:
        return exact_path_matches[0]

    options = ", ".join(sorted(str(entry["entry_path"]) for entry in matches))
    raise KeyError(
        f"Entry '{name}' is ambiguous. Use a more specific entry path. Matches: {options}"
    )


def resolve_entry(
    name: str,
    *,
    entry_kind: str | None = None,
    overrides: dict | None = None,
) -> dict[str, object]:
    key = name.strip()
    normalized_key = key.casefold()
    normalized_registry = _normalized_registry()
    incoming_overrides = dict(overrides or {})

    if entry_kind is not None and entry_kind not in ENTRY_KIND_VALUES:
        valid = ", ".join(ENTRY_KIND_VALUES)
        raise KeyError(f"Unknown entry_kind '{entry_kind}'. Valid values: {valid}")

    tutorials_root = Path(incoming_overrides.get("tutorials_root", tutorials_root_default()))

    if entry_kind in {None, "registered_tutorial"} and normalized_key in normalized_registry:
        return {
            "resolution": "registered",
            "requested_name": key,
            "requested_entry_kind": entry_kind,
            "resolved_name": key,
            "factory": normalized_registry[normalized_key],
            "factory_overrides": incoming_overrides,
            "entry_name": key,
            "entry_kind": "registered_tutorial",
            "entry_path": _registered_tutorial_entry(key, tutorials_root)["entry_path"],
            "is_runnable": True,
            "source_type": "spec_factory",
            "workflow_family": None,
        }

    matched_entry = _match_entry(key, entry_kind, tutorials_root)
    if matched_entry is not None:
        generic_overrides = dict(incoming_overrides)
        generic_overrides.setdefault("case_dir_name", str(matched_entry["entry_path"]))
        return {
            "resolution": "case_folder",
            "requested_name": key,
            "requested_entry_kind": entry_kind,
            "resolved_name": str(matched_entry["entry_name"]),
            "factory": make_generic_case_spec,
            "factory_overrides": generic_overrides,
            **matched_entry,
        }

    if entry_kind in {None, "case_folder"} and normalized_key in {"genericcase", "randomcase"}:
        if "case_dir_name" not in incoming_overrides:
            raise KeyError(
                f"Entry '{name}' requires a 'case_dir_name' override to select a case folder."
            )
        matched_case_dir = str(incoming_overrides["case_dir_name"])
        return {
            "resolution": "generic_alias",
            "requested_name": key,
            "requested_entry_kind": entry_kind,
            "resolved_name": matched_case_dir,
            "factory": make_generic_case_spec,
            "factory_overrides": incoming_overrides,
            "entry_name": Path(matched_case_dir).name,
            "entry_kind": "case_folder",
            "entry_path": matched_case_dir,
            "is_runnable": False,
            "source_type": "generic_alias",
            "workflow_family": None,
        }

    valid = ", ".join(list_tutorials())
    raise KeyError(
        f"Unknown entry '{name}'. Valid registered tutorials: {valid}. "
        "You can also pass any existing tutorial case folder or workflow entry path."
    )


def resolve_tutorial(name: str, overrides: dict | None = None) -> dict[str, object]:
    return resolve_entry(name, overrides=overrides)
