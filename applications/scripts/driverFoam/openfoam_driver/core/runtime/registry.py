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

import os
from dataclasses import replace
from pathlib import Path
from typing import TYPE_CHECKING, Callable

from .models import TutorialSpec
from .generic_case import make_generic_case_spec
from openfoam_driver.core.specs.common import tutorials_root_default

if TYPE_CHECKING:
    from ..plugin_interface import DriverContext















SpecFactory = Callable[..., TutorialSpec]



ENTRY_KIND_VALUES = (
    "registered_tutorial",
    "case_folder",
)

def _is_case_directory(
    path: Path,
    driver_context: "DriverContext | None" = None,
) -> bool:
    if not path.is_dir() or path.name.startswith(".") or path.name == "__pycache__":
        return False
    from ..compatibility import resolve_public_driver_context
    from ..plugin_capabilities import CaseCompatibilityRequest

    driver_context = resolve_public_driver_context(driver_context)
    return (
        driver_context.capabilities.case_compatibility.has_case_marker(
            CaseCompatibilityRequest(path),
        )
        or (path / "Allrun").is_file()
    )


def _case_is_runnable(
    case_root: Path,
    *,
    driver_context: "DriverContext | None" = None,
) -> bool:
    if (case_root / "Allrun").is_file():
        return True

    from ..compatibility import resolve_public_driver_context
    from ..plugin_capabilities import CaseCompatibilityRequest

    driver_context = resolve_public_driver_context(driver_context)
    return driver_context.capabilities.case_compatibility.is_runnable_without_workflow(
        CaseCompatibilityRequest(case_root),
    )


def _iter_case_directories_recursive(
    tutorials_root: Path,
    driver_context: "DriverContext | None" = None,
) -> list[Path]:
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
        if _is_case_directory(path, driver_context):
            discovered.append(path)
            dirnames[:] = []
    return discovered


def _registered_tutorial_entry(
    tutorial: str,
    tutorials_root: Path,
    driver_context: "DriverContext | None" = None,
) -> dict[str, object]:
    factory = _normalized_registry(driver_context)[tutorial.casefold()]
    try:
        spec = factory(tutorials_root=tutorials_root)
    except Exception:
        # Cataloging is best-effort: describe_entry() derives tutorials_root
        # from the *queried* entry's own case_root parent (see introspection.
        # describe_entry), which does not necessarily hold every other
        # registered tutorial's case directory -- e.g. singleCell nests one
        # level deeper than the manufactured-solution tutorials. A factory
        # that can't build its spec under this particular root (missing
        # case files, wrong nesting, ...) is simply not runnable from here;
        # that must not crash the listing for every other tutorial.
        return {
            "entry_name": tutorial,
            "entry_kind": "registered_tutorial",
            "entry_path": tutorial,
            "is_runnable": False,
            "source_type": "spec_factory",
            "workflow_family": None,
        }
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


def _classify_case_entry(
    case_root: Path,
    tutorials_root: Path,
    driver_context: "DriverContext | None" = None,
) -> dict[str, object]:
    relative_path = str(case_root.relative_to(tutorials_root))
    return {
        "entry_name": case_root.name,
        "entry_kind": "case_folder",
        "entry_path": relative_path,
        "is_runnable": _case_is_runnable(case_root, driver_context=driver_context),
        "source_type": "filesystem_case",
        "workflow_family": None,
    }


def _entry_catalog_for_root(
    tutorials_root: Path,
    driver_context: "DriverContext | None" = None,
) -> list[dict[str, object]]:
    from ..compatibility import resolve_public_driver_context

    driver_context = resolve_public_driver_context(driver_context)
    entries: list[dict[str, object]] = [
        _registered_tutorial_entry(tutorial, tutorials_root, driver_context)
        for tutorial in list_tutorials(driver_context)
    ]
    known_registered = {tutorial.casefold() for tutorial in list_tutorials(driver_context)}
    for case_root in _iter_case_directories_recursive(tutorials_root, driver_context):
        classified = _classify_case_entry(case_root, tutorials_root, driver_context)
        entries.append(classified)

    return sorted(
        entries,
        key=lambda entry: (
            ENTRY_KIND_VALUES.index(str(entry["entry_kind"])),
            str(entry["entry_name"]).casefold(),
            str(entry["entry_path"]).casefold(),
        ),
    )





def list_case_directories(
    tutorials_root: Path | None = None,
    *,
    driver_context: "DriverContext | None" = None,
) -> list[str]:
    resolved_root = Path(tutorials_root) if tutorials_root is not None else tutorials_root_default()
    if not resolved_root.exists():
        return []
    return sorted(
        child.name
        for child in resolved_root.iterdir()
        if _is_case_directory(child, driver_context)
    )


def list_available_tutorials(
    tutorials_root: Path | None = None,
    *,
    driver_context: "DriverContext | None" = None,
) -> list[str]:
    available = list_tutorials(driver_context)
    known = {name.casefold() for name in available}
    for case_dir in list_case_directories(
        tutorials_root, driver_context=driver_context,
    ):
        if case_dir.casefold() in known:
            continue
        available.append(case_dir)
        known.add(case_dir.casefold())
    return available


def list_entries(
    tutorials_root: Path | None = None,
    *,
    driver_context: "DriverContext | None" = None,
) -> list[dict[str, object]]:
    resolved_root = Path(tutorials_root) if tutorials_root is not None else tutorials_root_default()
    return _entry_catalog_for_root(resolved_root, driver_context)







def load_tutorial_spec(
    name: str,
    overrides: dict | None = None,
    *,
    driver_context: "DriverContext | None" = None,
) -> TutorialSpec:
    resolution = resolve_tutorial(name, overrides=overrides, driver_context=driver_context)
    spec = resolution["factory"](**resolution["factory_overrides"])
    return _with_entry_metadata(spec, resolution)


def load_entry_spec(
    name: str,
    *,
    entry_kind: str | None = None,
    overrides: dict | None = None,
    driver_context: "DriverContext | None" = None,
) -> TutorialSpec:
    resolution = resolve_entry(
        name,
        entry_kind=entry_kind,
        overrides=overrides,
        driver_context=driver_context,
    )
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
    # Plain case folders are owned by their on-disk Allrun. If a discovered
    # folder has no Allrun, do not preserve the generic-spec placeholder DAG.
    if (
        resolution["resolution"] == "case_folder"
        and not (Path(spec.case_root) / "Allrun").is_file()
    ):
        metadata["workflow_dag"] = None
    return replace(spec, metadata=metadata)


def _match_entry(
    name: str,
    entry_kind: str | None,
    tutorials_root: Path,
    driver_context: "DriverContext | None" = None,
) -> dict[str, object] | None:
    normalized_name = name.strip().casefold()
    matches = [
        entry
        for entry in list_entries(tutorials_root, driver_context=driver_context)
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
    driver_context: "DriverContext | None" = None,
) -> dict[str, object]:
    from ..compatibility import resolve_public_driver_context
    from ..plugin_capabilities import CaseCompatibilityRequest

    driver_context = resolve_public_driver_context(driver_context)
    key = name.strip()
    normalized_key = key.casefold()
    normalized_registry = _normalized_registry(driver_context)
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
            "entry_path": _registered_tutorial_entry(key, tutorials_root, driver_context)["entry_path"],
            "is_runnable": True,
            "source_type": "spec_factory",
            "workflow_family": None,
        }

    matched_entry = _match_entry(key, entry_kind, tutorials_root, driver_context)
    if matched_entry is not None:
        generic_overrides = dict(incoming_overrides)
        generic_overrides.setdefault("case_dir_name", str(matched_entry["entry_path"]))
        matched_case_root = tutorials_root / str(matched_entry["entry_path"])
        # Keep existing cardiac case-folder semantics while moving truly
        # solver-neutral folders to the core implementation.  The marker is
        # deliberately narrow: an electroProperties file belongs to the
        # cardiac plugin; its absence must not prevent generic OpenFOAM use.
        generic_factory = (
            _get_plugin_tutorials(driver_context).get("make_generic_case_spec")
            if driver_context.capabilities.case_compatibility.has_case_marker(
                CaseCompatibilityRequest(matched_case_root),
            )
            else make_generic_case_spec
        )
        return {
            "resolution": "case_folder",
            "requested_name": key,
            "requested_entry_kind": entry_kind,
            "resolved_name": str(matched_entry["entry_name"]),
            "factory": generic_factory,
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

    valid = ", ".join(list_tutorials(driver_context))
    raise KeyError(
        f"Unknown entry '{name}'. Valid registered tutorials: {valid}. "
        "You can also pass any existing tutorial case folder or workflow entry path."
    )


def resolve_tutorial(
    name: str,
    overrides: dict | None = None,
    *,
    driver_context: "DriverContext | None" = None,
) -> dict[str, object]:
    return resolve_entry(name, overrides=overrides, driver_context=driver_context)


def _get_plugin_tutorials(driver_context: "DriverContext | None" = None):
    from ..compatibility import resolve_public_driver_context

    driver_context = resolve_public_driver_context(driver_context)
    return driver_context.capabilities.tutorials.catalog()

def _normalized_registry(driver_context: "DriverContext | None" = None) -> dict[str, object]:
    spec_factories = _get_plugin_tutorials(driver_context).get("spec_factories", {})
    return {name.casefold(): factory for name, factory in spec_factories.items()}

def list_tutorials(driver_context: "DriverContext | None" = None) -> list[str]:
    registered = _get_plugin_tutorials(driver_context).get("registered_tutorials", ())
    return list(registered)
