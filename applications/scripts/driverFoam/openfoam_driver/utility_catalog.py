"""
Utility Manifest Catalog

A static, eagerly-loaded catalog of cardiacFoam utilities discovered from
``utility.manifest.toml`` sidecar files placed next to each utility source
directory under ``applications/utilities/<name>/``.

Schema
------
Each ``utility.manifest.toml`` must contain the following required fields:

    name        (str)   – utility name; must match the parent directory name.
    description (str)   – one-line purpose summary.
    category    (str)   – one of ALLOWED_CATEGORIES.

The following fields are optional:

    purpose       (str)                – 2–4 sentence extended description.
    inputs        (list[str])          – case-relative paths this utility reads.
    outputs       (list[str])          – case-relative paths this utility writes.
    requires_mesh (bool, default True) – False for 0-D / file-only tools.
    example       (str)                – representative command-line invocation.

    [[flags]]                          – zero or more [[flags]] tables, each
        name        (str)              – flag name, e.g. "-noScale".
        description (str)              – what the flag does.
        takes_value (bool)             – whether the flag accepts an argument.

Public API
----------
    ALLOWED_CATEGORIES   – frozenset of valid category strings.
    UtilityFlag          – frozen dataclass for a single CLI flag.
    UtilityManifest      – frozen dataclass for a parsed manifest.
    load_utility_manifests(utilities_root) -> dict[str, UtilityManifest]
    UTILITY_CATALOG      – module-level dict populated at import time.
"""

from __future__ import annotations

import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Final

if sys.version_info >= (3, 11):
    import tomllib
else:
    try:
        import tomllib  # type: ignore[no-reattr]
    except ImportError:
        try:
            import tomli as tomllib  # type: ignore[no-reattr]
        except ImportError as exc:
            raise ImportError(
                "Python < 3.11 requires the 'tomli' package: pip install tomli"
            ) from exc

MANIFEST_FILENAME: Final[str] = "utility.manifest.toml"

ALLOWED_CATEGORIES: Final[frozenset[str]] = frozenset(
    {
        "mesh",
        "field-setup",
        "post-processing",
        "io-conversion",
        "verification",
        "parametric-sweep",
    }
)

_REQUIRED_FIELDS: Final[frozenset[str]] = frozenset({"name", "description", "category"})

_KNOWN_FIELDS: Final[frozenset[str]] = frozenset(
    {
        "name",
        "description",
        "purpose",
        "inputs",
        "outputs",
        "requires_mesh",
        "flags",
        "example",
        "category",
    }
)

_KNOWN_FLAG_FIELDS: Final[frozenset[str]] = frozenset(
    {"name", "description", "takes_value"}
)


@dataclass(frozen=True)
class UtilityFlag:
    """Metadata for a single CLI flag exposed by a utility."""

    name: str
    """Flag name, e.g. '-noScale'."""

    description: str
    """What the flag does."""

    takes_value: bool = False
    """Whether the flag accepts a follow-on argument."""


@dataclass(frozen=True)
class UtilityManifest:
    """Parsed metadata for a single cardiacFoam utility."""

    name: str
    """Utility name; matches the parent directory name."""

    description: str
    """One-line purpose summary."""

    purpose: str
    """Extended description (2–4 sentences). May be empty."""

    inputs: tuple[str, ...]
    """Case-relative paths the utility reads. May be empty."""

    outputs: tuple[str, ...]
    """Case-relative paths the utility writes. May be empty."""

    requires_mesh: bool
    """True when the utility needs a meshed OpenFOAM case."""

    flags: tuple[UtilityFlag, ...]
    """CLI flags the utility registers beyond the OpenFOAM defaults."""

    example: str
    """Representative command-line invocation. May be empty."""

    category: str
    """Functional category; one of ALLOWED_CATEGORIES."""

    source_path: Path
    """Absolute path to the manifest file; populated by the loader."""


def _parse_flag(raw: object, manifest_path: Path) -> UtilityFlag:
    if not isinstance(raw, dict):
        raise ValueError(
            f"{manifest_path}: each [[flags]] entry must be a TOML table, "
            f"got {type(raw).__name__!r}"
        )
    unknown = set(raw) - _KNOWN_FLAG_FIELDS
    if unknown:
        raise ValueError(
            f"{manifest_path}: unknown field(s) in [[flags]]: {sorted(unknown)}"
        )
    if "name" not in raw:
        raise ValueError(f"{manifest_path}: [[flags]] entry is missing 'name'")
    if "description" not in raw:
        raise ValueError(
            f"{manifest_path}: [[flags]] entry {raw['name']!r} is missing 'description'"
        )
    return UtilityFlag(
        name=raw["name"],
        description=raw["description"],
        takes_value=bool(raw.get("takes_value", False)),
    )


def _parse_manifest(toml_path: Path) -> UtilityManifest:
    """Parse and validate a single ``utility.manifest.toml`` file."""
    with toml_path.open("rb") as fh:
        raw = tomllib.load(fh)

    unknown = set(raw) - _KNOWN_FIELDS
    if unknown:
        raise ValueError(
            f"{toml_path}: unknown TOML field(s): {sorted(unknown)}"
        )

    missing = _REQUIRED_FIELDS - set(raw)
    if missing:
        raise ValueError(
            f"{toml_path}: missing required field(s): {sorted(missing)}"
        )

    name: str = raw["name"]
    dir_name: str = toml_path.parent.name

    if name != dir_name:
        raise ValueError(
            f"{toml_path}: 'name' field {name!r} does not match "
            f"directory name {dir_name!r}"
        )

    category: str = raw["category"]
    if category not in ALLOWED_CATEGORIES:
        raise ValueError(
            f"{toml_path}: 'category' {category!r} is not in "
            f"ALLOWED_CATEGORIES: {sorted(ALLOWED_CATEGORIES)}"
        )

    raw_flags = raw.get("flags", [])
    if not isinstance(raw_flags, list):
        raise ValueError(f"{toml_path}: 'flags' must be a TOML array of tables")

    flags = tuple(_parse_flag(f, toml_path) for f in raw_flags)

    return UtilityManifest(
        name=name,
        description=raw["description"],
        purpose=raw.get("purpose", ""),
        inputs=tuple(raw.get("inputs", [])),
        outputs=tuple(raw.get("outputs", [])),
        requires_mesh=bool(raw.get("requires_mesh", True)),
        flags=flags,
        example=raw.get("example", ""),
        category=category,
        source_path=toml_path,
    )


def load_utility_manifests(utilities_root: Path) -> dict[str, UtilityManifest]:
    """Walk ``utilities_root/*/utility.manifest.toml`` and return a name → manifest map.

    Loader guarantees:
    - Non-directory entries under ``utilities_root`` are silently skipped.
    - A directory with no manifest is not loaded (but the contract test will
      flag it as a missing manifest).
    - Raises ``ValueError`` for: name/directory mismatch, unknown TOML fields,
      missing required fields, unknown category, or duplicate utility names.
    """
    catalog: dict[str, UtilityManifest] = {}

    if not utilities_root.is_dir():
        raise ValueError(f"utilities_root is not a directory: {utilities_root}")

    for child in sorted(utilities_root.iterdir()):
        if not child.is_dir():
            continue
        manifest_path = child / MANIFEST_FILENAME
        if not manifest_path.exists():
            continue

        manifest = _parse_manifest(manifest_path)

        if manifest.name in catalog:
            raise ValueError(
                f"Duplicate utility name {manifest.name!r}: "
                f"found in both {catalog[manifest.name].source_path} "
                f"and {manifest_path}"
            )

        catalog[manifest.name] = manifest

    return catalog


_UTILITIES_ROOT: Final[Path] = (
    Path(__file__).resolve().parents[3] / "utilities"
)

UTILITY_CATALOG: Final[dict[str, UtilityManifest]] = load_utility_manifests(
    _UTILITIES_ROOT
)
