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
#     utility_catalog
#
# Description
#     Manages specification manifests for executable utilities.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

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
                                         Deprecated: prefer ``produces``. A
                                         UserWarning is emitted if ``outputs``
                                         and ``produces`` disagree.
    requires_mesh (bool, default True) – False for 0-D / file-only tools.
    example       (str)                – representative command-line invocation.

    [[flags]]                          – zero or more [[flags]] tables, each
        name          (str)              – flag name, e.g. "-noScale".
        description   (str)              – what the flag does.
        takes_value   (bool)             – whether the flag accepts an argument.
        argument_kind (str, optional)    – one of ALLOWED_ARGUMENT_KINDS:
                                           scalar|label|path|word|word_list|switch
        required      (bool, default False) – whether the flag is mandatory.
        default       (str, optional)    – default value as a string.

    positional_args = [...]            – ordered list of positional argument
        tables, each with:
        name          (str)
        argument_kind (str)  – one of ALLOWED_ARGUMENT_KINDS
        description   (str)

    produces = [...]                   – structured output declarations that
        supersede ``outputs`` over time. Each entry:
        artifact_id   (str)   – stable identifier.
        path_pattern  (str)   – case-relative path; may contain {case_id} or
                                 {time} placeholders (validated at load time
                                 via models._validate_path_pattern — Gap B).
        format        (str)   – one of the ArtifactFormat values from
                                 core.runtime.models.
        description   (str, optional)
        produced_by   (str, optional)  – utility/solver name.
        variables     (list[str], optional)
        optional      (bool, default False)
        time_indexed  (bool, default False)

Public API
----------
    ALLOWED_CATEGORIES     – frozenset of valid category strings.
    ALLOWED_ARGUMENT_KINDS – frozenset of valid argument_kind strings.
    ALLOWED_ARTIFACT_FORMATS – frozenset of valid produces[*].format strings.
    PositionalArg          – frozen dataclass for a positional argument.
    UtilityFlag            – frozen dataclass for a single CLI flag.
    ProducesEntry          – frozen dataclass for a produces entry.
    UtilityManifest        – frozen dataclass for a parsed manifest.
    load_utility_manifests(utilities_root) -> dict[str, UtilityManifest]
    UTILITY_CATALOG        – module-level dict populated at import time.
"""

from __future__ import annotations

import sys
import warnings
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

from .runtime.models import _validate_path_pattern

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

ALLOWED_ARGUMENT_KINDS: Final[frozenset[str]] = frozenset(
    {"scalar", "label", "path", "word", "word_list", "switch"}
)

ALLOWED_ARTIFACT_FORMATS: Final[frozenset[str]] = frozenset(
    {
        "csv_probe",
        "csv_sweep",
        "vtk_sequence",
        "openfoam_time_dirs",
        "openfoam_log",
        "json_summary",
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
        "positional_args",
        "produces",
    }
)

_KNOWN_FLAG_FIELDS: Final[frozenset[str]] = frozenset(
    {"name", "description", "takes_value", "argument_kind", "required", "default"}
)

_KNOWN_POSITIONAL_ARG_FIELDS: Final[frozenset[str]] = frozenset(
    {"name", "argument_kind", "description"}
)

_KNOWN_PRODUCES_FIELDS: Final[frozenset[str]] = frozenset(
    {
        "artifact_id",
        "path_pattern",
        "format",
        "description",
        "produced_by",
        "variables",
        "optional",
        "time_indexed",
    }
)


@dataclass(frozen=True)
class PositionalArg:
    """Metadata for a single positional argument of a utility."""

    name: str
    """Argument name, e.g. 'vtk_file'."""

    argument_kind: str
    """One of ALLOWED_ARGUMENT_KINDS."""

    description: str
    """What the argument represents."""


@dataclass(frozen=True)
class UtilityFlag:
    """Metadata for a single CLI flag exposed by a utility."""

    name: str
    """Flag name, e.g. '-noScale'."""

    description: str
    """What the flag does."""

    takes_value: bool = False
    """Whether the flag accepts a follow-on argument."""

    argument_kind: str = ""
    """One of ALLOWED_ARGUMENT_KINDS. Empty string means not specified."""

    required: bool = False
    """Whether the flag is mandatory."""

    default: str | None = None
    """Default value as a string, or None if not specified."""


@dataclass(frozen=True)
class ProducesEntry:
    """Structured output declaration for a utility."""

    artifact_id: str
    """Stable identifier within the manifest."""

    path_pattern: str
    """Case-relative path; may contain {case_id} / {time} placeholders."""

    format: str
    """One of ALLOWED_ARTIFACT_FORMATS."""

    description: str = ""
    """Human-readable purpose. May be empty."""

    produced_by: str = ""
    """Utility/solver name that writes this artifact. Empty means implicit."""

    variables: tuple[str, ...] = ()
    """Per-variable structure inside the file. () means no per-variable structure."""

    optional: bool = False
    """True when the artifact appears only under specific configurations."""

    time_indexed: bool = False
    """True for time-directory style outputs."""


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
    """Case-relative paths the utility writes. May be empty.
    Deprecated: prefer ``produces``. Kept for backward compatibility."""

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

    positional_args: tuple[PositionalArg, ...] = ()
    """Ordered list of positional arguments. May be empty."""

    produces: tuple[ProducesEntry, ...] = ()
    """Structured output declarations. May be empty."""


def _parse_positional_arg(raw: object, manifest_path: Path) -> PositionalArg:
    if not isinstance(raw, dict):
        raise ValueError(
            f"{manifest_path}: each positional_args entry must be a TOML inline table, "
            f"got {type(raw).__name__!r}"
        )
    unknown = set(raw) - _KNOWN_POSITIONAL_ARG_FIELDS
    if unknown:
        raise ValueError(
            f"{manifest_path}: unknown field(s) in positional_args entry: {sorted(unknown)}"
        )
    for required_field in ("name", "argument_kind", "description"):
        if required_field not in raw:
            raise ValueError(
                f"{manifest_path}: positional_args entry is missing {required_field!r}"
            )
    kind: str = raw["argument_kind"]
    if kind not in ALLOWED_ARGUMENT_KINDS:
        raise ValueError(
            f"{manifest_path}: positional_args entry {raw['name']!r} has unknown "
            f"argument_kind {kind!r}; allowed: {sorted(ALLOWED_ARGUMENT_KINDS)}"
        )
    return PositionalArg(
        name=raw["name"],
        argument_kind=kind,
        description=raw["description"],
    )


def _parse_produces_entry(raw: object, manifest_path: Path) -> ProducesEntry:
    if not isinstance(raw, dict):
        raise ValueError(
            f"{manifest_path}: each produces entry must be a TOML inline table or table, "
            f"got {type(raw).__name__!r}"
        )
    unknown = set(raw) - _KNOWN_PRODUCES_FIELDS
    if unknown:
        raise ValueError(
            f"{manifest_path}: unknown field(s) in produces entry: {sorted(unknown)}"
        )
    for required_field in ("artifact_id", "path_pattern", "format"):
        if required_field not in raw:
            raise ValueError(
                f"{manifest_path}: produces entry is missing {required_field!r}"
            )
    fmt: str = raw["format"]
    if fmt not in ALLOWED_ARTIFACT_FORMATS:
        raise ValueError(
            f"{manifest_path}: produces entry {raw['artifact_id']!r} has unknown "
            f"format {fmt!r}; allowed: {sorted(ALLOWED_ARTIFACT_FORMATS)}"
        )
    pattern: str = raw["path_pattern"]
    # Gap B: validate placeholders at TOML-load time reusing the same validator
    # used by DataArtifact.__post_init__ so utility-manifest authors cannot
    # introduce a third placeholder convention.
    try:
        _validate_path_pattern(pattern)
    except ValueError as exc:
        raise ValueError(
            f"{manifest_path}: produces entry {raw['artifact_id']!r}: {exc}"
        ) from exc

    raw_variables = raw.get("variables", [])
    if not isinstance(raw_variables, list):
        raise ValueError(
            f"{manifest_path}: produces entry {raw['artifact_id']!r}: "
            f"'variables' must be a list, got {type(raw_variables).__name__!r}"
        )

    return ProducesEntry(
        artifact_id=raw["artifact_id"],
        path_pattern=pattern,
        format=fmt,
        description=raw.get("description", ""),
        produced_by=raw.get("produced_by", ""),
        variables=tuple(raw_variables),
        optional=bool(raw.get("optional", False)),
        time_indexed=bool(raw.get("time_indexed", False)),
    )


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
    argument_kind: str = raw.get("argument_kind", "")
    if argument_kind and argument_kind not in ALLOWED_ARGUMENT_KINDS:
        raise ValueError(
            f"{manifest_path}: [[flags]] entry {raw['name']!r} has unknown "
            f"argument_kind {argument_kind!r}; allowed: {sorted(ALLOWED_ARGUMENT_KINDS)}"
        )
    raw_default = raw.get("default", None)
    if raw_default is not None and not isinstance(raw_default, str):
        raise ValueError(
            f"{manifest_path}: [[flags]] entry {raw['name']!r}: "
            f"'default' must be a string, got {type(raw_default).__name__!r}"
        )
    return UtilityFlag(
        name=raw["name"],
        description=raw["description"],
        takes_value=bool(raw.get("takes_value", False)),
        argument_kind=argument_kind,
        required=bool(raw.get("required", False)),
        default=raw_default,
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

    raw_positional = raw.get("positional_args", [])
    if not isinstance(raw_positional, list):
        raise ValueError(f"{toml_path}: 'positional_args' must be a TOML array")

    positional_args = tuple(_parse_positional_arg(a, toml_path) for a in raw_positional)

    raw_produces = raw.get("produces", [])
    if not isinstance(raw_produces, list):
        raise ValueError(f"{toml_path}: 'produces' must be a TOML array")

    produces = tuple(_parse_produces_entry(e, toml_path) for e in raw_produces)

    outputs: tuple[str, ...] = tuple(raw.get("outputs", []))

    # Emit UserWarning if outputs and produces disagree (outputs has paths not
    # mirrored in any produces entry). This is a deprecation signal — outputs
    # is kept for backward compatibility but produces should be the authority.
    if outputs and produces:
        produces_paths = {e.path_pattern for e in produces}
        disagreeing = [p for p in outputs if p not in produces_paths]
        if disagreeing:
            warnings.warn(
                f"{toml_path}: 'outputs' contains path(s) {disagreeing!r} not "
                f"mirrored in any 'produces' entry. Prefer 'produces' as the "
                f"authoritative output declaration.",
                UserWarning,
                stacklevel=2,
            )

    return UtilityManifest(
        name=name,
        description=raw["description"],
        purpose=raw.get("purpose", ""),
        inputs=tuple(raw.get("inputs", [])),
        outputs=outputs,
        requires_mesh=bool(raw.get("requires_mesh", True)),
        flags=flags,
        example=raw.get("example", ""),
        category=category,
        source_path=toml_path,
        positional_args=positional_args,
        produces=produces,
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
        return catalog

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


#: Repo root holding the ``utility.manifest.toml`` sidecars. Public so that a
#: plugin declaring the same root derives it from here rather than recomputing
#: its own ``Path(__file__).parents[N]`` arithmetic, which can silently drift.
UTILITIES_ROOT: Final[Path] = (
    Path(__file__).resolve().parents[4] / "utilities"
)

UTILITY_CATALOG: Final[dict[str, UtilityManifest]] = load_utility_manifests(
    UTILITIES_ROOT
)
