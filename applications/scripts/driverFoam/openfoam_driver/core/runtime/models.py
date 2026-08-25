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
#     models
#
# Description
#     Defines runtime data models and configuration specification classes.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable, Final, Literal


ArtifactFormat = Literal[
    "csv_probe",
    "csv_sweep",
    "vtk_sequence",
    "openfoam_time_dirs",
    "openfoam_log",
    "json_summary",
]
"""Closed enum of artifact output formats. Extending this set is a contract
change — every consumer that branches on ``DataArtifact.format`` must be
updated alongside the addition, and the plan (section 2.1) must be amended."""


@dataclass(frozen=True)
class DataArtifact:
    """Declarative description of a raw data output produced by a run or utility.

    Shared vocabulary between the engine and ``utility.manifest.toml`` ``produces``
    entries (which declare what a utility writes). Agents consume both through
    the same shape.

    ``path_pattern`` is case-relative. The only recognised placeholders are
    ``{case_id}`` (substituted with the sweep case identifier) and ``{time}``
    (substituted with an OpenFOAM time directory name). Anything else is a
    literal path component.
    """

    artifact_id: str
    """Stable identifier within the manifest. Predictor merges static + derived
    artifacts by ``artifact_id`` (static wins on collision)."""

    path_pattern: str
    """Case-relative path; may contain ``{case_id}`` / ``{time}`` placeholders."""

    format: ArtifactFormat
    """One of the documented :data:`ArtifactFormat` values."""

    variables: tuple[str, ...] = ()
    """Per-variable structure inside the file. ``()`` means "no per-variable
    structure" (e.g., an opaque log). Never ``None`` — ambiguity-free merging
    in the predictor depends on this."""

    description: str = ""
    """Human-readable purpose. May be empty."""

    produced_by: str = ""
    """Solver or utility name that writes this artifact. Empty means
    engine-implicit (e.g., OpenFOAM time directories produced by cardiacFoam)."""

    optional: bool = False
    """True when the artifact appears only under specific configurations
    (e.g., a probe CSV that requires probes to be enabled in ``controlDict``)."""

    time_indexed: bool = False
    """True for OpenFOAM time-directory style outputs that produce one file
    per write interval. ``path_pattern`` will typically contain ``{time}``."""

    def __post_init__(self) -> None:
        # Catch typos like {caseId} or {run_id} at construction so they never
        # be finalized. Expansion-time validation alone is
        # not enough: an agent may read path_pattern literally without ever
        # calling expand_path_pattern.
        _validate_path_pattern(self.path_pattern)


_PATH_PATTERN_PLACEHOLDER: Final[re.Pattern[str]] = re.compile(r"\{([a-zA-Z_][a-zA-Z0-9_]*)\}")
_KNOWN_PATH_PLACEHOLDERS: Final[frozenset[str]] = frozenset({"case_id", "time"})


def _validate_path_pattern(pattern: str) -> None:
    """Raise ``ValueError`` if ``pattern`` contains any placeholder not in
    :data:`_KNOWN_PATH_PLACEHOLDERS`. Does not require values — this is shape
    validation, not expansion. Called from :class:`DataArtifact.__post_init__`
    so typos surface at construction rather than at expand-time."""
    for match in _PATH_PATTERN_PLACEHOLDER.finditer(pattern):
        name = match.group(1)
        if name not in _KNOWN_PATH_PLACEHOLDERS:
            raise ValueError(
                f"unknown placeholder {{{name}}} in path pattern {pattern!r}; "
                f"recognised placeholders are {sorted(_KNOWN_PATH_PLACEHOLDERS)}"
            )


def expand_path_pattern(
    pattern: str,
    *,
    case_id: str | None = None,
    time: str | None = None,
) -> str:
    """Substitute ``{case_id}`` and ``{time}`` placeholders in a path pattern.

    The set of recognised placeholders is closed (see
    :data:`_KNOWN_PATH_PLACEHOLDERS`). Encountering an unknown ``{foo}`` token
    raises ``ValueError`` so authors cannot silently introduce a third
    placeholder convention.

    Passing an unused keyword (e.g., ``time=`` when the pattern has no
    ``{time}``) is tolerated — callers compose artifacts uniformly and should
    not have to inspect every pattern before invoking the helper.
    """
    _validate_path_pattern(pattern)
    values = {"case_id": case_id, "time": time}

    def _resolve(match: re.Match[str]) -> str:
        name = match.group(1)
        value = values[name]
        if value is None:
            raise ValueError(
                f"path pattern {pattern!r} references {{{name}}} but no "
                f"{name}= was supplied"
            )
        return value

    return _PATH_PATTERN_PLACEHOLDER.sub(_resolve, pattern)


def data_artifact_from_json(data: dict[str, Any]) -> DataArtifact:
    """Reconstruct a :class:`DataArtifact` from its JSON dict form.

    Inverse of ``dataclasses.asdict(artifact)``. ``variables`` is coerced
    back to a tuple of strings; absent optional keys fall back to the
    dataclass defaults. ``DataArtifact.__post_init__`` still runs, so a
    malformed ``path_pattern`` (unknown placeholder) raises ``ValueError``
    here rather than reaching the executor.
    """
    return DataArtifact(
        artifact_id=str(data["artifact_id"]),
        path_pattern=str(data["path_pattern"]),
        format=str(data["format"]),
        variables=tuple(str(v) for v in data.get("variables", ())),
        description=str(data.get("description", "")),
        produced_by=str(data.get("produced_by", "")),
        optional=bool(data.get("optional", False)),
        time_indexed=bool(data.get("time_indexed", False)),
    )


@dataclass(frozen=True)
class CaseConfig:
    """A single simulation configuration inside a tutorial sweep."""

    case_id: str
    params: dict[str, Any]


BuildCasesFn = Callable[[], list[CaseConfig]]
ApplyCaseFn = Callable[[Path, CaseConfig], None]


@dataclass(frozen=True)
class TutorialSpec:
    """Full tutorial definition consumed by the driver engine."""

    name: str
    case_root: Path
    setup_root: Path
    output_dir: Path
    build_cases: BuildCasesFn
    apply_case: ApplyCaseFn
    metadata: dict[str, Any] = field(default_factory=dict)
