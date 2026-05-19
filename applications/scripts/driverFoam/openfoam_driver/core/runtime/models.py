from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable, Literal


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

    Shared vocabulary between the engine (which writes ``artifacts_manifest.json``
    listing what a run actually produced) and ``utility.manifest.toml`` ``produces``
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


@dataclass(frozen=True)
class CaseConfig:
    """A single simulation configuration inside a tutorial sweep."""

    case_id: str
    params: dict[str, Any]


BuildCasesFn = Callable[[], list[CaseConfig]]
ApplyCaseFn = Callable[[Path, CaseConfig], None]
RunCaseFn = Callable[[Path, Path, CaseConfig], None]
CollectOutputsFn = Callable[[Path, Path], None]
PostprocessFn = Callable[[Path, Path], None]


@dataclass(frozen=True)
class TutorialSpec:
    """Full tutorial definition consumed by the driver engine."""

    name: str
    case_root: Path
    setup_root: Path
    output_dir: Path
    build_cases: BuildCasesFn
    apply_case: ApplyCaseFn
    run_case: RunCaseFn
    collect_outputs: CollectOutputsFn | None = None
    postprocess: PostprocessFn | None = None
    metadata: dict[str, Any] = field(default_factory=dict)
