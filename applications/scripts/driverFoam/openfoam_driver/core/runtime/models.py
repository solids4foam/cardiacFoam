from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable


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
