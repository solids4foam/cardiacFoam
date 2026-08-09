"""Solver-neutral fallback for an existing OpenFOAM case folder.

This deliberately models only the execution contract common to driverFOAM:
one on-disk case directory and an ``Allrun`` workflow fallback.  Solver
plugins may provide richer materialisation, dictionary mutation, or
post-processing, but filesystem discovery must never need those services.
"""

from __future__ import annotations

from pathlib import Path

from ...specs.common import tutorials_root_default
from .models import CaseConfig, TutorialSpec


def _no_case_mutation(_case_root: Path, _case: CaseConfig) -> None:
    """Generic case folders are not mutated unless a caller requests it explicitly."""


def _run_is_owned_by_workflow(
    _case_root: Path,
    _setup_root: Path,
    _case: CaseConfig,
) -> None:
    """Compatibility callback; strict execution invokes the normalized DAG."""


def make_generic_case_spec(
    *,
    tutorials_root: Path | None = None,
    case_dir_name: str,
    output_dir_name: str = "postProcessing",
) -> TutorialSpec:
    """Create the minimal spec for an existing case folder.

    ``case_dir_name`` is kept relative to ``tutorials_root`` so the registry
    controls containment.  The generated DAG is intentionally a compact form;
    normalisation supplies argv/cwd/retry metadata before execution.
    """

    if not str(case_dir_name).strip():
        raise ValueError("case_dir_name cannot be empty")
    root = Path(tutorials_root) if tutorials_root is not None else tutorials_root_default()
    case_root = root / case_dir_name
    return TutorialSpec(
        name=str(case_dir_name),
        case_root=case_root,
        setup_root=case_root,
        output_dir=case_root / output_dir_name,
        build_cases=lambda: [CaseConfig(case_id="default", params={})],
        apply_case=_no_case_mutation,
        run_case=_run_is_owned_by_workflow,
        metadata={
            "notes": "Core generic case-folder runner.",
            "workflow_dag": {"steps": [{"id": "run", "command": "Allrun", "depends_on": []}]},
            "generic_case": True,
        },
    )
