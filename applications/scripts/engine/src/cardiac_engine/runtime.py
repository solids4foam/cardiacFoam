"""Shared high-level engine facade for step-based preprocessing."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

from cardiac_preproc.pipeline import register_default_steps
from cardiac_preproc.pipeline.context import StepContext, StepResult
from cardiac_preproc.pipeline.runner import run_registered_step


@dataclass
class CardiacPreprocEngine:
    """Facade for running registered preprocessing steps."""

    project_root: Path
    output_dir: Path | None = None

    def __post_init__(self) -> None:
        register_default_steps()

    def run_step(
        self,
        step_name: str,
        *,
        current_mesh: str | None,
        options: dict[str, Any],
    ) -> StepResult:
        return run_registered_step(
            step_name,
            StepContext(
                project_root=self.project_root,
                current_mesh=current_mesh,
                output_dir=self.output_dir,
                options=options,
            ),
        )

    def run_pipeline(
        self,
        *,
        current_mesh: str | None,
        steps: list[tuple[str, dict[str, Any]]],
    ) -> StepResult:
        last = StepResult(output_mesh=current_mesh)
        mesh = current_mesh
        for step_name, options in steps:
            last = self.run_step(step_name, current_mesh=mesh, options=options)
            mesh = last.output_mesh or mesh
        return last
