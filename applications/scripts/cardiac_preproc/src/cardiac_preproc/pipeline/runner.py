"""Pipeline step runner."""

from __future__ import annotations

from cardiac_preproc.pipeline.context import StepContext, StepResult
from cardiac_preproc.pipeline.registry import get_step


def run_registered_step(name: str, context: StepContext) -> StepResult:
    runner = get_step(name)
    return runner(context)
