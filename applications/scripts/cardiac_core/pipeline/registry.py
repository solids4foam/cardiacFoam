"""Registry for internal pipeline steps."""

from __future__ import annotations

from collections.abc import Callable

from cardiac_core.pipeline.context import StepContext, StepResult

StepRunner = Callable[[StepContext], StepResult]

_REGISTRY: dict[str, StepRunner] = {}


def register_step(name: str, runner: StepRunner) -> None:
    _REGISTRY[name] = runner


def get_step(name: str) -> StepRunner:
    if name not in _REGISTRY:
        raise KeyError(f"Step '{name}' is not registered.")
    return _REGISTRY[name]


def list_steps() -> tuple[str, ...]:
    return tuple(sorted(_REGISTRY.keys()))
