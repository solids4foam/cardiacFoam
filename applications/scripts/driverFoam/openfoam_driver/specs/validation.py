"""Run-document validator.

The validator reports three kinds of issue:

1. **Required-field omissions.** Each ``DictEntry`` flagged ``required=True``
   must have a value in the Run document slice owned by its *primary* phase.
2. **Enum violations.** Entries with ``value_kind="enum"`` whose value is not
   one of the declared ``enum_values``.
3. **Cross-field constraints.** A small set of v1 sanity checks (currently:
   ``eikonalSolver`` is incompatible with an explicit ``ionicModel``).

The *primary* phase is the first phase in workflow order
(``anatomy → physics → stimulus → solver``) that the entry claims.
Multi-phase entries are validated there; the other phases do not duplicate
validation errors.
"""

from __future__ import annotations

from dataclasses import dataclass

from openfoam_driver.dict_entries import (
    ELECTRO_PROPERTY_ENTRY_GROUPS,
    PHYSICS_PROPERTY_ENTRIES,
    Phase,
)

_PHASE_ORDER: tuple[Phase, ...] = (
    "anatomy", "physics", "stimulus", "solver",
)


@dataclass(frozen=True)
class ValidationError:
    phase: Phase
    field: str
    message: str
    level: str  # "error" | "warning"


def _all_entries():
    yield from PHYSICS_PROPERTY_ENTRIES
    for group in ELECTRO_PROPERTY_ENTRY_GROUPS.values():
        yield from group


def primary_phase(entry) -> Phase | None:
    """Return the editing phase for a (possibly multi-phase) entry.

    Walks ``_PHASE_ORDER`` and returns the first phase the entry claims;
    every other declared phase is a read-only mirror.
    """
    for ph in _PHASE_ORDER:
        if ph in entry.phases:
            return ph
    return None


def _slice_value(run, phase: Phase, driver_path: str):
    """Look up a leaf-name key inside a phase slice of the Run config."""
    slice_ = run.config.get(phase, {}) or {}
    key = driver_path.split(".")[-1]
    return slice_.get(key)


def validate_run(run) -> list[ValidationError]:
    errors: list[ValidationError] = []

    # 1) Required-field checks.
    for e in _all_entries():
        if not e.required:
            continue
        ph = primary_phase(e)
        if ph is None:
            continue
        val = _slice_value(run, ph, e.driver_path)
        if val in (None, ""):
            errors.append(ValidationError(
                phase=ph,
                field=e.driver_path,
                message=f"{e.driver_path} is required.",
                level="error",
            ))

    # 2) Enum checks. A missing value is already covered by (1) when required;
    #    skip otherwise so optional enums don't fire when unset.
    for e in _all_entries():
        if e.value_kind != "enum" or not e.enum_values:
            continue
        ph = primary_phase(e)
        if ph is None:
            continue
        val = _slice_value(run, ph, e.driver_path)
        if val is None or val == "":
            continue
        if val not in e.enum_values:
            errors.append(ValidationError(
                phase=ph,
                field=e.driver_path,
                message=f"{val!r} is not one of {list(e.enum_values)}.",
                level="error",
            ))

    # 3) Cross-field constraints (conservative — only well-known ones for v1).
    physics = run.config.get("physics", {}) or {}
    if (
        physics.get("myocardiumSolver") == "eikonalSolver"
        and physics.get("ionicModel")
    ):
        errors.append(ValidationError(
            phase="physics",
            field="ionicModel",
            message=(
                "ionicModel is not applicable when "
                "myocardiumSolver=eikonalSolver."
            ),
            level="error",
        ))

    return errors
