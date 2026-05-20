"""Run-document validator.

The validator reports three kinds of issue:

1. **Required-field omissions.** Each ``DictEntry`` flagged ``required=True``
   must have a value in the Run document slice owned by its *primary* phase.
2. **Enum violations.** Entries with ``value_kind="enum"`` whose value is not
   one of the declared ``enum_values``.
3. **Structured constraints (P5c/P5b).** Each ``DictEntry`` may declare
   ``applicable_when``, ``forbidden_when``, ``required_when``, and
   ``mutually_exclusive_with``. The validator evaluates these against a
   flattened view of the run config; rule families and semantics are
   documented in plan §5.1. The legacy hardcoded ``eikonalSolver``/
   ``ionicModel`` cross-field check was removed in P5b — the
   ``ionicModel`` entry now carries ``forbidden_when={"myocardiumSolver":
   "eikonalSolver"}`` which is evaluated programmatically here.

The *primary* phase is the first phase in workflow order
(``anatomy → physics → stimulus → solver``) that the entry claims.
Multi-phase entries are validated there; the other phases do not duplicate
validation errors.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Iterable

from openfoam_driver.dict_entries import (
    DictEntry,
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


_COEFFS_PREFIX = "$ELECTRO_MODEL_COEFFS."


def slot_key(driver_path: str) -> str:
    """Map a driver_path to its slot key inside a phase slice.

    Strips the ``$ELECTRO_MODEL_COEFFS.`` prefix when present; otherwise
    returns the path as-is. Multi-segment unprefixed paths are kept intact
    so that nested-group leaves don't collide with top-level keys of the
    same name (e.g. ``$ELECTRO_MODEL_COEFFS.bathPotentialDomain.phiEReferenceValue`` must not overwrite the
    top-level ``type`` entry inside the physics slice).
    """
    if driver_path.startswith(_COEFFS_PREFIX):
        return driver_path[len(_COEFFS_PREFIX):]
    return driver_path


def _slice_value(run, phase: Phase, driver_path: str):
    """Look up the slot value for a driver_path inside a phase slice."""
    slice_ = run.config.get(phase, {}) or {}
    return slice_.get(slot_key(driver_path))


def _flatten_context(run) -> dict[str, Any]:
    """Build a flat predicate-key → value view across every phase slice.

    Structured-constraint predicates reference slot-keys (the post-
    ``slot_key`` form: ``myocardiumSolver``, ``ionicModel``,
    ``singleCellStimulus.stim_amplitude``, ...). Multiple phases never
    write the same slot-key in practice; if they do, the last wins —
    document the convention rather than silently merging.
    """
    context: dict[str, Any] = {}
    for slice_ in run.config.values():
        if not slice_:
            continue
        for key, val in slice_.items():
            if val not in (None, ""):
                context[key] = val
    return context


def _predicate_matches(
    context: dict[str, Any],
    key: str,
    expected: str | tuple[str, ...],
) -> bool:
    """Return True iff ``context[key]`` matches ``expected``.

    Scalar ``expected`` → equality. Tuple ``expected`` → membership.
    A missing key is treated as not-matching (the predicate's
    precondition is absent).
    """
    if key not in context:
        return False
    actual = context[key]
    if isinstance(expected, tuple):
        return actual in expected
    return actual == expected


def _entry_is_applicable(entry: DictEntry, context: dict[str, Any]) -> bool:
    """Evaluate ``applicable_when`` — every predicate must match (AND)."""
    if not entry.applicable_when:
        return True
    return all(
        _predicate_matches(context, key, expected)
        for key, expected in entry.applicable_when.items()
    )


def _entry_value_present(entry: DictEntry, context: dict[str, Any]) -> bool:
    """Is the entry's own slot set in the flattened context?"""
    key = slot_key(entry.driver_path)
    return key in context and context[key] not in (None, "")


def _format_predicate(predicate: dict[str, Any]) -> str:
    parts: list[str] = []
    for key, expected in predicate.items():
        if isinstance(expected, tuple):
            parts.append(f"{key} ∈ {{{', '.join(expected)}}}")
        else:
            parts.append(f"{key}={expected}")
    return " and ".join(parts)


def _all_entries_list():
    return list(_all_entries())


def validate_run(
    run,
    *,
    entries: Iterable[DictEntry] | None = None,
) -> list[ValidationError]:
    """Validate ``run`` against the dict-entry catalog.

    ``entries`` overrides the live catalog for testability and for callers
    that want to validate against a curated subset (e.g., dict_builder).
    When omitted, the full live catalog is used.
    """
    entry_list: list[DictEntry] = (
        list(entries) if entries is not None else _all_entries_list()
    )
    context = _flatten_context(run)
    errors: list[ValidationError] = []

    # 1) Required-field checks. Skip entries whose applicable_when fails —
    #    requiredness is conditional on applicability.
    for e in entry_list:
        if not e.required:
            continue
        if not _entry_is_applicable(e, context):
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

    # 2) Enum checks. Skip inapplicable entries for the same reason.
    for e in entry_list:
        if e.value_kind != "enum" or not e.enum_values:
            continue
        if not _entry_is_applicable(e, context):
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

    # 3) Structured constraints (P5c/P5b).
    # (Legacy hardcoded eikonalSolver+ionicModel check removed in P5b — the
    # ionicModel entry now carries forbidden_when={"myocardiumSolver": "eikonalSolver"}
    # which the section below evaluates programmatically.)
    errors.extend(_evaluate_structured(entry_list, context))

    return errors


def _evaluate_structured(
    entries: list[DictEntry],
    context: dict[str, Any],
) -> list[ValidationError]:
    """Evaluate the four structured-constraint families per entry."""
    errors: list[ValidationError] = []
    paths_set = {slot_key(e.driver_path) for e in entries
                 if _entry_value_present(e, context)}

    for e in entries:
        # Skip entries whose applicable_when precondition fails entirely.
        if not _entry_is_applicable(e, context):
            continue
        ph = primary_phase(e) or "physics"

        # forbidden_when: fires when ANY predicate matches AND the entry's
        # own slot has a value. Each matching predicate emits its own
        # ValidationError so the reason text stays specific.
        if _entry_value_present(e, context):
            for key, expected in e.forbidden_when.items():
                if _predicate_matches(context, key, expected):
                    errors.append(ValidationError(
                        phase=ph,
                        field=e.driver_path,
                        message=(
                            f"{e.driver_path} is forbidden when "
                            f"{_format_predicate({key: expected})}."
                        ),
                        level="error",
                    ))

        # required_when: fires when ANY predicate matches AND the entry's
        # slot is unset.
        if not _entry_value_present(e, context):
            for key, expected in e.required_when.items():
                if _predicate_matches(context, key, expected):
                    errors.append(ValidationError(
                        phase=ph,
                        field=e.driver_path,
                        message=(
                            f"{e.driver_path} is required when "
                            f"{_format_predicate({key: expected})}."
                        ),
                        level="error",
                    ))

        # mutually_exclusive_with: fires when BOTH this entry's slot is set
        # AND any of the listed sibling slots is also set. To avoid
        # double-reporting symmetric relations we only flag the side that
        # *declares* the relation.
        if _entry_value_present(e, context):
            for sibling_path in e.mutually_exclusive_with:
                sibling_slot = slot_key(sibling_path)
                if sibling_slot in paths_set:
                    errors.append(ValidationError(
                        phase=ph,
                        field=e.driver_path,
                        message=(
                            f"{e.driver_path} is mutually exclusive with "
                            f"{sibling_path}."
                        ),
                        level="error",
                    ))

    return errors
